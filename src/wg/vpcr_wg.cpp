// BRANCH v0.5 — Phase 3 Stage 2 implementation.
//
// design()  : per CandidateRegion, find up to N primer pairs in the
//             low-entropy windows that span the region. Tag the
//             amplicon with the region's OR'd flag set.
//
// count()   : k-mer-seed inverted index. Build {seed_bits → list of
//             (amp_id, primer_offset, side)} once. Per read: slide a
//             canonical seed window, on each hit verify the full
//             primer (≤ max_mm mismatches), accumulate fwd/rev hit
//             positions per amplicon. Emit a count when both fwd and
//             rev hits land in PCR-able orientation + distance.
//
// Cost: O(reads × read_len) for the scan pass, independent of amplicon
// count. OpenMP-parallel over reads with thread-local hit accumulators
// merged at the end. Identical algorithm runs on GPU via
// wg_stage2_kernel.cu when CUDA is enabled.

#include "wg/vpcr_wg.hpp"
#include "vpcr/vpcr.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace branch::wg {

namespace {

inline std::uint8_t encode_base(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

inline std::uint64_t rc_kmer(std::uint64_t kmer, std::size_t k) noexcept {
    std::uint64_t rc = 0;
    for (std::size_t i = 0; i < k; ++i) {
        rc = (rc << 2) | ((~kmer) & 3ULL);
        kmer >>= 2;
    }
    return rc;
}

inline std::uint64_t canonical_bits(std::uint64_t kmer, std::size_t k) noexcept {
    auto rc = rc_kmer(kmer, k);
    return kmer < rc ? kmer : rc;
}

double dimer_entropy(std::string_view window) {
    if (window.size() < 2) return 0.0;
    std::array<std::uint32_t, 16> counts{};
    std::uint32_t total = 0;
    for (std::size_t i = 0; i + 1 < window.size(); ++i) {
        std::uint8_t a = encode_base(window[i]);
        std::uint8_t b = encode_base(window[i + 1]);
        if (a > 3 || b > 3) continue;
        ++counts[a * 4 + b];
        ++total;
    }
    if (total == 0) return 0.0;
    double H = 0.0;
    for (std::uint32_t c : counts) {
        if (c == 0) continue;
        double p = static_cast<double>(c) / total;
        H -= p * std::log2(p);
    }
    return H / 4.0;
}

std::pair<double, double> wilson_ci(std::uint32_t count, std::uint32_t expected_single_copy) {
    if (expected_single_copy == 0) return {0.0, 0.0};
    double n = static_cast<double>(count);
    double e = static_cast<double>(expected_single_copy);
    double mean = n / e;
    double sigma = std::sqrt(std::max(n, 1.0)) / e;
    return { std::max(0.0, mean - 1.96 * sigma), mean + 1.96 * sigma };
}

inline int count_mm(const char* text, std::size_t text_len, std::size_t pos,
                    const std::string& pattern, int max_mm) {
    if (pos + pattern.size() > text_len) return max_mm + 1;
    int mm = 0;
    for (std::size_t i = 0; i < pattern.size(); ++i) {
        char a = text[pos + i];
        char b = pattern[i];
        if (a != b && (std::toupper(static_cast<unsigned char>(a))
                    != std::toupper(static_cast<unsigned char>(b)))) {
            if (++mm > max_mm) return mm;
        }
    }
    return mm;
}

}  // namespace

VpcrAutoDesigner::VpcrAutoDesigner(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void VpcrAutoDesigner::design(
    const std::vector<std::string>& hap_seqs,
    const std::vector<CandidateRegion>& regions,
    std::vector<VpcrAmplicon>& out,
    std::uint32_t n_per_region) const {
    out.clear();
    if (regions.empty() || hap_seqs.empty()) return;

    const int wmin = profile_.primer_min_window_bp;
    const int wmax = profile_.primer_max_window_bp;
    const double max_ent = profile_.primer_max_entropy;
    const int amp_min = 200;
    const int amp_max = (profile_.name && profile_.name[0] == 'h') ? 1500 : 3000;
    const int wuse = (wmin + wmax) / 2;

    std::uint64_t aidx = 0;
    for (const auto& r : regions) {
        if (r.hap_idx >= hap_seqs.size()) continue;
        const auto& seq = hap_seqs[r.hap_idx];
        if (r.end_bp > seq.size()) continue;

        // Mark which positions in this region have low-entropy windows.
        std::size_t span = r.end_bp - r.start_bp;
        if (span < static_cast<std::size_t>(amp_min + 2 * wuse)) continue;
        std::vector<std::uint8_t> ok(span + 1, 0);
        for (std::size_t i = 0; i + wuse <= span; ++i) {
            double H = dimer_entropy(std::string_view(&seq[r.start_bp + i], wuse));
            if (H <= max_ent) ok[i] = 1;
        }

        std::uint32_t emitted = 0;
        for (std::size_t i = 0; i + wuse + amp_min <= span && emitted < n_per_region; ++i) {
            if (!ok[i]) continue;
            std::size_t lo = i + wuse + amp_min;
            std::size_t hi = std::min<std::size_t>(span - wuse,
                                                   i + wuse + static_cast<std::size_t>(amp_max));
            for (std::size_t j = lo; j < hi; ++j) {
                if (!ok[j]) continue;
                VpcrAmplicon a{};
                std::ostringstream name;
                name << "amp_h" << (r.hap_idx + 1)
                     << "_r" << r.start_bp << "_" << (++aidx);
                a.amplicon_name = name.str();
                a.fwd_primer.assign(seq, r.start_bp + i, wuse);
                std::string rev_window(seq, r.start_bp + j, wuse);
                a.rev_primer = branch::vpcr::reverse_complement(rev_window);
                a.hap_idx = r.hap_idx;
                a.start_bp = static_cast<std::uint32_t>(r.start_bp + i);
                a.expected_size_bp = static_cast<std::uint32_t>(j + wuse - i);
                a.region_flags = r.flags;
                out.push_back(std::move(a));
                ++emitted;
                i += wuse - 1;
                break;
            }
        }
    }
}

void VpcrAutoDesigner::count(
    std::vector<VpcrAmplicon>& amplicons,
    const std::vector<std::pair<std::string, std::string>>& reads,
    int expected_single_copy_count) const {
    if (amplicons.empty() || reads.empty()) return;

    constexpr std::size_t kSeed = 12;
    const int max_mm = 2;
    const int amp_min = 200;
    const int amp_max = (profile_.name && profile_.name[0] == 'h') ? 1500 : 3000;

    struct SeedHit { std::uint32_t amp_id; std::uint16_t off; };
    std::unordered_map<std::uint64_t, std::vector<SeedHit>> fwd_idx, rev_idx;
    fwd_idx.reserve(amplicons.size() * 2);
    rev_idx.reserve(amplicons.size() * 2);

    auto add_seeds = [&](const std::string& primer, std::uint32_t aid,
                         std::unordered_map<std::uint64_t, std::vector<SeedHit>>& idx) {
        if (primer.size() < kSeed) return;
        const std::uint64_t mask = (1ULL << (2 * kSeed)) - 1ULL;
        std::uint64_t k = 0;
        std::size_t filled = 0;
        for (std::size_t i = 0; i < primer.size(); ++i) {
            std::uint8_t b = encode_base(primer[i]);
            if (b > 3) { filled = 0; k = 0; continue; }
            k = ((k << 2) | b) & mask;
            ++filled;
            if (filled < kSeed) continue;
            std::uint64_t cb = canonical_bits(k, kSeed);
            std::uint16_t off = static_cast<std::uint16_t>(i + 1 - kSeed);
            idx[cb].push_back(SeedHit{aid, off});
        }
    };
    for (std::uint32_t a = 0; a < amplicons.size(); ++a) {
        add_seeds(amplicons[a].fwd_primer, a, fwd_idx);
        add_seeds(amplicons[a].rev_primer, a, rev_idx);
    }

    constexpr std::size_t kSeedBucketCap = 64;
    auto prune = [](std::unordered_map<std::uint64_t, std::vector<SeedHit>>& idx) {
        for (auto it = idx.begin(); it != idx.end(); ) {
            if (it->second.size() > kSeedBucketCap) it = idx.erase(it);
            else ++it;
        }
    };
    prune(fwd_idx);
    prune(rev_idx);

    std::vector<std::uint32_t> hit_count(amplicons.size(), 0);

    auto scan_one_read = [&](const std::string& seq,
                             std::vector<std::uint32_t>& local_count) {
        if (seq.size() < kSeed) return;
        const std::uint64_t mask = (1ULL << (2 * kSeed)) - 1ULL;
        struct HR {
            std::uint32_t fwd = std::numeric_limits<std::uint32_t>::max();
            std::uint32_t rev = std::numeric_limits<std::uint32_t>::max();
        };
        std::unordered_map<std::uint32_t, HR> per_amp;

        std::uint64_t k = 0;
        std::size_t filled = 0;
        for (std::size_t i = 0; i < seq.size(); ++i) {
            std::uint8_t b = encode_base(seq[i]);
            if (b > 3) { filled = 0; k = 0; continue; }
            k = ((k << 2) | b) & mask;
            ++filled;
            if (filled < kSeed) continue;
            std::uint64_t cb = canonical_bits(k, kSeed);
            std::size_t kpos = i + 1 - kSeed;

            auto fit = fwd_idx.find(cb);
            if (fit != fwd_idx.end()) {
                for (const auto& sh : fit->second) {
                    const auto& primer = amplicons[sh.amp_id].fwd_primer;
                    if (kpos < sh.off) continue;
                    std::size_t pos = kpos - sh.off;
                    if (count_mm(seq.data(), seq.size(), pos, primer, max_mm) <= max_mm) {
                        auto& hr = per_amp[sh.amp_id];
                        if (pos < hr.fwd) hr.fwd = static_cast<std::uint32_t>(pos);
                    }
                }
            }
            auto rit = rev_idx.find(cb);
            if (rit != rev_idx.end()) {
                for (const auto& sh : rit->second) {
                    const auto& primer = amplicons[sh.amp_id].rev_primer;
                    if (kpos < sh.off) continue;
                    std::size_t pos = kpos - sh.off;
                    if (count_mm(seq.data(), seq.size(), pos, primer, max_mm) <= max_mm) {
                        auto& hr = per_amp[sh.amp_id];
                        if (pos > hr.rev ||
                            hr.rev == std::numeric_limits<std::uint32_t>::max())
                            hr.rev = static_cast<std::uint32_t>(pos);
                    }
                }
            }
        }
        for (const auto& [amp_id, hr] : per_amp) {
            if (hr.fwd == std::numeric_limits<std::uint32_t>::max()) continue;
            if (hr.rev == std::numeric_limits<std::uint32_t>::max()) continue;
            if (hr.rev <= hr.fwd) continue;
            std::uint32_t span = hr.rev
                + static_cast<std::uint32_t>(amplicons[amp_id].rev_primer.size())
                - hr.fwd;
            if (span < static_cast<std::uint32_t>(amp_min)) continue;
            if (span > static_cast<std::uint32_t>(amp_max)) continue;
            ++local_count[amp_id];
        }
    };

#ifdef _OPENMP
    int n_threads = omp_get_max_threads();
#else
    int n_threads = 1;
#endif
    std::vector<std::vector<std::uint32_t>> per_thread(n_threads,
        std::vector<std::uint32_t>(amplicons.size(), 0));

    // Chunked outer loop so we get a heartbeat between chunks. The
    // per-chunk barrier costs ~1% throughput but gives live progress
    // visibility — vital for whole-genome runs that take hours.
    const std::size_t kChunk = 200'000;
    const auto t_start = std::chrono::steady_clock::now();
    for (std::size_t off = 0; off < reads.size(); off += kChunk) {
        const std::size_t end = std::min(off + kChunk, reads.size());
        #pragma omp parallel for schedule(dynamic, 64)
        for (std::size_t i = off; i < end; ++i) {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            scan_one_read(reads[i].second, per_thread[tid]);
        }
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();
        const double frac = static_cast<double>(end)
                          / static_cast<double>(reads.size());
        const double eta = el / std::max(frac, 1e-6) - el;
        std::fprintf(stderr,
            "[wg/vpcr] count %zu/%zu reads (%.1f%%) "
            "elapsed=%.0fs ETA=%.0fs\n",
            end, reads.size(), frac * 100.0, el, eta);
    }

    for (auto& tc : per_thread) {
        for (std::size_t a = 0; a < amplicons.size(); ++a) hit_count[a] += tc[a];
    }

    for (std::size_t a = 0; a < amplicons.size(); ++a) {
        amplicons[a].read_count = hit_count[a];
        if (expected_single_copy_count > 0) {
            amplicons[a].cn_estimate =
                static_cast<double>(hit_count[a]) / expected_single_copy_count;
            auto ci = wilson_ci(hit_count[a],
                static_cast<std::uint32_t>(expected_single_copy_count));
            amplicons[a].ci_low = ci.first;
            amplicons[a].ci_high = ci.second;
        }
    }
}

}  // namespace branch::wg
