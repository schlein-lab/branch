// BRANCH v0.5 — Phase 1.1 implementation.
//
// Pipeline:
//   1. load (canonical_bits, count) table from Phase 0 dump
//   2. find_het_pairs(): for each k-mer K with count in the het window,
//      enumerate the 3*k 1-bp-substitution variants, look up their
//      canonical forms in the counter, and emit a het-pair when both
//      counts plus their sum fit the expected diploid pattern.
//      Allele labels (0/1) are assigned deterministically — the
//      lex-smaller canonical bits get label 0.
//   3. set_anchor_from_reads(): pick the read with the highest het-pair
//      hit-count from a bounded scan; that read's local allele pattern
//      is stamped onto hap1_pattern_.
//   4. route_read(): match k-mers against the het-pair table; count
//      hap1-pattern matches vs mismatches; assign Hap1/Hap2/Ambig.

#include "wg/haplotype_router.hpp"
#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

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

}  // namespace

double ReadRouting::confidence() const noexcept {
    auto m = std::max(votes_hap1, votes_hap2);
    auto total = votes_hap1 + votes_hap2;
    if (total == 0) return 0.0;
    return static_cast<double>(m) / static_cast<double>(total + 1);
}

HaplotypeRouter::HaplotypeRouter(const ::branch::graph::TechProfile& profile,
                                 const std::string& kmer_hist_path,
                                 int expected_coverage,
                                 bool lazy_load_counter)
    : profile_(profile), expected_coverage_(expected_coverage),
      counter_path_(kmer_hist_path) {
    if (!lazy_load_counter) {
        load_counter_(kmer_hist_path);
        counter_loaded_ = true;
    }
}

void HaplotypeRouter::load_counter_(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open Phase 0 output: " + path);

    char magic[8];
    f.read(magic, 8);
    if (std::memcmp(magic, "KHIST05", 7) != 0)
        throw std::runtime_error("bad magic in Phase 0 file: " + path);

    std::uint64_t n_recurrent = 0, n_reads = 0, n_bases = 0;
    std::int32_t cov = 0, het = 0;
    double z = 0;
    std::uint32_t k = 0;
    f.read(reinterpret_cast<char*>(&n_recurrent), sizeof(n_recurrent));
    f.read(reinterpret_cast<char*>(&n_reads), sizeof(n_reads));
    f.read(reinterpret_cast<char*>(&n_bases), sizeof(n_bases));
    f.read(reinterpret_cast<char*>(&cov), sizeof(cov));
    f.read(reinterpret_cast<char*>(&het), sizeof(het));
    f.read(reinterpret_cast<char*>(&z), sizeof(z));
    f.read(reinterpret_cast<char*>(&k), sizeof(k));

    std::uint32_t hist_n = 0;
    f.read(reinterpret_cast<char*>(&hist_n), sizeof(hist_n));
    f.seekg(sizeof(std::uint64_t) * hist_n, std::ios::cur);

    counter_.reserve(static_cast<std::size_t>(n_recurrent));
    for (std::uint64_t i = 0; i < n_recurrent; ++i) {
        std::uint64_t key = 0;
        std::uint32_t cnt = 0;
        f.read(reinterpret_cast<char*>(&key), sizeof(key));
        f.read(reinterpret_cast<char*>(&cnt), sizeof(cnt));
        counter_.emplace(key, cnt);
    }
}

void HaplotypeRouter::find_het_pairs(double low_frac, double high_frac,
                                     double sum_frac_lo, double sum_frac_hi) {
    het_.clear();
    n_het_pairs_ = 0;
    // Lazy-load the counter on first CPU het-pair invocation. GPU path
    // never reaches this code (results piped via install_het_pairs_from_pairs),
    // so the 95 GB unordered_map blow-up is avoided.
    if (!counter_loaded_ && !counter_path_.empty()) {
        load_counter_(counter_path_);
        counter_loaded_ = true;
    }
    if (expected_coverage_ <= 0 || counter_.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::uint32_t lo = static_cast<std::uint32_t>(expected_coverage_ * low_frac);
    const std::uint32_t hi = static_cast<std::uint32_t>(expected_coverage_ * high_frac);
    const std::uint32_t sum_lo = static_cast<std::uint32_t>(expected_coverage_ * sum_frac_lo);
    const std::uint32_t sum_hi = static_cast<std::uint32_t>(expected_coverage_ * sum_frac_hi);
    if (lo == 0 || hi == 0) return;

    std::uint32_t pair_id = 0;
    for (const auto& [bits_a, count_a] : counter_) {
        if (count_a < lo || count_a > hi) continue;
        if (het_.find(bits_a) != het_.end()) continue;

        // Try all 3*k 1-bp-substitution variants; record the FIRST that
        // satisfies the het-frequency + diploid-sum constraints.
        bool paired = false;
        for (std::size_t pos = 0; pos < k && !paired; ++pos) {
            std::uint64_t old_base = (bits_a >> (2 * pos)) & 3ULL;
            std::uint64_t cleared = bits_a & ~(3ULL << (2 * pos));
            for (std::uint64_t alt = 0; alt < 4 && !paired; ++alt) {
                if (alt == old_base) continue;
                std::uint64_t fwd_bits = cleared | (alt << (2 * pos));
                std::uint64_t alt_canon = canonical_kmer_bits(fwd_bits, k);
                if (alt_canon == bits_a) continue;
                auto it = counter_.find(alt_canon);
                if (it == counter_.end()) continue;
                std::uint32_t count_b = it->second;
                if (count_b < lo || count_b > hi) continue;
                std::uint32_t sum = count_a + count_b;
                if (sum < sum_lo || sum > sum_hi) continue;
                if (het_.find(alt_canon) != het_.end()) continue;

                // Deterministic allele label: lex-smaller bits = 0.
                std::uint8_t label_a = (bits_a < alt_canon) ? 0 : 1;
                std::uint8_t label_b = label_a ^ 1;
                het_[bits_a]    = HetEntry{pair_id, label_a};
                het_[alt_canon] = HetEntry{pair_id, label_b};
                ++pair_id;
                paired = true;
            }
        }
    }
    n_het_pairs_ = pair_id;
}

void HaplotypeRouter::install_het_pairs_from_pairs(
    const std::vector<std::uint64_t>& a,
    const std::vector<std::uint64_t>& b) {
    het_.clear();
    n_het_pairs_ = 0;
    if (a.size() != b.size()) return;
    het_.reserve(a.size() * 2);
    for (std::uint32_t i = 0; i < a.size(); ++i) {
        std::uint64_t bits_a = a[i];
        std::uint64_t bits_b = b[i];
        if (bits_a == bits_b) continue;
        // Deterministic allele label: lex-smaller = 0.
        std::uint8_t la = bits_a < bits_b ? 0 : 1;
        std::uint8_t lb = la ^ 1;
        het_[bits_a] = HetEntry{i, la};
        het_[bits_b] = HetEntry{i, lb};
    }
    n_het_pairs_ = static_cast<std::size_t>(a.size());
}

void HaplotypeRouter::set_anchor_from_reads(
    const std::vector<std::pair<std::string, std::string>>& reads,
    std::size_t max_scan) {
    hap1_pattern_.clear();
    if (het_.empty() || reads.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);

    // Pass 1: find the read with the most het-pair hits within max_scan.
    std::size_t best_idx = 0;
    std::size_t best_hits = 0;
    std::size_t scanned = std::min(reads.size(), max_scan);
    auto t_anchor = std::chrono::steady_clock::now();
    std::fprintf(stderr,
        "[wg/anchor] scanning %zu reads for anchor candidate "
        "(het_ table size %zu)\n",
        scanned, het_.size());
    std::fflush(stderr);
    for (std::size_t r = 0; r < scanned; ++r) {
        if ((r + 1) % 10'000 == 0 || r + 1 == scanned) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_anchor).count();
            const double frac = static_cast<double>(r + 1)
                              / static_cast<double>(scanned);
            const double eta = el / std::max(frac, 1e-6) - el;
            std::fprintf(stderr,
                "[wg/anchor] %zu/%zu (%.1f%%) best_hits=%zu "
                "elapsed=%.0fs ETA=%.0fs\n",
                r + 1, scanned, frac * 100.0, best_hits, el, eta);
            std::fflush(stderr);
        }
        const auto& seq = reads[r].second;
        if (seq.size() < k) continue;
        std::uint64_t kmer = 0;
        std::size_t filled = 0;
        std::size_t hits = 0;
        for (char c : seq) {
            std::uint8_t b = encode_base(c);
            if (b > 3) { filled = 0; kmer = 0; continue; }
            kmer = ((kmer << 2) | b) & mask;
            ++filled;
            if (filled < k) continue;
            std::uint64_t cb = canonical_kmer_bits(kmer, k);
            if (het_.find(cb) != het_.end()) ++hits;
        }
        if (hits > best_hits) { best_hits = hits; best_idx = r; }
    }
    if (best_hits == 0) return;

    // Pass 2: stamp the anchor read's local allele pattern onto hap1.
    // For each het-pair hit, record the observed allele as the hap1 label.
    const auto& seq = reads[best_idx].second;
    std::uint64_t kmer = 0;
    std::size_t filled = 0;
    for (char c : seq) {
        std::uint8_t b = encode_base(c);
        if (b > 3) { filled = 0; kmer = 0; continue; }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;
        std::uint64_t cb = canonical_kmer_bits(kmer, k);
        auto it = het_.find(cb);
        if (it == het_.end()) continue;
        // Only set the first observation per pair (avoid noise from
        // sequencing errors flipping interpretation later in the read).
        hap1_pattern_.try_emplace(it->second.pair_id, it->second.allele);
    }
}

ReadRouting HaplotypeRouter::route_read(std::string_view seq) const noexcept {
    ReadRouting r;
    if (het_.empty() || hap1_pattern_.empty()) return r;
    const std::size_t k = profile_.minimizer_k;
    if (seq.size() < k) return r;

    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t kmer = 0;
    std::size_t filled = 0;
    std::uint32_t match = 0, mismatch = 0;

    for (std::size_t i = 0; i < seq.size(); ++i) {
        std::uint8_t b = encode_base(seq[i]);
        if (b > 3) { filled = 0; kmer = 0; continue; }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;
        std::uint64_t cb = canonical_kmer_bits(kmer, k);
        auto eit = het_.find(cb);
        if (eit == het_.end()) continue;
        auto pit = hap1_pattern_.find(eit->second.pair_id);
        if (pit == hap1_pattern_.end()) continue;
        if (pit->second == eit->second.allele) ++match;
        else                                   ++mismatch;
    }

    r.votes_hap1 = match;
    r.votes_hap2 = mismatch;
    if (match == 0 && mismatch == 0) {
        r.hap = Haplotype::Ambig;
    } else if (match > mismatch) {
        r.hap = Haplotype::Hap1;
    } else if (mismatch > match) {
        r.hap = Haplotype::Hap2;
    } else {
        r.hap = Haplotype::Ambig;
    }
    return r;
}

void HaplotypeRouter::write_assignments_tsv(
    const std::vector<std::pair<std::string, std::string>>& reads,
    const std::string& out_path,
    double min_confidence) const {
    std::ofstream f(out_path);
    if (!f) throw std::runtime_error("cannot open: " + out_path);
    f << "read_id\thap\tvotes_hap1\tvotes_hap2\tconfidence\n";
    for (const auto& [id, seq] : reads) {
        auto r = route_read(seq);
        const char* hap = "ambig";
        if (r.confidence() >= min_confidence) {
            switch (r.hap) {
                case Haplotype::Hap1: hap = "hap1"; break;
                case Haplotype::Hap2: hap = "hap2"; break;
                default: hap = "ambig"; break;
            }
        }
        f << id << '\t' << hap << '\t'
          << r.votes_hap1 << '\t' << r.votes_hap2 << '\t'
          << r.confidence() << '\n';
    }
}

std::size_t HaplotypeRouter::bytes() const noexcept {
    return counter_.size() * (sizeof(std::uint64_t) + sizeof(std::uint32_t))
         + het_.size() * (sizeof(std::uint64_t) + sizeof(HetEntry))
         + hap1_pattern_.size() * (sizeof(std::uint32_t) + sizeof(std::uint8_t));
}

std::vector<std::pair<std::uint64_t, std::uint64_t>>
HaplotypeRouter::export_het_pair_kmers() const {
    std::vector<std::pair<std::uint64_t, std::uint64_t>> out(n_het_pairs_,
        {0ULL, 0ULL});
    for (const auto& [bits, ent] : het_) {
        if (ent.pair_id >= n_het_pairs_) continue;
        if (ent.allele == 0) out[ent.pair_id].first  = bits;
        else                 out[ent.pair_id].second = bits;
    }
    // Drop any pair where one side is still 0 (defensive — shouldn't
    // happen since install_/find_ both set both sides).
    out.erase(std::remove_if(out.begin(), out.end(),
        [](const auto& p) { return p.first == 0 || p.second == 0; }),
        out.end());
    return out;
}

}  // namespace branch::wg
