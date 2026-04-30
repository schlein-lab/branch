// BRANCH v0.5 — Phase 3 implementation.
//
// Auto-primer design:
//   - For each haplotype, slide a primer-window (k_primer ∈
//     [primer_min_window_bp, primer_max_window_bp]) and compute
//     2-mer entropy. Windows below primer_max_entropy are "low-entropy"
//     (= conserved) and qualify as primer candidates.
//   - Pair two primers at PCR-able distance (200-3000 bp ONT,
//     200-1500 bp HiFi) with the second on the reverse complement strand.
//
// Counting: delegate to branch::vpcr::scan_amplicons (existing module
// in src/vpcr/), then estimate CN = read_count / expected_single_copy
// with a Wilson 95% CI.

#include "wg/vpcr_wg.hpp"
#include "vpcr/vpcr.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <sstream>

namespace branch::wg {

namespace {

double dimer_entropy(std::string_view window) {
    if (window.size() < 2) return 0.0;
    std::array<std::uint32_t, 16> counts{};
    auto code = [](char c) -> int {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return -1;
        }
    };
    std::uint32_t total = 0;
    for (std::size_t i = 0; i + 1 < window.size(); ++i) {
        int a = code(window[i]);
        int b = code(window[i + 1]);
        if (a < 0 || b < 0) continue;
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
    return H / 4.0;  // normalize to [0,1] (max H = 4 for uniform 16 dimers)
}

// Wilson 95% CI for a Poisson-ish count. Returns {lo, hi} as a fraction of mean.
std::pair<double, double> wilson_ci(std::uint32_t count, std::uint32_t expected_single_copy) {
    if (expected_single_copy == 0) return {0.0, 0.0};
    double n = static_cast<double>(count);
    double e = static_cast<double>(expected_single_copy);
    // Approximate CI on count/e using sqrt(count) as Poisson sigma.
    double mean = n / e;
    double sigma = std::sqrt(std::max(n, 1.0)) / e;
    return { std::max(0.0, mean - 1.96 * sigma), mean + 1.96 * sigma };
}

}  // namespace

VpcrAutoDesigner::VpcrAutoDesigner(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void VpcrAutoDesigner::design(
    const std::vector<std::string>& hap_seqs,
    std::vector<VpcrAmplicon>& out) const {
    out.clear();
    const int wmin = profile_.primer_min_window_bp;
    const int wmax = profile_.primer_max_window_bp;
    const double max_ent = profile_.primer_max_entropy;
    // PCR-able distances by tech.
    const int amp_min = 200;
    const int amp_max = (profile_.name && profile_.name[0] == 'h') ? 1500 : 3000;
    const int wuse = (wmin + wmax) / 2;

    std::uint32_t aidx = 0;
    for (std::size_t h = 0; h < hap_seqs.size(); ++h) {
        const auto& seq = hap_seqs[h];
        if (seq.size() < static_cast<std::size_t>(amp_min + 2 * wuse)) continue;
        // Mark every window's entropy.
        std::vector<std::uint8_t> ok(seq.size(), 0);
        for (std::size_t i = 0; i + wuse <= seq.size(); ++i) {
            double H = dimer_entropy(std::string_view(&seq[i], wuse));
            if (H <= max_ent) ok[i] = 1;
        }
        // For each fwd primer candidate, look for a rev primer at right distance.
        // Bound: at most 50 amplicons per haplotype — vPCR counting is O(reads
        // * primer_len) per amplicon, so 100 total amplicons * 250k reads at
        // ~1500 bp/read is the runtime budget on a SLURM std-partition node.
        std::size_t emitted = 0;
        for (std::size_t i = 0; i + wuse + amp_min <= seq.size() && emitted < 50; ++i) {
            if (!ok[i]) continue;
            std::size_t lo = i + wuse + amp_min;
            std::size_t hi = std::min(seq.size() - wuse,
                                      i + wuse + static_cast<std::size_t>(amp_max));
            for (std::size_t j = lo; j < hi; ++j) {
                if (!ok[j]) continue;
                VpcrAmplicon a{};
                std::ostringstream name;
                name << "wg_amp_h" << (h + 1) << "_" << (++aidx);
                a.amplicon_name = name.str();
                a.fwd_primer.assign(seq, i, wuse);
                std::string rev_window(seq, j, wuse);
                a.rev_primer = branch::vpcr::reverse_complement(rev_window);
                a.expected_size_bp = static_cast<std::uint32_t>(j + wuse - i);
                out.push_back(std::move(a));
                ++emitted;
                // Skip a primer-window worth of positions to avoid near-dupes.
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
    // Allow up to 2 mismatches per primer (matches branch::vpcr default).
    const int max_mm = 2;

    // For very large read sets (whole-genome runs), cap the per-amplicon
    // scan at 250 k reads. Above that the count is a sample — we record
    // an effective expected_single_copy_count = expected * (sample / total)
    // so the CN estimate stays unbiased. Below that, the cap is a no-op.
    constexpr std::size_t kMaxReadsPerAmplicon = 250'000;
    std::size_t scan_n = std::min(reads.size(), kMaxReadsPerAmplicon);
    int eff_single_copy = expected_single_copy_count;
    if (scan_n < reads.size() && eff_single_copy > 0) {
        double scale = static_cast<double>(scan_n)
                       / static_cast<double>(reads.size());
        eff_single_copy = static_cast<int>(
            static_cast<double>(eff_single_copy) * scale);
        if (eff_single_copy < 1) eff_single_copy = 1;
    }
    for (auto& a : amplicons) {
        std::uint32_t hits = 0;
        branch::vpcr::PrimerPair pp{a.amplicon_name, a.fwd_primer, a.rev_primer};
        for (std::size_t i = 0; i < scan_n; ++i) {
            const auto& [id, seq] = reads[i];
            auto h = branch::vpcr::scan_amplicons(pp, id, seq, max_mm);
            if (!h.empty()) ++hits;
        }
        a.read_count = hits;
        if (eff_single_copy > 0) {
            a.cn_estimate = static_cast<double>(hits) / eff_single_copy;
            auto ci = wilson_ci(hits, static_cast<std::uint32_t>(eff_single_copy));
            a.ci_low = ci.first;
            a.ci_high = ci.second;
        }
    }
}

}  // namespace branch::wg
