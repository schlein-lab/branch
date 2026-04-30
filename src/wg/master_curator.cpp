// BRANCH v0.5 — Phase 1.5 implementation.
//
// Anchor pile-up algorithm:
//   1. compute the haplotype's minimizer index (canonical_bits → tile_pos)
//   2. for each read, compute its minimizers; collect shared hits as
//      (read_pos, hap_pos) anchor pairs
//   3. estimate the read→hap offset = median(hap_pos - read_pos) across
//      anchors; if too few anchors, drop the read
//   4. walk the read at that offset; for every base, increment the
//      per-position observation counter (master-base or alt-base)
//   5. final pass over positions:
//        if alt-base count >= min_coverage_branch:
//          emit BranchCandidate with alt_seq = single-base alt + Q-evidence
//        if master-base agreement < min_base_agreement AND majority is alt:
//          emit CurationEvent (master gets corrected to majority)
//
// Memory: 4-counter array (A,C,G,T) per haplotype position = 16 bytes/bp.
// For a 100 Mbp chr14 haplotype that's 1.6 GB per hap — within budget.

#include "wg/master_curator.hpp"
#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

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

inline char decode_base(std::uint8_t b) noexcept {
    static const char kBases[4] = {'A', 'C', 'G', 'T'};
    return b < 4 ? kBases[b] : 'N';
}

struct Anchor { std::uint32_t hap_pos; std::uint32_t read_pos; };

void minimize_seq(std::string_view seq, std::size_t k, std::size_t w,
                  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>& out) {
    if (seq.size() < k) return;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::vector<std::pair<std::uint64_t, std::uint32_t>> wbuf;
    wbuf.reserve(w);
    std::uint64_t kmer = 0;
    std::size_t filled = 0;
    std::uint64_t prev = ~0ULL;
    for (std::size_t i = 0; i < seq.size(); ++i) {
        std::uint8_t b = encode_base(seq[i]);
        if (b > 3) { filled = 0; kmer = 0; wbuf.clear(); prev = ~0ULL; continue; }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;
        std::uint64_t cb = canonical_kmer_bits(kmer, k);
        std::uint32_t pos = static_cast<std::uint32_t>(i - k + 1);
        wbuf.emplace_back(cb, pos);
        if (wbuf.size() > w) wbuf.erase(wbuf.begin());
        if (wbuf.size() == w) {
            auto mn = wbuf[0];
            for (std::size_t j = 1; j < wbuf.size(); ++j)
                if (wbuf[j].first < mn.first) mn = wbuf[j];
            if (mn.first != prev) {
                out[mn.first].push_back(mn.second);
                prev = mn.first;
            }
        }
    }
}

// Sigmoid Q-weight: ramp from 0 at Q≈8 to 1 at Q≈16, centered on Q=12.
inline double q_weight(int q) noexcept {
    return 1.0 / (1.0 + std::exp(-(q - 12) * 0.7));
}

}  // namespace

MasterCurator::MasterCurator(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void MasterCurator::curate(
    std::uint32_t hap_idx,
    const std::vector<MasterTile>& tiles,
    const std::string& hap_seq,
    const std::vector<std::pair<std::string, std::string>>& reads,
    std::vector<CurationEvent>& out_events,
    std::vector<BranchCandidate>& out_branches) const {
    out_events.clear();
    out_branches.clear();
    if (hap_seq.empty() || reads.empty() || tiles.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::size_t w = std::max<std::size_t>(profile_.minimizer_w, 5);

    // 1. Index the haplotype.
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> hap_minz_idx;
    minimize_seq(hap_seq, k, w, hap_minz_idx);

    // 2. Per-position 4-counter.
    std::vector<std::array<std::uint16_t, 4>> counts(hap_seq.size(), {0, 0, 0, 0});
    // 3. Q-weighted alt-base evidence per (pos, alt_base).
    //    Sparse: only non-zero entries kept.
    std::unordered_map<std::uint64_t, double> alt_q_sum;  // (pos<<2 | alt_base) → q-evidence
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> alt_supporters;

    // ONT R10.4.1 SUP mean Q ≈ 18; HiFi ≈ 30. Use profile name as proxy.
    int default_q = (profile_.name && profile_.name[0] == 'h') ? 30 : 18;
    double qw = q_weight(default_q);

    // 4. Walk every read, anchor it on the haplotype, walk bases.
    for (std::size_t r = 0; r < reads.size(); ++r) {
        const auto& seq = reads[r].second;
        if (seq.size() < k) continue;

        // Per-read minimizers as inline scan (avoid building a temp vector).
        const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
        std::vector<std::pair<std::uint64_t, std::uint32_t>> wbuf;
        wbuf.reserve(w);
        std::uint64_t kmer = 0;
        std::size_t filled = 0;
        std::uint64_t prev = ~0ULL;

        // Collect anchors via shared minimizers.
        std::vector<std::int64_t> offsets;  // hap_pos - read_pos
        offsets.reserve(64);
        for (std::size_t i = 0; i < seq.size(); ++i) {
            std::uint8_t b = encode_base(seq[i]);
            if (b > 3) { filled = 0; kmer = 0; wbuf.clear(); prev = ~0ULL; continue; }
            kmer = ((kmer << 2) | b) & mask;
            ++filled;
            if (filled < k) continue;
            std::uint64_t cb = canonical_kmer_bits(kmer, k);
            std::uint32_t pos = static_cast<std::uint32_t>(i - k + 1);
            wbuf.emplace_back(cb, pos);
            if (wbuf.size() > w) wbuf.erase(wbuf.begin());
            if (wbuf.size() < w) continue;
            auto mn = wbuf[0];
            for (std::size_t j = 1; j < wbuf.size(); ++j)
                if (wbuf[j].first < mn.first) mn = wbuf[j];
            if (mn.first == prev) continue;
            prev = mn.first;
            auto it = hap_minz_idx.find(mn.first);
            if (it == hap_minz_idx.end()) continue;
            // Take the first hap position; multiple positions = repetitive,
            // skip to avoid mis-anchoring.
            if (it->second.size() == 1) {
                std::int64_t off = static_cast<std::int64_t>(it->second[0]) -
                                   static_cast<std::int64_t>(mn.second);
                offsets.push_back(off);
            }
        }

        if (offsets.size() < 5) continue;
        std::nth_element(offsets.begin(),
                         offsets.begin() + offsets.size() / 2,
                         offsets.end());
        std::int64_t off = offsets[offsets.size() / 2];

        // Pile-up the read against the haplotype at offset `off`.
        for (std::size_t i = 0; i < seq.size(); ++i) {
            std::int64_t hp = static_cast<std::int64_t>(i) + off;
            if (hp < 0 || hp >= static_cast<std::int64_t>(hap_seq.size())) continue;
            std::uint8_t b = encode_base(seq[i]);
            if (b > 3) continue;
            auto& cnt = counts[hp];
            if (cnt[b] < std::numeric_limits<std::uint16_t>::max())
                ++cnt[b];

            // Track Q-weighted alt-evidence and supporters when read-base
            // disagrees with master-base.
            std::uint8_t mb = encode_base(hap_seq[hp]);
            if (mb < 4 && b != mb) {
                std::uint64_t key = (static_cast<std::uint64_t>(hp) << 2) | b;
                alt_q_sum[key] += qw;
                auto& sup = alt_supporters[key];
                if (sup.size() < 64) sup.push_back(static_cast<std::uint32_t>(r));
            }
        }
    }

    // 5. Emit events + branches based on the per-position pile-up.
    const int min_cov = profile_.min_coverage_branch;
    const double min_agreement = profile_.min_base_agreement;
    for (std::size_t p = 0; p < counts.size(); ++p) {
        const auto& c = counts[p];
        std::uint32_t total = c[0] + c[1] + c[2] + c[3];
        if (total < static_cast<std::uint32_t>(min_cov * 2)) continue;
        std::uint8_t mb = encode_base(hap_seq[p]);
        if (mb >= 4) continue;
        std::uint32_t mb_count = c[mb];

        // CurationEvent: master is in the minority and an alt has > 50%.
        std::uint8_t majority = mb;
        std::uint32_t maj_count = mb_count;
        for (int b = 0; b < 4; ++b) {
            if (c[b] > maj_count) { maj_count = c[b]; majority = static_cast<std::uint8_t>(b); }
        }
        double agreement = total > 0 ? (double)mb_count / total : 1.0;
        if (agreement < min_agreement && majority != mb) {
            CurationEvent ev;
            ev.hap_idx = hap_idx;
            ev.pos_in_hap = static_cast<std::uint32_t>(p);
            ev.original_base = decode_base(mb);
            ev.curated_base = decode_base(majority);
            ev.n_supporting_reads = maj_count;
            out_events.push_back(ev);
        }

        // BranchCandidate: any alt with count >= min_cov.
        for (int b = 0; b < 4; ++b) {
            if (b == mb) continue;
            if (c[b] < min_cov) continue;
            std::uint64_t key = (static_cast<std::uint64_t>(p) << 2) | b;
            BranchCandidate bc;
            bc.hap_idx = hap_idx;
            bc.pos_in_hap = static_cast<std::uint32_t>(p);
            bc.alt_seq.assign(1, decode_base(static_cast<std::uint8_t>(b)));
            auto qit = alt_q_sum.find(key);
            bc.q_weighted_evidence = qit != alt_q_sum.end() ? qit->second : 0.0;
            auto sit = alt_supporters.find(key);
            if (sit != alt_supporters.end()) bc.supporting_read_ids = sit->second;
            out_branches.push_back(std::move(bc));
        }
    }
}

}  // namespace branch::wg
