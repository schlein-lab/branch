#include "refiner.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <unordered_map>

namespace branch::wg::phaser {

namespace {

inline std::uint64_t hash64(std::uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

inline std::uint8_t base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;
    }
}

void minimizers_of(const std::string& seq, int k, int w,
                   std::vector<std::uint64_t>& out) {
    out.clear();
    if (static_cast<int>(seq.size()) < k) return;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t fwd = 0, rev = 0;
    int valid = 0;
    std::vector<std::pair<std::uint64_t, int>> window;
    window.reserve(static_cast<std::size_t>(w + 1));
    auto push_w = [&](std::uint64_t h, int pos) {
        while (!window.empty() && window.back().first > h) window.pop_back();
        window.emplace_back(h, pos);
        while (!window.empty() && window.front().second <= pos - w)
            window.erase(window.begin());
    };
    out.reserve(seq.size() / static_cast<std::size_t>(w) + 16);
    std::uint64_t last = ~0ULL;
    for (int i = 0; i < static_cast<int>(seq.size()); ++i) {
        const std::uint8_t b = base_code(seq[i]);
        if (b == 4) { fwd = 0; rev = 0; valid = 0; window.clear(); continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << (2 * (k - 1)));
        ++valid;
        if (valid >= k) {
            const std::uint64_t canon = std::min(fwd, rev);
            push_w(hash64(canon), i);
            if (valid >= k + w - 1) {
                if (window.front().first != last) {
                    out.push_back(window.front().first);
                    last = window.front().first;
                }
            }
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

float jaccard_sorted(const std::vector<std::uint64_t>& a,
                     const std::vector<std::uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0f;
    std::size_t i = 0, j = 0, inter = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) { ++inter; ++i; ++j; }
        else if (a[i] < b[j]) ++i; else ++j;
    }
    const std::size_t uni = a.size() + b.size() - inter;
    return uni == 0 ? 0.0f : static_cast<float>(inter) / static_cast<float>(uni);
}

}  // namespace

RefinerProgress Refiner::refine(
    std::vector<UnitigNode>& unitigs,
    std::vector<Bubble>& bubbles,
    std::vector<ReadAssignment>& assignments,
    FastqIndex& reads,
    const RefinerOpts& opts)
{
    (void)bubbles;

    RefinerProgress p;
    p.reads_total = assignments.size();

    // Refiner skips at WG scale.
    //
    // The previous implementation rebuilds unitig sketches and re-scores
    // every assignment against every unitig: O(N_reads × N_unitigs × |sk|)
    // per pass. That is fine for the unit tests (a handful of unitigs)
    // but completely infeasible at WG scale (4.3M × 4M ≈ 1.7e13 ops per
    // pass) and would also need ~50 GB of cached per-read minimizers on
    // top of the ~60 GB raw reads. A scalable v2 refiner needs to score
    // only against neighbouring nodes in the overlap graph, not all
    // unitigs — that is a redesign, not a parameter tweak.
    //
    // For v10g we skip refinement at WG scale and rely on the build-pass
    // assignments produced by phase_engine. This keeps the unit tests
    // green for small inputs while letting WG runs complete.
    constexpr std::size_t kRefineMaxNodes = 50'000;
    if (unitigs.size() > kRefineMaxNodes) {
        std::cerr << "[phaser/refine] skipping at WG scale (n_unitigs="
                  << unitigs.size() << " > " << kRefineMaxNodes
                  << "); v2 will rescore only adjacent nodes\n";
        p.iter_n = 0;
        return p;
    }

    // Pre-compute unitig sketches once per refinement pass.
    const int k = 21, w = 11;
    auto rebuild_unitig_sketches = [&](std::vector<std::vector<std::uint64_t>>& out_sk) {
        out_sk.assign(unitigs.size(), {});
        for (std::size_t i = 0; i < unitigs.size(); ++i) {
            minimizers_of(unitigs[i].seq, k, w, out_sk[i]);
        }
    };

    std::vector<std::vector<std::uint64_t>> u_sk;

    std::string read_seq_buf;
    for (int it = 0; it < opts.max_iters; ++it) {
        rebuild_unitig_sketches(u_sk);
        std::size_t changed_this_iter = 0;
        std::vector<std::uint64_t> rd_sk;

        for (auto& a : assignments) {
            // Skip uncertain force-assigned reads — moving them isn't
            // signal-driven, just noise.
            if (a.tag == ReadTag::UNCERTAIN) continue;
            if (a.tag == ReadTag::UNTAGGED) continue;

            reads.read_seq(a.read_id, read_seq_buf);
            if (read_seq_buf.empty()) continue;
            minimizers_of(read_seq_buf, k, w, rd_sk);

            // Score current home + every other unitig of compatible kind.
            const NodeKind cur_kind = (a.home_node < unitigs.size())
                                       ? unitigs[a.home_node].kind
                                       : NodeKind::UNCLASSIFIED;
            float best_alt_score = 0.0f;
            NodeId best_alt = ~0u;
            float cur_score = (a.home_node < unitigs.size())
                              ? jaccard_sorted(rd_sk, u_sk[a.home_node])
                              : 0.0f;
            for (std::size_t ni = 0; ni < unitigs.size(); ++ni) {
                if (static_cast<NodeId>(ni) == a.home_node) continue;
                // Don't move H1_ONLY reads to BUBBLE_H2 unitigs etc. without
                // first being free of a hard tag — once tagged, only allow
                // moves between SHARED <-> same-kind bubble.
                const auto kn = unitigs[ni].kind;
                if (cur_kind == NodeKind::BUBBLE_H1 && kn == NodeKind::BUBBLE_H2) continue;
                if (cur_kind == NodeKind::BUBBLE_H2 && kn == NodeKind::BUBBLE_H1) continue;
                const float s = jaccard_sorted(rd_sk, u_sk[ni]);
                if (s > best_alt_score) {
                    best_alt_score = s;
                    best_alt = static_cast<NodeId>(ni);
                }
            }
            if (best_alt != ~0u && best_alt_score > cur_score + opts.min_score_margin) {
                // Move the read.
                if (a.home_node < unitigs.size()) {
                    auto& mr = unitigs[a.home_node].member_reads;
                    mr.erase(std::remove(mr.begin(), mr.end(), a.read_id), mr.end());
                }
                a.home_node = best_alt;
                a.confidence = std::min(1.0f, best_alt_score);
                unitigs[best_alt].member_reads.push_back(a.read_id);
                // Re-derive tag from new home.
                switch (unitigs[best_alt].kind) {
                    case NodeKind::SHARED:    a.tag = ReadTag::SHARED;  break;
                    case NodeKind::BUBBLE_H1: a.tag = ReadTag::H1_ONLY; break;
                    case NodeKind::BUBBLE_H2: a.tag = ReadTag::H2_ONLY; break;
                    case NodeKind::BRANCH:    a.tag = ReadTag::BRANCH;  break;
                    default: break;
                }
                ++changed_this_iter;
            }
        }

        p.iter_n = it + 1;
        p.reads_remapped += changed_this_iter;
        p.per_iter_changes.push_back(changed_this_iter);
        std::cerr << "[phaser/refine] iter " << (it + 1)
                  << " changed=" << changed_this_iter << "\n";
        if (changed_this_iter <
            static_cast<std::size_t>(opts.change_floor * static_cast<double>(p.reads_total))) {
            p.converged = true;
            break;
        }
    }
    return p;
}

}  // namespace branch::wg::phaser
