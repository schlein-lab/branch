#include "bubble_classifier.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <unordered_map>

namespace branch::wg::phaser {

namespace {

// Walk forward from each alt_paths[a] until reaching sink_anchor, marking
// each interior node as belonging to alt index `a`. Returns a flat
// node_id -> alt_idx map for this single bubble.
std::unordered_map<NodeId, std::uint8_t> bubble_node_to_alt(
    const Bubble& b, const std::vector<UnitigNode>& unitigs)
{
    std::unordered_map<NodeId, std::uint8_t> out;
    out.reserve(b.alt_paths.size() * 4);
    for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
        NodeId cur = b.alt_paths[a];
        int guard = 0;
        while (cur != b.sink_anchor && cur < unitigs.size() && guard++ < 64) {
            out[cur] = static_cast<std::uint8_t>(a);
            if (unitigs[cur].next_nodes.empty()) break;
            cur = unitigs[cur].next_nodes[0];
        }
    }
    return out;
}

// Patterns to match, in priority order. Each entry is the canonical VAF
// vector (sorted descending) plus the type and topological branch_depth.
//
// Depth-K interpretation: at depth K, the youngest sister-sub-clones have
// VAF = 2^-(K+1). The canonical pattern is:
//   depth-0 (HET):    [1/2, 1/2]                              founder pair
//   depth-1 (BRANCH): [1/2, 1/4, 1/4]                         founder + 2 sister-clones
//   depth-2 (BIB):    [1/2, 1/4, 1/8, 1/8]                    one of the depth-1 clones sub-branched
//   depth-3:          [1/2, 1/4, 1/8, 1/16, 1/16]             ...
// In percent that's e.g. depth-2 = [50, 25, 12.5, 12.5] — the
// "12.5/12.5/25/50" pattern Christian described.
struct Pattern {
    std::vector<float> vafs;
    BubbleType         type;
    std::uint8_t       branch_depth;
};

const std::vector<Pattern>& canonical_patterns() {
    static const std::vector<Pattern> p = {
        // Most specific (deeper) first so longer patterns are tried first.
        {{0.5f, 0.25f, 0.125f, 0.0625f, 0.0625f}, BubbleType::BRANCH_IN_BRANCH, 3},
        {{0.5f, 0.25f, 0.125f, 0.125f},           BubbleType::BRANCH_IN_BRANCH, 2},
        {{0.5f, 0.25f, 0.25f},                    BubbleType::BRANCH,           1},
        {{1.0f / 3, 1.0f / 3, 1.0f / 3},          BubbleType::TRIALLELIC,       0},
        {{0.5f, 0.5f},                            BubbleType::HET,              0},
        {{1.0f},                                  BubbleType::HOMOZ,            0},
    };
    return p;
}

// Result of classifying a VAF vector.
struct ClassResult {
    BubbleType   type;
    std::uint8_t branch_depth;
};

bool vec_within_tol(const std::vector<float>& a,
                    const std::vector<float>& b, float tol) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > tol) return false;
    }
    return true;
}

}  // namespace

BubbleType classify_vaf_vector(const std::vector<float>& vafs, float tol) {
    if (vafs.empty()) return BubbleType::UNKNOWN;
    for (const auto& p : canonical_patterns()) {
        if (vec_within_tol(vafs, p.vafs, tol)) return p.type;
    }
    return BubbleType::NOISY;
}

// Same as classify_vaf_vector but also returns branch_depth.
ClassResult classify_vaf_vector_full(const std::vector<float>& vafs, float tol) {
    if (vafs.empty()) return {BubbleType::UNKNOWN, 0};
    for (const auto& p : canonical_patterns()) {
        if (vec_within_tol(vafs, p.vafs, tol)) return {p.type, p.branch_depth};
    }
    return {BubbleType::NOISY, 0};
}

void BubbleClassifier::classify(
    const std::vector<UnitigNode>& unitigs,
    const std::vector<ReadAssignment>& assignments,
    std::vector<Bubble>& bubbles,
    const BubbleClassifierOpts& opts)
{
    const auto t0 = std::chrono::steady_clock::now();

    // For each bubble, count reads whose home_node sits on each alt-path.
    for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
        if (opts.log_every > 0 && bi % opts.log_every == 0 && bi > 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            std::fprintf(stderr,
                "[phaser/bc] classify %zu/%zu (%.1f%%) elapsed=%.0fs\n",
                bi, bubbles.size(),
                100.0 * static_cast<double>(bi) /
                static_cast<double>(bubbles.size()),
                el);
            std::fflush(stderr);
        }
        auto& b = bubbles[bi];
        const std::size_t n_alts = b.alt_paths.size();
        if (n_alts == 0) continue;

        auto n2a = bubble_node_to_alt(b, unitigs);

        std::vector<std::uint32_t> counts(n_alts, 0);
        // Iterate assignments — at WG scale ~5M; small constant per bubble.
        // For every read whose home_node is on this bubble's alt path,
        // bump the corresponding counter.
        for (const auto& a : assignments) {
            if (a.home_node >= unitigs.size()) continue;
            auto it = n2a.find(a.home_node);
            if (it == n2a.end()) continue;
            ++counts[it->second];
        }
        std::uint32_t total = 0;
        for (auto c : counts) total += c;

        b.read_counts = counts;

        if (total < opts.min_total_reads) {
            b.type = BubbleType::UNKNOWN;
            b.vafs.assign(n_alts, 0.0f);
            continue;
        }

        std::vector<float> vafs(n_alts, 0.0f);
        for (std::size_t i = 0; i < n_alts; ++i) {
            vafs[i] = static_cast<float>(counts[i]) /
                      static_cast<float>(total);
        }
        b.vafs = vafs;

        // For pattern matching we sort a copy descending; b.vafs keeps
        // the original alt-order so downstream code can map back to
        // alt_paths[i].
        std::vector<float> sorted_vafs = vafs;
        std::sort(sorted_vafs.begin(), sorted_vafs.end(),
                  std::greater<float>());
        auto cls = classify_vaf_vector_full(sorted_vafs, opts.vaf_tol);
        b.type = cls.type;
        b.branch_depth = cls.branch_depth;

        // For BRANCH / BRANCH_IN_BRANCH, identify the parent alt
        // (the alt with the largest VAF, ~0.5 in the canonical pattern).
        if (b.type == BubbleType::BRANCH ||
            b.type == BubbleType::BRANCH_IN_BRANCH) {
            std::size_t parent = 0;
            for (std::size_t i = 1; i < n_alts; ++i) {
                if (vafs[i] > vafs[parent]) parent = i;
            }
            b.parent_alt = static_cast<std::uint32_t>(parent);
        }
    }

    // Tally for stderr summary.
    std::size_t het = 0, branch = 0, bib = 0, tri = 0, hom = 0,
                noisy = 0, unk = 0;
    for (const auto& b : bubbles) {
        switch (b.type) {
            case BubbleType::HET:              ++het;   break;
            case BubbleType::HOMOZ:            ++hom;   break;
            case BubbleType::BRANCH:           ++branch; break;
            case BubbleType::BRANCH_IN_BRANCH: ++bib;   break;
            case BubbleType::TRIALLELIC:       ++tri;   break;
            case BubbleType::NOISY:            ++noisy; break;
            case BubbleType::EXTINCT_INTERMEDIATE:      ++noisy; break;
            case BubbleType::UNKNOWN:          ++unk;   break;
        }
    }
    std::cerr << "[phaser/bc] classified " << bubbles.size()
              << " bubbles: het=" << het
              << " branch=" << branch
              << " branch_in_branch=" << bib
              << " triallelic=" << tri
              << " homoz=" << hom
              << " noisy=" << noisy
              << " unknown=" << unk << "\n";
    std::fflush(stderr);
}

void BubbleClassifier::classify_from_vafs(
    std::vector<Bubble>& bubbles,
    const BubbleClassifierOpts& opts)
{
    std::size_t het = 0, branch = 0, bib = 0, tri = 0, hom = 0,
                noisy = 0, unk = 0;
    for (auto& b : bubbles) {
        if (b.vafs.empty()) {
            b.type = BubbleType::UNKNOWN;
            ++unk; continue;
        }
        std::uint32_t total = 0;
        for (auto c : b.read_counts) total += c;
        if (total < opts.min_total_reads) {
            b.type = BubbleType::UNKNOWN;
            ++unk; continue;
        }
        std::vector<float> sorted_vafs = b.vafs;
        std::sort(sorted_vafs.begin(), sorted_vafs.end(),
                  std::greater<float>());
        auto cls = classify_vaf_vector_full(sorted_vafs, opts.vaf_tol);
        b.type = cls.type;
        b.branch_depth = cls.branch_depth;
        if (b.type == BubbleType::BRANCH ||
            b.type == BubbleType::BRANCH_IN_BRANCH) {
            std::size_t parent = 0;
            for (std::size_t i = 1; i < b.vafs.size(); ++i) {
                if (b.vafs[i] > b.vafs[parent]) parent = i;
            }
            b.parent_alt = static_cast<std::uint32_t>(parent);
        }
        switch (b.type) {
            case BubbleType::HET:              ++het;    break;
            case BubbleType::HOMOZ:            ++hom;    break;
            case BubbleType::BRANCH:           ++branch; break;
            case BubbleType::BRANCH_IN_BRANCH: ++bib;    break;
            case BubbleType::TRIALLELIC:       ++tri;    break;
            case BubbleType::NOISY:            ++noisy;  break;
            case BubbleType::EXTINCT_INTERMEDIATE:      ++noisy;  break;
            case BubbleType::UNKNOWN:          ++unk;    break;
        }
    }
    std::cerr << "[phaser/bc] classified-from-vafs " << bubbles.size()
              << " bubbles (clonally branched haplotype lineage): "
              << "het(D0)=" << het
              << " branch(D1)=" << branch
              << " bib(D2)=" << bib
              << " triallelic=" << tri
              << " homoz=" << hom
              << " noisy=" << noisy
              << " unknown=" << unk << "\n";
    std::fflush(stderr);
}

}  // namespace branch::wg::phaser
