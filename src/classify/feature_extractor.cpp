// BRANCH v0.2 — Feature extraction implementation.
//
// See feature_extractor.hpp for scope. This file populates six of the
// twelve FeatureVector slots from the fields on BubbleCandidate and the
// LosslessGraph's per-node / per-edge read support. The remaining six
// slots (indices 3-7, 10) stay at 0.0f and are a v0.3 deliverable.
//
// Sources of truth per feature:
//   FlankJaccardK31      — k=31 Jaccard similarity of left flank (from
//                          entry_node consensus) vs right flank (from
//                          exit_node consensus). Uses 155bp flanks (5×K).
//                          Returns 0.0f if consensus fields are empty.
//   DepthRatioDiploid    — (branch + alt) / diploid_baseline, where
//                          diploid_baseline is the mean per-edge read
//                          support across the whole graph. This matches
//                          the "elevated depth ⇒ Duplication" signal
//                          cascade Stage 2 consumes.
//   ReadSpanCoverageRatio — fraction of entry-node reads that actually
//                          traverse the bubble: (branch + alt) /
//                          max(entry.read_support, 1). Cascade Stage 3
//                          uses this as the branch-confirmation signal.
//   SegdupAnnotationFlag  — 0/1 from candidate.segdup_flag.
//   RepeatAnnotationFlag  — 0/1 from candidate.repeat_flag.
//   BubbleLengthBp        — candidate.bubble_length_bp (longest alt).

#include "classify/feature_extractor.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>

namespace branch::classify {

namespace {

// K-mer parameters for FlankJaccardK31.
constexpr std::size_t kFlankK = 31;
constexpr std::size_t kFlankLength = 5 * kFlankK;  // 155bp

// Extract k-mers from a sequence into a set. Returns empty set if
// sequence is shorter than K.
std::unordered_set<std::string_view> extract_kmers(
    std::string_view seq, std::size_t k) {
    std::unordered_set<std::string_view> kmers;
    if (seq.size() < k) {
        return kmers;
    }
    kmers.reserve(seq.size() - k + 1);
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        kmers.insert(seq.substr(i, k));
    }
    return kmers;
}

// Compute Jaccard similarity: |A ∩ B| / |A ∪ B|.
// Returns 0.0f if both sets are empty.
float jaccard_similarity(
    const std::unordered_set<std::string_view>& a,
    const std::unordered_set<std::string_view>& b) {
    if (a.empty() && b.empty()) {
        return 0.0f;
    }
    std::size_t intersection = 0;
    for (const auto& kmer : a) {
        if (b.count(kmer) > 0) {
            ++intersection;
        }
    }
    const std::size_t union_size = a.size() + b.size() - intersection;
    if (union_size == 0) {
        return 0.0f;
    }
    return static_cast<float>(intersection) / static_cast<float>(union_size);
}

// Compute FlankJaccardK31: Jaccard of k=31 kmers from left flank (right
// side of entry_node consensus) vs right flank (left side of exit_node
// consensus). Uses 155bp flanks. Returns 0.0f if consensus is empty.
float compute_flank_jaccard_k31(
    const branch::graph::LosslessGraph& graph,
    std::uint32_t entry_node,
    std::uint32_t exit_node) {
    // Bounds check
    if (entry_node >= graph.node_count() || exit_node >= graph.node_count()) {
        return 0.0f;  // TODO: invalid node IDs
    }

    const auto& entry_consensus = graph.node(entry_node).consensus;
    const auto& exit_consensus = graph.node(exit_node).consensus;

    // Fallback if consensus not populated
    if (entry_consensus.empty() || exit_consensus.empty()) {
        return 0.0f;  // TODO: consensus not available
    }

    // Left flank: rightmost 155bp of entry_node consensus
    std::string_view left_flank;
    if (entry_consensus.size() >= kFlankLength) {
        left_flank = std::string_view(entry_consensus).substr(
            entry_consensus.size() - kFlankLength, kFlankLength);
    } else {
        left_flank = entry_consensus;
    }

    // Right flank: leftmost 155bp of exit_node consensus
    std::string_view right_flank;
    if (exit_consensus.size() >= kFlankLength) {
        right_flank = std::string_view(exit_consensus).substr(0, kFlankLength);
    } else {
        right_flank = exit_consensus;
    }

    // Extract k-mers and compute Jaccard
    auto left_kmers = extract_kmers(left_flank, kFlankK);
    auto right_kmers = extract_kmers(right_flank, kFlankK);

    return jaccard_similarity(left_kmers, right_kmers);
}

// Mean read_support per edge across the whole graph. Clamped to 1.0 so
// downstream division stays finite on an empty / zero-support graph.
float diploid_baseline(const branch::graph::LosslessGraph& graph) {
    const auto& edges = graph.edges();
    if (edges.empty()) {
        return 1.0f;
    }
    std::uint64_t sum = 0;
    for (const auto& e : edges) {
        sum += e.read_support;
    }
    const float mean = static_cast<float>(sum) / static_cast<float>(edges.size());
    return std::max(mean, 1.0f);
}

}  // namespace

FeatureVector extract_features(
    const BubbleCandidate& candidate,
    const branch::graph::LosslessGraph& graph) {

    FeatureVector f{};  // zero-initialise all 12 slots

    // --- FlankJaccardK31 (index 0) ----------------------------------
    // k=31 Jaccard of entry_node right-flank vs exit_node left-flank.
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] =
        compute_flank_jaccard_k31(graph, candidate.entry_node, candidate.exit_node);

    // --- DepthRatioDiploid (index 1) --------------------------------
    const std::uint32_t total_support =
        candidate.read_support_branch + candidate.read_support_alt;
    const float baseline = diploid_baseline(graph);
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] =
        static_cast<float>(total_support) / baseline;

    // --- ReadSpanCoverageRatio (index 2) ----------------------------
    float span_ratio = 0.0f;
    if (candidate.entry_node < graph.node_count()) {
        const std::uint32_t entry_reads =
            graph.node(candidate.entry_node).read_support;
        if (entry_reads > 0) {
            span_ratio = static_cast<float>(total_support) /
                         static_cast<float>(entry_reads);
            if (span_ratio > 1.0f) span_ratio = 1.0f;
        }
    }
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = span_ratio;

    // --- Annotation flags (indices 8, 9) ----------------------------
    f[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] =
        candidate.segdup_flag ? 1.0f : 0.0f;
    f[static_cast<std::size_t>(Feature::RepeatAnnotationFlag)] =
        candidate.repeat_flag ? 1.0f : 0.0f;

    // --- BubbleLengthBp (index 11) ----------------------------------
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] =
        static_cast<float>(candidate.bubble_length_bp);

    return f;
}

}  // namespace branch::classify
