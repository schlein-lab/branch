// BRANCH v0.2 — Feature extraction implementation.
//
// See feature_extractor.hpp for scope. This file populates six of the
// twelve FeatureVector slots from the fields on BubbleCandidate and the
// LosslessGraph's per-node / per-edge read support. The remaining six
// slots (indices 3-7, 10) stay at 0.0f and are a v0.3 deliverable.
//
// Sources of truth per feature:
//   FlankJaccardK31      — stub 0.0f until sequence storage lands on
//                          graph::Node. v0.3 reads k=31 sketches of the
//                          left+right flanks and emits real Jaccard.
//                          Kept at 0.0f so cascade falls through to
//                          DepthRatio / ReadSpan which already work.
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

namespace branch::classify {

namespace {

// Mean read_support per edge across the whole graph. Clamped to 1.0 so
// downstream division stays finite on an empty / zero-support graph.
// Using the mean (not median) is sufficient for v0.2 — Stage 2's 2.0x
// threshold is loose enough that outlier-driven mean bias is absorbed.
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

    // Indices 3, 4, 5, 6, 7, 10 intentionally stay at 0.0f. See
    // feature_extractor.hpp / docs/CODE_REVIEW_classify.md (v0.3 item).

    // --- FlankJaccardK31 (index 0) ----------------------------------
    // Stub 0.0f: sequence storage is not on graph::Node yet. Cascade
    // Stage 1 needs >= 0.99 to fire so 0.0 cleanly falls through.
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.0f;

    // --- DepthRatioDiploid (index 1) --------------------------------
    const std::uint32_t total_support =
        candidate.read_support_branch + candidate.read_support_alt;
    const float baseline = diploid_baseline(graph);
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] =
        static_cast<float>(total_support) / baseline;

    // --- ReadSpanCoverageRatio (index 2) ----------------------------
    // Fraction of reads at the entry node that traverse the bubble.
    // Guard out-of-range entry_node: for a default-constructed
    // candidate with no graph context we emit 0.0 and let cascade
    // decide on other signals.
    float span_ratio = 0.0f;
    if (candidate.entry_node < graph.node_count()) {
        const std::uint32_t entry_reads =
            graph.node(candidate.entry_node).read_support;
        if (entry_reads > 0) {
            span_ratio = static_cast<float>(total_support) /
                         static_cast<float>(entry_reads);
            // Clamp to [0, 1]; ratios > 1 happen on multi-allelic
            // bubbles where alt supports sum above the entry count and
            // would falsely trip Stage 3 on spurious cases.
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
