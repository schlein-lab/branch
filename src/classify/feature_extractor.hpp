#pragma once

// BRANCH v0.2 — Bubble feature extraction.
//
// Bridges detect::Bubble (structural candidates) → classify::BubbleCandidate
// (feature-tagged candidates consumed by cascade). Reads the raw graph
// topology + read-support counts stored on BubbleCandidate and fills the
// 12-slot FeatureVector that cascade consumes.
//
// Scope (v0.2): the six features cascade.hpp currently uses —
//   FlankJaccardK31, DepthRatioDiploid, ReadSpanCoverageRatio,
//   SegdupAnnotationFlag, RepeatAnnotationFlag, BubbleLengthBp.
// The other six (ReadSpanCoverageIQR, RefAlignScore{L,R},
// PathLengthDifference, GcContentDivergence, AlleleFrequencyEstimate)
// stay at 0.0; they are a v0.3 concern once ref-alignment + per-position
// depth tracks are wired. See docs/CODE_REVIEW_classify.md.

#include "classify/features.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::classify {

// Extract feature vector from a BubbleCandidate + its graph context.
// The candidate must have been populated from a detect::Bubble — i.e.
// entry_node/exit_node/read_support_{branch,alt}/bubble_length_bp set.
//
// Returns a new FeatureVector; does not mutate `candidate`. Caller is
// expected to assign the result into candidate.features before invoking
// cascade::classify_one.
[[nodiscard]] FeatureVector extract_features(
    const BubbleCandidate& candidate,
    const branch::graph::LosslessGraph& graph);

}  // namespace branch::classify
