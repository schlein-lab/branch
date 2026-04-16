#pragma once

// BRANCH v0.1 - Bubble Classification Cascade.
//
// Per algorithm + ML expert consultation: sequential cascade with
// early-exit per stage outperforms parallel scoring. Each stage is a
// model (eventually LightGBM per v0.3), here represented as a plain
// function returning a probability in [0, 1]. If the probability
// exceeds early_exit_threshold, the cascade stops and returns the
// current class.
//
// Stage order (fixed, v0.3 will tune thresholds via trained models):
//   1. FlankJaccard     -> branch-prior
//   2. DepthRatio       -> duplication-prior
//   3. ReadSpanRatio    -> branch-confirmation
//   4. AnnotationPrior  -> tie-breaker from UCSC tracks
//
// Recursion termination (5 guarantees per algorithm expert):
//   (1) min_bubble_length_bp = 500
//   (2) min_coverage = 3.0
//   (3) max_recursion_depth = 8
//   (4) cycle-detection via parent flank-hash comparison
//   (5) BIC/AIC info-theoretic stop

#include <array>
#include <cstdint>

#include "backend/backend_vtable.hpp"
#include "classify/features.hpp"

namespace branch::classify {

struct CascadeConfig {
    float flank_jaccard_early_exit{0.99f};   // stage 1: >0.99 -> Branch
    float depth_ratio_dup_threshold{1.8f};   // stage 2: >=1.8 -> Duplication
    float read_span_branch_threshold{0.5f};  // stage 3: >=0.5 -> Branch-confirm
    float min_confidence_to_emit{0.5f};      // below -> NonSeparable

    // Termination guarantees for recursive Mixed decomposition
    std::uint32_t min_bubble_length_bp{500};
    float min_coverage{3.0f};
    std::uint8_t max_recursion_depth{8};
};

// Stage function signature. Returns probability in [0, 1] that the
// bubble belongs to the stages class, or 0 if the stage has no signal.
using StageFn = float (*)(const FeatureVector&);

struct StageResult {
    branch::backend::BubbleClass label;
    float confidence;
    std::uint8_t stage_index;   // which stage decided (0..3)
};

// Static reference implementations (regel-based v0.1 placeholders).
// These are replaced by LightGBM inference in v0.3.
inline float stage_flank_jaccard(const FeatureVector& f) {
    return f[static_cast<std::size_t>(Feature::FlankJaccardK31)];
}
inline float stage_depth_ratio(const FeatureVector& f) {
    auto r = f[static_cast<std::size_t>(Feature::DepthRatioDiploid)];
    return r > 0.0f ? r : 0.0f;
}
inline float stage_read_span(const FeatureVector& f) {
    return f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)];
}
inline float stage_annotation_prior(const FeatureVector& f) {
    return f[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)];
}

// Single-candidate classification. For batch use, call inside a loop
// from a Backend implementation (enables GPU kernels to replace it).
inline StageResult classify_one(const FeatureVector& features,
                                const CascadeConfig& cfg = {}) {
    // Stage 1: Flank Jaccard
    float p1 = stage_flank_jaccard(features);
    if (p1 >= cfg.flank_jaccard_early_exit) {
        return {branch::backend::BubbleClass::Branch, p1, 0};
    }

    // Stage 2: Depth-Sum
    float p2 = stage_depth_ratio(features);
    if (p2 >= cfg.depth_ratio_dup_threshold) {
        return {branch::backend::BubbleClass::Duplication, 1.0f, 1};
    }

    // Stage 3: Read-Span Ratio
    float p3 = stage_read_span(features);
    if (p3 >= cfg.read_span_branch_threshold) {
        return {branch::backend::BubbleClass::Branch, p3, 2};
    }

    // Stage 4: Annotation prior tie-break
    float p4 = stage_annotation_prior(features);
    if (p4 > 0.5f) {
        return {branch::backend::BubbleClass::Duplication, p4, 3};
    }

    // No stage decisively fired -> non-separable
    return {branch::backend::BubbleClass::NonSeparable,
            cfg.min_confidence_to_emit, 3};
}

}  // namespace branch::classify
