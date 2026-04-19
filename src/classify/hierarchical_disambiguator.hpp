#pragma once

// BRANCH v0.2 - Hierarchical Branch-vs-Duplication-vs-Mixed disambiguator.
//
// Task P1.2. Given a bubble that the cascade did not resolve cleanly, or
// that landed in the "widersprueliche Diskriminatoren" regime of the
// cascade, run a dedicated 4-signal classifier that returns one of
//   { Branch, Duplication, Mixed, NonSeparable }
// together with a calibrated per-bubble confidence in [0, 1].
//
// Four discriminators (all pulled from classify::FeatureVector, so
// nothing new is required from the feature extractor):
//
//   1. Flank similarity        - Feature::FlankJaccardK31
//                                High Jaccard -> the two flanks agree ->
//                                genuine allelic branch point.
//   2. Depth-sum ratio         - Feature::DepthRatioDiploid
//                                Observed depth relative to expected
//                                diploid. ~1 -> Branch, >=2 -> Dup.
//   3. Read-span coverage      - Feature::ReadSpanCoverageRatio
//                                Fraction of reads that span the full
//                                bubble. High -> Branch-confirmation.
//   4. Prior                   - configurable per-bubble prior (e.g. from
//                                an annotation track). Neutral default
//                                is 0.5, i.e. no information.
//
// On Mixed, we invoke mixed_decomposer::decompose_mixed(). Each resulting
// sub-bubble is re-classified with the same four discriminators, but the
// recursion depth is bounded by DisambiguatorConfig::max_recursion_depth
// (hard-capped at 1 by this task: depth-0 is the top-level call, depth-1
// re-classifies sub-bubbles, depth-2 is rejected and returned as the
// sub-bubble's top-level verdict without further decomposition).
//
// -- Confidence curve -------------------------------------------------
//
// Each discriminator's deviation from its ambiguity midpoint is mapped
// through a smooth logistic of the form
//
//     s_i = 1 / (1 + exp(-k * (x_i - mid_i)))
//
// which lives in [0, 1]. The per-bubble confidence is the fused signal
//
//     conf = w_flank * s_flank
//          + w_depth * s_depth
//          + w_span  * s_span
//          + w_prior * s_prior
//
// with weights that sum to 1.0 (DisambiguatorConfig::weight_*). So:
//
//   - A clean Branch (Jaccard ~= 1.0, depth ~= 1.0, span ~= 1.0)
//     -> conf ~= 1.0.
//   - A clean Duplication (Jaccard ~= 0.5, depth ~= 2.0, span ~= 0.5)
//     -> conf ~= 1.0 (for Duplication's verdict).
//   - A truly ambiguous bubble (all signals mid-range)
//     -> conf ~= 0.5.
//
// Confidence is monotone in signal strength and always reported relative
// to the emitted label, so downstream filters can threshold it directly.

#include <cstdint>
#include <optional>

#include "backend/backend_vtable.hpp"
#include "classify/cascade.hpp"
#include "classify/features.hpp"

namespace branch::classify {

// Configuration for the hierarchical disambiguator. All thresholds,
// weights, and recursion bounds live here -- no magic numbers in the
// implementation.
struct DisambiguatorConfig {
    // ----- Branch vs. Duplication midpoints (s_i = 0.5 at midpoint) ---
    float flank_jaccard_midpoint{0.80f};   // Jaccard: 0.8 ambiguous, 1.0 clean
    float depth_ratio_branch{1.0f};        // Diploid expectation for Branch
    float depth_ratio_duplication{2.0f};   // Diploid expectation for Duplication
    float read_span_midpoint{0.50f};       // 0.5 span ratio ambiguous
    float prior_neutral{0.50f};            // Neutral prior contribution

    // Logistic steepness (shared across signals). Higher = sharper.
    float logistic_slope{12.0f};

    // ----- Fusion weights (must sum to ~1.0) --------------------------
    float weight_flank{0.35f};
    float weight_depth{0.35f};
    float weight_span{0.20f};
    float weight_prior{0.10f};

    // ----- Label-decision thresholds ----------------------------------
    // Branch: flank agrees AND depth near 1x AND span high.
    float branch_flank_floor{0.90f};       // clean Branch floor
    float branch_span_floor{0.60f};        // clean Branch span floor
    float branch_depth_tolerance{0.40f};   // |depth - 1.0| < tol => Branch regime

    // Duplication: depth near 2x, flank need not be super high.
    float dup_depth_floor{1.50f};          // minimum elevated depth

    // Mixed: flank and depth disagree (flank says Branch, depth says Dup).
    float mixed_flank_floor{0.70f};        // moderate flank support
    float mixed_depth_floor{1.30f};        // elevated (not yet Dup) depth

    // Minimum confidence to commit to a decisive (non-NonSeparable) label.
    float min_emit_confidence{0.50f};

    // ----- Hierarchical recursion bound -------------------------------
    // 0 = top-level only, 1 = one level of decomposition (task default).
    std::uint8_t max_recursion_depth{1};

    // Minimum bubble length to even attempt classification.
    std::uint32_t min_bubble_length_bp{500};
};

// Result of a single disambiguator call. For Mixed with successful
// decomposition, sub_bubble_results holds the re-classified sub-paths.
struct DisambiguatorResult {
    branch::backend::BubbleClass label{branch::backend::BubbleClass::Unknown};
    float confidence{0.0f};           // [0, 1]
    std::uint8_t recursion_depth{0};  // 0 at top level, >=1 after decompose
    bool decomposed{false};           // true if mixed_decomposer was run
    // Per-signal diagnostics (kept for BED + regression tests).
    float score_flank{0.0f};
    float score_depth{0.0f};
    float score_span{0.0f};
    float score_prior{0.0f};
};

// Classify a single bubble via the four-signal disambiguator.
// `prior` is the per-bubble prior in [0, 1]; pass std::nullopt to use
// cfg.prior_neutral.
[[nodiscard]] DisambiguatorResult disambiguate(
    const FeatureVector& features,
    const DisambiguatorConfig& cfg = {},
    std::optional<float> prior = std::nullopt) noexcept;

// Hierarchical version: on Mixed, invokes mixed_decomposer and
// re-classifies each sub-bubble with the four-signal classifier. The
// aggregate result's confidence is the mean of sub-bubble confidences
// weighted by sub-bubble read-support. Sub-bubble results are returned
// via `sub_results` (may be nullptr if the caller only needs the top
// verdict).
//
// `per_read_snp_vectors` is the input to mixed_decomposer. Pass an empty
// vector to disable decomposition entirely (the call degrades to
// disambiguate()).
[[nodiscard]] DisambiguatorResult disambiguate_hierarchical(
    const FeatureVector& features,
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    const DisambiguatorConfig& cfg = {},
    std::optional<float> prior = std::nullopt,
    std::vector<DisambiguatorResult>* sub_results = nullptr) noexcept;

}  // namespace branch::classify
