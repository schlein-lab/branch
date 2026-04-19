#include "classify/hierarchical_disambiguator.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "classify/mixed_decomposer.hpp"

namespace branch::classify {

namespace {

// Logistic squashing function. Returns a value in (0, 1).
inline float logistic(float x, float midpoint, float slope) noexcept {
    return 1.0f / (1.0f + std::exp(-slope * (x - midpoint)));
}

// Extract the discriminator scores for a feature vector. Each score is
// a [0, 1] "confidence this bubble is NOT purely random noise on this
// signal" measure, relative to the respective ambiguity midpoint.
struct SignalScores {
    float flank;    // 0 -> dissimilar flanks (Dup-ish), 1 -> identical flanks (Branch-ish)
    float depth;    // 0 -> near 1x diploid (Branch), 1 -> near 2x (Duplication)
    float span;     // 0 -> few spanning reads, 1 -> fully spanned (Branch-confirm)
    float prior;    // 0 -> prior says Duplication, 1 -> prior says Branch

    float flank_raw;
    float depth_raw;
    float span_raw;
};

SignalScores compute_scores(const FeatureVector& f,
                            const DisambiguatorConfig& cfg,
                            float prior_value) noexcept {
    const float flank_raw = f[static_cast<std::size_t>(Feature::FlankJaccardK31)];
    const float depth_raw = f[static_cast<std::size_t>(Feature::DepthRatioDiploid)];
    const float span_raw  = f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)];

    SignalScores s{};
    s.flank_raw = flank_raw;
    s.depth_raw = depth_raw;
    s.span_raw  = span_raw;

    // Flank: high Jaccard -> Branch. Logistic mid at cfg.flank_jaccard_midpoint.
    s.flank = logistic(flank_raw, cfg.flank_jaccard_midpoint, cfg.logistic_slope);

    // Depth: we want "how close to 2x" as Duplication score. Midpoint is
    // the arithmetic mean of Branch and Dup expectations.
    const float depth_mid = 0.5f * (cfg.depth_ratio_branch + cfg.depth_ratio_duplication);
    s.depth = logistic(depth_raw, depth_mid, cfg.logistic_slope);

    // Span: higher ratio -> Branch-confirm.
    s.span = logistic(span_raw, cfg.read_span_midpoint, cfg.logistic_slope);

    // Prior: treat as already-in-[0,1]. "1" = prior says Branch.
    s.prior = std::clamp(prior_value, 0.0f, 1.0f);

    return s;
}

// Fuse four signals into a single scalar confidence. For label
// `Branch` we want signals (flank high, depth low, span high, prior
// high). For label `Duplication` we invert the depth term: a high raw
// depth score contributes positively to the Dup verdict.
//
// Both cases are written symmetrically in the spirit of softmax-over-
// labels, but since we only need a per-label monotone confidence the
// simpler "label-aligned" formulation suffices.
float fuse_branch_confidence(const SignalScores& s,
                             const DisambiguatorConfig& cfg) noexcept {
    // For Branch: flank high, depth low (= 1 - s.depth), span high, prior high.
    return cfg.weight_flank * s.flank
         + cfg.weight_depth * (1.0f - s.depth)
         + cfg.weight_span  * s.span
         + cfg.weight_prior * s.prior;
}

float fuse_dup_confidence(const SignalScores& s,
                          const DisambiguatorConfig& cfg) noexcept {
    // For Duplication: flank low (= 1 - s.flank), depth high, span low, prior low.
    return cfg.weight_flank * (1.0f - s.flank)
         + cfg.weight_depth * s.depth
         + cfg.weight_span  * (1.0f - s.span)
         + cfg.weight_prior * (1.0f - s.prior);
}

float fuse_mixed_confidence(const SignalScores& s,
                            const DisambiguatorConfig& cfg) noexcept {
    // Mixed fires when flank AND depth both point strongly to different
    // labels. Its confidence is the geometric-mean-ish of "is my flank
    // Branch-ish?" and "is my depth Dup-ish?", modulated by span (low
    // span partially supports ambiguity) and by a neutral prior term.
    const float flank_term = s.flank;           // high -> flank says Branch
    const float depth_term = s.depth;           // high -> depth says Dup
    const float disagree = std::min(flank_term, depth_term);
    const float span_penalty = 1.0f - std::fabs(s.span - 0.5f) * 2.0f;  // maxed at s.span=0.5
    return cfg.weight_flank * disagree
         + cfg.weight_depth * disagree
         + cfg.weight_span  * std::max(0.0f, span_penalty)
         + cfg.weight_prior * (1.0f - std::fabs(s.prior - 0.5f) * 2.0f);
}

// Decide the label + confidence given the raw feature values and the
// derived signal scores. Confidence is always aligned to the emitted
// label and clamped to [0, 1].
DisambiguatorResult decide(const SignalScores& s,
                           const DisambiguatorConfig& cfg) noexcept {
    DisambiguatorResult r{};
    r.score_flank = s.flank;
    r.score_depth = s.depth;
    r.score_span  = s.span;
    r.score_prior = s.prior;

    // Pre-guard: if depth is essentially zero, no support at all.
    if (s.depth_raw <= 0.0f && s.span_raw <= 0.0f && s.flank_raw <= 0.0f) {
        r.label = branch::backend::BubbleClass::NonSeparable;
        r.confidence = 0.0f;
        return r;
    }

    // Rule order: clear Branch -> clear Duplication -> Mixed -> NonSeparable.
    const bool branch_regime =
        (s.flank_raw >= cfg.branch_flank_floor) &&
        (s.span_raw  >= cfg.branch_span_floor) &&
        (std::fabs(s.depth_raw - cfg.depth_ratio_branch) <= cfg.branch_depth_tolerance);

    const bool dup_regime = (s.depth_raw >= cfg.dup_depth_floor) &&
                            (s.flank_raw < cfg.branch_flank_floor);

    const bool mixed_regime = (s.flank_raw >= cfg.mixed_flank_floor) &&
                              (s.depth_raw >= cfg.mixed_depth_floor);

    const float conf_branch = std::clamp(fuse_branch_confidence(s, cfg), 0.0f, 1.0f);
    const float conf_dup    = std::clamp(fuse_dup_confidence(s, cfg),    0.0f, 1.0f);
    const float conf_mixed  = std::clamp(fuse_mixed_confidence(s, cfg),  0.0f, 1.0f);

    if (branch_regime && conf_branch >= cfg.min_emit_confidence) {
        r.label = branch::backend::BubbleClass::Branch;
        r.confidence = conf_branch;
        return r;
    }
    // Mixed is checked BEFORE Duplication: Mixed's predicate is a strict
    // conjunction (flank AND depth both elevated) and must win over the
    // looser single-signal Duplication regime when both match.
    if (mixed_regime && conf_mixed >= cfg.min_emit_confidence) {
        r.label = branch::backend::BubbleClass::Mixed;
        r.confidence = conf_mixed;
        return r;
    }
    if (dup_regime && conf_dup >= cfg.min_emit_confidence) {
        r.label = branch::backend::BubbleClass::Duplication;
        r.confidence = conf_dup;
        return r;
    }

    // Otherwise, pick the label with the highest confidence but mark it
    // NonSeparable if everything is below the emit threshold.
    float best = conf_branch;
    branch::backend::BubbleClass best_label = branch::backend::BubbleClass::Branch;
    if (conf_dup > best) { best = conf_dup; best_label = branch::backend::BubbleClass::Duplication; }
    if (conf_mixed > best) { best = conf_mixed; best_label = branch::backend::BubbleClass::Mixed; }

    if (best < cfg.min_emit_confidence) {
        r.label = branch::backend::BubbleClass::NonSeparable;
        r.confidence = best;
        return r;
    }

    r.label = best_label;
    r.confidence = best;
    return r;
}

}  // namespace

DisambiguatorResult disambiguate(const FeatureVector& features,
                                 const DisambiguatorConfig& cfg,
                                 std::optional<float> prior) noexcept {
    // Length guard. If the bubble is too short, report NonSeparable with
    // zero confidence (same contract as cascade::classify_one).
    const float bubble_length = features[static_cast<std::size_t>(Feature::BubbleLengthBp)];
    if (bubble_length < static_cast<float>(cfg.min_bubble_length_bp)) {
        DisambiguatorResult r{};
        r.label = branch::backend::BubbleClass::NonSeparable;
        r.confidence = 0.0f;
        return r;
    }

    const float prior_value = prior.value_or(cfg.prior_neutral);
    const SignalScores s = compute_scores(features, cfg, prior_value);
    DisambiguatorResult r = decide(s, cfg);
    r.recursion_depth = 0;
    r.decomposed = false;
    return r;
}

DisambiguatorResult disambiguate_hierarchical(
    const FeatureVector& features,
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    const DisambiguatorConfig& cfg,
    std::optional<float> prior,
    std::vector<DisambiguatorResult>* sub_results) noexcept {

    DisambiguatorResult top = disambiguate(features, cfg, prior);

    if (top.label != branch::backend::BubbleClass::Mixed) {
        return top;
    }
    if (cfg.max_recursion_depth == 0) {
        // Decomposition disabled by config.
        return top;
    }
    if (per_read_snp_vectors.empty()) {
        return top;  // nothing to decompose
    }

    // Decompose into sub-bubbles. The decomposer returns at most two
    // clusters in the current implementation; treat each as a sub-bubble
    // inheriting the parent's feature vector but with its depth ratio
    // rescaled by the cluster's estimated VAF.
    auto sub_bubbles = decompose_mixed(per_read_snp_vectors);
    if (sub_bubbles.empty()) {
        // Decomposition failed -> keep Mixed verdict at depth 0.
        return top;
    }

    top.decomposed = true;

    // Depth-bound enforcement: we ONLY allow one level. Construct sub
    // feature vectors and re-classify with a disambiguator copy whose
    // remaining recursion budget is zero, so that if a sub-bubble is
    // ALSO mixed we stop there (depth-1 bound).
    DisambiguatorConfig sub_cfg = cfg;
    sub_cfg.max_recursion_depth = 0;

    float total_weight = 0.0f;
    float weighted_conf_sum = 0.0f;

    for (const auto& sb : sub_bubbles) {
        FeatureVector sub_features = features;
        // Scale depth by 2 * vaf so a perfectly balanced Mixed decomposes
        // into two "1x diploid" Branch-ish sub-bubbles.
        const float depth_scale = 2.0f * static_cast<float>(sb.estimated_vaf);
        sub_features[static_cast<std::size_t>(Feature::DepthRatioDiploid)] =
            features[static_cast<std::size_t>(Feature::DepthRatioDiploid)] * depth_scale;
        // Span ratio for sub-bubble is the fraction of reads this cluster
        // owns, floored at the parent's span ratio to avoid double-penalty.
        sub_features[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] =
            std::max(features[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)],
                     static_cast<float>(sb.estimated_vaf));

        DisambiguatorResult sub = disambiguate(sub_features, sub_cfg, prior);
        sub.recursion_depth = 1;
        sub.decomposed = false;

        const float w = static_cast<float>(sb.read_ids.size());
        total_weight += w;
        weighted_conf_sum += w * sub.confidence;

        if (sub_results) sub_results->push_back(sub);
    }

    if (total_weight > 0.0f) {
        // Aggregate confidence = read-support-weighted mean of sub-bubble
        // confidences. Keep the Mixed label at the top level (the caller
        // can inspect `sub_results` for the decomposition verdicts).
        top.confidence = std::clamp(weighted_conf_sum / total_weight, 0.0f, 1.0f);
    }
    top.recursion_depth = 0;  // top-level verdict is still at depth 0
    return top;
}

}  // namespace branch::classify
