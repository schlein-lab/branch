#include <gtest/gtest.h>

#include "classify/cascade.hpp"

using branch::backend::BubbleClass;
using branch::classify::CascadeConfig;
using branch::classify::classify_one;
using branch::classify::Feature;
using branch::classify::FeatureVector;

namespace {
// Cascade has early guards: BubbleLengthBp >= 500 and DepthRatio >= 3.0
// (= coverage proxy). All tests must satisfy these to reach the stages.
FeatureVector make_features() {
    FeatureVector f{};
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 1000.0f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 5.0f;
    return f;
}
}  // namespace

TEST(CascadeTest, Stage1_fires_on_high_flank_jaccard) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.995f;
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::Branch);
    EXPECT_EQ(r.stage_index, 0u);
    EXPECT_GE(r.confidence, 0.995f);
}

TEST(CascadeTest, Stage2_fires_on_high_depth_ratio) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.7f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.0f;
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::Duplication);
    EXPECT_EQ(r.stage_index, 1u);
}

TEST(CascadeTest, Stage3_fires_on_read_span_ratio) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.5f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.0f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.7f;
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::Branch);
    EXPECT_EQ(r.stage_index, 2u);
}

TEST(CascadeTest, Stage4_tiebreak_uses_segdup_annotation) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.5f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.0f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.1f;
    f[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] = 1.0f;
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::Duplication);
    EXPECT_EQ(r.stage_index, 3u);
}

TEST(CascadeTest, No_stage_fires_returns_NonSeparable) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.3f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.0f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.1f;
    f[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] = 0.0f;
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::NonSeparable);
}

TEST(CascadeTest, Config_threshold_changes_affect_decision) {
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.90f;
    // Default threshold 0.99 -> no fire
    auto r1 = classify_one(f);
    EXPECT_NE(r1.label, BubbleClass::Branch);
    // Relaxed threshold 0.85 -> fires
    CascadeConfig cfg{.flank_jaccard_early_exit = 0.85f};
    auto r2 = classify_one(f, cfg);
    EXPECT_EQ(r2.label, BubbleClass::Branch);
}

// ============ Mixed Class / Ambiguous Cases ============

TEST(CascadeTest, Mixed_high_flank_but_high_depth_ratio) {
    // Edge case: FlankJaccard suggests Branch (high), but DepthRatio suggests Dup
    // Stage 1 checks FlankJaccard first (>= 0.99), so high enough flank wins
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.995f;  // High -> Branch
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.5f;  // Also high (Dup signal)
    auto r = classify_one(f);
    // Stage 1 fires first due to cascade order
    EXPECT_EQ(r.label, BubbleClass::Branch);
    EXPECT_EQ(r.stage_index, 0u);
}

TEST(CascadeTest, Mixed_moderate_flank_high_depth_goes_to_dup) {
    // FlankJaccard not high enough for Stage 1, but DepthRatio >= 2.0 fires Stage 2
    auto f = make_features();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.8f;    // Not enough for Stage 1
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.1f;  // >= 2.0 -> Stage 2 Dup
    auto r = classify_one(f);
    EXPECT_EQ(r.label, BubbleClass::Duplication);
    EXPECT_EQ(r.stage_index, 1u);
}

// ============ Guard Tests (min_bubble_length, min_coverage) ============

TEST(CascadeTest, Guard_min_bubble_length_rejects_short_bubble) {
    FeatureVector f{};
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 200.0f;   // Below 500bp guard
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 5.0f;  // OK
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.999f;  // Would fire Stage 1
    auto r = classify_one(f);
    // Guard rejects before any stage fires
    EXPECT_EQ(r.label, BubbleClass::NonSeparable);
}

// NOTE: the min_coverage guard was removed from classify_one() because
// it used the same DepthRatioDiploid feature that Stage 2 needs for
// Duplication detection, making Stage 3 unreachable. Low coverage is
// handled by per-stage confidence damping rather than a hard reject.
// When that damping lands, re-add a test here that verifies reduced
// confidence at low coverage rather than a NonSeparable return.
