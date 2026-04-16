#include <gtest/gtest.h>

#include "classify/cascade.hpp"

using branch::backend::BubbleClass;
using branch::classify::CascadeConfig;
using branch::classify::classify_one;
using branch::classify::Feature;
using branch::classify::FeatureVector;

namespace {
FeatureVector make_features() {
    FeatureVector f{};
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
