#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "classify/cascade.hpp"
#include "classify/features.hpp"
#include "classify/hierarchical_disambiguator.hpp"

using branch::backend::BubbleClass;
using branch::classify::DisambiguatorConfig;
using branch::classify::DisambiguatorResult;
using branch::classify::disambiguate;
using branch::classify::disambiguate_hierarchical;
using branch::classify::Feature;
using branch::classify::FeatureVector;

namespace {

// Build a feature vector with sensible baseline values that pass the
// length guard. Tests then override only the signals under examination.
FeatureVector baseline() {
    FeatureVector f{};
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 2000.0f;
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.50f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.0f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.50f;
    return f;
}

}  // namespace

// ---------------------------------------------------------------------
// 1. Clean Branch -> confidence close to 1.0.
// ---------------------------------------------------------------------
TEST(Disambiguator, CleanBranchHighConfidence) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.99f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.00f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.95f;

    auto r = disambiguate(f);
    EXPECT_EQ(r.label, BubbleClass::Branch);
    EXPECT_GE(r.confidence, 0.90f);
    EXPECT_LE(r.confidence, 1.0f);
    EXPECT_FALSE(r.decomposed);
    EXPECT_EQ(r.recursion_depth, 0u);
}

// ---------------------------------------------------------------------
// 2. Clean Duplication -> label = Duplication, high confidence.
// ---------------------------------------------------------------------
TEST(Disambiguator, CleanDuplicationHighConfidence) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.40f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.10f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.30f;

    auto r = disambiguate(f);
    EXPECT_EQ(r.label, BubbleClass::Duplication);
    EXPECT_GE(r.confidence, 0.85f);
}

// ---------------------------------------------------------------------
// 3. Mixed regime (flank high AND depth high) -> label = Mixed.
// ---------------------------------------------------------------------
TEST(Disambiguator, AmbiguousMixedLabel) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.85f;   // Branch-ish
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.60f; // elevated but not Dup
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.50f;

    auto r = disambiguate(f);
    EXPECT_EQ(r.label, BubbleClass::Mixed);
    EXPECT_GE(r.confidence, 0.40f);
    EXPECT_LE(r.confidence, 0.85f);
}

// ---------------------------------------------------------------------
// 4. Truly ambiguous (all signals mid) -> confidence near 0.5
//    (NonSeparable or low-confidence label).
// ---------------------------------------------------------------------
TEST(Disambiguator, NonSeparableLowConfidence) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.50f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.40f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.40f;

    auto r = disambiguate(f);
    EXPECT_LE(r.confidence, 0.75f);
}

// ---------------------------------------------------------------------
// 5. Length guard: short bubble -> NonSeparable with zero confidence.
// ---------------------------------------------------------------------
TEST(Disambiguator, LengthGuardRejectsShortBubble) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 100.0f;
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.99f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.0f;

    auto r = disambiguate(f);
    EXPECT_EQ(r.label, BubbleClass::NonSeparable);
    EXPECT_FLOAT_EQ(r.confidence, 0.0f);
}

// ---------------------------------------------------------------------
// 6. Prior can push a borderline call toward Branch or Duplication.
// ---------------------------------------------------------------------
TEST(Disambiguator, PriorBiasesBorderlineCall) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.60f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.30f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.45f;

    auto r_branch_prior = disambiguate(f, {}, /*prior=*/0.95f);
    auto r_dup_prior    = disambiguate(f, {}, /*prior=*/0.05f);
    // Branch prior pushes confidence for Branch higher than Dup prior does.
    EXPECT_GT(r_branch_prior.score_prior, r_dup_prior.score_prior);
}

// ---------------------------------------------------------------------
// 7. Mixed -> decompose -> reclassify.
// ---------------------------------------------------------------------
TEST(Disambiguator, MixedDecomposesAndReclassifies) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.85f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.70f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.55f;

    // Two-haplotype read clusters -> decompose_mixed returns 2 subs of 3.
    std::vector<std::vector<char>> snp_vectors = {
        {'A'}, {'A'}, {'A'},
        {'T'}, {'T'}, {'T'},
    };

    std::vector<DisambiguatorResult> subs;
    auto r = disambiguate_hierarchical(f, snp_vectors, {}, std::nullopt, &subs);

    EXPECT_EQ(r.label, BubbleClass::Mixed);
    EXPECT_TRUE(r.decomposed);
    ASSERT_EQ(subs.size(), 2u);
    // Each sub must have recursion_depth 1, per the depth-1 bound.
    for (const auto& s : subs) {
        EXPECT_EQ(s.recursion_depth, 1u);
        EXPECT_FALSE(s.decomposed);
    }
}

// ---------------------------------------------------------------------
// 8. Depth-1 bound enforcement: setting max_recursion_depth=0 disables
//    decomposition even on clearly-Mixed input.
// ---------------------------------------------------------------------
TEST(Disambiguator, DepthBoundDisablesDecomposition) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.85f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.70f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.55f;

    std::vector<std::vector<char>> snp_vectors = {
        {'A'}, {'A'}, {'A'},
        {'T'}, {'T'}, {'T'},
    };

    DisambiguatorConfig cfg{};
    cfg.max_recursion_depth = 0;

    std::vector<DisambiguatorResult> subs;
    auto r = disambiguate_hierarchical(f, snp_vectors, cfg, std::nullopt, &subs);

    EXPECT_EQ(r.label, BubbleClass::Mixed);
    EXPECT_FALSE(r.decomposed);
    EXPECT_TRUE(subs.empty());
}

// ---------------------------------------------------------------------
// 9. Depth bound: sub-bubbles can NEVER themselves decompose further,
//    even if their features would suggest Mixed again.
// ---------------------------------------------------------------------
TEST(Disambiguator, SubBubblesDoNotRecurseFurther) {
    auto f = baseline();
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.85f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.70f;
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = 0.55f;

    std::vector<std::vector<char>> snp_vectors = {
        {'A', 'A'}, {'A', 'C'}, {'A', 'A'},
        {'T', 'G'}, {'T', 'T'}, {'T', 'G'},
    };

    std::vector<DisambiguatorResult> subs;
    DisambiguatorConfig cfg{};
    cfg.max_recursion_depth = 1;  // explicit
    auto r = disambiguate_hierarchical(f, snp_vectors, cfg, std::nullopt, &subs);

    EXPECT_EQ(r.label, BubbleClass::Mixed);
    ASSERT_FALSE(subs.empty());
    for (const auto& s : subs) {
        EXPECT_EQ(s.recursion_depth, 1u);
        EXPECT_FALSE(s.decomposed);  // CRITICAL: no depth-2 recursion.
    }
}

// ---------------------------------------------------------------------
// 10. Cascade Mixed thresholds are parameterised (no magic 0.7/1.5).
// ---------------------------------------------------------------------
TEST(Disambiguator, CascadeMixedThresholdsParameterised) {
    branch::classify::CascadeConfig cfg{};
    EXPECT_FLOAT_EQ(cfg.mixed_flank_lower, 0.70f);
    EXPECT_FLOAT_EQ(cfg.mixed_depth_lower, 1.50f);

    // Tighten flank threshold so a 0.72 flank no longer counts as Mixed.
    cfg.mixed_flank_lower = 0.80f;
    cfg.mixed_depth_lower = 1.50f;

    FeatureVector f{};
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 1000.0f;
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.72f;
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 1.55f;

    auto r = branch::classify::classify_one(f, cfg);
    EXPECT_NE(r.label, BubbleClass::Mixed);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
