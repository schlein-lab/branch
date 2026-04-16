#include <gtest/gtest.h>

#include "classify/features.hpp"

using branch::classify::BubbleCandidate;
using branch::classify::Feature;
using branch::classify::FeatureVector;
using branch::classify::kNumFeatures;

TEST(FeaturesTest, kNumFeatures_is_12) {
    EXPECT_EQ(kNumFeatures, 12u);
}

TEST(FeaturesTest, FeatureVector_size_matches_kNumFeatures) {
    FeatureVector v{};
    EXPECT_EQ(v.size(), kNumFeatures);
    EXPECT_EQ(sizeof(v), kNumFeatures * sizeof(float));
}

TEST(FeaturesTest, FeatureEnum_is_contiguous_and_zero_based) {
    EXPECT_EQ(static_cast<std::uint8_t>(Feature::FlankJaccardK31), 0u);
    EXPECT_EQ(static_cast<std::uint8_t>(Feature::BubbleLengthBp),
              kNumFeatures - 1);
}

TEST(FeaturesTest, BubbleCandidate_features_indexable_by_enum) {
    BubbleCandidate bc{};
    bc.features[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.0F;
    bc.features[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] = 1.0F;

    EXPECT_FLOAT_EQ(bc.features[1], 2.0F);
    EXPECT_FLOAT_EQ(bc.features[8], 1.0F);
}
