#include <gtest/gtest.h>
#include "classify/feature_extractor.hpp"
#include "classify/features.hpp"
#include "graph/lossless_graph.hpp"

using namespace branch::classify;
using namespace branch::graph;

class FlankJaccardTest : public ::testing::Test {
protected:
    LosslessGraph make_graph_with_consensus(
        const std::string& cons_entry,
        const std::string& cons_exit) {
        LosslessGraph g;
        auto n0 = g.add_node(static_cast<uint32_t>(cons_entry.size()));
        g.node(n0).consensus = cons_entry;
        auto n1 = g.add_node(static_cast<uint32_t>(cons_exit.size()));
        g.node(n1).consensus = cons_exit;
        g.add_edge(n0, n1, 10);
        return g;
    }

    BubbleCandidate make_candidate(uint32_t entry, uint32_t exit) {
        BubbleCandidate c{};
        c.entry_node = entry;
        c.exit_node = exit;
        c.read_support_branch = 5;
        c.read_support_alt = 5;
        c.bubble_length_bp = 100;
        return c;
    }
};

TEST_F(FlankJaccardTest, IdenticalFlanksReturnNearOne) {
    std::string seq;
    seq.reserve(200);
    for (int i = 0; i < 50; ++i) seq += "ACGT";

    auto g = make_graph_with_consensus(seq, seq);
    auto c = make_candidate(0, 1);

    auto fv = extract_features(c, g);
    float jaccard = fv[static_cast<size_t>(Feature::FlankJaccardK31)];

    EXPECT_NEAR(jaccard, 1.0f, 0.01f);
}

TEST_F(FlankJaccardTest, DisjointFlanksReturnZero) {
    std::string seq1(200, 'A');
    std::string seq2(200, 'C');

    auto g = make_graph_with_consensus(seq1, seq2);
    auto c = make_candidate(0, 1);

    auto fv = extract_features(c, g);
    float jaccard = fv[static_cast<size_t>(Feature::FlankJaccardK31)];

    EXPECT_NEAR(jaccard, 0.0f, 0.01f);
}

TEST_F(FlankJaccardTest, EmptyConsensusFallbackZero) {
    auto g = make_graph_with_consensus("", "");
    auto c = make_candidate(0, 1);

    auto fv = extract_features(c, g);
    float jaccard = fv[static_cast<size_t>(Feature::FlankJaccardK31)];

    EXPECT_EQ(jaccard, 0.0f);
}

// Original tests from test_features.cpp
TEST(FeaturesTest, kNumFeatures_is_12) {
    EXPECT_EQ(kNumFeatures, 12u);
}

TEST(FeaturesTest, FeatureVector_size_matches_kNumFeatures) {
    FeatureVector v{};
    EXPECT_EQ(v.size(), kNumFeatures);
}

TEST(FeaturesTest, BubbleCandidate_features_indexable_by_enum) {
    BubbleCandidate bc{};
    bc.features[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.0f;
    bc.features[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] = 1.0f;
    EXPECT_FLOAT_EQ(bc.features[1], 2.0f);
    EXPECT_FLOAT_EQ(bc.features[8], 1.0f);
}
