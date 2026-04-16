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

// ============== ReadSpanCoverageIQR Tests ==============

class ReadSpanIQRTest : public ::testing::Test {
protected:
    LosslessGraph make_empty_graph() {
        LosslessGraph g;
        g.add_node(100);
        g.add_node(100);
        g.add_edge(0, 1, 10);
        return g;
    }

    BubbleCandidate make_candidate_with_spans(std::vector<uint32_t> spans) {
        BubbleCandidate c{};
        c.entry_node = 0;
        c.exit_node = 1;
        c.read_support_branch = 5;
        c.read_support_alt = 5;
        c.bubble_length_bp = 100;
        c.read_spans = std::move(spans);
        return c;
    }
};

TEST_F(ReadSpanIQRTest, FourElementsComputesIQR) {
    // spans = {100, 200, 300, 400}
    // sorted: 100, 200, 300, 400
    // Q1 at index 0.75 = 100*0.25 + 200*0.75 = 175
    // Q3 at index 2.25 = 300*0.75 + 400*0.25 = 325  
    // IQR = 325 - 175 = 150
    auto g = make_empty_graph();
    auto c = make_candidate_with_spans({100, 200, 300, 400});
    auto fv = extract_features(c, g);
    float iqr = fv[static_cast<size_t>(Feature::ReadSpanCoverageIQR)];
    EXPECT_NEAR(iqr, 150.0f, 1.0f);
}

TEST_F(ReadSpanIQRTest, LessThanFourReturnsZero) {
    auto g = make_empty_graph();
    auto c = make_candidate_with_spans({100, 200, 300});
    auto fv = extract_features(c, g);
    float iqr = fv[static_cast<size_t>(Feature::ReadSpanCoverageIQR)];
    EXPECT_EQ(iqr, 0.0f);
}

TEST_F(ReadSpanIQRTest, EmptySpansReturnsZero) {
    auto g = make_empty_graph();
    auto c = make_candidate_with_spans({});
    auto fv = extract_features(c, g);
    float iqr = fv[static_cast<size_t>(Feature::ReadSpanCoverageIQR)];
    EXPECT_EQ(iqr, 0.0f);
}

// ============== GcContentDivergence Tests ==============

class GcContentTest : public ::testing::Test {
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

    BubbleCandidate make_candidate() {
        BubbleCandidate c{};
        c.entry_node = 0;
        c.exit_node = 1;
        c.read_support_branch = 5;
        c.read_support_alt = 5;
        c.bubble_length_bp = 100;
        return c;
    }
};

TEST_F(GcContentTest, IdenticalGCReturnsZero) {
    // Both 50% GC
    auto g = make_graph_with_consensus("ACGT", "ACGT");
    auto c = make_candidate();
    auto fv = extract_features(c, g);
    float gc_div = fv[static_cast<size_t>(Feature::GcContentDivergence)];
    EXPECT_NEAR(gc_div, 0.0f, 0.01f);
}

TEST_F(GcContentTest, DifferentGCReturnsDivergence) {
    // Entry: GGGG = 100% GC, Exit: AAAA = 0% GC -> divergence = 1.0
    auto g = make_graph_with_consensus("GGGG", "AAAA");
    auto c = make_candidate();
    auto fv = extract_features(c, g);
    float gc_div = fv[static_cast<size_t>(Feature::GcContentDivergence)];
    EXPECT_NEAR(gc_div, 1.0f, 0.01f);
}

TEST_F(GcContentTest, EmptyConsensusReturnsZero) {
    auto g = make_graph_with_consensus("", "");
    auto c = make_candidate();
    auto fv = extract_features(c, g);
    float gc_div = fv[static_cast<size_t>(Feature::GcContentDivergence)];
    EXPECT_EQ(gc_div, 0.0f);
}
