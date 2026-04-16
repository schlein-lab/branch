#include <gtest/gtest.h>

#include <sstream>

#include "graph/graph_io.hpp"

using branch::graph::LosslessGraph;
using branch::graph::NodeId;
using branch::graph::write_gfa;
using branch::graph::read_gfa;

TEST(GraphIOTest, WriteGfa_produces_header_and_S_lines) {
    LosslessGraph g;
    g.add_node(1000, 2);
    g.add_node(500);
    g.add_edge(0, 1, 30);

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));

    std::string gfa = oss.str();
    EXPECT_NE(gfa.find("H\tVN:Z:1.2"), std::string::npos);
    EXPECT_NE(gfa.find("S\t0\t*\tLN:i:1000"), std::string::npos);
    EXPECT_NE(gfa.find("CN:i:2"), std::string::npos);
    EXPECT_NE(gfa.find("L\t0\t+\t1\t+\t0M"), std::string::npos);
    EXPECT_NE(gfa.find("VC:i:30"), std::string::npos);
}

TEST(GraphIOTest, Roundtrip_preserves_graph_structure) {
    LosslessGraph g;
    g.add_node(2000, 3);
    g.add_node(1500);
    g.add_edge(0, 1, 42);
    g.set_copy_count(0, 3, 0.95f);
    g.set_edge_vaf(0, 0.07f, 0.88f);

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));

    LosslessGraph g2;
    std::istringstream iss(oss.str());
    ASSERT_TRUE(read_gfa(g2, iss));

    EXPECT_EQ(g2.node_count(), 2u);
    EXPECT_EQ(g2.edge_count(), 1u);
    EXPECT_EQ(g2.node(0).length_bp, 2000u);
    EXPECT_EQ(g2.node(0).copy_count, 3u);
    EXPECT_NEAR(g2.node(0).copy_count_confidence, 0.95f, 0.01f);
    EXPECT_EQ(g2.edges()[0].from, 0u);
    EXPECT_EQ(g2.edges()[0].to, 1u);
    EXPECT_EQ(g2.edges()[0].read_support, 42u);
    EXPECT_NEAR(g2.edges()[0].vaf, 0.07f, 0.001f);
    EXPECT_NEAR(g2.edges()[0].vaf_confidence, 0.88f, 0.01f);
}

TEST(GraphIOTest, ReadGfa_skips_unknown_lines_gracefully) {
    std::string gfa = "H\tVN:Z:1.2\n"
                      "# this is a comment\n"
                      "S\t0\t*\tLN:i:100\n"
                      "W\tunknown\twalk\tline\n"
                      "S\t1\t*\tLN:i:200\n";
    LosslessGraph g;
    std::istringstream iss(gfa);
    ASSERT_TRUE(read_gfa(g, iss));
    EXPECT_EQ(g.node_count(), 2u);
}
