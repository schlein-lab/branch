#include <gtest/gtest.h>

#include "graph/lossless_graph.hpp"

using branch::graph::Edge;
using branch::graph::LosslessGraph;
using branch::graph::Node;
using branch::graph::NodeId;

TEST(LosslessGraphTest, Edge_packs_to_24_bytes) {
    EXPECT_EQ(sizeof(Edge), 24u);
}

TEST(LosslessGraphTest, Empty_graph_has_zero_counts) {
    LosslessGraph g;
    EXPECT_EQ(g.node_count(), 0u);
    EXPECT_EQ(g.edge_count(), 0u);
}

TEST(LosslessGraphTest, add_node_returns_sequential_ids) {
    LosslessGraph g;
    NodeId a = g.add_node(1000);
    NodeId b = g.add_node(2000);
    NodeId c = g.add_node(3000, 3);
    EXPECT_EQ(a, 0u);
    EXPECT_EQ(b, 1u);
    EXPECT_EQ(c, 2u);
    EXPECT_EQ(g.node_count(), 3u);
    EXPECT_EQ(g.node(c).copy_count, 3u);
    EXPECT_EQ(g.node(c).length_bp, 3000u);
}

TEST(LosslessGraphTest, add_edge_stores_read_support) {
    LosslessGraph g;
    NodeId a = g.add_node(500);
    NodeId b = g.add_node(500);
    g.add_edge(a, b, 30);
    ASSERT_EQ(g.edge_count(), 1u);
    EXPECT_EQ(g.edges()[0].from, a);
    EXPECT_EQ(g.edges()[0].to, b);
    EXPECT_EQ(g.edges()[0].read_support, 30u);
    EXPECT_FLOAT_EQ(g.edges()[0].vaf, 1.0f);
}

TEST(LosslessGraphTest, set_copy_count_and_vaf_update_in_place) {
    LosslessGraph g;
    NodeId a = g.add_node(1000);
    NodeId b = g.add_node(1000);
    g.add_edge(a, b);

    g.set_copy_count(a, 4, 0.85f);
    g.set_edge_vaf(0, 0.05f, 0.9f);

    EXPECT_EQ(g.node(a).copy_count, 4u);
    EXPECT_FLOAT_EQ(g.node(a).copy_count_confidence, 0.85f);
    EXPECT_FLOAT_EQ(g.edges()[0].vaf, 0.05f);
    EXPECT_FLOAT_EQ(g.edges()[0].vaf_confidence, 0.9f);
}
