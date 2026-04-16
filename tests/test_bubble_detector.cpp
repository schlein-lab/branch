#include <gtest/gtest.h>

#include "detect/bubble_detector.hpp"
#include "graph/lossless_graph.hpp"

using branch::detect::Bubble;
using branch::detect::BubbleDetectorConfig;
using branch::detect::detect_bubbles;
using branch::graph::LosslessGraph;
using branch::graph::NodeId;

TEST(BubbleDetectorTest, Empty_graph_yields_no_bubbles) {
    LosslessGraph g;
    auto bubbles = detect_bubbles(g);
    EXPECT_EQ(bubbles.size(), 0u);
}

TEST(BubbleDetectorTest, Linear_chain_yields_no_bubbles) {
    LosslessGraph g;
    NodeId a = g.add_node(100);
    NodeId b = g.add_node(100);
    NodeId c = g.add_node(100);
    g.add_edge(a, b, 10);
    g.add_edge(b, c, 10);
    auto bubbles = detect_bubbles(g);
    EXPECT_EQ(bubbles.size(), 0u);
}

TEST(BubbleDetectorTest, Simple_binary_bubble_detected) {
    // Topology:
    //   entry --> mid1 --> exit
    //         `-> mid2 --'
    LosslessGraph g;
    NodeId entry = g.add_node(100);
    NodeId mid1  = g.add_node(100);
    NodeId mid2  = g.add_node(100);
    NodeId exit  = g.add_node(100);
    g.add_edge(entry, mid1, 10);
    g.add_edge(entry, mid2, 10);
    g.add_edge(mid1, exit, 10);
    g.add_edge(mid2, exit, 10);

    auto bubbles = detect_bubbles(g);
    ASSERT_EQ(bubbles.size(), 1u);
    EXPECT_EQ(bubbles[0].entry, entry);
    EXPECT_EQ(bubbles[0].exit, exit);
    EXPECT_EQ(bubbles[0].alts.size(), 2u);
}

TEST(BubbleDetectorTest, Direct_both_edges_to_exit_detected) {
    // entry --> exit with two parallel direct edges.
    // The reachable-set approach places exit itself in the reachable
    // set of each fan-out target, so this shape is picked up.
    LosslessGraph g;
    NodeId entry = g.add_node(100);
    NodeId mid1  = g.add_node(100);
    NodeId mid2  = g.add_node(100);
    NodeId exit  = g.add_node(100);
    g.add_edge(entry, mid1);
    g.add_edge(entry, mid2);
    g.add_edge(mid1, exit);
    g.add_edge(mid2, exit);
    auto bubbles = detect_bubbles(g);
    EXPECT_GE(bubbles.size(), 1u);
}

TEST(BubbleDetectorTest, Max_alt_path_length_caps_detection) {
    // A very long parallel path must be filtered by max_alt_path_length.
    LosslessGraph g;
    NodeId entry = g.add_node(100);
    NodeId prev1 = g.add_node(100);
    NodeId prev2 = g.add_node(100);
    g.add_edge(entry, prev1);
    g.add_edge(entry, prev2);
    // Extend each side with 20 nodes (beyond default cap of 8).
    for (int i = 0; i < 20; ++i) {
        NodeId n1 = g.add_node(100);
        NodeId n2 = g.add_node(100);
        g.add_edge(prev1, n1);
        g.add_edge(prev2, n2);
        prev1 = n1;
        prev2 = n2;
    }
    NodeId exit = g.add_node(100);
    g.add_edge(prev1, exit);
    g.add_edge(prev2, exit);

    BubbleDetectorConfig cfg{.max_alt_path_length = 4};  // too small
    auto bubbles = detect_bubbles(g, cfg);
    EXPECT_EQ(bubbles.size(), 0u);
}
