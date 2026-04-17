#include <gtest/gtest.h>
#include "branch/graph/graph_builder.hpp"
#include "branch/overlapper/overlap_pair.hpp"

using namespace branch::graph;
using namespace branch::overlapper;

TEST(GraphBuilderTest, BuildFromOverlaps) {
    std::vector<OverlapPair> overlaps = {
        {0, 1, 100, 200, 100, 200, true},
        {1, 2, 100, 200, 100, 200, true}
    };
    
    OverlapGraph graph = GraphBuilder::build_from_overlaps(overlaps, 3);
    
    EXPECT_EQ(graph.num_nodes(), 3);
    EXPECT_EQ(graph.num_edges(), 2);
}

TEST(GraphBuilderTest, NodesHaveReadSupportOne) {
    std::vector<OverlapPair> overlaps = {
        {0, 1, 100, 200, 100, 200, true}
    };
    
    OverlapGraph graph = GraphBuilder::build_from_overlaps(overlaps, 2);
    
    // Each node from a read should have read_support = 1
    const auto& nodes = graph.nodes();
    EXPECT_EQ(nodes.size(), 2);
    for (const auto& node : nodes) {
        EXPECT_EQ(node.read_support, 1);
    }
}

TEST(GraphBuilderTest, EdgeAggregation) {
    // Two overlaps between same pair of reads (0,1)
    std::vector<OverlapPair> overlaps = {
        {0, 1, 100, 200, 100, 200, true},
        {0, 1, 150, 250, 150, 250, true}  // same from/to, different positions
    };
    
    OverlapGraph graph = GraphBuilder::build_from_overlaps(overlaps, 2);
    
    // Should have 1 edge with read_support=2, not 2 edges with read_support=1
    EXPECT_EQ(graph.num_nodes(), 2);
    EXPECT_EQ(graph.num_edges(), 1);
    
    const auto& edges = graph.edges();
    EXPECT_EQ(edges[0].read_support, 2);
}
