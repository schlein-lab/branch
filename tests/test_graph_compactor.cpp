#include <gtest/gtest.h>
#include "branch/graph/graph_compactor.hpp"
#include "branch/graph/overlap_graph.hpp"

using namespace branch::graph;

TEST(GraphCompactorTest, FindUnitigChains) {
    OverlapGraph graph;
    
    // Linear chain: 0 -> 1 -> 2
    graph.add_node(100, 1);
    graph.add_node(100, 1);
    graph.add_node(100, 1);
    graph.add_edge(0, 1, 1);
    graph.add_edge(1, 2, 1);
    
    auto chains = GraphCompactor::find_unitig_chains(graph);
    
    // Should find one chain
    EXPECT_GE(chains.size(), 1);
}

TEST(GraphCompactorTest, CompactLinearChain) {
    OverlapGraph graph;
    
    // Linear chain: 0 -> 1 -> 2
    graph.add_node(100, 1);
    graph.add_node(200, 1);
    graph.add_node(150, 1);
    graph.add_edge(0, 1, 1);
    graph.add_edge(1, 2, 1);
    
    std::vector<std::string> seqs = {"AAA", "BBB", "CCC"};
    auto compacted = GraphCompactor::compact_unitigs_with_sequences(graph, seqs);
    
    // Should compact to 1 node
    EXPECT_EQ(compacted.num_nodes(), 1);
    EXPECT_EQ(compacted.num_edges(), 0);
    
    // Total length should be sum: 100 + 200 + 150 = 450
    EXPECT_EQ(compacted.nodes()[0].length_bp, 450);
}

TEST(GraphCompactorTest, CompactedNodeSumsReadSupport) {
    OverlapGraph graph;
    
    // Linear chain: 0 -> 1 -> 2, each with read_support=1
    graph.add_node(100, 1);  // RC=1
    graph.add_node(100, 1);  // RC=1
    graph.add_node(100, 1);  // RC=1
    graph.add_edge(0, 1, 1);
    graph.add_edge(1, 2, 1);
    
    std::vector<std::string> seqs = {"AAA", "BBB", "CCC"};
    auto compacted = GraphCompactor::compact_unitigs_with_sequences(graph, seqs);
    
    // Should compact to 1 node with read_support = 3 (sum of chain)
    EXPECT_EQ(compacted.num_nodes(), 1);
    EXPECT_EQ(compacted.nodes()[0].read_support, 3);
}

TEST(GraphCompactorTest, BranchingPreservesNodes) {
    OverlapGraph graph;
    
    // Branching: 0 -> 1, 0 -> 2
    graph.add_node(100, 1);
    graph.add_node(100, 1);
    graph.add_node(100, 1);
    graph.add_edge(0, 1, 1);
    graph.add_edge(0, 2, 1);
    
    std::vector<std::string> seqs = {"AAA", "BBB", "CCC"};
    auto compacted = GraphCompactor::compact_unitigs_with_sequences(graph, seqs);
    
    // No compaction possible with branching
    EXPECT_EQ(compacted.num_nodes(), 3);
    EXPECT_EQ(compacted.num_edges(), 2);
}
