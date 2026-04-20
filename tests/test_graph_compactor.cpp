#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

#include "graph/graph_compactor.hpp"
#include "graph/lossless_graph.hpp"

using branch::graph::compact_unitigs;
using branch::graph::compact_unitigs_with_sequences;
using branch::graph::CompactionResult;
using branch::graph::LosslessGraph;
using branch::graph::NodeId;

TEST(GraphCompactorTest, Empty_graph_compacts_to_empty_graph) {
    LosslessGraph g;
    CompactionResult r = compact_unitigs(g);
    EXPECT_EQ(r.compacted.node_count(), 0u);
    EXPECT_EQ(r.compacted.edge_count(), 0u);
    EXPECT_TRUE(r.old_to_new.empty());
}

TEST(GraphCompactorTest, Linear_chain_collapses_to_single_unitig) {
    // A(100) -> B(200) -> C(300) -> D(400), copy_count of A is 2
    LosslessGraph g;
    NodeId a = g.add_node(100, 2);
    NodeId b = g.add_node(200, 1);
    NodeId c = g.add_node(300, 1);
    NodeId d = g.add_node(400, 1);
    g.add_edge(a, b, 10);
    g.add_edge(b, c, 20);
    g.add_edge(c, d, 30);

    CompactionResult r = compact_unitigs(g);

    // 1 unitig, no edges (all intra-chain).
    ASSERT_EQ(r.compacted.node_count(), 1u);
    EXPECT_EQ(r.compacted.edge_count(), 0u);

    // Length-sum and copy_count inherited from chain start A.
    EXPECT_EQ(r.compacted.node(0).length_bp, 100u + 200u + 300u + 400u);
    EXPECT_EQ(r.compacted.node(0).copy_count, 2u);

    // All old nodes map to the same unitig 0.
    ASSERT_EQ(r.old_to_new.size(), 4u);
    EXPECT_EQ(r.old_to_new[a], 0u);
    EXPECT_EQ(r.old_to_new[b], 0u);
    EXPECT_EQ(r.old_to_new[c], 0u);
    EXPECT_EQ(r.old_to_new[d], 0u);
}

TEST(GraphCompactorTest, Chain_then_branch_yields_two_unitigs) {
    // A(10) -> B(20) -> C(30) -> { D1(40), D2(50) }
    LosslessGraph g;
    NodeId a  = g.add_node(10, 3);
    NodeId b  = g.add_node(20, 1);
    NodeId c  = g.add_node(30, 1);
    NodeId d1 = g.add_node(40, 1);
    NodeId d2 = g.add_node(50, 1);
    g.add_edge(a, b, 1);
    g.add_edge(b, c, 2);
    g.add_edge(c, d1, 7);
    g.add_edge(c, d2, 9);

    CompactionResult r = compact_unitigs(g);

    ASSERT_EQ(r.compacted.node_count(), 3u);
    ASSERT_EQ(r.compacted.edge_count(), 2u);

    // A, B, C share one unitig id.
    NodeId u_abc = r.old_to_new[a];
    EXPECT_EQ(r.old_to_new[b], u_abc);
    EXPECT_EQ(r.old_to_new[c], u_abc);

    // D1 and D2 each have their own (distinct) unitig.
    NodeId u_d1 = r.old_to_new[d1];
    NodeId u_d2 = r.old_to_new[d2];
    EXPECT_NE(u_d1, u_abc);
    EXPECT_NE(u_d2, u_abc);
    EXPECT_NE(u_d1, u_d2);

    // ABC unitig has length 10+20+30 = 60 and copy_count from A = 3.
    EXPECT_EQ(r.compacted.node(u_abc).length_bp, 60u);
    EXPECT_EQ(r.compacted.node(u_abc).copy_count, 3u);
    // D1 / D2 kept as singletons.
    EXPECT_EQ(r.compacted.node(u_d1).length_bp, 40u);
    EXPECT_EQ(r.compacted.node(u_d2).length_bp, 50u);

    // Two inter-unitig edges, both starting from the ABC unitig.
    bool saw_d1 = false, saw_d2 = false;
    for (const auto& e : r.compacted.edges()) {
        EXPECT_EQ(e.from, u_abc);
        if (e.to == u_d1) {
            EXPECT_EQ(e.read_support, 7u);
            saw_d1 = true;
        } else if (e.to == u_d2) {
            EXPECT_EQ(e.read_support, 9u);
            saw_d2 = true;
        } else {
            FAIL() << "unexpected edge target " << e.to;
        }
    }
    EXPECT_TRUE(saw_d1);
    EXPECT_TRUE(saw_d2);
}

TEST(GraphCompactorTest, Consensus_from_three_reads_with_snp) {
    // Test: 3 reads with 1 SNP position (majority voting)
    // Read 1: ACGTACGT
    // Read 2: ACGAACGT  (T->A at position 3)
    // Read 3: ACGAACGT  (T->A at position 3)
    // Expected consensus: ACGAACGT (A wins 2:1 at position 3)
    
    LosslessGraph g;
    NodeId n0 = g.add_node(8, 1);  // length 8bp
    NodeId n1 = g.add_node(8, 1);
    NodeId n2 = g.add_node(8, 1);
    
    // Linear chain: n0 -> n1 -> n2
    g.add_edge(n0, n1, 1);
    g.add_edge(n1, n2, 1);
    
    // Provide sequences for each node
    std::vector<std::string> sequences = {
        "ACGTACGT",  // node 0
        "ACGAACGT",  // node 1 (SNP at pos 3: T->A)
        "ACGAACGT"   // node 2 (SNP at pos 3: T->A)
    };
    
    CompactionResult r = compact_unitigs_with_sequences(g, sequences);
    
    // Should collapse to 1 unitig
    ASSERT_EQ(r.compacted.node_count(), 1u);
    
    // Consensus should be ACGAACGT (majority wins)
    const std::string& cons = r.compacted.node(0).consensus;
    EXPECT_FALSE(cons.empty()) << "Consensus should not be empty";
    
    // At position 3, 'A' should win (2 votes) over 'T' (1 vote)
    if (cons.size() >= 4) {
        EXPECT_EQ(cons[3], 'A') << "Position 3 should be 'A' (majority)";
    }
    
    // Full consensus check
    EXPECT_EQ(cons, "ACGAACGT") << "Full consensus mismatch";
}

TEST(GraphCompactorTest, Single_read_consensus_direct_copy) {
    // Single read: consensus is the read itself
    LosslessGraph g;
    (void)g.add_node(5, 1);
    
    std::vector<std::string> sequences = { "ATCGA" };
    
    CompactionResult r = compact_unitigs_with_sequences(g, sequences);
    
    ASSERT_EQ(r.compacted.node_count(), 1u);
    EXPECT_EQ(r.compacted.node(0).consensus, "ATCGA");
}

TEST(GraphCompactorTest, Unitig_read_support_sums_across_members) {
    // Linear chain A->B->C->D with per-node read_support {3, 5, 1, 4}.
    // The compactor should collapse them into one unitig whose
    // read_support is 3+5+1+4 = 13, matching the RC:i that will be
    // emitted on the S-line.
    LosslessGraph g;
    NodeId a = g.add_node(10, 1);
    NodeId b = g.add_node(20, 1);
    NodeId c = g.add_node(30, 1);
    NodeId d = g.add_node(40, 1);
    g.node(a).read_support = 3;
    g.node(b).read_support = 5;
    g.node(c).read_support = 1;
    g.node(d).read_support = 4;
    g.add_edge(a, b, 1);
    g.add_edge(b, c, 1);
    g.add_edge(c, d, 1);

    CompactionResult r = compact_unitigs(g);
    ASSERT_EQ(r.compacted.node_count(), 1u);
    EXPECT_EQ(r.compacted.node(0).read_support, 13u);
}

TEST(GraphCompactorTest, Unitig_read_support_sums_in_with_sequences_overload) {
    // Same propagation invariant on the parallel overload used by the
    // CLI. Three-node chain with read_support {2, 2, 2} -> 6.
    LosslessGraph g;
    NodeId n0 = g.add_node(5, 1);
    NodeId n1 = g.add_node(5, 1);
    NodeId n2 = g.add_node(5, 1);
    g.node(n0).read_support = 2;
    g.node(n1).read_support = 2;
    g.node(n2).read_support = 2;
    g.add_edge(n0, n1, 1);
    g.add_edge(n1, n2, 1);

    std::vector<std::string> seqs = {"AAAAA", "AAAAA", "AAAAA"};
    CompactionResult r = compact_unitigs_with_sequences(g, seqs);
    ASSERT_EQ(r.compacted.node_count(), 1u);
    EXPECT_EQ(r.compacted.node(0).read_support, 6u);
}

TEST(GraphCompactorTest, Two_reads_consensus_majority) {
    // Two reads with 1 difference
    // Read 1: AAAA
    // Read 2: AATA
    // Expected: AATA or AAAA (tie at position 2, either valid)
    
    LosslessGraph g;
    NodeId n0 = g.add_node(4, 1);
    NodeId n1 = g.add_node(4, 1);
    g.add_edge(n0, n1, 1);
    
    std::vector<std::string> sequences = { "AAAA", "AATA" };
    
    CompactionResult r = compact_unitigs_with_sequences(g, sequences);
    
    ASSERT_EQ(r.compacted.node_count(), 1u);
    const std::string& cons = r.compacted.node(0).consensus;
    EXPECT_FALSE(cons.empty());
    EXPECT_EQ(cons.size(), 4u);
    // First, second, fourth positions should be 'A'
    EXPECT_EQ(cons[0], 'A');
    EXPECT_EQ(cons[1], 'A');
    EXPECT_EQ(cons[3], 'A');
    // Position 2 could be 'A' or 'T' (tie)
    EXPECT_TRUE(cons[2] == 'A' || cons[2] == 'T');
}
