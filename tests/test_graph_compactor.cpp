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

TEST(GraphCompactorTest, Tournament_clique_of_equal_length_nodes_merges_to_one_unitig) {
    // Regression test for P2.3: the synthetic 4-allele fixture (and real
    // read-overlap graphs generated from equal-length reads of a single
    // allele) produces a tournament-shaped subgraph where every pair of
    // nodes has a directed edge a→b (with a < b, the canonical ordering
    // emitted by cpu_backend). In a tournament every node has
    // in_degree>=1 and out_degree>=1, so the linear-chain pass cannot
    // touch them — each remains its own unitig. Before this fix the
    // compactor left the clique unchanged (N nodes → N unitigs); now
    // the closed-twin pass collapses them to a single unitig.
    //
    // Guarding criterion: closed neighborhoods (N(v) ∪ {v}) agree across
    // every member of the clique, AND the members share length_bp (the
    // "different-allele firewall" that stops a 2-copy and a 6-copy read
    // from collapsing into the same unitig when they happen to share
    // neighbors).
    LosslessGraph g;
    constexpr std::uint32_t kLen = 19000;
    constexpr int kCliqueSize = 6;
    for (int i = 0; i < kCliqueSize; ++i) g.add_node(kLen, 1);
    // Emit canonical-ordered directed edges for every pair.
    for (int a = 0; a < kCliqueSize; ++a) {
        for (int b = a + 1; b < kCliqueSize; ++b) {
            g.add_edge(static_cast<NodeId>(a), static_cast<NodeId>(b), 1);
        }
    }

    CompactionResult r = compact_unitigs(g);
    EXPECT_EQ(r.compacted.node_count(), 1u)
        << "closed-twin merge should collapse a tournament of equal-length "
           "nodes into a single unitig";
    // No inter-unitig edges: all 15 clique edges are now intra-unitig.
    EXPECT_EQ(r.compacted.edge_count(), 0u);
    // The merged unitig takes the representative (shared) length.
    EXPECT_EQ(r.compacted.node(0).length_bp, kLen);
    for (int i = 0; i < kCliqueSize; ++i) {
        EXPECT_EQ(r.old_to_new[i], 0u)
            << "clique member " << i << " should map to the merged unitig";
    }
}

TEST(GraphCompactorTest, Closed_twins_of_different_length_do_not_merge) {
    // Different-allele firewall. Two nodes that share every neighbor
    // (closed-twin topology) but have distinct length_bp must remain
    // separate — otherwise a 9kb k=2 read and a 19kb k=6 read would
    // collapse into one unitig, destroying the allele distinction.
    LosslessGraph g;
    NodeId a = g.add_node(9000, 1);    // different length
    NodeId b = g.add_node(19000, 1);
    NodeId c = g.add_node(19000, 1);
    NodeId d = g.add_node(19000, 1);
    // a, b, c, d form a complete tournament on 4 equal-length nodes
    // EXCEPT a has a different length. b,c,d are closed-twins; a is not.
    g.add_edge(a, b, 1); g.add_edge(a, c, 1); g.add_edge(a, d, 1);
    g.add_edge(b, c, 1); g.add_edge(b, d, 1);
    g.add_edge(c, d, 1);

    CompactionResult r = compact_unitigs(g);
    EXPECT_EQ(r.compacted.node_count(), 2u)
        << "b,c,d should collapse; a stays separate due to length mismatch";
    EXPECT_NE(r.old_to_new[a], r.old_to_new[b]);
    EXPECT_EQ(r.old_to_new[b], r.old_to_new[c]);
    EXPECT_EQ(r.old_to_new[b], r.old_to_new[d]);
}

TEST(GraphCompactorTest, Isolated_nodes_merge_by_length_into_orphan_bundles) {
    // Regression test for P2.3: graph_filter's containment drop
    // removes every edge touching a contained node but cannot renumber
    // NodeIds, so the node entries linger with in_degree=0 AND
    // out_degree=0. Before this fix each isolated node remained its own
    // singleton unitig, bloating the compacted graph 1:1 with dropped
    // reads (73 extra unitigs on the 80-read fixture). Now they fold
    // into one orphan unitig per distinct length_bp — preserving the
    // coverage signal that length-bucket VAF aggregation reads out,
    // without the 1:1 node blow-up.
    LosslessGraph g;
    // Three reads of length 4000 (k=0 family), five of length 9000
    // (k=2 family), no edges among them. All eight end up isolated.
    for (int i = 0; i < 3; ++i) g.add_node(4000, 1);
    for (int i = 0; i < 5; ++i) g.add_node(9000, 1);

    CompactionResult r = compact_unitigs(g);
    EXPECT_EQ(r.compacted.node_count(), 2u)
        << "two distinct isolated-read lengths should yield two orphan bundles";

    // The two orphan bundles carry the representative length of their
    // members, not the sum.
    std::vector<std::uint32_t> observed;
    for (NodeId u = 0; u < r.compacted.node_count(); ++u) {
        observed.push_back(r.compacted.node(u).length_bp);
    }
    std::sort(observed.begin(), observed.end());
    EXPECT_EQ(observed, (std::vector<std::uint32_t>{4000, 9000}));

    // Members with matching length_bp share a unitig id.
    EXPECT_EQ(r.old_to_new[0], r.old_to_new[1]);
    EXPECT_EQ(r.old_to_new[0], r.old_to_new[2]);
    EXPECT_EQ(r.old_to_new[3], r.old_to_new[4]);
    EXPECT_EQ(r.old_to_new[3], r.old_to_new[7]);
    EXPECT_NE(r.old_to_new[0], r.old_to_new[3]);
}

TEST(GraphCompactorTest, Fourway_synthetic_fixture_collapses_to_small_graph) {
    // Integration-shaped regression: simulate the topology the
    // assemble pipeline feeds to the compactor after graph_filter on
    // the synthetic 4-allele fixture (copies {0,2,4,6}, 80 reads,
    // cassette 2500bp, read_len 19000):
    //   - 7 surviving nodes after containment drop, forming a
    //     tournament: 1 k=2 node (len 9000) with outgoing edges to all
    //     6 k=6 nodes (len 19000), plus the canonical a<b directed
    //     edges among the k=6 peers (tournament of 6).
    //   - 73 isolated nodes (reads dropped by containment but still
    //     present in the node vector). Split by length to mirror the
    //     real fixture: 12 × 4000 (k=0), 44 × 9000 (k=2 minus node 0),
    //     17 × 14000 (k=4).
    //
    // Before P2.3 fix: 80 raw nodes → 80 compacted nodes (no-op).
    // After: the tournament collapses (closed-twin merge, length gates
    // node-0 out); the isolated nodes fold to 3 orphan bundles by
    // length. Expected: 5 compacted nodes.
    LosslessGraph g;

    // Active survivor set: indices [0..6]. Node 0 is k=2 (9000); 1..6 k=6 (19000).
    const NodeId node0 = g.add_node(9000, 1);
    std::vector<NodeId> k6;
    for (int i = 0; i < 6; ++i) k6.push_back(g.add_node(19000, 1));
    for (NodeId v : k6) g.add_edge(node0, v, 1);
    for (std::size_t a = 0; a < k6.size(); ++a) {
        for (std::size_t b = a + 1; b < k6.size(); ++b) {
            g.add_edge(k6[a], k6[b], 1);
        }
    }

    // Isolated pool mirroring the real fixture's containment-dropped set.
    for (int i = 0; i < 12; ++i) g.add_node(4000, 1);
    for (int i = 0; i < 44; ++i) g.add_node(9000, 1);
    for (int i = 0; i < 17; ++i) g.add_node(14000, 1);

    CompactionResult r = compact_unitigs(g);
    EXPECT_LE(r.compacted.node_count(), 10u)
        << "fourway fixture should compact to ≤10 nodes (P2.3 deliverable)";
    EXPECT_EQ(r.compacted.node_count(), 5u)
        << "expected node-0 + merged-k6-cluster + 3 per-length orphan bundles";
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
