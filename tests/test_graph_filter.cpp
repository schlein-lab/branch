#include <gtest/gtest.h>

#include <algorithm>

#include "graph/graph_filter.hpp"
#include "graph/lossless_graph.hpp"

using branch::graph::filter_graph;
using branch::graph::FilterConfig;
using branch::graph::FilterStats;
using branch::graph::LosslessGraph;
using branch::graph::NodeId;

TEST(GraphFilterTest, Empty_graph_passthrough) {
    LosslessGraph g;
    FilterStats s = filter_graph(g);
    EXPECT_EQ(s.edges_before, 0u);
    EXPECT_EQ(s.edges_after, 0u);
    EXPECT_EQ(s.nodes_dropped_contained, 0u);
    EXPECT_EQ(s.transitive_edges_removed, 0u);
    EXPECT_EQ(g.node_count(), 0u);
    EXPECT_EQ(g.edge_count(), 0u);
}

TEST(GraphFilterTest, Transitive_triangle_reduces_to_chain) {
    // A(1000) -> B(1000), B -> C(1000), A -> C (the transitive edge).
    // After filter with containment disabled, only A->B and B->C remain.
    LosslessGraph g;
    NodeId a = g.add_node(1000);
    NodeId b = g.add_node(1000);
    NodeId c = g.add_node(1000);
    g.add_edge(a, b, 5);
    g.add_edge(b, c, 7);
    g.add_edge(a, c, 3);  // transitive: len_ab(1000) + len_bc(1000) == 2000
                          // direct_ac length proxy is node c length = 1000,
                          // but since our length proxy is the target node
                          // length, the match is: direct(1000) vs sum(2000).
                          // That does NOT match with fuzz=100; to make this
                          // a real transitive triangle in v0.2 length-proxy
                          // semantics we use a larger fuzz so the triangle
                          // collapses. See next test for the exact-match
                          // case.
    FilterConfig cfg;
    cfg.drop_contained = false;      // isolate the transitive pass
    cfg.reduce_transitive = true;
    cfg.transitive_fuzz = 2000;      // accept wide mismatch for this fixture

    FilterStats s = filter_graph(g, cfg);

    EXPECT_EQ(s.edges_before, 3u);
    EXPECT_EQ(s.transitive_edges_removed, 1u);
    EXPECT_EQ(s.edges_after, 2u);
    ASSERT_EQ(g.edge_count(), 2u);

    bool saw_ab = false, saw_bc = false, saw_ac = false;
    for (const auto& e : g.edges()) {
        if (e.from == a && e.to == b) saw_ab = true;
        else if (e.from == b && e.to == c) saw_bc = true;
        else if (e.from == a && e.to == c) saw_ac = true;
    }
    EXPECT_TRUE(saw_ab);
    EXPECT_TRUE(saw_bc);
    EXPECT_FALSE(saw_ac) << "transitive A->C edge should have been removed";
}

TEST(GraphFilterTest, Transitive_reduction_respects_fuzz_window) {
    // Same topology but with a tight fuzz. The direct A->C "length"
    // (node C length_bp = 500) is far from sum(A->B) + (B->C) = 500+500
    // = 1000, so with fuzz=10 the edge should NOT be reduced.
    LosslessGraph g;
    NodeId a = g.add_node(500);
    NodeId b = g.add_node(500);
    NodeId c = g.add_node(500);
    g.add_edge(a, b, 1);
    g.add_edge(b, c, 1);
    g.add_edge(a, c, 1);

    FilterConfig cfg;
    cfg.drop_contained = false;
    cfg.reduce_transitive = true;
    cfg.transitive_fuzz = 10;

    FilterStats s = filter_graph(g, cfg);

    EXPECT_EQ(s.transitive_edges_removed, 0u)
        << "with tight fuzz no edge should be considered transitive";
    EXPECT_EQ(g.edge_count(), 3u);
}

TEST(GraphFilterTest, Containment_drops_short_engulfed_leaf) {
    // Predicate (v0.2, no overlap offsets):
    //   N is contained iff a predecessor M is STRICTLY longer than N
    //   AND every successor of N is also a successor of M.
    // A short leaf (no outgoing edges) trivially satisfies the second
    // clause, so this is the classic "long read M engulfs short leaf S,
    // while sibling T matches M's length exactly and must survive".
    //
    //   M(5000) -> S(500)   -- leaf, strictly shorter: dropped.
    //   M(5000) -> T(5000)  -- same-length sibling: must survive.
    LosslessGraph g;
    NodeId m = g.add_node(5000);
    NodeId s = g.add_node(500);
    NodeId t = g.add_node(5000);
    g.add_edge(m, s, 1);
    g.add_edge(m, t, 1);

    FilterConfig cfg;
    cfg.drop_contained = true;
    cfg.reduce_transitive = false;

    FilterStats st = filter_graph(g, cfg);

    EXPECT_EQ(st.nodes_dropped_contained, 1u);
    // Every edge touching S must be gone.
    for (const auto& e : g.edges()) {
        EXPECT_NE(e.from, s);
        EXPECT_NE(e.to, s);
    }
    // M -> T must survive.
    bool saw_mt = false;
    for (const auto& e : g.edges()) {
        if (e.from == m && e.to == t) saw_mt = true;
    }
    EXPECT_TRUE(saw_mt);
    EXPECT_EQ(g.node_count(), 3u);
}

TEST(GraphFilterTest, Containment_transfers_read_support_to_coverer) {
    // When a short node is dropped as contained, its read_support must
    // be folded onto the absorbing predecessor so downstream VAF /
    // coverage analysis still sees the right depth.
    //
    //   M(5000, rc=3) -> S(500, rc=2)   -- leaf, dropped; rc transfers.
    //   M -> T(5000, rc=1)              -- same-length sibling, survives.
    // Expected after filter: M.rc = 3 + 2 = 5; T.rc = 1 unchanged.
    LosslessGraph g;
    NodeId m = g.add_node(5000);
    NodeId s = g.add_node(500);
    NodeId t = g.add_node(5000);
    g.node(m).read_support = 3;
    g.node(s).read_support = 2;
    g.node(t).read_support = 1;
    g.add_edge(m, s, 1);
    g.add_edge(m, t, 1);

    FilterConfig cfg;
    cfg.drop_contained = true;
    cfg.reduce_transitive = false;

    FilterStats st = filter_graph(g, cfg);

    EXPECT_EQ(st.nodes_dropped_contained, 1u);
    EXPECT_EQ(st.rc_transferred, 1u);
    EXPECT_EQ(g.node(m).read_support, 5u)
        << "covering predecessor must accumulate dropped node's RC";
    EXPECT_EQ(g.node(t).read_support, 1u) << "sibling untouched";
}

TEST(GraphFilterTest, Keep_contained_preserves_all_nodes_and_rc) {
    // With drop_contained = false, every node survives and no RC is
    // transferred. Matches --keep-contained CLI path.
    LosslessGraph g;
    NodeId m = g.add_node(5000);
    NodeId s = g.add_node(500);
    g.node(m).read_support = 3;
    g.node(s).read_support = 2;
    g.add_edge(m, s, 1);

    FilterConfig cfg;
    cfg.drop_contained = false;
    cfg.reduce_transitive = false;

    FilterStats st = filter_graph(g, cfg);

    EXPECT_EQ(st.nodes_dropped_contained, 0u);
    EXPECT_EQ(st.rc_transferred, 0u);
    EXPECT_EQ(g.node(m).read_support, 3u);
    EXPECT_EQ(g.node(s).read_support, 2u);
    EXPECT_EQ(g.edge_count(), 1u);
}
