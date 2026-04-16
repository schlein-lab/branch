#include <gtest/gtest.h>

#include <array>

#include "backend/backend_vtable.hpp"
#include "graph/graph_builder.hpp"

using branch::backend::OverlapPair;
using branch::graph::build_graph;
using branch::graph::BuildResult;
using branch::graph::ReadMeta;

TEST(GraphBuilderTest, Empty_inputs_yield_empty_graph) {
    auto r = build_graph({}, {});
    EXPECT_EQ(r.graph.node_count(), 0u);
    EXPECT_EQ(r.graph.edge_count(), 0u);
}

TEST(GraphBuilderTest, One_read_one_node_no_edges) {
    std::array<ReadMeta, 1> reads{{{.read_id = 0, .length_bp = 15000}}};
    auto r = build_graph(reads, {});
    EXPECT_EQ(r.graph.node_count(), 1u);
    EXPECT_EQ(r.graph.edge_count(), 0u);
    EXPECT_EQ(r.graph.node(0).length_bp, 15000u);
    EXPECT_EQ(r.read_to_node[0], 0u);
}

TEST(GraphBuilderTest, Overlap_pair_creates_directed_edge) {
    std::array<ReadMeta, 2> reads{{
        {.read_id = 0, .length_bp = 10000},
        {.read_id = 1, .length_bp = 12000},
    }};
    std::array<OverlapPair, 1> overlaps{{
        {.read_a = 0, .read_b = 1, .offset_a = 0, .offset_b = 3000,
         .overlap_len = 7000, .diff_count = 10, .strand = 0, ._pad = 0},
    }};
    auto r = build_graph(reads, overlaps);
    EXPECT_EQ(r.graph.node_count(), 2u);
    ASSERT_EQ(r.graph.edge_count(), 1u);
    EXPECT_EQ(r.graph.edges()[0].from, 0u);
    EXPECT_EQ(r.graph.edges()[0].to, 1u);
    EXPECT_EQ(r.graph.edges()[0].read_support, 1u);
}

TEST(GraphBuilderTest, Sparse_read_ids_get_packed_node_ids) {
    // Read IDs 7 and 42 -> node IDs 0 and 1 (insertion order).
    std::array<ReadMeta, 2> reads{{
        {.read_id = 7,  .length_bp = 1000},
        {.read_id = 42, .length_bp = 2000},
    }};
    std::array<OverlapPair, 1> overlaps{{
        {.read_a = 7, .read_b = 42, .offset_a = 0, .offset_b = 0,
         .overlap_len = 500, .diff_count = 0, .strand = 0, ._pad = 0},
    }};
    auto r = build_graph(reads, overlaps);
    EXPECT_EQ(r.graph.node_count(), 2u);
    EXPECT_EQ(r.read_to_node[7], 0u);
    EXPECT_EQ(r.read_to_node[42], 1u);
    ASSERT_EQ(r.graph.edge_count(), 1u);
    EXPECT_EQ(r.graph.edges()[0].from, 0u);
    EXPECT_EQ(r.graph.edges()[0].to, 1u);
}

TEST(GraphBuilderTest, Overlap_to_unknown_read_is_silently_dropped) {
    std::array<ReadMeta, 1> reads{{{.read_id = 0, .length_bp = 1000}}};
    // Overlap references read_id 99 which wasn't in reads[].
    std::array<OverlapPair, 1> overlaps{{
        {.read_a = 0, .read_b = 99, .offset_a = 0, .offset_b = 0,
         .overlap_len = 500, .diff_count = 0, .strand = 0, ._pad = 0},
    }};
    auto r = build_graph(reads, overlaps);
    EXPECT_EQ(r.graph.node_count(), 1u);
    EXPECT_EQ(r.graph.edge_count(), 0u);
}
