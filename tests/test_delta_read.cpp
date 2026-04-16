#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "graph/delta_read.hpp"
#include "graph/graph_builder.hpp"

using branch::backend::OverlapPair;
using branch::graph::build_graph;
using branch::graph::NodeId;
using branch::graph::PositionalDelta;
using branch::graph::ReadMeta;
using branch::graph::ReadPath;
using branch::graph::ReadPathView;
using branch::graph::populate_read_paths;

TEST(DeltaReadTest, PositionalDeltaPacks_to_8_bytes) {
    EXPECT_EQ(sizeof(PositionalDelta), 8u);
}

TEST(DeltaReadTest, EmptyReadPath_has_no_nodes_or_deltas) {
    ReadPath rp{};
    EXPECT_TRUE(rp.path.empty());
    EXPECT_TRUE(rp.deltas.empty());
    EXPECT_EQ(rp.read_length, 0u);
}

TEST(DeltaReadTest, ReadPath_stores_node_traversal) {
    ReadPath rp{
        .read_id = 42,
        .read_length = 15000,
        .path = {1, 2, 3, 4},
        .deltas = {},
    };
    EXPECT_EQ(rp.read_id, 42u);
    EXPECT_EQ(rp.path.size(), 4u);
    EXPECT_EQ(rp.path.back(), 4u);
}

TEST(DeltaReadTest, ReadPathView_borrows_without_copy) {
    ReadPath rp{
        .read_id = 1,
        .read_length = 100,
        .path = {10, 20},
        .deltas = {{.offset = 5, .length = 1, .alt_base = 0, .alt_bases_index = 0}},
    };
    ReadPathView v{
        .read_id = rp.read_id,
        .read_length = rp.read_length,
        .path = rp.path,
        .deltas = rp.deltas,
    };
    EXPECT_EQ(v.path.size(), 2u);
    EXPECT_EQ(v.deltas.size(), 1u);
    EXPECT_EQ(v.deltas[0].offset, 5u);
}

TEST(PopulateReadPathsTest, Populate_read_paths_simple_case) {
    // Three reads, chain R0-R1-R2 (two overlaps).
    std::array<ReadMeta, 3> reads{{
        {.read_id = 0, .length_bp = 10000},
        {.read_id = 1, .length_bp = 10000},
        {.read_id = 2, .length_bp = 10000},
    }};
    std::array<OverlapPair, 2> overlaps{{
        {.read_a = 0, .read_b = 1, .offset_a = 5000, .offset_b = 0,
         .overlap_len = 5000, .diff_count = 0, .strand = 0, ._pad = 0},
        {.read_a = 1, .read_b = 2, .offset_a = 5000, .offset_b = 0,
         .overlap_len = 5000, .diff_count = 0, .strand = 0, ._pad = 0},
    }};
    auto built = build_graph(reads, overlaps);

    std::vector<std::uint32_t> lengths{10000, 10000, 10000};
    std::vector<ReadPath> paths;
    populate_read_paths(3, built.read_to_node, lengths, overlaps, paths);

    ASSERT_EQ(paths.size(), 3u);
    // Every read must have at least one node (its own), deltas empty in v0.2.
    for (std::size_t i = 0; i < paths.size(); ++i) {
        EXPECT_EQ(paths[i].read_id, static_cast<branch::graph::ReadId>(i));
        EXPECT_GE(paths[i].path.size(), 1u)
            << "read " << i << " should map to at least its own node";
        EXPECT_TRUE(paths[i].deltas.empty()) << "deltas remain empty until v0.3";
    }
    // R0 overlaps R1, R1 overlaps R0 and R2, R2 overlaps R1.
    EXPECT_EQ(paths[0].path.size(), 2u);  // self + R1
    EXPECT_EQ(paths[1].path.size(), 3u);  // self + R0 + R2
    EXPECT_EQ(paths[2].path.size(), 2u);  // self + R1
}

TEST(PopulateReadPathsTest, Populate_read_paths_empty_overlap) {
    // Two reads, no overlaps at all.
    std::array<ReadMeta, 2> reads{{
        {.read_id = 0, .length_bp = 1000},
        {.read_id = 1, .length_bp = 1000},
    }};
    auto built = build_graph(reads, {});

    std::vector<std::uint32_t> lengths{1000, 1000};
    std::vector<ReadPath> paths;
    populate_read_paths(2, built.read_to_node, lengths, {}, paths);

    ASSERT_EQ(paths.size(), 2u);
    // Each read still carries its own self-node; no partners, no crash.
    for (const auto& p : paths) {
        EXPECT_EQ(p.path.size(), 1u);
        EXPECT_TRUE(p.deltas.empty());
    }
}

TEST(PopulateReadPathsTest, Populate_read_paths_path_nodes_are_ordered) {
    // Read 1 overlaps with R0 at offset 1000 and with R2 at offset 6000.
    // The resulting ReadPath for R1 must list nodes in overlap-order:
    //   [self_R1, node_of_R0, node_of_R2]
    // since the self-node is emitted first (offset 0) and partners are
    // appended sorted by offset along R1.
    std::array<ReadMeta, 3> reads{{
        {.read_id = 0, .length_bp = 8000},
        {.read_id = 1, .length_bp = 15000},
        {.read_id = 2, .length_bp = 8000},
    }};
    std::array<OverlapPair, 2> overlaps{{
        // R0 overlaps R1 at R1-offset=1000
        {.read_a = 0, .read_b = 1, .offset_a = 0, .offset_b = 1000,
         .overlap_len = 4000, .diff_count = 0, .strand = 0, ._pad = 0},
        // R1 overlaps R2 at R1-offset=6000
        {.read_a = 1, .read_b = 2, .offset_a = 6000, .offset_b = 0,
         .overlap_len = 4000, .diff_count = 0, .strand = 0, ._pad = 0},
    }};
    auto built = build_graph(reads, overlaps);

    std::vector<std::uint32_t> lengths{8000, 15000, 8000};
    std::vector<ReadPath> paths;
    populate_read_paths(3, built.read_to_node, lengths, overlaps, paths);

    const NodeId n0 = built.read_to_node[0];
    const NodeId n1 = built.read_to_node[1];
    const NodeId n2 = built.read_to_node[2];

    ASSERT_EQ(paths.size(), 3u);
    // R1: self first, then R0 (offset 1000), then R2 (offset 6000).
    ASSERT_EQ(paths[1].path.size(), 3u);
    EXPECT_EQ(paths[1].path[0], n1);
    EXPECT_EQ(paths[1].path[1], n0);
    EXPECT_EQ(paths[1].path[2], n2);
    // Sanity: the ordered partner subsequence must be monotone by offset.
    // Since only one partner each appears after self for R0 and R2, trivially ordered.
    EXPECT_EQ(paths[0].path.front(), n0);
    EXPECT_EQ(paths[2].path.front(), n2);
    // read_length carried through.
    EXPECT_EQ(paths[1].read_length, 15000u);
}
