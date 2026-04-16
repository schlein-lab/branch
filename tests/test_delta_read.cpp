#include <gtest/gtest.h>

#include "graph/delta_read.hpp"

using branch::graph::PositionalDelta;
using branch::graph::ReadPath;
using branch::graph::ReadPathView;

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
