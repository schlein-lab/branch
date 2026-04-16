#include <gtest/gtest.h>

#include "graph/reverse_index.hpp"

using branch::graph::NodeId;
using branch::graph::ReadId;
using branch::graph::ReverseIndex;

TEST(ReverseIndexTest, Empty_index_reports_zero) {
    ReverseIndex ri;
    EXPECT_EQ(ri.node_count(), 0u);
    EXPECT_EQ(ri.total_entries(), 0u);
    EXPECT_FALSE(ri.is_frozen());
}

TEST(ReverseIndexTest, add_grows_buckets_on_demand) {
    ReverseIndex ri;
    ri.add(5u, 100u);
    ri.add(5u, 101u);
    ri.add(2u, 42u);
    EXPECT_GE(ri.node_count(), 6u);
    EXPECT_EQ(ri.total_entries(), 3u);
}

TEST(ReverseIndexTest, freeze_collapses_to_CSR_and_preserves_entries) {
    ReverseIndex ri;
    ri.reserve_nodes(10);
    ri.add(3u, 10u);
    ri.add(3u, 11u);
    ri.add(3u, 12u);
    ri.add(7u, 99u);

    ri.freeze();

    EXPECT_TRUE(ri.is_frozen());
    EXPECT_EQ(ri.total_entries(), 4u);

    auto n3 = ri.reads_for(3u);
    ASSERT_EQ(n3.size(), 3u);
    EXPECT_EQ(n3[0], 10u);
    EXPECT_EQ(n3[1], 11u);
    EXPECT_EQ(n3[2], 12u);

    auto n7 = ri.reads_for(7u);
    ASSERT_EQ(n7.size(), 1u);
    EXPECT_EQ(n7[0], 99u);

    auto n0 = ri.reads_for(0u);
    EXPECT_EQ(n0.size(), 0u);
}

TEST(ReverseIndexTest, reads_for_returns_span_pointing_into_csr) {
    ReverseIndex ri;
    ri.reserve_nodes(3);
    ri.add(0u, 1u);
    ri.add(0u, 2u);
    ri.add(1u, 3u);
    ri.freeze();

    auto span0 = ri.reads_for(0u);
    ASSERT_EQ(span0.size(), 2u);
    // Contiguous layout -> data pointers differ by sizeof(ReadId)
    EXPECT_EQ(&span0[1] - &span0[0], 1);
}

// ============ Death Tests ============

TEST(ReverseIndexDeathTest, add_after_freeze_throws) {
    ReverseIndex ri;
    ri.add(0u, 1u);
    ri.freeze();
    EXPECT_TRUE(ri.is_frozen());
    // Adding after freeze must throw
    EXPECT_THROW(ri.add(0u, 2u), std::logic_error);
    EXPECT_THROW(ri.add(5u, 99u), std::logic_error);
}

// ============ Edge Cases ============

TEST(ReverseIndexTest, freeze_with_zero_entries) {
    ReverseIndex ri;
    // No entries added
    EXPECT_EQ(ri.total_entries(), 0u);
    EXPECT_FALSE(ri.is_frozen());
    // Freeze empty index must succeed
    ri.freeze();
    EXPECT_TRUE(ri.is_frozen());
    EXPECT_EQ(ri.total_entries(), 0u);
    // reads_for on any node returns empty span
    auto span = ri.reads_for(0u);
    EXPECT_EQ(span.size(), 0u);
}
