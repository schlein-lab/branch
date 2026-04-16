#include <gtest/gtest.h>

#include "graph/kmer_sketch.hpp"

using branch::graph::kMinimizerK;
using branch::graph::kMinimizerW;
using branch::graph::MinimizerHit;
using branch::graph::MinimizerHitsView;
using branch::graph::ReadSketch;

TEST(KmerSketchTest, MinimizerHit_packs_to_16_bytes) {
    EXPECT_EQ(sizeof(MinimizerHit), 16u);
}

TEST(KmerSketchTest, Parameters_are_HiFi_appropriate) {
    EXPECT_EQ(kMinimizerK, 21u);
    EXPECT_EQ(kMinimizerW, 19u);
}

TEST(KmerSketchTest, ReadSketch_holds_hits_in_order) {
    ReadSketch s{};
    s.read_id = 7;
    s.read_length = 15000;
    s.hits.push_back({.hash = 0xDEADBEEF, .read_id = 7, .pos = 100});
    s.hits.push_back({.hash = 0xCAFEBABE, .read_id = 7, .pos = 120});
    ASSERT_EQ(s.hits.size(), 2u);
    EXPECT_EQ(s.hits[0].pos, 100u);
    EXPECT_EQ(s.hits[1].pos, 120u);
}

TEST(KmerSketchTest, MinimizerHitsView_is_non_owning_span) {
    std::vector<MinimizerHit> storage = {
        {.hash = 0x1, .read_id = 0, .pos = 0},
        {.hash = 0x2, .read_id = 0, .pos = 50},
    };
    MinimizerHitsView v{.hits = storage};
    EXPECT_EQ(v.hits.size(), 2u);
    EXPECT_EQ(v.hits[1].hash, 0x2u);
}
