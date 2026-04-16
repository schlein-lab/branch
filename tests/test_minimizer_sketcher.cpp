#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include "graph/minimizer_sketcher.hpp"

using branch::graph::canonical_hash;
using branch::graph::encode_base;
using branch::graph::kMinimizerK;
using branch::graph::kMinimizerW;
using branch::graph::MinimizerHit;
using branch::graph::rc_kmer;
using branch::graph::sketch_read;
using branch::graph::splitmix64;

TEST(MinimizerSketcherTest, encode_base_covers_ACGTN_case_insensitive) {
    EXPECT_EQ(encode_base('A'), 0);
    EXPECT_EQ(encode_base('a'), 0);
    EXPECT_EQ(encode_base('C'), 1);
    EXPECT_EQ(encode_base('c'), 1);
    EXPECT_EQ(encode_base('G'), 2);
    EXPECT_EQ(encode_base('T'), 3);
    EXPECT_EQ(encode_base('N'), 4);
    EXPECT_EQ(encode_base('X'), 4);
    EXPECT_EQ(encode_base('-'), 4);
}

TEST(MinimizerSketcherTest, splitmix64_is_deterministic_and_bijective_on_small_set) {
    EXPECT_EQ(splitmix64(0x1234ULL), splitmix64(0x1234ULL));
    EXPECT_NE(splitmix64(0x1234ULL), splitmix64(0x1235ULL));
}

TEST(MinimizerSketcherTest, rc_kmer_inverts_correctly) {
    // k=4: AACG = 0b 00 00 01 10 = 0x06
    // rc  : CGTT = 0b 01 10 11 11 = 0x6F
    std::uint64_t fwd = 0x06ULL;
    std::uint64_t rc = rc_kmer(fwd, 4);
    // rc(rc(x)) == x
    EXPECT_EQ(rc_kmer(rc, 4), fwd);
}

TEST(MinimizerSketcherTest, canonical_hash_is_strand_symmetric) {
    // Hash of a kmer and its reverse complement must be identical.
    std::uint64_t fwd = 0x12345ULL;
    auto rc = rc_kmer(fwd, 21);
    EXPECT_EQ(canonical_hash(fwd, 21), canonical_hash(rc, 21));
}

TEST(MinimizerSketcherTest, sketch_short_read_emits_nothing) {
    std::vector<MinimizerHit> out;
    std::string short_read(10, 'A');  // below k=21
    sketch_read(short_read, 1, out);
    EXPECT_EQ(out.size(), 0u);
}

TEST(MinimizerSketcherTest, sketch_random_read_emits_reasonable_number_of_minimizers) {
    // Produce ~1 minimizer per w bases on non-pathological sequence.
    std::string s;
    s.reserve(1000);
    std::uint64_t state = 0xdeadbeefULL;
    for (int i = 0; i < 1000; ++i) {
        state = splitmix64(state);
        s.push_back("ACGT"[state & 3ULL]);
    }
    std::vector<MinimizerHit> out;
    sketch_read(s, 7, out);
    // Expect roughly len / w minimizers; allow a wide range for consecutive-duplicate elision.
    const std::size_t expected_min = 1000 / (kMinimizerW * 4);
    const std::size_t expected_max = 1000 / 1;
    EXPECT_GE(out.size(), expected_min);
    EXPECT_LE(out.size(), expected_max);
    // All hits must carry the provided read_id.
    for (const auto& h : out) {
        EXPECT_EQ(h.read_id, 7u);
    }
}

TEST(MinimizerSketcherTest, sketch_resets_on_ambiguous_base) {
    // Two valid windows separated by Ns; sketch continues after reset.
    std::string s;
    for (int i = 0; i < 100; ++i) s.push_back("ACGT"[i % 4]);
    s += std::string(25, 'N');
    for (int i = 0; i < 100; ++i) s.push_back("ACGT"[i % 4]);

    std::vector<MinimizerHit> out;
    sketch_read(s, 99, out);
    EXPECT_GT(out.size(), 0u);
    // No minimizer should point inside the N-run.
    for (const auto& h : out) {
        EXPECT_TRUE(h.pos < 100 || h.pos >= 125);
    }
}

TEST(MinimizerSketcherTest, positions_are_non_decreasing_within_a_read) {
    std::string s(500, 'A');
    for (std::size_t i = 0; i < s.size(); i += 4) s[i] = "ACGT"[i / 4 % 4];
    std::vector<MinimizerHit> out;
    sketch_read(s, 0, out);
    for (std::size_t i = 1; i < out.size(); ++i) {
        EXPECT_GE(out[i].pos, out[i - 1].pos);
    }
}
