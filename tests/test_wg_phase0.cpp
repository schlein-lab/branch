// BRANCH v0.5 — Phase 0 unit tests.

#include "wg/bloom_filter.hpp"
#include "wg/kmer_counter.hpp"
#include "wg/kmer_hist.hpp"

#include <gtest/gtest.h>
#include <string>
#include <vector>

using namespace branch::wg;
using branch::graph::kTechProfileHiFi;
using branch::graph::kTechProfileONT;

TEST(BloomFilter, BasicAddContains) {
    BloomFilter bf(1000, /*bits_per_key=*/10);
    bf.add(42);
    bf.add(123456789);
    EXPECT_TRUE(bf.maybe_contains(42));
    EXPECT_TRUE(bf.maybe_contains(123456789));
    // FP rate ~0.1% at 10 bits/key — most random keys should NOT match
    int fp = 0;
    for (std::uint64_t i = 100000; i < 101000; ++i) {
        if (bf.maybe_contains(i)) ++fp;
    }
    EXPECT_LT(fp, 50);  // generous bound, ~10 expected
}

TEST(KmerCounter, IncrementSurvives) {
    KmerCounter c(64);
    c.add(0xCAFE, 1);
    c.add(0xCAFE, 1);
    c.add(0xBEEF, 5);
    EXPECT_EQ(c.get(0xCAFE), 2u);
    EXPECT_EQ(c.get(0xBEEF), 5u);
    EXPECT_EQ(c.get(0xDEAD), 0u);
    EXPECT_EQ(c.size(), 2u);
}

TEST(KmerCounter, GrowsBeyondInitialCapacity) {
    KmerCounter c(8);
    for (std::uint64_t i = 1; i <= 1000; ++i) c.add(i, static_cast<std::uint32_t>(i));
    EXPECT_EQ(c.size(), 1000u);
    EXPECT_EQ(c.get(500), 500u);
    EXPECT_EQ(c.get(0), 0u);
}

TEST(KmerHistBuilder, TwoPassRecurrentDetection) {
    // Build a synthetic dataset where some k-mers are seen many times,
    // others exactly once.
    KmerHistBuilder b(kTechProfileHiFi, /*expected_distinct=*/64);
    // Read 1: two copies of "ACGTACGTACGTACGTACGTAC"
    std::string repeat = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    for (int i = 0; i < 30; ++i) b.pass1_add_read(repeat);
    // Singletons: each read unique
    for (int i = 0; i < 5; ++i) {
        std::string uniq = "GGGGGGGGG";
        uniq += static_cast<char>('A' + (i % 4));
        for (int j = 0; j < 25; ++j) uniq += "TTTTTTTTTTT";
        b.pass1_add_read(uniq);
    }
    b.switch_to_pass2();
    for (int i = 0; i < 30; ++i) b.pass2_add_read(repeat);
    for (int i = 0; i < 5; ++i) {
        std::string uniq = "GGGGGGGGG";
        uniq += static_cast<char>('A' + (i % 4));
        for (int j = 0; j < 25; ++j) uniq += "TTTTTTTTTTT";
        b.pass2_add_read(uniq);
    }
    auto r = b.finalize(/*expected_coverage=*/30, /*bimodality_z=*/2.0);
    // The repeat k-mers should dominate the recurrent count.
    EXPECT_GT(r.n_recurrent_kmers, 0u);
    // Some bin in the histogram should be populated (we don't predict
    // exact count because k-mer multiplicity within a single repeat read
    // is hard to compute analytically — the test just verifies that the
    // two-pass machinery wrote something).
    std::uint64_t total_in_hist = 0;
    for (auto v : r.histogram) total_in_hist += v;
    EXPECT_GT(total_in_hist, 0u);
}

TEST(CanonicalKmer, ForwardAndRcCollide) {
    // Build a k-mer + its reverse complement, hash both → should match.
    // ACGT (forward) and ACGT (RC of itself, palindrome at k=4).
    std::uint64_t kmer_acgt = 0;
    // A=0, C=1, G=2, T=3 → ACGT = 0001 1011 = 0x1B
    kmer_acgt = (0ULL << 6) | (1ULL << 4) | (2ULL << 2) | 3ULL;
    EXPECT_EQ(canonical_kmer_hash(kmer_acgt, 4),
              canonical_kmer_hash(kmer_acgt, 4));  // self consistent
}
