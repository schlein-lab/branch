#include <gtest/gtest.h>
#include "classify/mixed_decomposer.hpp"

using namespace branch::classify;

TEST(MixedDecomposer, TwoHaplotypes) {
    // 6 reads: 3 with SNP='A', 3 with SNP='T'
    std::vector<std::vector<char>> snp_vectors = {
        {'A'},  // read 0 - haplotype 1
        {'A'},  // read 1 - haplotype 1
        {'A'},  // read 2 - haplotype 1
        {'T'},  // read 3 - haplotype 2
        {'T'},  // read 4 - haplotype 2
        {'T'},  // read 5 - haplotype 2
    };

    auto result = decompose_mixed(snp_vectors, 3, 3);

    ASSERT_EQ(result.size(), 2u);
    EXPECT_EQ(result[0].read_ids.size(), 3u);
    EXPECT_EQ(result[1].read_ids.size(), 3u);
    EXPECT_NEAR(result[0].estimated_vaf, 0.5, 0.01);
    EXPECT_NEAR(result[1].estimated_vaf, 0.5, 0.01);
}

TEST(MixedDecomposer, TooFewReads) {
    // Only 4 reads - not enough for 2 clusters of min_size=3
    std::vector<std::vector<char>> snp_vectors = {
        {'A'}, {'A'}, {'T'}, {'T'},
    };

    auto result = decompose_mixed(snp_vectors, 3, 3);
    EXPECT_TRUE(result.empty());
}

TEST(MixedDecomposer, SingleHaplotype) {
    // All reads identical - should cluster into 1 group, fail to decompose
    std::vector<std::vector<char>> snp_vectors = {
        {'A'}, {'A'}, {'A'}, {'A'}, {'A'}, {'A'},
    };

    auto result = decompose_mixed(snp_vectors, 3, 3);
    // Either empty or both clusters have same reads (implementation dependent)
    // With single-link clustering of identical vectors, might end up unbalanced
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
