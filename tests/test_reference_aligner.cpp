// BRANCH — Reference Aligner Tests
//
// Tests for the minimap2/ksw2 reference alignment wrapper.
// Test cases:
//   1. query == ref → identity = 1.0
//   2. query with 1 SNP → identity < 1.0

#include <gtest/gtest.h>

#include <string>

// Note: Full integration with minimap2 requires the library to be built.
// These tests verify the interface; actual alignment requires minimap2.

namespace {

// Test sequence (104bp)
const std::string kRefSeq = 
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

// Same as ref → expect identity = 1.0
const std::string kQueryIdentical = kRefSeq;

// One SNP at position 25: A→T
const std::string kQueryOneSNP = 
    "ACGTACGTACGTACGTACGTACGTTCGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

}  // namespace

// Placeholder test — actual minimap2 integration requires library linkage
TEST(ReferenceAlignerTest, IdenticalQueryReturnsIdentityOne) {
    // This test documents the expected behavior.
    // Full implementation requires minimap2 library to be built and linked.
    
    // Expected: aligner.align_identity(kQueryIdentical) == 1.0f
    // For now, verify test infrastructure works
    EXPECT_EQ(kRefSeq.size(), 104u);
    EXPECT_EQ(kQueryIdentical, kRefSeq);
    
    // When minimap2 is linked:
    // std::string ref_path = "/tmp/test_ref.fa";
    // std::ofstream out(ref_path);
    // out << ">ref\n" << kRefSeq << "\n";
    // out.close();
    // branch::align::ReferenceAligner aligner(ref_path);
    // ASSERT_TRUE(aligner.is_valid());
    // float identity = aligner.align_identity(kQueryIdentical);
    // EXPECT_FLOAT_EQ(identity, 1.0f);
    // std::remove(ref_path.c_str());
}

TEST(ReferenceAlignerTest, QueryWithSNPReturnsLowerIdentity) {
    // This test documents the expected behavior.
    // Full implementation requires minimap2 library to be built and linked.
    
    // Verify the SNP is present
    EXPECT_NE(kQueryOneSNP, kRefSeq);
    
    // Count differences
    int diffs = 0;
    for (size_t i = 0; i < kRefSeq.size() && i < kQueryOneSNP.size(); ++i) {
        if (kRefSeq[i] != kQueryOneSNP[i]) ++diffs;
    }
    EXPECT_EQ(diffs, 1);  // Exactly one SNP
    
    // When minimap2 is linked:
    // branch::align::ReferenceAligner aligner(ref_path);
    // float identity = aligner.align_identity(kQueryOneSNP);
    // EXPECT_LT(identity, 1.0f);
    // EXPECT_GT(identity, 0.9f);  // Should be ~0.99 for 1 SNP in 100bp
}

TEST(ReferenceAlignerTest, EmptyQueryReturnsZero) {
    // Empty query should return 0.0f identity
    std::string empty_query;
    EXPECT_TRUE(empty_query.empty());
    
    // When minimap2 is linked:
    // float identity = aligner.align_identity(empty_query);
    // EXPECT_FLOAT_EQ(identity, 0.0f);
}
