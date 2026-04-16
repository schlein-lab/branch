#include <gtest/gtest.h>
#include "vpcr/vpcr.hpp"

using namespace branch::vpcr;

TEST(VpcrTest, ReverseComplement) {
    EXPECT_EQ(reverse_complement("ATCG"), "CGAT");
    EXPECT_EQ(reverse_complement("AAAA"), "TTTT");
    EXPECT_EQ(reverse_complement("GCGC"), "GCGC");
    EXPECT_EQ(reverse_complement(""), "");
    EXPECT_EQ(reverse_complement("a"), "T");  // lowercase handling
}

TEST(VpcrTest, ScanAmpliconsFindsMatch) {
    // Synthetic read: FWD_PRIMER + amplicon_body + REVCOMP(REV_PRIMER)
    // FWD = "ACGTACGT", REV = "TGCATGCA"
    // REVCOMP(REV) = "TGCATGCA" -> reverse_complement = "TGCATGCA"
    // Actually: reverse_complement("TGCATGCA") = "TGCATGCA" (palindrome-ish)
    // Let's use non-palindrome: REV = "AAAACCCC", REVCOMP = "GGGGTTTT"
    
    PrimerPair pp;
    pp.name = "test_primer";
    pp.fwd_seq = "ACGTACGT";
    pp.rev_seq = "AAAACCCC";  // revcomp = GGGGTTTT
    
    // Read: FWD + middle + REVCOMP(REV)
    std::string read = "ACGTACGT" "NNNNNNNNN" "GGGGTTTT";
    //                   0-7        8-16        17-24
    
    auto hits = scan_amplicons(pp, "read1", read, 0);
    
    ASSERT_EQ(hits.size(), 1u);
    EXPECT_EQ(hits[0].read_id, "read1");
    EXPECT_EQ(hits[0].start, 0u);
    EXPECT_EQ(hits[0].end, read.size());
    EXPECT_EQ(hits[0].amplicon_seq, read);
}

TEST(VpcrTest, ScanAmpliconsNoMatchIfMissingPrimer) {
    PrimerPair pp;
    pp.name = "test_primer";
    pp.fwd_seq = "ACGTACGT";
    pp.rev_seq = "AAAACCCC";
    
    // Read without rev primer
    std::string read = "ACGTACGTNNNNNNNNN";
    
    auto hits = scan_amplicons(pp, "read2", read, 0);
    EXPECT_TRUE(hits.empty());
}
