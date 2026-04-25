#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include "analysis/annotation_bed.hpp"

using namespace branch::analysis;

TEST(AnnotationBed, LoadBed) {
    const char* tmp = "/tmp/test_annot.bed";
    {
        std::ofstream f(tmp);
        f << "chr1\t100\t200\tentry1\n";
        f << "chr1\t300\t400\tentry2\n";
        f << "chr2\t500\t600\tentry3\n";
    }
    auto entries = load_bed(tmp);
    ASSERT_EQ(entries.size(), 3u);
    EXPECT_EQ(entries[0].chrom, "chr1");
    EXPECT_EQ(entries[0].start, 100u);
    EXPECT_EQ(entries[0].end, 200u);
    EXPECT_EQ(entries[0].name, "entry1");
    EXPECT_EQ(entries[2].chrom, "chr2");
    std::remove(tmp);
}

TEST(AnnotationBed, IntervalIndexQuery) {
    std::vector<BedEntry> entries = {
        {"chr1", "seg1", 100, 200, std::nullopt, {}, {}},
        {"chr1", "seg2", 300, 400, std::nullopt, {}, {}},
        {"chr2", "seg3", 500, 600, std::nullopt, {}, {}}
    };
    IntervalIndex idx(entries);
    
    auto hits = idx.query("chr1", 150, 250);
    ASSERT_EQ(hits.size(), 1u);
    EXPECT_EQ(hits[0]->name, "seg1");
    
    hits = idx.query("chr1", 200, 300);
    EXPECT_EQ(hits.size(), 0u);
    
    hits = idx.query("chr1", 50, 350);
    ASSERT_EQ(hits.size(), 2u);
    
    hits = idx.query("chr2", 550, 650);
    ASSERT_EQ(hits.size(), 1u);
    EXPECT_EQ(hits[0]->name, "seg3");
    
    hits = idx.query("chr3", 0, 1000);
    EXPECT_EQ(hits.size(), 0u);
}

TEST(AnnotationBed, EmptyIndex) {
    std::vector<BedEntry> entries;
    IntervalIndex idx(entries);
    EXPECT_TRUE(idx.empty());
    auto hits = idx.query("chr1", 0, 1000);
    EXPECT_EQ(hits.size(), 0u);
}

// P1.2: confidence column round-trips via write_bed_entry + load_bed.
TEST(AnnotationBed, ConfidenceColumnRoundTrip) {
    const char* tmp = "/tmp/test_annot_conf.bed";
    {
        std::ofstream f(tmp);
        BedEntry with_conf{"chr1", "bubble_1", 100, 200, 0.873f, {}, {}};
        BedEntry without_conf{"chr1", "bubble_2", 300, 400, std::nullopt, {}, {}};
        write_bed_entry(f, with_conf);
        write_bed_entry(f, without_conf);
    }
    auto entries = load_bed(tmp);
    ASSERT_EQ(entries.size(), 2u);
    ASSERT_TRUE(entries[0].confidence.has_value());
    EXPECT_NEAR(*entries[0].confidence, 0.873f, 1e-3f);
    EXPECT_FALSE(entries[1].confidence.has_value());
    std::remove(tmp);
}

// P2.4: per-alt read-support column round-trips and min_vaf / total are
// emitted for >=2 alts with non-zero sum.
TEST(AnnotationBed, AltReadSupportColumnRoundTrip) {
    const char* tmp = "/tmp/test_annot_vaf.bed";
    {
        std::ofstream f(tmp);
        BedEntry two_alts{"graph", "bubble_1", 0, 50000, 0.9f, {45u, 3u}, {}};
        BedEntry three_alts{"graph", "bubble_2", 0, 12000, 0.5f, {100u, 10u, 2u}, {}};
        BedEntry no_alts{"graph", "bubble_3", 0, 100, 0.1f, {}, {}};
        write_bed_entry(f, two_alts);
        write_bed_entry(f, three_alts);
        write_bed_entry(f, no_alts);
    }
    auto entries = load_bed(tmp);
    ASSERT_EQ(entries.size(), 3u);
    ASSERT_EQ(entries[0].alt_read_supports.size(), 2u);
    EXPECT_EQ(entries[0].alt_read_supports[0], 45u);
    EXPECT_EQ(entries[0].alt_read_supports[1], 3u);
    ASSERT_EQ(entries[1].alt_read_supports.size(), 3u);
    EXPECT_EQ(entries[1].alt_read_supports[2], 2u);
    EXPECT_TRUE(entries[2].alt_read_supports.empty());
    std::remove(tmp);
}

// P2.4: serialised output columns are exactly 8 with alts, 5 without.
TEST(AnnotationBed, AltSupportSerialisedColumnCounts) {
    std::ostringstream os_with, os_without;
    BedEntry with_alts{"graph", "b1", 0, 1000, 0.8f, {40u, 10u}, {}};
    BedEntry without_alts{"chr1", "b2", 100, 200, 0.5f, {}, {}};
    write_bed_entry(os_with, with_alts);
    write_bed_entry(os_without, without_alts);

    auto count_cols = [](const std::string& line) {
        std::size_t n = 1;
        for (char c : line) if (c == '\t') ++n;
        return n;
    };
    // Strip trailing newline.
    std::string l1 = os_with.str();   l1.pop_back();
    std::string l2 = os_without.str(); l2.pop_back();
    EXPECT_EQ(count_cols(l1), 8u);
    EXPECT_EQ(count_cols(l2), 5u);

    // min_vaf = 10/(40+10) = 0.2000, total = 50
    EXPECT_NE(l1.find("\t0.2000\t"), std::string::npos);
    EXPECT_NE(l1.rfind("\t50"),     std::string::npos);
}

// P2.4: min_vaf is "." when fewer than 2 alts.
TEST(AnnotationBed, AltSupportSingleAltHasNoVaf) {
    std::ostringstream os;
    BedEntry one_alt{"graph", "b1", 0, 100, std::nullopt, {17u}, {}};
    write_bed_entry(os, one_alt);
    // Expected tail: "\t17\t.\t17\n"
    EXPECT_NE(os.str().find("\t17\t.\t17\n"), std::string::npos);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
