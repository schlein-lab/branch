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
        {"chr1", "seg1", 100, 200, std::nullopt},
        {"chr1", "seg2", 300, 400, std::nullopt},
        {"chr2", "seg3", 500, 600, std::nullopt}
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
        BedEntry with_conf{"chr1", "bubble_1", 100, 200, 0.873f};
        BedEntry without_conf{"chr1", "bubble_2", 300, 400, std::nullopt};
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

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
