// Tests for branch::io::read_mosdepth_regions.
//
// Covers:
//   - Plain-text BED input (5-column and 4-column variants)
//   - Gzip-compressed input (.bed.gz) — auto-detected by suffix
//   - Skipping of comment and empty lines
//   - Windows \r\n line endings

#include "io/mosdepth_reader.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

#include <gtest/gtest.h>
#include <zlib.h>

namespace {

std::string make_tmp_path(const std::string& suffix) {
    auto base = std::filesystem::temp_directory_path();
    auto p = base / ("branch_mosdepth_test_" + std::to_string(::getpid()) + "_" +
                     std::to_string(std::rand()) + suffix);
    return p.string();
}

void write_plain(const std::string& path, const std::string& content) {
    std::ofstream ofs(path);
    ofs << content;
}

void write_gz(const std::string& path, const std::string& content) {
    gzFile gz = gzopen(path.c_str(), "wb");
    ASSERT_NE(gz, nullptr);
    gzwrite(gz, content.data(), static_cast<unsigned>(content.size()));
    gzclose(gz);
}

}  // namespace

TEST(MosdepthReader, ReadsPlain5ColBed) {
    auto path = make_tmp_path(".bed");
    write_plain(path,
                "chr1\t100\t200\tgeneA\t42.5\n"
                "chr2\t300\t400\tgeneB\t17.0\n");
    auto rows = branch::io::read_mosdepth_regions(path);
    ASSERT_EQ(rows.size(), 2u);
    EXPECT_EQ(rows[0].chrom, "chr1");
    EXPECT_EQ(rows[0].start, 100u);
    EXPECT_EQ(rows[0].end, 200u);
    EXPECT_EQ(rows[0].name, "geneA");
    EXPECT_FLOAT_EQ(rows[0].mean_coverage, 42.5f);
    EXPECT_EQ(rows[1].name, "geneB");
    EXPECT_FLOAT_EQ(rows[1].mean_coverage, 17.0f);
    std::filesystem::remove(path);
}

TEST(MosdepthReader, ReadsPlain4ColBed) {
    auto path = make_tmp_path(".bed");
    write_plain(path,
                "chr1\t100\t200\t42.5\n"
                "chr2\t300\t400\t17.0\n");
    auto rows = branch::io::read_mosdepth_regions(path);
    ASSERT_EQ(rows.size(), 2u);
    EXPECT_TRUE(rows[0].name.empty());
    EXPECT_FLOAT_EQ(rows[0].mean_coverage, 42.5f);
    std::filesystem::remove(path);
}

TEST(MosdepthReader, ReadsGzippedBed) {
    auto path = make_tmp_path(".bed.gz");
    write_gz(path,
             "chr1\t100\t200\tREF_RPL11\t31.39\n"
             "chr14\t105620506\t105626066\tIGHG4\t67.85\n"
             "chr19\t49470000\t49475000\tREF_EPOR\t38.70\n");
    auto rows = branch::io::read_mosdepth_regions(path);
    ASSERT_EQ(rows.size(), 3u);
    EXPECT_EQ(rows[0].name, "REF_RPL11");
    EXPECT_FLOAT_EQ(rows[0].mean_coverage, 31.39f);
    EXPECT_EQ(rows[1].name, "IGHG4");
    EXPECT_EQ(rows[1].start, 105620506u);
    EXPECT_EQ(rows[1].end, 105626066u);
    EXPECT_FLOAT_EQ(rows[1].mean_coverage, 67.85f);
    EXPECT_EQ(rows[2].name, "REF_EPOR");
    std::filesystem::remove(path);
}

TEST(MosdepthReader, SkipsCommentAndEmptyLines) {
    auto path = make_tmp_path(".bed");
    write_plain(path,
                "# this is a comment\n"
                "\n"
                "chr1\t100\t200\tgeneA\t42.5\n"
                "# another comment\n"
                "chr2\t300\t400\tgeneB\t17.0\n");
    auto rows = branch::io::read_mosdepth_regions(path);
    EXPECT_EQ(rows.size(), 2u);
    std::filesystem::remove(path);
}

TEST(MosdepthReader, HandlesCRLFLineEndings) {
    auto path = make_tmp_path(".bed");
    write_plain(path,
                "chr1\t100\t200\tgeneA\t42.5\r\n"
                "chr2\t300\t400\tgeneB\t17.0\r\n");
    auto rows = branch::io::read_mosdepth_regions(path);
    ASSERT_EQ(rows.size(), 2u);
    EXPECT_EQ(rows[0].name, "geneA");
    EXPECT_EQ(rows[1].name, "geneB");
    std::filesystem::remove(path);
}

TEST(MosdepthReader, ReturnsEmptyOnMissingFile) {
    auto rows = branch::io::read_mosdepth_regions("/nonexistent/path/xyz.bed");
    EXPECT_TRUE(rows.empty());
}

TEST(MosdepthReader, ReturnsEmptyOnMissingGzFile) {
    auto rows = branch::io::read_mosdepth_regions("/nonexistent/path/xyz.bed.gz");
    EXPECT_TRUE(rows.empty());
}
