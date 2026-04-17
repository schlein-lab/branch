// Unit tests for write_branch_report_json — verifies JSON structure,
// null fields for unknown VAF/coverage, and per-branch mapping arrays.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "project/project_report.hpp"

namespace {

std::string read_file(const std::string& path) {
    std::ifstream f(path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(ProjectReport, EmptyBranchListIsValid) {
    auto tmp = std::filesystem::temp_directory_path() /
        "project_report_empty.json";
    std::string err;
    std::vector<branch::project::BranchEntry> empty;

    ASSERT_TRUE(branch::project::write_branch_report_json(
        tmp.string(), "sample1", empty, &err));
    EXPECT_TRUE(err.empty());

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find("\"sample\": \"sample1\""), std::string::npos);
    EXPECT_NE(content.find("\"branch_count\": 0"), std::string::npos);
    EXPECT_NE(content.find("\"unannotated_count\": 0"), std::string::npos);
    EXPECT_NE(content.find("\"branches\": []"), std::string::npos);
}

TEST(ProjectReport, UnknownVafAndCoverageWrittenAsNull) {
    branch::project::BranchEntry e;
    e.branch_id = "branch_42";
    e.length_bp = 1500;
    e.vaf = -1.0;
    e.coverage = -1.0;
    e.unannotated = true;

    auto tmp = std::filesystem::temp_directory_path() /
        "project_report_null.json";
    std::string err;
    std::vector<branch::project::BranchEntry> b{e};
    ASSERT_TRUE(branch::project::write_branch_report_json(
        tmp.string(), "sample_null", b, &err));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find("\"vaf\": null"), std::string::npos);
    EXPECT_NE(content.find("\"coverage\": null"), std::string::npos);
    EXPECT_NE(content.find("\"unannotated\": true"), std::string::npos);
    EXPECT_NE(content.find("\"unannotated_count\": 1"), std::string::npos);
    EXPECT_NE(content.find("\"branch_id\": \"branch_42\""), std::string::npos);
}

TEST(ProjectReport, KnownVafSerialisedAsNumber) {
    branch::project::BranchEntry e;
    e.branch_id = "b";
    e.length_bp = 800;
    e.vaf = 0.25;
    e.coverage = 35.5;
    e.unannotated = false;

    auto tmp = std::filesystem::temp_directory_path() /
        "project_report_values.json";
    std::string err;
    std::vector<branch::project::BranchEntry> b{e};
    ASSERT_TRUE(branch::project::write_branch_report_json(
        tmp.string(), "s", b, &err));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_EQ(content.find("\"vaf\": null"), std::string::npos);
    EXPECT_EQ(content.find("\"coverage\": null"), std::string::npos);
    EXPECT_NE(content.find("\"unannotated\": false"), std::string::npos);
}

TEST(ProjectReport, LinearMappingsArraySerialised) {
    branch::project::BranchEntry e;
    e.branch_id = "b1";
    e.length_bp = 1200;
    e.linear_mappings = {
        {"b1", "CHM13", "chr1", 1000, 2000, 60, '+', 1200, 0, 1200},
        {"b1", "GRCh38", "chr1", 1100, 2100, 55, '-', 1200, 0, 1200},
    };

    auto tmp = std::filesystem::temp_directory_path() /
        "project_report_mappings.json";
    std::string err;
    std::vector<branch::project::BranchEntry> b{e};
    ASSERT_TRUE(branch::project::write_branch_report_json(
        tmp.string(), "s", b, &err));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find("\"ref\": \"CHM13\""), std::string::npos);
    EXPECT_NE(content.find("\"ref\": \"GRCh38\""), std::string::npos);
    EXPECT_NE(content.find("\"mapq\": 60"), std::string::npos);
    EXPECT_NE(content.find("\"strand\": \"-\""), std::string::npos);
}

TEST(ProjectReport, SpecialCharactersInSampleNameEscaped) {
    auto tmp = std::filesystem::temp_directory_path() /
        "project_report_escape.json";
    std::string err;
    std::vector<branch::project::BranchEntry> b;

    ASSERT_TRUE(branch::project::write_branch_report_json(
        tmp.string(), "quoted\"name", b, &err));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find("\"sample\": \"quoted\\\"name\""),
              std::string::npos);
}
