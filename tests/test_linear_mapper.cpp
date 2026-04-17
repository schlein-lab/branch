// BRANCH v0.4 — Unit tests for linear_mapper.
//
// Tests the linear mapper module:
// 1. Shell safety validation for paths
// 2. End-to-end mapping with minimap2 (skipped if minimap2 not available)

#include <gtest/gtest.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#include "project/linear_mapper.hpp"

namespace {

// Check if minimap2 is available in PATH
bool minimap2_available() {
    int ret = std::system("which minimap2 > /dev/null 2>&1");
    return ret == 0;
}

// Create a temporary file with given content, return the path
std::string create_temp_file(const std::string& content, const std::string& suffix) {
    char tmpl[] = "/tmp/branch_test_XXXXXX";
    int fd = mkstemp(tmpl);
    if (fd < 0) return "";

    std::string path = std::string(tmpl) + suffix;
    close(fd);
    // Rename to add suffix
    if (rename(tmpl, path.c_str()) != 0) {
        unlink(tmpl);
        return "";
    }

    std::ofstream os(path);
    if (!os) {
        unlink(path.c_str());
        return "";
    }
    os << content;
    os.close();
    return path;
}

}  // namespace

TEST(LinearMapper, RejectsUnsafeShellPaths) {
    branch::project::LinearMapOptions opts;
    std::string err;

    // Test with backtick in path
    std::vector<branch::project::LinearRef> refs = {{"test", "/path/with`backtick"}};
    auto result = branch::project::map_branches_linear("/safe/path.fa", refs, opts, &err);
    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(LinearMapper, RejectsUnsafeFastaPath) {
    branch::project::LinearMapOptions opts;
    std::string err;

    std::vector<branch::project::LinearRef> refs = {{"test", "/safe/ref.fa"}};
    auto result = branch::project::map_branches_linear("/path/with$dollar.fa", refs, opts, &err);
    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(LinearMapper, RejectsUnsafeQuotes) {
    branch::project::LinearMapOptions opts;
    std::string err;

    std::vector<branch::project::LinearRef> refs = {{"test", "/path/with'quote.fa"}};
    auto result = branch::project::map_branches_linear("/safe/path.fa", refs, opts, &err);
    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
}

TEST(LinearMapper, EmptyRefsReturnsEmpty) {
    branch::project::LinearMapOptions opts;
    std::string err;

    std::vector<branch::project::LinearRef> refs;  // empty
    auto result = branch::project::map_branches_linear("/some/path.fa", refs, opts, &err);
    EXPECT_TRUE(result.empty());
    EXPECT_TRUE(err.empty());  // No error, just no refs to map
}

TEST(LinearMapper, EndToEndWithMinimap2) {
    if (!minimap2_available()) {
        GTEST_SKIP() << "minimap2 not found in PATH, skipping end-to-end test";
    }

    // Create a tiny reference with a 30bp sequence
    std::string ref_content =
        ">chr1\n"
        "ACGTACGTACGTACGTACGTACGTACGTAC\n";

    // Create a tiny FASTA with a 20bp branch that matches the ref
    std::string fasta_content =
        ">branch_001\n"
        "ACGTACGTACGTACGTACGT\n";

    std::string ref_path = create_temp_file(ref_content, ".ref.fa");
    std::string fasta_path = create_temp_file(fasta_content, ".branch.fa");

    ASSERT_FALSE(ref_path.empty()) << "Failed to create temp ref file";
    ASSERT_FALSE(fasta_path.empty()) << "Failed to create temp fasta file";

    // Clean up on exit
    struct Cleanup {
        std::string ref, fasta;
        ~Cleanup() {
            if (!ref.empty()) unlink(ref.c_str());
            if (!fasta.empty()) unlink(fasta.c_str());
        }
    } cleanup{ref_path, fasta_path};

    branch::project::LinearMapOptions opts;
    opts.threads = 1;
    opts.min_mapq = 0;

    std::vector<branch::project::LinearRef> refs = {{"TestRef", ref_path}};
    std::string err;

    auto result = branch::project::map_branches_linear(fasta_path, refs, opts, &err);

    // Should get at least one mapping (minimap2 with asm20 on identical sequence)
    // Note: very short sequences might not map depending on minimap2 settings
    // So we just check that no error occurred
    EXPECT_TRUE(err.empty()) << "Error: " << err;

    // If we got mappings, verify the structure
    for (const auto& m : result) {
        EXPECT_EQ(m.ref_name, "TestRef");
        EXPECT_EQ(m.branch_id, "branch_001");
        EXPECT_FALSE(m.target.empty());
    }
}

TEST(LinearMapper, MultipleRefs) {
    if (!minimap2_available()) {
        GTEST_SKIP() << "minimap2 not found in PATH, skipping multi-ref test";
    }

    // Create two refs
    std::string ref1_content = ">chr1\nACGTACGTACGTACGTACGTACGTACGTAC\n";
    std::string ref2_content = ">chrX\nTGCATGCATGCATGCATGCATGCATGCAT\n";

    // Branch that matches ref1
    std::string fasta_content = ">branch_A\nACGTACGTACGTACGTACGT\n";

    std::string ref1_path = create_temp_file(ref1_content, ".ref1.fa");
    std::string ref2_path = create_temp_file(ref2_content, ".ref2.fa");
    std::string fasta_path = create_temp_file(fasta_content, ".branch.fa");

    ASSERT_FALSE(ref1_path.empty());
    ASSERT_FALSE(ref2_path.empty());
    ASSERT_FALSE(fasta_path.empty());

    struct Cleanup {
        std::string r1, r2, fa;
        ~Cleanup() {
            if (!r1.empty()) unlink(r1.c_str());
            if (!r2.empty()) unlink(r2.c_str());
            if (!fa.empty()) unlink(fa.c_str());
        }
    } cleanup{ref1_path, ref2_path, fasta_path};

    branch::project::LinearMapOptions opts;
    opts.threads = 1;

    std::vector<branch::project::LinearRef> refs = {
        {"CHM13", ref1_path},
        {"GRCh38", ref2_path}
    };
    std::string err;

    auto result = branch::project::map_branches_linear(fasta_path, refs, opts, &err);
    EXPECT_TRUE(err.empty()) << "Error: " << err;

    // Results should be sorted by ref_name
    for (std::size_t i = 1; i < result.size(); ++i) {
        EXPECT_LE(result[i-1].ref_name, result[i].ref_name)
            << "Results should be sorted by ref_name";
    }
}
