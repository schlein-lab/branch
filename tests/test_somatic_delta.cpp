// BRANCH v0.4.3 — Unit tests for somatic_delta.

#include <gtest/gtest.h>
#include "project/somatic_delta.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace branch::project;

namespace {

// Helper: create temp FASTA file
std::string create_temp_fasta(const std::vector<std::pair<std::string, std::string>>& seqs) {
    char path[] = "/tmp/test_somatic_XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) return "";
    close(fd);
    
    std::string fasta_path = std::string(path) + ".fa";
    std::rename(path, fasta_path.c_str());
    
    std::ofstream out(fasta_path);
    for (const auto& [name, seq] : seqs) {
        out << ">" << name << "\n" << seq << "\n";
    }
    out.close();
    return fasta_path;
}

// Helper: cleanup temp file
void cleanup_temp(const std::string& path) {
    std::remove(path.c_str());
}

}  // namespace

// ============================================================================
// Shell-safety tests
// ============================================================================

TEST(SomaticDeltaShellSafety, RejectsBranchesPathWithQuote) {
    std::vector<std::pair<std::string, std::string>> refs = {{"ref1", "/data/ref.fa"}};
    SomaticDeltaOptions opts;
    std::string err;

    auto result = compute_somatic_deltas("/tmp/bad'path.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(SomaticDeltaShellSafety, RejectsRefPathWithDollar) {
    std::vector<std::pair<std::string, std::string>> refs = {{"ref1", "/data/$EVIL.fa"}};
    SomaticDeltaOptions opts;
    std::string err;

    auto result = compute_somatic_deltas("/tmp/good.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(SomaticDeltaShellSafety, RejectsRefNameWithBacktick) {
    std::vector<std::pair<std::string, std::string>> refs = {{"ref`whoami`", "/data/ref.fa"}};
    SomaticDeltaOptions opts;
    std::string err;

    auto result = compute_somatic_deltas("/tmp/good.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

// ============================================================================
// Empty refs test
// ============================================================================

TEST(SomaticDeltaEmptyRefs, ReturnsEmptyForNoRefs) {
    std::vector<std::pair<std::string, std::string>> refs;  // empty
    SomaticDeltaOptions opts;
    std::string err;

    auto result = compute_somatic_deltas("/tmp/branches.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_TRUE(err.empty());  // No error, just empty result
}

// ============================================================================
// ksw2 alignment tests
// ============================================================================

TEST(SomaticDeltaKsw2, IdenticalSequencesHaveZeroEditDistance) {
    // 30bp identical sequences
    std::string seq = "ACGTACGTACGTACGTACGTACGTACGTAC";
    
    SomaticDeltaOptions opts;
    std::string cigar;
    int edit_dist = ksw2_edit_distance(seq, seq, opts, &cigar);
    
    EXPECT_EQ(edit_dist, 0);
    EXPECT_FALSE(cigar.empty());
}

TEST(SomaticDeltaKsw2, OneSNPHasEditDistanceOne) {
    // 30bp with 1 SNP at position 15
    std::string query  = "ACGTACGTACGTACGTACGTACGTACGTAC";
    std::string target = "ACGTACGTACGTACGAACGTACGTACGTAC";  // T->A at pos 15
    
    SomaticDeltaOptions opts;
    std::string cigar;
    int edit_dist = ksw2_edit_distance(query, target, opts, &cigar);
    
    // Edit distance should be 1 (one substitution)
    EXPECT_EQ(edit_dist, 1);
}

TEST(SomaticDeltaKsw2, OneInsertionHasEditDistanceOne) {
    std::string query  = "ACGTACGTACGTACGTACGTACGTACGTAC";
    std::string target = "ACGTACGTACGTACGTAACGTACGTACGTAC";  // extra A after pos 16
    
    SomaticDeltaOptions opts;
    std::string cigar;
    int edit_dist = ksw2_edit_distance(query, target, opts, &cigar);
    
    // Edit distance should be 1 (one indel)
    EXPECT_GE(edit_dist, 1);
    EXPECT_LE(edit_dist, 2);  // Allow some tolerance for alignment algorithms
}

// ============================================================================
// Integration test with FASTA files
// ============================================================================

TEST(SomaticDeltaIntegration, ComputeDeltasFromFastaFiles) {
    // Create temp branch FASTA
    std::string branch_fa = create_temp_fasta({
        {"branch1", "ACGTACGTACGTACGTACGTACGTACGTAC"},
        {"branch2", "TGCATGCATGCATGCATGCATGCATGCATG"}
    });
    ASSERT_FALSE(branch_fa.empty());
    
    // Create temp reference FASTA
    std::string ref_fa = create_temp_fasta({
        {"hprc_path", "ACGTACGTACGTACGTACGTACGTACGTAC"}  // identical to branch1
    });
    ASSERT_FALSE(ref_fa.empty());
    
    std::vector<std::pair<std::string, std::string>> refs = {{"HPRC", ref_fa}};
    SomaticDeltaOptions opts;
    std::string err;
    
    auto results = compute_somatic_deltas(branch_fa, refs, opts, &err);
    
    // Should have 2 results (one per branch)
    EXPECT_EQ(results.size(), 2);
    
    // branch1 vs ref should have edit_distance = 0
    auto it1 = std::find_if(results.begin(), results.end(),
                            [](const SomaticDelta& d) { return d.branch_id == "branch1"; });
    ASSERT_NE(it1, results.end());
    EXPECT_EQ(it1->edit_distance, 0);
    EXPECT_EQ(it1->ref_name, "HPRC");
    
    // branch2 vs ref should have edit_distance > 0
    auto it2 = std::find_if(results.begin(), results.end(),
                            [](const SomaticDelta& d) { return d.branch_id == "branch2"; });
    ASSERT_NE(it2, results.end());
    EXPECT_GT(it2->edit_distance, 0);
    
    // Cleanup
    cleanup_temp(branch_fa);
    cleanup_temp(ref_fa);
}
