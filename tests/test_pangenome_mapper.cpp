// BRANCH v0.4 — Unit tests for pangenome_mapper.

#include <gtest/gtest.h>
#include "project/pangenome_mapper.hpp"

#include <cstdlib>
#include <string>
#include <vector>

using namespace branch::project;

namespace {

// Helper: check if GraphAligner is available in PATH
bool graphaligner_in_path() {
    int ret = std::system("which GraphAligner >/dev/null 2>&1");
    return ret == 0;
}

}  // namespace

// ============================================================================
// Shell-safety tests
// ============================================================================

TEST(PangenomeMapperShellSafety, RejectsFastaPathWithQuote) {
    std::vector<PangenomeRef> refs = {{"HPRC", "/data/hprc.gbz"}};
    PangenomeMapOptions opts;
    std::string err;

    auto result = map_branches_pangenome("/tmp/bad'path.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(PangenomeMapperShellSafety, RejectsToolPathWithDollar) {
    std::vector<PangenomeRef> refs = {{"HPRC", "/data/hprc.gbz"}};
    PangenomeMapOptions opts;
    opts.tool_path = "/usr/bin/$EVIL";
    std::string err;

    auto result = map_branches_pangenome("/tmp/good.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

TEST(PangenomeMapperShellSafety, RejectsRefPathWithBacktick) {
    std::vector<PangenomeRef> refs = {{"HPRC", "/data/hprc`whoami`.gbz"}};
    PangenomeMapOptions opts;
    std::string err;

    auto result = map_branches_pangenome("/tmp/good.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_FALSE(err.empty());
    EXPECT_NE(err.find("unsafe"), std::string::npos);
}

// ============================================================================
// Empty refs test
// ============================================================================

TEST(PangenomeMapperEmptyRefs, ReturnsEmptyForNoRefs) {
    std::vector<PangenomeRef> refs;  // empty
    PangenomeMapOptions opts;
    std::string err;

    auto result = map_branches_pangenome("/tmp/branches.fa", refs, opts, &err);

    EXPECT_TRUE(result.empty());
    EXPECT_TRUE(err.empty());  // No error, just empty result
}

// ============================================================================
// Integration test (requires GraphAligner in PATH)
// ============================================================================

TEST(PangenomeMapperIntegration, SkipIfNoGraphAligner) {
    if (!graphaligner_in_path()) {
        GTEST_SKIP() << "GraphAligner not found in PATH — skipping integration test";
    }

    // If GraphAligner is available but no real GBZ, we expect it to fail gracefully
    std::vector<PangenomeRef> refs = {{"test", "/nonexistent/graph.gfa"}};
    PangenomeMapOptions opts;
    std::string err;

    auto result = map_branches_pangenome("/nonexistent/query.fa", refs, opts, &err);

    // Should return empty (tool ran but found nothing / failed on missing files)
    EXPECT_TRUE(result.empty());
}
