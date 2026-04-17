// BRANCH v0.4 — Smoke tests for `branch project` subcommand.
//
// Tests:
// 1. `branch project --help` exits 0 with usage
// 2. `branch project` (no args) exits 2 with usage
// 3. Missing required args exits 2
// 4. Nonexistent files with all required args exits with error (not 0)
// 5. E2E: real FASTA + real (toy) reference -> JSON report written
//    with expected structure.

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#ifndef BRANCH_BINARY
#error "BRANCH_BINARY must be defined (path to branch executable)"
#endif

namespace {

int run_branch(const std::string& args) {
    std::string cmd = std::string(BRANCH_BINARY) + " " + args + " > /dev/null 2>&1";
    return WEXITSTATUS(std::system(cmd.c_str()));
}

std::string run_branch_capture(const std::string& args) {
    std::string cmd = std::string(BRANCH_BINARY) + " " + args + " 2>&1";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "";

    std::string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        output += buffer;
    }
    pclose(pipe);
    return output;
}

}  // namespace

TEST(ProjectSmoke, HelpExitsZero) {
    EXPECT_EQ(0, run_branch("project --help"));
}

TEST(ProjectSmoke, NoArgsExitsTwo) {
    EXPECT_EQ(2, run_branch("project"));
}

TEST(ProjectSmoke, MissingFastaExitsTwo) {
    // Has --ref-linear and --out-prefix but missing --fasta
    EXPECT_EQ(2, run_branch("project --ref-linear CHM13=/some/path --out-prefix /tmp/test"));
}

TEST(ProjectSmoke, MissingRefLinearExitsTwo) {
    // Has --fasta and --out-prefix but missing --ref-linear
    EXPECT_EQ(2, run_branch("project --fasta /some/path.fa --out-prefix /tmp/test"));
}

TEST(ProjectSmoke, MissingOutPrefixExitsTwo) {
    // Has --fasta and --ref-linear but missing --out-prefix
    EXPECT_EQ(2, run_branch("project --fasta /some/path.fa --ref-linear CHM13=/some/ref"));
}

TEST(ProjectSmoke, NonexistentFilesExitNonzero) {
    // All required args present but files don't exist — should fail (exit 3)
    int code = run_branch("project --fasta /nonexistent/path.fa "
                          "--ref-linear CHM13=/nonexistent/ref.fa "
                          "--out-prefix /tmp/branch_test_output");
    EXPECT_NE(0, code) << "Should fail when FASTA file doesn't exist";
}

TEST(ProjectSmoke, HelpWithDashH) {
    EXPECT_EQ(0, run_branch("project -h"));
}

TEST(ProjectSmoke, HelpShowsRequiredArgs) {
    std::string output = run_branch_capture("project --help");
    EXPECT_NE(output.find("--fasta"), std::string::npos)
        << "Help should mention --fasta";
    EXPECT_NE(output.find("--ref-linear"), std::string::npos)
        << "Help should mention --ref-linear";
    EXPECT_NE(output.find("--out-prefix"), std::string::npos)
        << "Help should mention --out-prefix";
}

TEST(ProjectSmoke, HelpShowsVersion) {
    std::string output = run_branch_capture("project --help");
    EXPECT_NE(output.find("0.4"), std::string::npos)
        << "Help should mention version";
}

TEST(ProjectSmoke, E2E_RealFasta_ProducesJsonReport) {
    // Skip gracefully if minimap2 is not in PATH — the linear_mapper
    // shells out to it and returns an error otherwise.
    if (std::system("which minimap2 > /dev/null 2>&1") != 0) {
        GTEST_SKIP() << "minimap2 not in PATH; skipping real-FASTA e2e";
    }

    namespace fs = std::filesystem;
    auto tmp = fs::temp_directory_path();

    auto branches_fa = tmp / "branch_project_e2e_branches.fa";
    auto ref_fa      = tmp / "branch_project_e2e_ref.fa";
    auto out_prefix  = (tmp / "branch_project_e2e_out").string();
    auto report_json = tmp / "branch_project_e2e_out.branch-report.json";

    {
        std::ofstream b(branches_fa);
        // 2 branches, second uniquely long enough to map.
        b << ">branch_1\n"
          << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
          << ">branch_2\n"
          << "TTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGA\n";
    }
    {
        // Reference contains both branch sequences at known positions.
        std::ofstream r(ref_fa);
        r << ">chrTest\n"
          << std::string(500, 'A')
          << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
          << std::string(200, 'G')
          << "TTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGATTGA"
          << std::string(300, 'C')
          << "\n";
    }

    std::string cmd = std::string(BRANCH_BINARY) +
        " project --fasta " + branches_fa.string() +
        " --ref-linear TestRef=" + ref_fa.string() +
        " --out-prefix " + out_prefix +
        " 2>/dev/null";
    int code = WEXITSTATUS(std::system(cmd.c_str()));

    fs::remove(branches_fa);
    fs::remove(ref_fa);

    ASSERT_EQ(code, 0) << "branch project failed (exit " << code << ")";
    ASSERT_TRUE(fs::exists(report_json))
        << "expected JSON report at " << report_json;

    std::ifstream rep(report_json);
    std::ostringstream ss;
    ss << rep.rdbuf();
    const std::string content = ss.str();
    fs::remove(report_json);

    // Structural checks on the JSON.
    EXPECT_NE(content.find("\"sample\":"),           std::string::npos);
    EXPECT_NE(content.find("\"version\": \"0.4.1\""), std::string::npos);
    EXPECT_NE(content.find("\"branch_count\": 2"),    std::string::npos);
    EXPECT_NE(content.find("\"branches\":"),          std::string::npos);
    EXPECT_NE(content.find("\"branch_id\": \"branch_1\""),
              std::string::npos);
    EXPECT_NE(content.find("\"branch_id\": \"branch_2\""),
              std::string::npos);
    // At least one mapping against TestRef should have succeeded on
    // an identical substring embedded in the reference.
    EXPECT_NE(content.find("\"ref\": \"TestRef\""),   std::string::npos);
}

TEST(ProjectSmoke, MainHelpShowsProject) {
    // `branch --help` should list "project" subcommand
    std::string cmd = std::string(BRANCH_BINARY) + " --help 2>&1";
    FILE* pipe = popen(cmd.c_str(), "r");
    ASSERT_NE(pipe, nullptr);

    std::string output;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        output += buffer;
    }
    int status = pclose(pipe);
    EXPECT_EQ(0, WEXITSTATUS(status));
    EXPECT_NE(output.find("project"), std::string::npos)
        << "Main help should list 'project' subcommand";
}
