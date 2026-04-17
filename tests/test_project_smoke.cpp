// BRANCH v0.4 — Smoke tests for `branch project` subcommand.
//
// Tests:
// 1. `branch project --help` exits 0 with usage
// 2. `branch project` (no args) exits 2 with usage
// 3. Missing required args exits 2
// 4. Nonexistent files with all required args exits with error (not 0)

#include <gtest/gtest.h>

#include <cstdlib>
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
