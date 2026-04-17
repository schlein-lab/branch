// BRANCH v0.4 — Smoke tests for `branch project` subcommand.
//
// Tests:
// 1. `branch project --help` exits 0 with usage
// 2. `branch project` (no args) exits 2 with usage
// 3. `branch project --gfa x` (missing required args) exits 2

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

}  // namespace

TEST(ProjectSmoke, HelpExitsZero) {
    EXPECT_EQ(0, run_branch("project --help"));
}

TEST(ProjectSmoke, NoArgsExitsTwo) {
    EXPECT_EQ(2, run_branch("project"));
}

TEST(ProjectSmoke, MissingRequiredArgsExitsTwo) {
    // Only --gfa provided, missing --fasta, --ref-linear, --ref-pangenome, --out-prefix
    EXPECT_EQ(2, run_branch("project --gfa /nonexistent.gfa"));
}

TEST(ProjectSmoke, HelpWithDashH) {
    EXPECT_EQ(0, run_branch("project -h"));
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
