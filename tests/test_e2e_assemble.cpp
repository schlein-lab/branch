// BRANCH v0.2 — end-to-end integration test.
//
// Spawns the real `branch` binary via fork+execv and feeds it a
// synthetic FASTQ produced by `generate_toy_fastq`. Verifies that the
// full FASTQ -> sketcher -> CPU overlap -> graph -> GFA pipeline runs
// to completion and emits a structurally valid BRANCH GFA-1.2 file.
//
// This is the guardrail for the whole v0.2 happy path: if any stage
// regresses (reader crash, sketcher underflow, backend seg-fault, GFA
// writer format drift), this test turns red.
//
// The compile definitions BRANCH_BINARY and GEN_FASTQ_BINARY are
// injected by CMake (target_compile_definitions) and point at the
// built executables.

#include <gtest/gtest.h>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#ifndef BRANCH_BINARY
#error "BRANCH_BINARY must be defined by the build system"
#endif
#ifndef GEN_FASTQ_BINARY
#error "GEN_FASTQ_BINARY must be defined by the build system"
#endif

namespace fs = std::filesystem;

namespace {

struct RunResult {
    int exit_code{-1};
    bool signaled{false};
    int term_signal{0};
    std::string stdout_text;
    std::string stderr_text;
};

// Fork+execv with pipe capture. Returns a RunResult; exit_code < 0 on
// spawn failure. Argv must be null-terminated.
RunResult run_capture(const std::vector<std::string>& argv) {
    RunResult r;
    if (argv.empty()) { r.stderr_text = "empty argv"; return r; }

    int out_pipe[2]{-1, -1};
    int err_pipe[2]{-1, -1};
    if (pipe(out_pipe) != 0 || pipe(err_pipe) != 0) {
        r.stderr_text = std::string("pipe() failed: ") + std::strerror(errno);
        return r;
    }

    pid_t pid = fork();
    if (pid < 0) {
        r.stderr_text = std::string("fork() failed: ") + std::strerror(errno);
        ::close(out_pipe[0]); ::close(out_pipe[1]);
        ::close(err_pipe[0]); ::close(err_pipe[1]);
        return r;
    }

    if (pid == 0) {
        // Child: dup pipes onto stdout/stderr, close stdin.
        ::dup2(out_pipe[1], STDOUT_FILENO);
        ::dup2(err_pipe[1], STDERR_FILENO);
        ::close(out_pipe[0]); ::close(out_pipe[1]);
        ::close(err_pipe[0]); ::close(err_pipe[1]);
        int devnull = ::open("/dev/null", O_RDONLY);
        if (devnull >= 0) { ::dup2(devnull, STDIN_FILENO); ::close(devnull); }

        std::vector<char*> cargv;
        cargv.reserve(argv.size() + 1);
        for (const auto& s : argv) {
            cargv.push_back(const_cast<char*>(s.c_str()));
        }
        cargv.push_back(nullptr);
        ::execv(cargv[0], cargv.data());
        // execv only returns on failure.
        std::fprintf(stderr, "execv(%s) failed: %s\n",
                     cargv[0], std::strerror(errno));
        std::_Exit(127);
    }

    // Parent: read both pipes to EOF.
    ::close(out_pipe[1]);
    ::close(err_pipe[1]);

    auto drain = [](int fd, std::string& into) {
        char buf[4096];
        for (;;) {
            ssize_t n = ::read(fd, buf, sizeof(buf));
            if (n > 0) into.append(buf, buf + n);
            else if (n == 0) break;
            else if (errno == EINTR) continue;
            else break;
        }
    };
    // Naive serial drain is fine: both pipes have 64KB kernel buffers
    // and the outputs here are tiny (< a few KB).
    drain(out_pipe[0], r.stdout_text);
    drain(err_pipe[0], r.stderr_text);
    ::close(out_pipe[0]);
    ::close(err_pipe[0]);

    int status = 0;
    while (::waitpid(pid, &status, 0) < 0) {
        if (errno == EINTR) continue;
        r.stderr_text += std::string("\nwaitpid() failed: ") + std::strerror(errno);
        return r;
    }
    if (WIFEXITED(status)) {
        r.exit_code = WEXITSTATUS(status);
    } else if (WIFSIGNALED(status)) {
        r.signaled = true;
        r.term_signal = WTERMSIG(status);
        r.exit_code = 128 + r.term_signal;
    }
    return r;
}

std::string read_file(const fs::path& p) {
    std::ifstream ifs(p);
    std::stringstream ss;
    ss << ifs.rdbuf();
    return ss.str();
}

// Count lines in `text` that start with `prefix` (prefix is matched at
// the start of each line; the trailing tab is part of the prefix so a
// full GFA record-type match works).
std::size_t count_prefix(const std::string& text, const std::string& prefix) {
    std::size_t count = 0;
    std::size_t pos = 0;
    while (pos < text.size()) {
        std::size_t eol = text.find('\n', pos);
        std::size_t line_end = (eol == std::string::npos) ? text.size() : eol;
        std::string_view line(text.data() + pos, line_end - pos);
        if (line.size() >= prefix.size() &&
            line.compare(0, prefix.size(), prefix) == 0) {
            ++count;
        }
        if (eol == std::string::npos) break;
        pos = eol + 1;
    }
    return count;
}

class E2EAssemble : public ::testing::Test {
 protected:
    fs::path tmpdir_;

    void SetUp() override {
        // Unique tempdir per test instance: temp_directory_path /
        // branch_e2e_<pid>_<rand>.
        std::random_device rd;
        const auto tag = std::to_string(::getpid()) + "_" +
                         std::to_string(rd());
        tmpdir_ = fs::temp_directory_path() / (std::string("branch_e2e_") + tag);
        std::error_code ec;
        fs::create_directories(tmpdir_, ec);
        ASSERT_FALSE(ec) << "failed to create " << tmpdir_ << ": " << ec.message();

        // Sanity: the binaries we were compiled against must exist.
        ASSERT_TRUE(fs::exists(BRANCH_BINARY))
            << "branch binary missing at " << BRANCH_BINARY;
        ASSERT_TRUE(fs::exists(GEN_FASTQ_BINARY))
            << "generate_toy_fastq missing at " << GEN_FASTQ_BINARY;
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmpdir_, ec);
        // Do not fail the test on cleanup errors — just note them.
        if (ec) {
            std::fprintf(stderr,
                         "warning: could not remove %s: %s\n",
                         tmpdir_.c_str(), ec.message().c_str());
        }
    }

    // Build a toy FASTQ at tmpdir_/name. Deterministic by seed.
    void generate(const std::string& name, std::uint64_t seed,
                  std::size_t n_reads, std::size_t read_len,
                  std::size_t mean_overlap) {
        const auto path = tmpdir_ / name;
        std::vector<std::string> argv = {
            GEN_FASTQ_BINARY,
            path.string(),
            std::to_string(seed),
            std::to_string(n_reads),
            std::to_string(read_len),
            std::to_string(mean_overlap),
        };
        auto r = run_capture(argv);
        ASSERT_EQ(r.exit_code, 0)
            << "generate_toy_fastq failed (exit " << r.exit_code
            << ", signaled=" << r.signaled << " sig=" << r.term_signal
            << ")\nstdout:\n" << r.stdout_text
            << "\nstderr:\n" << r.stderr_text;
        ASSERT_TRUE(fs::exists(path));
        ASSERT_GT(fs::file_size(path), 0u);
    }

    // Run `branch assemble --fastq <in> --out <out>` and return result.
    RunResult run_assemble(const fs::path& in_fastq, const fs::path& out_gfa) {
        std::vector<std::string> argv = {
            BRANCH_BINARY,
            "assemble",
            "--fastq", in_fastq.string(),
            "--out", out_gfa.string(),
        };
        return run_capture(argv);
    }
};

TEST_F(E2EAssemble, E2E_toy_fastq_produces_nonempty_gfa) {
    // 10 reads x 1000 bp each, overlap 500 bp -> graph should be a
    // chain of ~9 adjacency edges over 10 segments.
    ASSERT_NO_FATAL_FAILURE(generate("toy.fastq", /*seed=*/42,
                                     /*n_reads=*/10,
                                     /*read_len=*/1000,
                                     /*mean_overlap=*/500));
    const auto in_fastq = tmpdir_ / "toy.fastq";
    const auto out_gfa = tmpdir_ / "toy.gfa";

    auto r = run_assemble(in_fastq, out_gfa);
    ASSERT_EQ(r.exit_code, 0)
        << "branch assemble failed (exit " << r.exit_code
        << ", signaled=" << r.signaled << " sig=" << r.term_signal
        << ")\nstdout:\n" << r.stdout_text
        << "\nstderr:\n" << r.stderr_text;

    ASSERT_TRUE(fs::exists(out_gfa)) << "no GFA at " << out_gfa;
    ASSERT_GT(fs::file_size(out_gfa), 0u);

    const std::string gfa = read_file(out_gfa);

    // GFA-1.2 header must be first record.
    EXPECT_NE(gfa.find("H\tVN:Z:1.2"), std::string::npos)
        << "missing GFA-1.2 header";

    // BRANCH custom tag signature — CN:i + CV:f on S-lines are
    // unique-to-BRANCH extensions required by the output contract.
    EXPECT_NE(gfa.find("CN:i:"), std::string::npos) << "missing BRANCH CN:i tag";
    EXPECT_NE(gfa.find("CV:f:"), std::string::npos) << "missing BRANCH CV:f tag";

    const std::size_t s_lines = count_prefix(gfa, "S\t");
    const std::size_t l_lines = count_prefix(gfa, "L\t");

    EXPECT_GE(s_lines, 1u) << "expected >=1 S-line, gfa:\n" << gfa;
    // With unitig compaction enabled (v0.3+), a linear 10-read chain collapses
    // to 1 unitig with 0 internal edges. Older behavior (no compaction) kept
    // 10 S-lines + 9 L-lines. Accept either as long as the output is coherent
    // — the compacted form is the better answer biologically.
    EXPECT_LE(s_lines, 10u)
        << "expected <=10 S (10 reads, maybe compacted), got " << s_lines;
    EXPECT_LE(l_lines, s_lines)
        << "edges must not exceed segments, got L=" << l_lines << " S=" << s_lines;
}

TEST_F(E2EAssemble, E2E_toy_fastq_respects_min_minimizers) {
    // 50 bp reads are shorter than the typical minimizer window —
    // expect the pipeline to still exit cleanly with zero (or very
    // few) overlaps and a degenerate graph. The assembler must not
    // crash on degenerate input.
    ASSERT_NO_FATAL_FAILURE(generate("short.fastq", /*seed=*/7,
                                     /*n_reads=*/6,
                                     /*read_len=*/50,
                                     /*mean_overlap=*/10));
    const auto in_fastq = tmpdir_ / "short.fastq";
    const auto out_gfa = tmpdir_ / "short.gfa";

    auto r = run_assemble(in_fastq, out_gfa);
    ASSERT_EQ(r.exit_code, 0)
        << "branch assemble failed on short reads (exit " << r.exit_code
        << ", signaled=" << r.signaled << " sig=" << r.term_signal
        << ")\nstderr:\n" << r.stderr_text;

    ASSERT_TRUE(fs::exists(out_gfa));
    const std::string gfa = read_file(out_gfa);
    EXPECT_NE(gfa.find("H\tVN:Z:1.2"), std::string::npos);

    // Degenerate input: reads too short for a reliable overlap, so
    // expect zero links. If the backend ever gets more sensitive this
    // can be relaxed to "few" — for now 0 is the contract.
    const std::size_t l_lines = count_prefix(gfa, "L\t");
    EXPECT_EQ(l_lines, 0u)
        << "expected 0 L-lines for degenerate 50bp fixture, got "
        << l_lines << "\nstderr:\n" << r.stderr_text;
}

TEST_F(E2EAssemble, E2E_missing_input_fails_cleanly) {
    // Path we are reasonably confident does not exist.
    const fs::path missing = tmpdir_ / "definitely_not_here.fastq";
    ASSERT_FALSE(fs::exists(missing));

    const fs::path out_gfa = tmpdir_ / "unused.gfa";
    auto r = run_assemble(missing, out_gfa);

    EXPECT_NE(r.exit_code, 0)
        << "expected non-zero exit on missing input, got 0\nstderr:\n"
        << r.stderr_text;
    EXPECT_FALSE(r.signaled) << "branch was killed by signal "
                             << r.term_signal
                             << " on missing-input path (should exit cleanly)";

    const std::string err_lc = [&] {
        std::string s = r.stderr_text;
        for (auto& c : s) c = static_cast<char>(std::tolower(c));
        return s;
    }();
    EXPECT_TRUE(err_lc.find("not found") != std::string::npos ||
                err_lc.find("cannot open") != std::string::npos ||
                err_lc.find("no such file") != std::string::npos)
        << "stderr should mention missing file, got:\n" << r.stderr_text;
}

}  // namespace
