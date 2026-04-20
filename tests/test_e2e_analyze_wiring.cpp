// BRANCH P2.1 — `branch assemble` -> `branch analyze` wiring test.
//
// Exercises the end-to-end hand-off:
//   1. `branch assemble --fastq toy.fastq --out toy.gfa` writes a GFA
//      that downstream stages can consume (graph persistence check).
//   2. A hand-crafted `LosslessGraph` containing a clean two-alt bubble
//      is serialised to GFA via write_gfa, then `branch analyze --graph`
//      is invoked on it. The emitted BED is inspected for a non-trivial
//      confidence (strictly between 0.0 and 1.0) on at least one bubble.
//
// Splitting the two stages keeps the wiring check fast and deterministic
// — the assembler's per-run graph shape depends on seed, overlap noise,
// and the unitig compactor's input; those regress separately via
// test_e2e_assemble. Here we only need to know that the GFA handoff
// format survives a round-trip and that analyze classifies every bubble.

#include <gtest/gtest.h>

#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "graph/graph_io.hpp"
#include "graph/lossless_graph.hpp"

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
        ::dup2(out_pipe[1], STDOUT_FILENO);
        ::dup2(err_pipe[1], STDERR_FILENO);
        ::close(out_pipe[0]); ::close(out_pipe[1]);
        ::close(err_pipe[0]); ::close(err_pipe[1]);
        int devnull = ::open("/dev/null", O_RDONLY);
        if (devnull >= 0) { ::dup2(devnull, STDIN_FILENO); ::close(devnull); }

        std::vector<char*> cargv;
        cargv.reserve(argv.size() + 1);
        for (const auto& s : argv) cargv.push_back(const_cast<char*>(s.c_str()));
        cargv.push_back(nullptr);
        ::execv(cargv[0], cargv.data());
        std::fprintf(stderr, "execv(%s) failed: %s\n",
                     cargv[0], std::strerror(errno));
        std::_Exit(127);
    }

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

// Build a compact deterministic bubble graph.
// Topology:
//
//   src --> entry --+--> alt1 --+--> exit --> sink
//                   |           |
//                   +--> alt2 --+
//
// Alt 1 carries a branch-like support split (15 reads), alt 2 a
// duplication-like support split (25 reads). Entry/exit carry long
// flanking consensus sequences so feature extraction produces a
// non-zero FlankJaccardK31 and the disambiguator emits a real
// (non-0, non-1) confidence.
branch::graph::LosslessGraph build_bubble_graph() {
    using branch::graph::LosslessGraph;
    using branch::graph::NodeId;

    LosslessGraph g;
    NodeId src   = g.add_node(1000);
    NodeId entry = g.add_node(2000);
    NodeId alt1  = g.add_node(2500);
    NodeId alt2  = g.add_node(2500);
    NodeId exit_ = g.add_node(2000);
    NodeId sink  = g.add_node(1000);

    // Entry/exit share the same flanking sequence so FlankJaccardK31 is
    // high (pulling the disambiguator toward the Branch regime). The
    // alt cassettes carry distinct sequences; depth ratio of ~1x means
    // neither clean Branch nor clean Duplication — confidence lands in
    // the informative mid-range.
    std::mt19937_64 rng(0xBEEF);
    auto random_dna = [&](std::size_t n) {
        static const char kAlpha[4] = {'A','C','G','T'};
        std::string s(n, 'A');
        for (auto& c : s) c = kAlpha[rng() & 3];
        return s;
    };
    const std::string shared_flank = random_dna(2000);
    g.node(src).consensus   = random_dna(1000);
    g.node(entry).consensus = shared_flank;
    g.node(alt1).consensus  = random_dna(2500);
    g.node(alt2).consensus  = random_dna(2500);
    g.node(exit_).consensus = shared_flank;
    g.node(sink).consensus  = random_dna(1000);

    g.add_edge(src, entry, 40);
    g.add_edge(entry, alt1, 15);
    g.add_edge(entry, alt2, 25);
    g.add_edge(alt1, exit_, 15);
    g.add_edge(alt2, exit_, 25);
    g.add_edge(exit_, sink, 40);

    g.node(src).read_support   = 40;
    g.node(entry).read_support = 40;
    g.node(alt1).read_support  = 15;
    g.node(alt2).read_support  = 25;
    g.node(exit_).read_support = 40;
    g.node(sink).read_support  = 40;
    return g;
}

class E2EAnalyzeWiring : public ::testing::Test {
 protected:
    fs::path tmpdir_;

    void SetUp() override {
        std::random_device rd;
        const auto tag = std::to_string(::getpid()) + "_" + std::to_string(rd());
        tmpdir_ = fs::temp_directory_path() /
                  (std::string("branch_e2e_analyze_") + tag);
        std::error_code ec;
        fs::create_directories(tmpdir_, ec);
        ASSERT_FALSE(ec) << "failed to create " << tmpdir_ << ": " << ec.message();
        ASSERT_TRUE(fs::exists(BRANCH_BINARY))
            << "branch binary missing at " << BRANCH_BINARY;
        ASSERT_TRUE(fs::exists(GEN_FASTQ_BINARY))
            << "generate_toy_fastq missing at " << GEN_FASTQ_BINARY;
    }

    void TearDown() override {
        std::error_code ec;
        fs::remove_all(tmpdir_, ec);
    }
};

// 1. branch assemble produces a GFA consumable by branch analyze.
TEST_F(E2EAnalyzeWiring, Assemble_emits_gfa_that_analyze_can_load) {
    // Generate a linear toy FASTQ; the assembler collapses it to a
    // single unitig — no bubbles expected, but the GFA must parse.
    const auto fastq = tmpdir_ / "toy.fastq";
    {
        std::vector<std::string> argv = {
            GEN_FASTQ_BINARY, fastq.string(), "42", "10", "1000", "500",
        };
        auto r = run_capture(argv);
        ASSERT_EQ(r.exit_code, 0) << r.stderr_text;
    }
    const auto gfa = tmpdir_ / "toy.gfa";
    {
        std::vector<std::string> argv = {
            BRANCH_BINARY, "assemble",
            "--fastq", fastq.string(),
            "--out", gfa.string(),
            "--threads", "2",
        };
        auto r = run_capture(argv);
        ASSERT_EQ(r.exit_code, 0)
            << "branch assemble failed: " << r.stderr_text;
    }
    ASSERT_TRUE(fs::exists(gfa));
    ASSERT_GT(fs::file_size(gfa), 0u);

    // Graph-mode analyze on the assembled GFA. No bubble expected in
    // the linear toy, so bed is empty — but the graph load + the zero-
    // bubble classification path must complete cleanly.
    const auto bed = tmpdir_ / "toy.bubbles.bed";
    std::vector<std::string> argv = {
        BRANCH_BINARY, "analyze",
        "--graph", gfa.string(),
        "--out-bed", bed.string(),
    };
    auto r = run_capture(argv);
    EXPECT_EQ(r.exit_code, 0)
        << "branch analyze --graph failed on assembled toy.gfa: "
        << r.stderr_text;
    EXPECT_NE(r.stdout_text.find("# bubbles_classified="), std::string::npos)
        << "analyze must emit bubbles_classified summary";
    EXPECT_NE(r.stdout_text.find("# conservation iterations="),
              std::string::npos)
        << "analyze must invoke conservation solver on the loaded graph";
}

// 2. A hand-crafted bubble graph produces a non-trivial confidence.
//    This is the wiring assertion: detect::detect_bubbles finds the
//    bubble, feature_extractor fills the vector from the loaded graph,
//    disambiguate_hierarchical returns a confidence that the BED writer
//    serialises into column 5.
TEST_F(E2EAnalyzeWiring, Graph_mode_produces_nontrivial_confidence) {
    const auto gfa = tmpdir_ / "bubble.gfa";
    {
        auto g = build_bubble_graph();
        ASSERT_TRUE(branch::graph::write_gfa(g, gfa.string()));
    }
    ASSERT_TRUE(fs::exists(gfa));
    ASSERT_GT(fs::file_size(gfa), 0u);

    const auto bed = tmpdir_ / "bubble.bed";
    std::vector<std::string> argv = {
        BRANCH_BINARY, "analyze",
        "--graph", gfa.string(),
        "--out-bed", bed.string(),
    };
    auto r = run_capture(argv);
    ASSERT_EQ(r.exit_code, 0)
        << "branch analyze failed: " << r.stderr_text;
    ASSERT_TRUE(fs::exists(bed)) << "BED not emitted";

    // Parse the BED and find at least one row whose column-5 confidence
    // is strictly between 0 and 1 (the classifier produced a real
    // non-trivial signal, not the 0.0-default or a 1.0 saturation).
    const std::string text = read_file(bed);
    ASSERT_FALSE(text.empty()) << "BED is empty — no bubble reached writer";

    std::size_t rows = 0;
    std::size_t nontrivial = 0;
    float seen_confidence = 0.0f;
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.empty() || line[0] == '#') continue;
        ++rows;
        // BedEntry::write_bed_entry emits tab-separated
        // chrom/start/end/name/confidence.
        std::vector<std::string> cols;
        std::size_t p = 0;
        while (p < line.size()) {
            auto q = line.find('\t', p);
            if (q == std::string::npos) q = line.size();
            cols.push_back(line.substr(p, q - p));
            p = (q == line.size()) ? q : q + 1;
        }
        ASSERT_GE(cols.size(), 5u)
            << "analyze BED row must have >=5 columns, got: " << line;
        const std::string& conf_tok = cols[4];
        if (conf_tok == "." || conf_tok.empty()) continue;
        char* end = nullptr;
        float c = std::strtof(conf_tok.c_str(), &end);
        if (end == conf_tok.c_str()) continue;
        if (c > 0.0f && c < 1.0f) {
            ++nontrivial;
            seen_confidence = c;
        }
    }
    EXPECT_GE(rows, 1u) << "expected >=1 bubble row, got 0";
    EXPECT_GE(nontrivial, 1u)
        << "expected at least one bubble with 0 < confidence < 1, "
           "but every row was 0.0 or 1.0. stdout:\n"
        << r.stdout_text << "\nstderr:\n" << r.stderr_text
        << "\nBED:\n" << text;
    // Sanity: confidence is a probability in [0, 1]. Caught by above
    // but worth asserting separately for the diagnostic path.
    EXPECT_GT(seen_confidence, 0.0f);
    EXPECT_LT(seen_confidence, 1.0f);
}

}  // namespace
