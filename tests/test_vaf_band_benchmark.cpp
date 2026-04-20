// BRANCH P1.3 — VAF-band benchmark harness.
//
// Full end-to-end guardrail:
//   1. Generate the 4-way bubble fixture (copy counts {0,2,4,6},
//      target VAFs {0.10, 0.50, 0.25, 0.15}) via the
//      generate_branch_fixture helper.
//   2. Run the real `branch` binary (assemble) on the fixture FASTQ.
//   3. Parse the emitted GFA's S-lines as "called bubbles": read CN:i
//      for copy-count and RC:i/total for VAF.
//   4. Parse the fixture meta TSV as truth.
//   5. Bin truth alleles by their target VAF; compute per-band recall,
//      precision, median |est_vaf - true_vaf|, and copy-count accuracy.
//   6. Write JSON report + a stdout table.
//   7. Gate: recall for the union of low bands (vaf <= 0.10) must be
//      >= kLowBandRecallFloor (0.5 — set permissively; the purpose of
//      P1.3 is to SHIP the harness, not hit a bar).
//
// JSON report schema (stable within schema_version):
// {
//   "schema_version": "p1.3-vaf-band-1",
//   "sample":         "<name>",
//   "n_truth":        <int>,
//   "n_called":       <int>,
//   "bands": [
//     {
//       "name": "0.00-0.05",
//       "lo": 0.0, "hi": 0.05,
//       "n_truth": <int>, "n_called": <int>,
//       "tp": <int>, "fp": <int>,
//       "recall": <float|null>, "precision": <float|null>,
//       "median_abs_vaf_residual": <float|null>,
//       "copy_count_accuracy":    <float|null>
//     },
//     ...
//   ]
// }
//
// BRANCH_BINARY and GEN_BRANCH_FIXTURE_BINARY are injected by CMake
// target_compile_definitions; see tests/CMakeLists.txt.

#include <gtest/gtest.h>

#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "analysis/vaf_band_stats.hpp"

#ifndef BRANCH_BINARY
#error "BRANCH_BINARY must be defined by the build system"
#endif
#ifndef GEN_BRANCH_FIXTURE_BINARY
#error "GEN_BRANCH_FIXTURE_BINARY must be defined by the build system"
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

class VafBandBenchmark : public ::testing::Test {
 protected:
    fs::path tmpdir_;

    void SetUp() override {
        std::random_device rd;
        const auto tag = std::to_string(::getpid()) + "_" +
                         std::to_string(rd());
        tmpdir_ = fs::temp_directory_path() /
                  (std::string("branch_vaf_band_") + tag);
        std::error_code ec;
        fs::create_directories(tmpdir_, ec);
        ASSERT_FALSE(ec) << "failed to create " << tmpdir_
                         << ": " << ec.message();

        ASSERT_TRUE(fs::exists(BRANCH_BINARY))
            << "branch binary missing at " << BRANCH_BINARY;
        ASSERT_TRUE(fs::exists(GEN_BRANCH_FIXTURE_BINARY))
            << "generate_branch_fixture missing at "
            << GEN_BRANCH_FIXTURE_BINARY;
    }

    void TearDown() override {
        // Preserve artefacts when a test fails so the JSON report and
        // GFA are available for postmortem.
        if (::testing::Test::HasFailure()) {
            std::fprintf(stderr,
                         "VafBandBenchmark: preserved artefacts at %s\n",
                         tmpdir_.c_str());
            return;
        }
        std::error_code ec;
        fs::remove_all(tmpdir_, ec);
    }
};

// Generate the fixture with canonical parameters (cassette=1000 so the
// 4-allele pool still fits inside auto-sized reads yet keeps compile +
// assemble wall time bounded). VAF mix matches the P1.3 spec.
TEST_F(VafBandBenchmark, Benchmark_end_to_end_fourway_bubble) {
    const fs::path fastq = tmpdir_ / "fx.fastq";
    const fs::path meta  = tmpdir_ / "fx.meta.tsv";
    const fs::path gfa   = tmpdir_ / "fx.gfa";
    const fs::path json  = tmpdir_ / "vaf_band_report.json";

    // Step 1. Fixture generation. Small read pool keeps the O(N^2)
    // all-vs-all overlap under ctest's default timeout even in Debug
    // / ASan builds: 80 reads split 0.10/0.50/0.25/0.15 gives
    // ~8 reads per 10% VAF — enough to exercise the low-band path.
    //
    // cassette=500 * 6 max copies + 2*2000 flank = ~7kb auto-sized
    // reads. Long enough to span the biggest bubble, short enough
    // that minimizer overlap is fast.
    std::vector<std::string> gen_argv = {
        GEN_BRANCH_FIXTURE_BINARY,
        fastq.string(),
        /*seed=*/"31",
        /*total_reads=*/"80",
        /*read_len=*/"5000",        // clamped up by --read-len-auto
        "--meta", meta.string(),
        "--cassette-bp", "500",
        "--vaf-mix", "0.10,0.50,0.25,0.15",
        "--read-len-auto",
    };
    {
        auto r = run_capture(gen_argv);
        ASSERT_EQ(r.exit_code, 0)
            << "fixture generator failed\nstderr:\n" << r.stderr_text;
        ASSERT_TRUE(fs::exists(fastq));
        ASSERT_GT(fs::file_size(fastq), 0u);
        ASSERT_TRUE(fs::exists(meta));
    }

    // Step 2. Run branch assemble.
    {
        std::vector<std::string> argv = {
            BRANCH_BINARY, "assemble",
            "--fastq", fastq.string(),
            "--out", gfa.string(),
            "--threads", "2",
        };
        auto r = run_capture(argv);
        ASSERT_EQ(r.exit_code, 0)
            << "branch assemble failed (exit " << r.exit_code << ")\n"
            << "stderr:\n" << r.stderr_text;
        ASSERT_TRUE(fs::exists(gfa));
        ASSERT_GT(fs::file_size(gfa), 0u);
    }

    // Step 3. Parse truth + calls.
    std::string err;
    auto truth = branch::analysis::parse_fixture_meta_tsv(meta.string(), &err);
    ASSERT_TRUE(err.empty()) << "meta parse error: " << err;
    ASSERT_EQ(truth.size(), 4u) << "expected 4 truth alleles";

    auto nodes = branch::analysis::parse_gfa_nodes(gfa.string(), &err);
    ASSERT_TRUE(err.empty()) << "GFA parse error: " << err;
    ASSERT_GT(nodes.size(), 0u) << "BRANCH produced no nodes";

    // P2.2+: every node carries real RC:i (read support) propagated
    // from graph_builder (1-per-raw-read) through graph_compactor
    // (summed across unitig members). build_calls_from_gfa prefers the
    // RC:i-driven path and falls back to length-bucket aggregation
    // only when sum(RC:i) == 0 (pre-P2.2 GFAs).
    branch::analysis::CalledBubbleSource call_source{};
    auto calls = branch::analysis::build_calls_from_gfa(
        nodes, /*bucket_size_bp=*/500, &call_source);
    ASSERT_GT(calls.size(), 0u) << "no calls after build_calls_from_gfa";
    EXPECT_EQ(call_source, branch::analysis::CalledBubbleSource::kReadSupport)
        << "RC:i should be populated post-P2.2 — benchmark took the "
           "length-bucket fallback, which means the compactor didn't "
           "propagate read_support onto the S-lines.";

    // Step 4. Compute per-band stats, write JSON + stdout table.
    auto report = branch::analysis::compute_band_stats(
        truth, calls, "synth_4way_bubble");
    ASSERT_TRUE(branch::analysis::write_report_json(
        report, json.string(), &err))
        << "failed to write JSON: " << err;
    ASSERT_TRUE(fs::exists(json));

    branch::analysis::write_report_table(report, std::cout);

    // Step 5. Gate — low-band (true vaf <= 0.10) recall must be
    //         >= kLowBandRecallFloor. The floor is intentionally
    //         permissive (0.5); tighten as the pipeline improves.
    //
    // The union of bands (0, 0.05] + (0.05, 0.10] covers all truth
    // alleles whose target VAF is in the low-VAF regime. We compute a
    // weighted recall across those bands.
    std::size_t low_tp = 0;
    std::size_t low_truth = 0;
    for (const auto& b : report.bands) {
        if (b.hi <= 0.10 + 1e-9) {
            low_tp    += b.n_true_positive;
            low_truth += b.n_truth;
        }
    }
    // Degenerate case: no low-band truth alleles -> gate is trivially
    // satisfied. That's the only way the assertion below can skip,
    // and it would indicate a fixture regression (test input changed).
    if (low_truth == 0) {
        GTEST_SKIP()
            << "no low-band (<=0.10) truth alleles — fixture regressed?";
    }
    const double low_recall =
        static_cast<double>(low_tp) / static_cast<double>(low_truth);

    EXPECT_GE(low_recall, branch::analysis::kLowBandRecallFloor)
        << "low-band recall " << low_recall
        << " < floor " << branch::analysis::kLowBandRecallFloor
        << " (TP=" << low_tp << "/truth=" << low_truth << ")."
        << " JSON: " << json;
}

// Unit coverage for assign_band edge cases — cheap and catches silly
// off-by-one regressions in the band table.
TEST(VafBandStats, assign_band_boundaries) {
    using branch::analysis::assign_band;
    using branch::analysis::kVafBands;
    EXPECT_EQ(assign_band(-0.1), kVafBands.size());
    EXPECT_EQ(assign_band(0.0),  kVafBands.size());   // (0, ...] excludes 0
    EXPECT_EQ(assign_band(1e-6), 0u);
    EXPECT_EQ(assign_band(0.05), 0u);                 // right-inclusive
    EXPECT_EQ(assign_band(0.0500001), 1u);
    EXPECT_EQ(assign_band(0.10), 1u);
    EXPECT_EQ(assign_band(0.10001), 2u);
    EXPECT_EQ(assign_band(0.50), 3u);
    EXPECT_EQ(assign_band(0.75), 4u);
    EXPECT_EQ(assign_band(1.0),  4u);
    EXPECT_EQ(assign_band(1.01), kVafBands.size());
}

// Unit coverage for the P2.2 RC:i preference: when any node has RC:i>0
// the benchmark must use per-node RC:i-driven VAF, not length buckets.
TEST(VafBandStats, build_calls_from_gfa_prefers_RC_when_populated) {
    using branch::analysis::CalledBubble;
    using branch::analysis::CalledBubbleSource;
    using branch::analysis::build_calls_from_gfa;

    std::vector<CalledBubble> nodes = {
        {.node_id = 0, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 10, .length_bp = 500},
        {.node_id = 1, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 30, .length_bp = 500},
        {.node_id = 2, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 60, .length_bp = 500},
    };

    CalledBubbleSource src{};
    auto calls = build_calls_from_gfa(nodes, 500, &src);
    EXPECT_EQ(src, CalledBubbleSource::kReadSupport);
    ASSERT_EQ(calls.size(), 3u);
    // sum = 100 -> VAFs 0.10, 0.30, 0.60.
    EXPECT_NEAR(calls[0].est_vaf, 0.10, 1e-9);
    EXPECT_NEAR(calls[1].est_vaf, 0.30, 1e-9);
    EXPECT_NEAR(calls[2].est_vaf, 0.60, 1e-9);
}

TEST(VafBandStats, build_calls_from_gfa_falls_back_to_length_bucket) {
    using branch::analysis::CalledBubble;
    using branch::analysis::CalledBubbleSource;
    using branch::analysis::build_calls_from_gfa;

    // All RC:i zero -> must hit the length-bucket fallback path.
    std::vector<CalledBubble> nodes = {
        {.node_id = 0, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 0, .length_bp = 500},
        {.node_id = 1, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 0, .length_bp = 500},
        {.node_id = 2, .est_copy_count = 1, .est_vaf = 0.0,
         .read_support = 0, .length_bp = 1000},
    };

    CalledBubbleSource src{};
    auto calls = build_calls_from_gfa(nodes, 500, &src);
    EXPECT_EQ(src, CalledBubbleSource::kLengthBucket);
    // Length-bucket aggregation: 2 buckets (500bp and 1000bp).
    EXPECT_EQ(calls.size(), 2u);
}

// Unit coverage for compute_band_stats core semantics. No subprocess,
// no fixture — just the pure function.
TEST(VafBandStats, compute_band_stats_basic) {
    // Four truth alleles spanning three bands:
    //   a0: vaf=0.10 -> band 1 "0.05-0.10"
    //   a6: vaf=0.15 -> band 2 "0.10-0.25"
    //   a4: vaf=0.25 -> band 2 "0.10-0.25"
    //   a2: vaf=0.50 -> band 3 "0.25-0.50"
    std::vector<branch::analysis::TruthAllele> truth = {
        {.allele_id = "a0", .copy_count = 0, .vaf_target = 0.10,
         .n_reads = 60, .cassette_len_bp = 1000},
        {.allele_id = "a2", .copy_count = 2, .vaf_target = 0.50,
         .n_reads = 300, .cassette_len_bp = 1000},
        {.allele_id = "a4", .copy_count = 4, .vaf_target = 0.25,
         .n_reads = 150, .cassette_len_bp = 1000},
        {.allele_id = "a6", .copy_count = 6, .vaf_target = 0.15,
         .n_reads = 90, .cassette_len_bp = 1000},
    };
    // Calls designed so that distance-based matching claims one call
    // per truth allele (all within kVafMatchTolerance = 0.10):
    //   call@0.11 -> a0 (dist 0.01)
    //   call@0.48 -> a2 (dist 0.02)
    //   call@0.26 -> a4 (dist 0.01)
    //   call@0.16 -> a6 (dist 0.01)
    // Extra call@0.03 has no nearby truth -> FP in band 0.
    std::vector<branch::analysis::CalledBubble> calls = {
        {.node_id = 1, .est_copy_count = 2, .est_vaf = 0.48,
         .read_support = 48, .length_bp = 3000},
        {.node_id = 2, .est_copy_count = 4, .est_vaf = 0.26,
         .read_support = 26, .length_bp = 5000},
        {.node_id = 3, .est_copy_count = 6, .est_vaf = 0.16,
         .read_support = 16, .length_bp = 7000},
        {.node_id = 4, .est_copy_count = 1, .est_vaf = 0.11,
         .read_support = 11, .length_bp = 500},
        {.node_id = 5, .est_copy_count = 3, .est_vaf = 0.03,
         .read_support = 3, .length_bp = 500},
    };

    auto rep = branch::analysis::compute_band_stats(truth, calls, "unit");
    EXPECT_EQ(rep.n_truth_total, 4u);
    EXPECT_EQ(rep.n_called_total, 5u);

    // Band 0 (0, 0.05]: no truth, 1 FP (the 0.03 spurious call).
    EXPECT_EQ(rep.bands[0].n_truth, 0u);
    EXPECT_EQ(rep.bands[0].n_false_positive, 1u);
    EXPECT_DOUBLE_EQ(rep.bands[0].precision(), 0.0);

    // Band 1 (0.05, 0.10]: truth a0 matched by call@0.11 -> 1 TP.
    EXPECT_EQ(rep.bands[1].n_truth, 1u);
    EXPECT_EQ(rep.bands[1].n_true_positive, 1u);
    EXPECT_DOUBLE_EQ(rep.bands[1].recall(), 1.0);

    // Band 2 (0.10, 0.25]: truth a4 + a6, both matched -> 2 TP.
    EXPECT_EQ(rep.bands[2].n_truth, 2u);
    EXPECT_EQ(rep.bands[2].n_true_positive, 2u);
    EXPECT_DOUBLE_EQ(rep.bands[2].recall(), 1.0);

    // Band 3 (0.25, 0.50]: truth a2 matched -> 1 TP.
    EXPECT_EQ(rep.bands[3].n_truth, 1u);
    EXPECT_EQ(rep.bands[3].n_true_positive, 1u);
}

}  // namespace
