#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "backend/cpu_backend.hpp"
#include "backend/vaf_stats.hpp"
#include "classify/cascade.hpp"
#include "classify/features.hpp"
#include "graph/lossless_graph.hpp"

using branch::backend::Backend;
using branch::backend::BubbleClass;
using branch::backend::ClassificationResult;
using branch::backend::cpu_classify_batch_with_graph;
using branch::backend::cpu_estimate_vaf_batch_candidates;
using branch::backend::kCpuBackendName;
using branch::backend::make_cpu_backend;
using branch::backend::OverlapPair;
using branch::backend::ReadBatch;
using branch::backend::VAFEstimate;
using branch::backend::wilson_ci;
using branch::classify::BubbleCandidate;
using branch::classify::Feature;

TEST(CpuBackendTest, Creates_non_empty_backend) {
    auto b = make_cpu_backend();
    EXPECT_FALSE(b.empty());
    EXPECT_STREQ(b.name(), kCpuBackendName);
}

TEST(CpuBackendTest, Overlap_stub_returns_zero_pairs) {
    auto b = make_cpu_backend();
    std::size_t count = 999;
    b.compute_overlaps(nullptr, {}, &count);
    EXPECT_EQ(count, 0u);
}

TEST(CpuBackendTest, Classify_default_candidate_returns_NonSeparable) {
    // v0.2 wiring: classify_batch drives the real cascade. A default-
    // initialised candidate has BubbleLengthBp = 0 which trips the
    // cascade's min_bubble_length guard and short-circuits to
    // NonSeparable. (Pre-v0.2 this path returned Unknown from a stub.)
    auto b = make_cpu_backend();
    std::array<branch::backend::BubbleCandidate, 3> cands{};
    std::array<ClassificationResult, 3> out{};
    b.classify_batch(cands, std::span<ClassificationResult>(out));
    for (auto& r : out) {
        EXPECT_EQ(r.label, BubbleClass::NonSeparable);
        EXPECT_NE(r.label, BubbleClass::Unknown)
            << "classify_batch must not hand back the v0.1 stub label";
    }
}

TEST(CpuBackendTest, Classify_batch_calls_cascade) {
    // Synthetic candidate that the cascade should label as Duplication:
    // Stage 1 (FlankJaccard) fails (< 0.99), Stage 2 (DepthRatio >= 2.0)
    // fires. BubbleLength satisfies the early guard.
    BubbleCandidate cand{};
    cand.features[static_cast<std::size_t>(Feature::BubbleLengthBp)] = 1200.0f;
    cand.features[static_cast<std::size_t>(Feature::FlankJaccardK31)] = 0.6f;
    cand.features[static_cast<std::size_t>(Feature::DepthRatioDiploid)] = 2.5f;

    auto b = make_cpu_backend();
    std::array<BubbleCandidate, 1> cands{cand};
    std::array<ClassificationResult, 1> out{};
    b.classify_batch(cands, std::span<ClassificationResult>(out));

    EXPECT_EQ(out[0].label, BubbleClass::Duplication);
    EXPECT_NE(out[0].label, BubbleClass::Unknown);
    EXPECT_EQ(out[0].stage, 1u);
    EXPECT_GT(out[0].confidence, 0.0f);
}

TEST(CpuBackendTest, Classify_batch_with_graph_runs_feature_extraction) {
    // End-to-end: candidate carries only raw graph fields (no features
    // populated), extractor fills them, cascade labels the result.
    // Build a tiny graph whose entry node has read_support = 40 and
    // one baseline edge so the diploid-baseline is sensible.
    branch::graph::LosslessGraph g;
    const auto entry = g.add_node(/*length_bp=*/100);
    const auto exit  = g.add_node(/*length_bp=*/100);
    g.node(entry).read_support = 40;
    g.add_edge(entry, exit, /*read_support=*/20);

    BubbleCandidate cand{};
    cand.entry_node = entry;
    cand.exit_node  = exit;
    cand.bubble_length_bp = 1500;
    cand.read_support_branch = 22;  // branch + alt = 44, baseline = 20
    cand.read_support_alt    = 22;  // DepthRatio = 44/20 = 2.2 >= 2.0

    std::array<BubbleCandidate, 1> cands{cand};
    std::array<ClassificationResult, 1> out{};
    cpu_classify_batch_with_graph(cands, g, out);

    EXPECT_EQ(out[0].label, BubbleClass::Duplication);
    EXPECT_EQ(out[0].stage, 1u);
}

TEST(CpuBackendTest, VAFEstimate_legacy_id_path_returns_uniform_prior) {
    // Legacy ID-based entry keeps returning the uninformative prior
    // until a branch-id -> read-support index is introduced.
    auto b = make_cpu_backend();
    std::array<std::uint32_t, 2> ids = {42, 99};
    std::array<VAFEstimate, 2> out{};
    b.estimate_vaf_batch(ids, out);
    for (auto& v : out) {
        EXPECT_FLOAT_EQ(v.point, 0.5f);
        EXPECT_FLOAT_EQ(v.ci_low, 0.0f);
        EXPECT_FLOAT_EQ(v.ci_high, 1.0f);
    }
}

TEST(CpuBackendTest, Vaf_wilson_matches_reference) {
    // Reference values for n=100, k=80 (two-sided 95% Wilson interval):
    //   p        = 0.80
    //   ci_low  ~= 0.7111
    //   ci_high ~= 0.8666
    const auto est = wilson_ci(/*successes=*/80, /*total=*/100);
    EXPECT_NEAR(est.point,   0.80f, 0.01f);
    EXPECT_NEAR(est.ci_low,  0.71f, 0.01f);
    EXPECT_NEAR(est.ci_high, 0.87f, 0.01f);
}

TEST(CpuBackendTest, Vaf_zero_total_returns_uninformative) {
    const auto est = wilson_ci(0, 0);
    EXPECT_FLOAT_EQ(est.point, 0.5f);
    EXPECT_FLOAT_EQ(est.ci_low, 0.0f);
    EXPECT_FLOAT_EQ(est.ci_high, 1.0f);
}

TEST(CpuBackendTest, Vaf_batch_candidates_uses_read_support) {
    // Two candidates: one with real counts, one with zero counts.
    BubbleCandidate c_real{};
    c_real.read_support_branch = 80;
    c_real.read_support_alt    = 20;

    BubbleCandidate c_empty{};
    c_empty.read_support_branch = 0;
    c_empty.read_support_alt    = 0;

    std::array<BubbleCandidate, 2> cands{c_real, c_empty};
    std::array<VAFEstimate, 2> out{};
    cpu_estimate_vaf_batch_candidates(cands, out);

    EXPECT_NEAR(out[0].point,   0.80f, 0.01f);
    EXPECT_NEAR(out[0].ci_low,  0.71f, 0.01f);
    EXPECT_NEAR(out[0].ci_high, 0.87f, 0.01f);

    EXPECT_FLOAT_EQ(out[1].point, 0.5f);
    EXPECT_FLOAT_EQ(out[1].ci_low, 0.0f);
    EXPECT_FLOAT_EQ(out[1].ci_high, 1.0f);
}

TEST(CpuBackendTest, Destructor_cleans_up_without_leak) {
    // No explicit leak-check here; run with BRANCH_USE_ASAN=ON in CI
    // to detect allocator mismatches.
    {
        auto b = make_cpu_backend();
        EXPECT_FALSE(b.empty());
    }
    // b destroyed — under ASan a leak or double-free would fire.
}

TEST(CpuBackendTest, compute_overlaps_offsets_are_dovetail_consistent) {
    // Construct a 7 kb dovetail: the last 7000 bp of read A equal the first
    // 7000 bp of read B.
    //
    // Layout (with prefix_overhang = 8000):
    //   A: [0 ............... 8000][8000 ............... 15000]
    //                               ^-- shared 7000 bp --------^
    //   B:                         [0 ............... 7000][7000 .... 15000]
    //                               ^-- shared 7000 bp --------^
    //
    // BRANCH convention is "offset_x = 0-based start of overlap within read
    // x", so the expected pair (in canonical read_id order A=0 < B=1) is:
    //   offset_a ~= 8000, offset_b == 0.
    //
    // We use a deterministic LCG generator so the sequence is non-periodic
    // (repeating a short motif would make every minimizer window collide at
    // every motif-period offset, hiding the true dovetail offset).
    auto make_seq = [](std::size_t len, unsigned seed) {
        const char bases[] = "ACGT";
        std::string s;
        s.reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            seed = seed * 1103515245 + 12345;
            s.push_back(bases[(seed >> 16) & 3]);
        }
        return s;
    };

    constexpr std::size_t kReadLen      = 15000;
    constexpr std::size_t kPrefixOver   = 8000;   // A-only head, B-only tail
    constexpr std::size_t kSharedLen    = kReadLen - kPrefixOver;  // 7000

    const std::string read_a = make_seq(kReadLen, /*seed=*/0xB1A5C0DEu);

    // Read B = last 7000 bp of A (shared region), then 8000 bp of fresh
    // non-overlapping sequence so the suffix does NOT produce minimizer
    // collisions with A.
    std::string read_b;
    read_b.reserve(kReadLen);
    read_b.append(read_a, kPrefixOver, kSharedLen);
    read_b.append(make_seq(kPrefixOver, /*seed=*/0xDEADBEEFu));
    ASSERT_EQ(read_b.size(), kReadLen);

    ReadBatch batch;
    batch.reads.push_back({std::string_view(read_a), 0});
    batch.reads.push_back({std::string_view(read_b), 1});

    auto b = make_cpu_backend();
    std::array<OverlapPair, 16> out_pairs{};
    std::size_t count = 0;
    b.compute_overlaps(&batch, out_pairs, &count);

    // We expect exactly one overlap between reads 0 and 1.
    ASSERT_GE(count, 1u)
        << "Expected at least one overlap for a 7 kb dovetail";

    // Locate the pair (0,1) -- canonical order places the smaller read_id
    // in read_a.
    const OverlapPair* hit = nullptr;
    for (std::size_t i = 0; i < count; ++i) {
        if (out_pairs[i].read_a == 0 && out_pairs[i].read_b == 1) {
            hit = &out_pairs[i];
            break;
        }
    }
    ASSERT_NE(hit, nullptr) << "No (0,1) overlap emitted";

    // Dovetail convention: A starts first, so A's overlap begins deep
    // inside A (offset_a > 0), while B's overlap begins at B's very start
    // (offset_b == 0).
    EXPECT_GT(hit->offset_a, 0u)
        << "offset_a should mark how far into A the overlap starts";
    EXPECT_EQ(hit->offset_b, 0u)
        << "offset_b must be zero when B joins the overlap at its start";

    // offset_a should approximate the prefix overhang within +/- 200 bp
    // (minimizers sample every ~w bases, so the first shared minimizer
    // can land a short distance into the shared region).
    const std::int64_t expected = static_cast<std::int64_t>(kPrefixOver);
    const std::int64_t actual   = static_cast<std::int64_t>(hit->offset_a);
    EXPECT_LE(std::abs(actual - expected), 200)
        << "offset_a=" << actual
        << " should be within +/- 200 bp of prefix overhang " << expected;
}

TEST(CpuBackendTest, compute_overlaps_overlap_len_is_actual_span) {
    // Regression guard: overlap_len must be the actual span (in bp) of the
    // overlapping region, not a proxy like `match_count * 10`.
    //
    // Construct two reads with a controlled 5000 bp overlap:
    //   read A = 2000 bp unique prefix + 5000 bp shared
    //   read B = 5000 bp shared        + 2000 bp unique suffix
    // So the last 5000 bp of A == the first 5000 bp of B.
    auto make_seq = [](std::size_t len, unsigned seed) {
        const char bases[] = "ACGT";
        std::string s;
        s.reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            seed = seed * 1103515245 + 12345;
            s.push_back(bases[(seed >> 16) & 3]);
        }
        return s;
    };

    const std::string prefix_a = make_seq(2000, 0xABCDu);
    const std::string shared   = make_seq(5000, 0x1234u);
    const std::string suffix_b = make_seq(2000, 0xBEEFu);
    const std::string read_a = prefix_a + shared;
    const std::string read_b = shared + suffix_b;

    ReadBatch batch;
    batch.reads.push_back({std::string_view(read_a), 0});
    batch.reads.push_back({std::string_view(read_b), 1});

    auto b = make_cpu_backend();
    std::array<OverlapPair, 16> out_pairs{};
    std::size_t count = 0;
    b.compute_overlaps(&batch, out_pairs, &count);

    ASSERT_GE(count, 1u) << "Expected >= 1 overlap for 5000 bp dovetail";

    const OverlapPair* hit = nullptr;
    for (std::size_t i = 0; i < count; ++i) {
        if ((out_pairs[i].read_a == 0 && out_pairs[i].read_b == 1) ||
            (out_pairs[i].read_a == 1 && out_pairs[i].read_b == 0)) {
            hit = &out_pairs[i];
            break;
        }
    }
    ASSERT_NE(hit, nullptr) << "No overlap reported for reads (0, 1)";

    // Minimizer positions only sample the overlap — the first and last
    // matching minimizer sit at most one (w+k) window from the overlap
    // boundaries, so allow ~one-window slack each side. k=21 is accounted
    // for in the implementation.
    EXPECT_GE(hit->overlap_len, 4500u)
        << "overlap_len=" << hit->overlap_len
        << " too small — appears to undercount the 5000 bp span";
    EXPECT_LE(hit->overlap_len, 5500u)
        << "overlap_len=" << hit->overlap_len
        << " too large — looks like a non-span proxy (e.g. match_count * 10)";
}

TEST(CpuBackendTest, Overlap_detects_shared_minimizers) {
    // Two reads sharing a 200bp identical substring should produce >= 1 overlap.
    // Read 1: 300bp with shared region at [50..250)
    // Read 2: 300bp with shared region at [100..300)
    // The shared region is identical -> minimizer hits should match.

    // Generate deterministic "random" sequence
    auto make_seq = [](std::size_t len, unsigned seed) {
        const char bases[] = "ACGT";
        std::string s;
        s.reserve(len);
        for (std::size_t i = 0; i < len; ++i) {
            seed = seed * 1103515245 + 12345;
            s.push_back(bases[(seed >> 16) & 3]);
        }
        return s;
    };

    // Shared 200bp region
    std::string shared = make_seq(200, 42);

    // Read 1: 50bp prefix + 200bp shared + 50bp suffix
    std::string read1 = make_seq(50, 100) + shared + make_seq(50, 200);

    // Read 2: 100bp prefix + 200bp shared
    std::string read2 = make_seq(100, 300) + shared;

    // Build ReadBatch
    ReadBatch batch;
    batch.reads.push_back({std::string_view(read1), 0});
    batch.reads.push_back({std::string_view(read2), 1});

    auto b = make_cpu_backend();
    std::array<OverlapPair, 10> out_pairs{};
    std::size_t count = 0;

    b.compute_overlaps(&batch, out_pairs, &count);

    // Should detect at least 1 overlap between reads 0 and 1
    EXPECT_GE(count, 1u) << "Expected >= 1 overlap for 200bp shared region";

    if (count > 0) {
        // Verify the pair is (0, 1) in canonical order
        bool found_pair = false;
        for (std::size_t i = 0; i < count; ++i) {
            if ((out_pairs[i].read_a == 0 && out_pairs[i].read_b == 1) ||
                (out_pairs[i].read_a == 1 && out_pairs[i].read_b == 0)) {
                found_pair = true;
                break;
            }
        }
        EXPECT_TRUE(found_pair) << "Overlap should be between reads 0 and 1";
    }
}
