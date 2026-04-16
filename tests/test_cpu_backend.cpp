#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "backend/cpu_backend.hpp"
#include "classify/features.hpp"

using branch::backend::Backend;
using branch::backend::BubbleClass;
using branch::backend::ClassificationResult;
using branch::backend::kCpuBackendName;
using branch::backend::make_cpu_backend;
using branch::backend::OverlapPair;
using branch::backend::ReadBatch;
using branch::backend::VAFEstimate;

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

TEST(CpuBackendTest, Classify_stub_returns_Unknown) {
    auto b = make_cpu_backend();
    // Must provide candidates so the stub loop iterates.
    std::array<branch::backend::BubbleCandidate, 3> cands{};
    std::array<ClassificationResult, 3> out{};
    b.classify_batch(cands, std::span<ClassificationResult>(out));
    for (auto& r : out) {
        EXPECT_EQ(r.label, BubbleClass::Unknown);
        EXPECT_FLOAT_EQ(r.confidence, 0.0f);
    }
}

TEST(CpuBackendTest, VAFEstimate_stub_returns_uniform_prior) {
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

TEST(CpuBackendTest, Destructor_cleans_up_without_leak) {
    // No explicit leak-check here; run with BRANCH_USE_ASAN=ON in CI
    // to detect allocator mismatches.
    {
        auto b = make_cpu_backend();
        EXPECT_FALSE(b.empty());
    }
    // b destroyed — under ASan a leak or double-free would fire.
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
