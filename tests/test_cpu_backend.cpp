#include <gtest/gtest.h>

#include "backend/cpu_backend.hpp"
#include "classify/features.hpp"

using branch::backend::Backend;
using branch::backend::BubbleClass;
using branch::backend::ClassificationResult;
using branch::backend::kCpuBackendName;
using branch::backend::make_cpu_backend;
using branch::backend::OverlapPair;
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
