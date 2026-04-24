// Tests for the backend factory + GPU backend wiring.
//
// Covers the make_backend(mode) dispatch + make_gpu_backend runtime
// probe. The GPU compute_overlaps smoke test only runs when a CUDA
// device is visible, so CPU-only CI stays green.

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "backend/backend_factory.hpp"
#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "backend/gpu_backend.hpp"

using branch::backend::Backend;
using branch::backend::BackendMode;
using branch::backend::OverlapPair;
using branch::backend::parse_backend_mode;
using branch::backend::ReadBatch;

TEST(BackendFactory, Parse_mode_accepts_auto_cpu_gpu_and_flags_unknown) {
    EXPECT_EQ(parse_backend_mode("auto"), BackendMode::Auto);
    EXPECT_EQ(parse_backend_mode("cpu"),  BackendMode::Cpu);
    EXPECT_EQ(parse_backend_mode("gpu"),  BackendMode::Gpu);
    EXPECT_EQ(parse_backend_mode(""),     BackendMode::Auto);

    std::string err;
    auto m = parse_backend_mode("cupy", &err);
    EXPECT_EQ(m, BackendMode::Auto);
    EXPECT_NE(err.find("unrecognised"), std::string::npos);
}

TEST(BackendFactory, Cpu_mode_always_returns_usable_backend) {
    Backend b = branch::backend::make_backend(BackendMode::Cpu);
    ASSERT_FALSE(b.empty());
    EXPECT_STREQ(b.name(), branch::backend::kCpuBackendName);
}

TEST(BackendFactory, Auto_mode_returns_usable_backend_even_on_cpu_only) {
    // On CPU-only builds or machines without a visible GPU the Auto
    // path must still yield a working Backend so callers never stall.
    Backend b = branch::backend::make_backend(BackendMode::Auto);
    ASSERT_FALSE(b.empty());
}

TEST(GpuBackendWiring, Empty_on_cpu_only_build_or_device_absent) {
    const bool compiled = branch::backend::gpu_backend_compiled_in();
    Backend b = branch::backend::make_gpu_backend();

    if (!compiled) {
        EXPECT_TRUE(b.empty())
            << "non-CUDA build must return empty from make_gpu_backend";
        return;
    }
    // Compiled in: empty iff no visible device. Either is self-consistent.
    if (b.empty()) {
        GTEST_SKIP() << "CUDA compiled in but no visible device at test time";
    }
    EXPECT_STREQ(b.name(), branch::backend::kGpuBackendName);
}

TEST(GpuBackendWiring, Compute_overlaps_smoke_when_device_available) {
    if (!branch::backend::gpu_backend_compiled_in()) {
        GTEST_SKIP() << "CUDA not compiled in";
    }
    Backend b = branch::backend::make_gpu_backend();
    if (b.empty()) {
        GTEST_SKIP() << "no CUDA device available at runtime";
    }

    // Two synthetic reads sharing a long common subsequence so the
    // sketcher + kernel have something to chew on. We do not assert a
    // specific pair count — the exact number depends on kernel
    // thresholds and sketcher parameters that can shift between
    // versions. What we assert: the call returns, does not crash,
    // writes out_count within bounds.
    std::string seq_a;
    seq_a.reserve(5000);
    for (int i = 0; i < 5000; ++i) seq_a.push_back("ACGT"[i % 4]);
    std::string seq_b = seq_a;

    ReadBatch batch;
    batch.reads.push_back({seq_a, 0});
    batch.reads.push_back({seq_b, 1});

    std::vector<OverlapPair> out(1024);
    std::size_t n = 0;
    b.compute_overlaps(&batch, out, &n);
    EXPECT_LE(n, out.size());
}
