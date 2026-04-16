#include <gtest/gtest.h>

#include "gpu/overlap_kernel.cuh"

using branch::gpu::MinimizerEntry;
using branch::gpu::OverlapKernelLaunchConfig;
using branch::gpu::launch_overlap_kernel;

TEST(GpuStubTest, MinimizerEntry_packs_to_16_bytes) {
    EXPECT_EQ(sizeof(MinimizerEntry), 16u);
}

TEST(GpuStubTest, Stub_launcher_returns_zero_pairs) {
    OverlapKernelLaunchConfig cfg{
        .n_minimizers = 0,
        .max_output_pairs = 0,
    };
    std::size_t pairs = launch_overlap_kernel(nullptr, nullptr, cfg);
    EXPECT_EQ(pairs, 0u);
}
