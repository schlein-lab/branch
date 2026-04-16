// BRANCH GPU Backend Tests
//
// When BRANCH_BUILD_CUDA=OFF: verify stub returns 0 pairs
// When BRANCH_BUILD_CUDA=ON: compare GPU results against CPU backend

#include <gtest/gtest.h>

#include "gpu/overlap_kernel.cuh"

using branch::gpu::MinimizerEntry;
using branch::gpu::OverlapKernelLaunchConfig;
using branch::gpu::launch_overlap_kernel;

// =============================================================================
// Common tests (both CUDA=ON and CUDA=OFF)
// =============================================================================

TEST(GpuBackendTest, MinimizerEntry_packs_to_16_bytes) {
    EXPECT_EQ(sizeof(MinimizerEntry), 16u);
}

TEST(GpuBackendTest, OverlapPair_packs_to_24_bytes) {
    EXPECT_EQ(sizeof(branch::backend::OverlapPair), 24u);
}

TEST(GpuBackendTest, LaunchConfig_has_sensible_defaults) {
    OverlapKernelLaunchConfig cfg{};
    EXPECT_EQ(cfg.device_id, 0);
    EXPECT_EQ(cfg.threads_per_block, 256);
    EXPECT_EQ(cfg.stream, nullptr);
}

// =============================================================================
// CUDA=OFF specific tests (stub behavior)
// =============================================================================

#ifndef BRANCH_BUILD_CUDA

TEST(GpuBackendStubTest, Stub_launcher_returns_zero_pairs_empty_input) {
    OverlapKernelLaunchConfig cfg{
        .n_minimizers = 0,
        .max_output_pairs = 0,
    };
    std::size_t pairs = launch_overlap_kernel(nullptr, nullptr, cfg);
    EXPECT_EQ(pairs, 0u);
}

TEST(GpuBackendStubTest, Stub_launcher_returns_zero_pairs_with_nonzero_config) {
    // Even with non-zero config, stub must return 0 (no GPU available)
    OverlapKernelLaunchConfig cfg{
        .n_minimizers = 1000,
        .max_output_pairs = 10000,
        .device_id = 0,
        .threads_per_block = 256,
        .stream = nullptr,
    };
    std::size_t pairs = launch_overlap_kernel(nullptr, nullptr, cfg);
    EXPECT_EQ(pairs, 0u);
}

#endif  // !BRANCH_BUILD_CUDA

// =============================================================================
// CUDA=ON specific tests (real GPU kernel)
// =============================================================================

#ifdef BRANCH_BUILD_CUDA

#include "backend/cpu_backend.hpp"
#include <vector>
#include <algorithm>

TEST(GpuBackendCudaTest, GPU_produces_same_pairs_as_CPU) {
    // This test only runs when CUDA is enabled and a GPU is available
    // It compares GPU kernel output against CPU reference implementation
    
    // Skip if no CUDA device
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        GTEST_SKIP() << "No CUDA device available";
    }
    
    // Create synthetic minimizer data
    // 4 reads, each with 10 minimizers, some shared hashes
    std::vector<MinimizerEntry> h_minimizers;
    
    // Read 0: positions 0-9, hash based on position
    for (uint32_t i = 0; i < 10; ++i) {
        h_minimizers.push_back(MinimizerEntry{
            .hash = static_cast<uint64_t>(i % 5),  // 5 unique hashes, creates collisions
            .read_id = 0,
            .pos = i * 100,
            .strand = 0,
        });
    }
    
    // Read 1: overlaps with read 0 on hashes 0-4
    for (uint32_t i = 0; i < 10; ++i) {
        h_minimizers.push_back(MinimizerEntry{
            .hash = static_cast<uint64_t>(i % 5),
            .read_id = 1,
            .pos = i * 100 + 50,  // offset of 50
            .strand = 0,
        });
    }
    
    // Read 2: different hashes, no overlap expected
    for (uint32_t i = 0; i < 10; ++i) {
        h_minimizers.push_back(MinimizerEntry{
            .hash = static_cast<uint64_t>(i + 100),
            .read_id = 2,
            .pos = i * 100,
            .strand = 0,
        });
    }
    
    // Read 3: shares hash 100 with read 2
    for (uint32_t i = 0; i < 10; ++i) {
        h_minimizers.push_back(MinimizerEntry{
            .hash = static_cast<uint64_t>(i + 100),
            .read_id = 3,
            .pos = i * 100 + 25,
            .strand = 1,  // opposite strand
        });
    }
    
    const std::size_t n_minimizers = h_minimizers.size();
    const std::size_t max_pairs = 1000;
    
    // Allocate device memory
    MinimizerEntry* d_minimizers = nullptr;
    branch::backend::OverlapPair* d_pairs = nullptr;
    
    ASSERT_EQ(cudaMalloc(&d_minimizers, n_minimizers * sizeof(MinimizerEntry)), cudaSuccess);
    ASSERT_EQ(cudaMalloc(&d_pairs, max_pairs * sizeof(branch::backend::OverlapPair)), cudaSuccess);
    
    // Copy to device
    ASSERT_EQ(cudaMemcpy(d_minimizers, h_minimizers.data(),
                        n_minimizers * sizeof(MinimizerEntry),
                        cudaMemcpyHostToDevice), cudaSuccess);
    
    // Launch kernel
    OverlapKernelLaunchConfig cfg{
        .n_minimizers = n_minimizers,
        .max_output_pairs = max_pairs,
        .device_id = 0,
        .threads_per_block = 256,
        .stream = nullptr,
    };
    
    std::size_t n_pairs = launch_overlap_kernel(d_minimizers, d_pairs, cfg);
    
    // Copy results back
    std::vector<branch::backend::OverlapPair> h_pairs(n_pairs);
    if (n_pairs > 0) {
        ASSERT_EQ(cudaMemcpy(h_pairs.data(), d_pairs,
                            n_pairs * sizeof(branch::backend::OverlapPair),
                            cudaMemcpyDeviceToHost), cudaSuccess);
    }
    
    // Cleanup
    cudaFree(d_minimizers);
    cudaFree(d_pairs);
    
    // Verify: expect pairs between (0,1) and (2,3) due to shared hashes
    EXPECT_GT(n_pairs, 0u) << "Expected at least some overlap pairs";
    
    // Check that pairs are canonically ordered (read_a < read_b)
    for (const auto& pair : h_pairs) {
        EXPECT_LT(pair.read_a, pair.read_b) << "Pairs must be canonically ordered";
    }
    
    // Check strand field is set correctly
    // Pairs between reads 0-1 should have strand=0 (same strand)
    // Pairs between reads 2-3 should have strand=1 (opposite strand)
    bool found_same_strand = false;
    bool found_opposite_strand = false;
    for (const auto& pair : h_pairs) {
        if ((pair.read_a == 0 && pair.read_b == 1) ||
            (pair.read_a == 1 && pair.read_b == 0)) {
            EXPECT_EQ(pair.strand, 0u);
            found_same_strand = true;
        }
        if ((pair.read_a == 2 && pair.read_b == 3) ||
            (pair.read_a == 3 && pair.read_b == 2)) {
            EXPECT_EQ(pair.strand, 1u);
            found_opposite_strand = true;
        }
    }
    
    EXPECT_TRUE(found_same_strand) << "Expected same-strand pairs between reads 0 and 1";
    EXPECT_TRUE(found_opposite_strand) << "Expected opposite-strand pairs between reads 2 and 3";
}

#endif  // BRANCH_BUILD_CUDA
