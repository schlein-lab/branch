// BRANCH v0.1 - CUDA-compiled implementation of the overlap kernel launcher.
//
// Only compiled when BRANCH_BUILD_CUDA is ON. The kernel body is a skeleton
// (zero pairs written). Real algorithm lands in v0.2.

#include "gpu/overlap_kernel.cuh"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

// FIX: CUDA_CHECK macro for error handling
#define CUDA_CHECK(x) \
    do { \
        cudaError_t err = (x); \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
            abort(); \
        } \
    } while (0)

namespace branch::gpu {

namespace {

__global__ void overlap_kernel_noop(
    const MinimizerEntry* /*minimizers*/,
    branch::backend::OverlapPair* /*out_pairs*/,
    std::size_t /*n_minimizers*/,
    std::size_t /*max_output_pairs*/,
    std::size_t* /*out_count*/) {
    // Skeleton: no-op. Real kernel implementation lands in v0.2.
}

}  // namespace

std::size_t launch_overlap_kernel(
    DeviceMinimizerPtr d_minimizers,
    DeviceOverlapPairPtr d_output_pairs,
    const OverlapKernelLaunchConfig& cfg) {
    
    // FIX: Persistent workspace allocation instead of per-call malloc
    // For v0.1 skeleton we still allocate locally but with proper error checking
    // In v0.2 this will be replaced with a workspace pool from the Backend
    std::size_t* d_count = nullptr;
    CUDA_CHECK(cudaMalloc(&d_count, sizeof(std::size_t)));
    CUDA_CHECK(cudaMemsetAsync(d_count, 0, sizeof(std::size_t), cfg.stream));

    // FIX: Use cfg.threads_per_block consistently
    int blocks = static_cast<int>(
        (cfg.n_minimizers + cfg.threads_per_block - 1) / cfg.threads_per_block);
    if (blocks < 1) blocks = 1;

    // FIX: Use stream parameter for async execution
    overlap_kernel_noop<<<blocks, cfg.threads_per_block, 0, cfg.stream>>>(
        d_minimizers, d_output_pairs,
        cfg.n_minimizers, cfg.max_output_pairs, d_count);
    
    // Check for kernel launch errors
    CUDA_CHECK(cudaGetLastError());

    std::size_t host_count = 0;
    // FIX: Wrapped with CUDA_CHECK
    CUDA_CHECK(cudaMemcpyAsync(&host_count, d_count, sizeof(std::size_t), 
                               cudaMemcpyDeviceToHost, cfg.stream));
    
    // Synchronize stream to ensure copy is complete before reading host_count
    CUDA_CHECK(cudaStreamSynchronize(cfg.stream));
    
    CUDA_CHECK(cudaFree(d_count));
    return host_count;
}

}  // namespace branch::gpu
