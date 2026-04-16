// BRANCH v0.1 - CUDA-compiled implementation of the overlap kernel launcher.
//
// Only compiled when BRANCH_BUILD_CUDA is ON. The kernel body is a skeleton
// (zero pairs written). Real algorithm lands in v0.2.

#include "gpu/overlap_kernel.cuh"

#include <cuda_runtime.h>

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
    std::size_t* d_count = nullptr;
    cudaMalloc(&d_count, sizeof(std::size_t));
    cudaMemset(d_count, 0, sizeof(std::size_t));

    int blocks = static_cast<int>(
        (cfg.n_minimizers + cfg.threads_per_block - 1) / cfg.threads_per_block);
    if (blocks < 1) blocks = 1;

    overlap_kernel_noop<<<blocks, cfg.threads_per_block>>>(
        d_minimizers, d_output_pairs,
        cfg.n_minimizers, cfg.max_output_pairs, d_count);

    std::size_t host_count = 0;
    cudaMemcpy(&host_count, d_count, sizeof(std::size_t), cudaMemcpyDeviceToHost);
    cudaFree(d_count);
    return host_count;
}

}  // namespace branch::gpu
