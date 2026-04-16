// BRANCH v0.1 - CPU stub implementation of the GPU overlap kernel launcher.
//
// Built when BRANCH_BUILD_CUDA is OFF. Allows non-CUDA builds to link a
// Backend that references branch::gpu::launch_overlap_kernel without
// forcing the entire GPU toolchain.

#include "gpu/overlap_kernel.cuh"

namespace branch::gpu {

std::size_t launch_overlap_kernel(
    DeviceMinimizerPtr /*d_minimizers*/,
    DeviceOverlapPairPtr /*d_output_pairs*/,
    const OverlapKernelLaunchConfig& /*cfg*/) {
    // No-op on CPU-only builds: no GPU device, no overlap computation here.
    // Note: cfg.stream is ignored in stub - signature matches header for ABI compatibility.
    return 0;
}

}  // namespace branch::gpu
