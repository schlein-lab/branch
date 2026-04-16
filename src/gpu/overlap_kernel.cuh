#pragma once

// BRANCH v0.1 - CUDA Overlap Kernel (skeleton).
//
// Per CUDA expert consultation: v0.1 GPU-Scope is the All-vs-All
// Minimizer Overlap kernel only. Banded alignment and beam-search
// classification stay on CPU until profiling justifies GPU for them.
//
// This is a skeleton: kernel body is intentionally a no-op stub so
// the build system and host-side dispatch can be wired end-to-end
// before the real algorithm lands in v0.2. Target archs are sm_80
// (A100) and sm_90 (H100) via the CMake fat-binary config.

#include <cstddef>
#include <cstdint>

#include "backend/backend_vtable.hpp"

namespace branch::gpu {

// Device-side minimizer sketch entry. 16 bytes, cache-line friendly
// when used in shared memory (48 KB per SM = 3072 entries).
struct MinimizerEntry {
    std::uint64_t hash;
    std::uint32_t read_id;
    std::uint32_t pos;
};

static_assert(sizeof(MinimizerEntry) == 16,
              "MinimizerEntry must pack to 16 bytes for shared-mem tiling");

// Host-side launcher. Takes device-resident buffers, launches the
// kernel, returns number of pairs written. For v0.1 this stub writes
// zero pairs and returns 0; the real implementation arrives in v0.2.
//
// The host/device split is deliberate: this header is CPU-compilable
// (no __global__ here) so the broader project can include it even
// when BRANCH_BUILD_CUDA is OFF.
struct OverlapKernelLaunchConfig {
    std::size_t n_minimizers;
    std::size_t max_output_pairs;
    int device_id{0};
    int threads_per_block{256};
};

// Opaque device pointers passed in by a Backend implementation that
// owns the CUDA resources (memory pool, stream, etc.).
using DeviceMinimizerPtr = const MinimizerEntry*;
using DeviceOverlapPairPtr = branch::backend::OverlapPair*;

// CPU-callable entry point. On non-CUDA builds, this symbol is
// defined by a stub .cpp that returns 0 and sets out_count to 0.
std::size_t launch_overlap_kernel(
    DeviceMinimizerPtr d_minimizers,
    DeviceOverlapPairPtr d_output_pairs,
    const OverlapKernelLaunchConfig& cfg);

}  // namespace branch::gpu
