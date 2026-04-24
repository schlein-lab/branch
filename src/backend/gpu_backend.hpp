#pragma once

// BRANCH — GPU backend.
//
// Returns a CUDA-backed Backend when:
//   (a) BRANCH was compiled with BRANCH_CUDA_ENABLED=ON (i.e. a CUDA
//       toolchain was available at build time), AND
//   (b) at least one CUDA-capable device is visible at runtime.
//
// When either condition fails the returned Backend is empty (see
// Backend::empty()) and callers should fall back to make_cpu_backend().
// Never throws — GPU errors at init time turn into "empty backend".
//
// Scope of the GPU backend today:
//   - compute_overlaps runs on the GPU (host-side minimizer sketch,
//     device-side all-vs-all bucket matching via overlap_kernel.cu,
//     zero-copy MinimizerHit ↔ MinimizerEntry layout).
//   - classify_batch and estimate_vaf_batch delegate to the CPU
//     implementations. A GPU port of those paths is follow-up work;
//     what matters for scaling to tumor WGS is the O(N²) overlap step.
//
// Separate CPU→GPU regression testing goes through test_gpu_backend.cpp
// (GTEST_SKIP when no device) and a smoke test in the end-to-end
// assemble pipeline.

#include "backend/backend_vtable.hpp"

namespace branch::backend {

// Identifier constants for logging + registry.
inline constexpr const char* kGpuBackendName = "cuda-overlap";
inline constexpr const char* kGpuBackendUnavailableName = "cuda-unavailable";

// Make a GPU backend. Returns an empty Backend (see Backend::empty())
// when CUDA is unavailable — either at build time (CPU-only build) or
// at runtime (no visible CUDA device, no driver, etc.). Never throws.
Backend make_gpu_backend() noexcept;

// True if the binary was compiled with CUDA support. Runtime device
// availability is a separate probe performed inside make_gpu_backend.
bool gpu_backend_compiled_in() noexcept;

}  // namespace branch::backend
