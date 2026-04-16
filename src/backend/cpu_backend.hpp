#pragma once

// BRANCH v0.1 — CPU reference backend.
//
// Concrete VTable implementation that runs every backend operation on
// the host CPU. Intended as:
//   - the default fallback when BRANCH_BUILD_CUDA is OFF,
//   - the ground-truth reference against which GPU kernels are
//     regression-tested (same inputs must yield the same outputs),
//   - the developmental testbed for new algorithms before they are
//     ported to CUDA.
//
// v0.1 ships with stubs (zero overlaps found, constant VAF estimates).
// Real implementations follow in v0.2 once graph-construction and
// overlap-computation modules exist.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

#include "backend/backend_vtable.hpp"

namespace branch::backend {

// Per-instance state owned by the CPU backend. Hidden behind a
// forward declaration so the header stays minimal; the full struct
// is opaque to callers outside of the CPU backend .cpp.
struct CpuBackendContext;

// Construct a CPU-backed Backend. Never returns an empty Backend;
// context is always allocated. Throws std::bad_alloc on OOM.
Backend make_cpu_backend();

// Identifier constant for logging + registry.
inline constexpr const char* kCpuBackendName = "cpu-reference";

}  // namespace branch::backend
