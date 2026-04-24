// BRANCH — GPU backend CPU-only stub.
//
// Compiled when BRANCH_CUDA_ENABLED=OFF — i.e. on builds without a CUDA
// toolchain. Provides the same symbols as gpu_backend.cpp but without
// pulling any CUDA runtime code, so callers get a well-defined
// "GPU not available" answer instead of a link error.

#include "backend/gpu_backend.hpp"

namespace branch::backend {

bool gpu_backend_compiled_in() noexcept {
    return false;
}

Backend make_gpu_backend() noexcept {
    // Empty Backend: Backend::empty() returns true, callers should
    // fall back to make_cpu_backend().
    return Backend{};
}

}  // namespace branch::backend
