# Backend + GPU Review (Part 3/4)

## Files
| File | LOC | Purpose |
|------|-----|----------|
| `src/backend/backend_vtable.hpp` | 143 | Type-erased VTable, OverlapPair 24B, RAII Backend class |
| `src/gpu/overlap_kernel.cuh` | 56 | Host-side CUDA decl, MinimizerEntry 16B, CPU-compilable |
| `src/gpu/overlap_kernel.cu` | 44 | CUDA kernel impl, skeleton no-op, launch wrapper |
| `src/gpu/overlap_kernel_stub.cpp` | 17 | CPU fallback, returns 0, links without CUDA |

## Critical Issues
1. **Per-call cudaMalloc** — `overlap_kernel.cu:24-34` allocates/frees `d_count` every call. For batch processing this kills perf. Use persistent memory pool or pre-allocated workspace.
2. **No CUDA error checking** — Missing `cudaGetLastError()` after kernel launch. Silent failures will corrupt results.
3. **No stream parameter** — Kernel launches on default stream (0). Blocks multi-GPU concurrency on Hummel-2 H100s.
4. **Hardcoded grid/block** — `launch_overlap_kernel` ignores `config.threads_per_block`. Uses 256 unconditionally.

## Design Feedback
- **VTable pattern: GOOD** — Type-erasure allows CPU↔GPU swap without recompile. Move-only semantics correct.
- **Alignment: CORRECT** — `OverlapPair` 24B, `MinimizerEntry` 16B fit shared mem constraints (3072 entries/48KB).
- **Header separation: CORRECT** — `.cuh` has no `__global__`, compiles with g++ for CPU builds.
- **Stub fallback: FUNCTIONAL** — Returns 0, prevents link errors when `BRANCH_BUILD_CUDA=OFF`.
- **Fat-binary ready: PARTIAL** — sm_80/sm_90 in CMake, but no runtime device selection in kernel.

## TODOs
- [ ] Add `CUDA_CHECK()` macro wrapping all CUDA calls
- [ ] Implement memory pool in `GpuBackend` class for `d_count` workspace
- [ ] Add `cudaStream_t` parameter to `launch_overlap_kernel`
- [ ] Use `config.threads_per_block` instead of hardcoded 256
- [ ] Add `cudaGetDeviceCount`/`cudaSetDevice` for multi-GPU dispatch
- [ ] Implement actual overlap logic in kernel (currently no-op)

---
*Review: 2025-01-XX | Scope: Backend abstraction, CUDA correctness, Fat-Binary, Stub*
