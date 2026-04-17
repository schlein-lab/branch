# Backend + GPU Review (Part 3/4)

Status audit: 2026-04-17 — markers added per finding.

## Files
| File | LOC | Purpose |
|------|-----|----------|
| `src/backend/backend_vtable.hpp` | 143 | Type-erased VTable, OverlapPair 24B, RAII Backend class |
| `src/gpu/overlap_kernel.cuh` | 56 | Host-side CUDA decl, MinimizerEntry 16B, CPU-compilable |
| `src/gpu/overlap_kernel.cu` | 44 | CUDA kernel impl, skeleton no-op, launch wrapper |
| `src/gpu/overlap_kernel_stub.cpp` | 17 | CPU fallback, returns 0, links without CUDA |

## Critical Issues
1. ⬜ Open — **Per-call cudaMalloc** — `overlap_kernel.cu:24-34` allocates/frees `d_count` every call. For batch processing this kills perf. Use persistent memory pool or pre-allocated workspace.
   - Why open: `launch_overlap_kernel` allokiert weiterhin pro Call (overlap_kernel.cu L192-195, L209 `cudaMalloc`/`cudaFree`). Kein GPU-Backend-owned Memory-Pool.

2. ✅ Resolved — **No CUDA error checking** — Missing `cudaGetLastError()` after kernel launch.
   - Why resolved: `CUDA_CHECK` Macro um alle CUDA-Calls (overlap_kernel.cu L18-26), `cudaGetLastError()` nach `find_bucket_boundaries` (L227) und `process_buckets_emit_pairs` (L258).

3. ✅ Resolved — **No stream parameter** — Kernel launches on default stream (0).
   - Why resolved: `OverlapKernelLaunchConfig::stream{nullptr}` (overlap_kernel.cuh L59), stream durchgereicht an Kernel-Launch (overlap_kernel.cu L182, L224, L251) + `cudaMemsetAsync`/`cudaMemcpyAsync`.

4. ✅ Resolved — **Hardcoded grid/block** — `launch_overlap_kernel` ignores `config.threads_per_block`.
   - Why resolved: beide Kernel-Launches nutzen `cfg.threads_per_block` (overlap_kernel.cu L220, L249).

## Design Feedback
- **VTable pattern: GOOD** — Type-erasure allows CPU↔GPU swap without recompile. Move-only semantics correct.
- **Alignment: CORRECT** — `OverlapPair` 24B, `MinimizerEntry` 16B fit shared mem constraints (3072 entries/48KB).
- **Header separation: CORRECT** — `.cuh` has no `__global__`, compiles with g++ for CPU builds.
- **Stub fallback: FUNCTIONAL** — Returns 0, prevents link errors when `BRANCH_BUILD_CUDA=OFF`.
- **Fat-binary ready: PARTIAL** — sm_80/sm_90 in CMake, but no runtime device selection in kernel.

## TODOs
- [x] ✅ Resolved: Add `CUDA_CHECK()` macro wrapping all CUDA calls (overlap_kernel.cu L18-26)
- [ ] ⬜ Open: Implement memory pool in `GpuBackend` class for `d_count` workspace
- [x] ✅ Resolved: Add `cudaStream_t` parameter to `launch_overlap_kernel` (overlap_kernel.cuh L59)
- [x] ✅ Resolved: Use `config.threads_per_block` instead of hardcoded 256 (overlap_kernel.cu L220, L249)
- [ ] ⬜ Open: Add `cudaGetDeviceCount`/`cudaSetDevice` for multi-GPU dispatch
- [x] ✅ Resolved: Implement actual overlap logic in kernel (currently no-op) — Sort-by-hash via CUB + bucket-enumeration + pair-emission jetzt implementiert (overlap_kernel.cu L58-169). Offset-clustering + strand-majority filter auf Host verschoben per Kommentar L150-152.

---
*Review: 2025-01-XX | Scope: Backend abstraction, CUDA correctness, Fat-Binary, Stub*
