# Backend + GPU Review

**Reviewer**: service-zyrkel
**Datum**: 2026-04-16
**Scope**: Backend VTable + GPU Overlap Kernel Layer

## Files (4 Dateien, LOC)

| File | LOC | Zweck |
|------|-----|-------|
| backend_vtable.hpp | 154 | Type-erased Backend mit Batch-APIs, RAII-Ownership |
| overlap_kernel.cuh | 59 | CPU-kompilierbarer CUDA-Header, MinimizerEntry 16B |
| overlap_kernel.cu | 48 | CUDA Kernel-Stub (no-op), synchroner Launch |
| overlap_kernel_stub.cpp | 20 | CPU-Fallback, gibt 0 zurueck |

## Critical Issues

| Datei:Zeile | Severity | Problem |
|-------------|----------|---------|
| overlap_kernel.cu:30-31 | **HIGH** | `cudaMalloc`/`cudaFree` pro Launch — synchrone Allokation blockiert GPU-Pipeline. Memory-Pool erforderlich. |
| overlap_kernel.cu:30-43 | **HIGH** | Keine CUDA-Error-Checks. `cudaMalloc` kann OOM returnen, wird ignoriert. |
| overlap_kernel.cu:42 | **MEDIUM** | `cudaMemcpy` D2H ist synchron — blockiert bis Kernel fertig. Stream-async fehlt. |
| overlap_kernel.cuh:39-44 | **LOW** | `OverlapKernelLaunchConfig` hat keinen `cudaStream_t` Parameter — Multi-Stream unmoeglich. |

## Design Feedback

**Positiv:**
- VTable-Ansatz vermeidet virtual-dispatch auf Hot-Paths — korrekt
- `OverlapPair` 24B SIMD-aligned, `MinimizerEntry` 16B Shared-Mem-freundlich
- Stub-Fallback erlaubt CPU-only-Builds ohne Link-Fehler
- Header CPU-kompilierbar (kein `__global__` im .cuh)

**Fehlend fuer Hummel-2 Multi-GPU:**
- Kein Device-Selection (`cudaSetDevice` fehlt trotz `device_id` in Config)
- Kein Multi-Stream-Support fuer Overlap von Compute+Transfer
- Fat-Binary sm_80+sm_90 nur in CMake definiert, nicht validiert

## TODOs

- [ ] **P0**: CUDA-Error-Macro `CUDA_CHECK(...)` mit `cudaGetLastError()` einfuehren
- [ ] **P0**: Memory-Pool (`cudaMallocAsync`/`cudaFreeAsync` oder vorallokierter Pool)
- [ ] **P1**: `cudaStream_t` in `OverlapKernelLaunchConfig` + async Memcpy
- [ ] **P1**: `cudaSetDevice(cfg.device_id)` vor Launch
- [ ] **P2**: CI-Test mit `nvcc -arch=sm_80,sm_90` auf Hummel-2 GPU-Node

---
*Review abgeschlossen. 4 Files, 281 LOC total. 2 HIGH, 1 MEDIUM, 1 LOW.*
