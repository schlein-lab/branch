// BRANCH — GPU backend implementation (CUDA build).
//
// This file is compiled ONLY when BRANCH_CUDA_ENABLED=ON — see
// src/backend/CMakeLists.txt. The CPU-only build picks up
// gpu_backend_stub.cpp instead, which returns an empty Backend from
// make_gpu_backend() without pulling in any CUDA symbols.
//
// Responsibilities:
//   - Own a cudaStream_t + device_id in a BackendContext.
//   - compute_overlaps: host-side minimizer sketching (shared with the
//     CPU backend via graph/minimizer_sketcher.hpp), cudaMalloc +
//     cudaMemcpyAsync upload, launch overlap_kernel.cu, cudaMemcpyAsync
//     download, stream sync. Emits OverlapPair records into the
//     caller-provided output span.
//   - classify_batch, estimate_vaf_batch: delegate to cpu_backend. A
//     GPU port of those paths is future work; we keep the same
//     backend instance usable end-to-end so callers don't have to juggle
//     two backend handles.

#include "backend/gpu_backend.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <new>
#include <span>
#include <vector>

#include <cuda_runtime.h>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "classify/features.hpp"
#include "graph/kmer_sketch.hpp"
#include "graph/minimizer_sketcher.hpp"
#include "gpu/overlap_kernel.cuh"

namespace branch::backend {

namespace {

struct GpuBackendContext {
    cudaStream_t stream{nullptr};
    int device_id{0};
    Backend cpu_delegate;  // used for classify + vaf until GPU versions land
};

// ---------- vtable entry points ----------

void gpu_destroy(BackendContext ctx) noexcept {
    auto* c = static_cast<GpuBackendContext*>(ctx);
    if (!c) return;
    if (c->stream) {
        // Ignore errors on teardown — nothing to recover.
        cudaStreamDestroy(c->stream);
    }
    delete c;
}

const char* gpu_name(BackendContext /*ctx*/) noexcept {
    return kGpuBackendName;
}

// Dispatches to the CPU delegate for classification — see header notes.
void gpu_classify_batch(BackendContext ctx,
                        std::span<const BubbleCandidate> candidates,
                        std::span<ClassificationResult> out_results) {
    auto* c = static_cast<GpuBackendContext*>(ctx);
    if (!c || c->cpu_delegate.empty()) {
        for (auto& r : out_results) r = ClassificationResult{};
        return;
    }
    c->cpu_delegate.classify_batch(candidates, out_results);
}

void gpu_estimate_vaf_batch(BackendContext ctx,
                            std::span<const std::uint32_t> branch_ids,
                            std::span<VAFEstimate> out_estimates) {
    auto* c = static_cast<GpuBackendContext*>(ctx);
    if (!c || c->cpu_delegate.empty()) {
        for (auto& v : out_estimates) v = VAFEstimate{0.5f, 0.0f, 1.0f};
        return;
    }
    c->cpu_delegate.estimate_vaf_batch(branch_ids, out_estimates);
}

void gpu_compute_overlaps(BackendContext ctx,
                          const ReadBatch* batch,
                          std::span<OverlapPair> out_pairs,
                          std::size_t* out_count) {
    *out_count = 0;
    auto* c = static_cast<GpuBackendContext*>(ctx);
    if (!c || !batch || batch->reads.empty() || out_pairs.empty()) return;

    // Step 1: Host-side minimizer sketch. MinimizerHit has the same
    // 16-byte layout as gpu::MinimizerEntry (static_assert enforced
    // in kmer_sketch.hpp), so the device upload is a plain memcpy.
    std::vector<branch::graph::MinimizerHit> host_hits;
    // Rough per-read bound: one minimizer per ~w bases at w=19. We
    // reserve generously to avoid reallocation mid-sketch.
    host_hits.reserve(batch->reads.size() * 32);
    for (const auto& r : batch->reads) {
        branch::graph::sketch_read(r.seq, r.id, host_hits);
    }
    if (host_hits.empty()) return;

    // Step 2: Allocate device buffers. We size the output buffer to
    // the caller's out_pairs span so the kernel can never overflow
    // past what the caller promised to read.
    const std::size_t n_hits = host_hits.size();
    const std::size_t max_pairs = out_pairs.size();

    branch::gpu::MinimizerEntry* d_minimizers = nullptr;
    OverlapPair* d_pairs = nullptr;

    if (cudaMallocAsync(&d_minimizers,
                         n_hits * sizeof(branch::gpu::MinimizerEntry),
                         c->stream) != cudaSuccess) {
        std::fprintf(stderr,
                     "[branch gpu-backend] cudaMalloc minimizers failed (%zu hits); falling back to CPU\n",
                     n_hits);
        if (c->cpu_delegate.empty()) return;
        c->cpu_delegate.compute_overlaps(batch, out_pairs, out_count);
        return;
    }
    if (cudaMallocAsync(&d_pairs, max_pairs * sizeof(OverlapPair),
                         c->stream) != cudaSuccess) {
        cudaFreeAsync(d_minimizers, c->stream);
        std::fprintf(stderr,
                     "[branch gpu-backend] cudaMalloc pair buffer failed (%zu pairs); falling back to CPU\n",
                     max_pairs);
        if (c->cpu_delegate.empty()) return;
        c->cpu_delegate.compute_overlaps(batch, out_pairs, out_count);
        return;
    }

    // Zero-copy layout: reinterpret host MinimizerHit bytes as
    // MinimizerEntry before upload. Static assert in kmer_sketch.hpp
    // guarantees the sizes match.
    static_assert(sizeof(branch::graph::MinimizerHit)
                      == sizeof(branch::gpu::MinimizerEntry),
                  "MinimizerHit / MinimizerEntry layout must match");
    cudaMemcpyAsync(
        d_minimizers,
        reinterpret_cast<const branch::gpu::MinimizerEntry*>(host_hits.data()),
        n_hits * sizeof(branch::gpu::MinimizerEntry),
        cudaMemcpyHostToDevice, c->stream);

    // Step 3: Launch kernel.
    branch::gpu::OverlapKernelLaunchConfig kcfg{};
    kcfg.n_minimizers = n_hits;
    kcfg.max_output_pairs = max_pairs;
    kcfg.device_id = c->device_id;
    kcfg.threads_per_block = 256;
    kcfg.stream = c->stream;

    const std::size_t n_pairs =
        branch::gpu::launch_overlap_kernel(d_minimizers, d_pairs, kcfg);

    // Step 4: Download results. n_pairs is clamped inside the kernel
    // to max_output_pairs so the memcpy can never run past the
    // caller's span.
    if (n_pairs > 0) {
        cudaMemcpyAsync(out_pairs.data(), d_pairs,
                        n_pairs * sizeof(OverlapPair),
                        cudaMemcpyDeviceToHost, c->stream);
    }
    cudaStreamSynchronize(c->stream);

    cudaFreeAsync(d_minimizers, c->stream);
    cudaFreeAsync(d_pairs, c->stream);

    *out_count = n_pairs;
}

const BackendVTable kGpuVTable = {
    .destroy             = gpu_destroy,
    .compute_overlaps    = gpu_compute_overlaps,
    .classify_batch      = gpu_classify_batch,
    .estimate_vaf_batch  = gpu_estimate_vaf_batch,
    .name                = gpu_name,
};

}  // namespace

bool gpu_backend_compiled_in() noexcept {
    return true;
}

Backend make_gpu_backend() noexcept {
    // Probe for a usable device at runtime. A CUDA-compiled binary
    // running on a CPU-only node (no driver, no hardware) must fail
    // gracefully rather than aborting.
    int n_devices = 0;
    cudaError_t err = cudaGetDeviceCount(&n_devices);
    if (err != cudaSuccess || n_devices <= 0) {
        return Backend{};
    }

    auto* ctx = new (std::nothrow) GpuBackendContext{};
    if (!ctx) return Backend{};

    ctx->device_id = 0;
    if (cudaSetDevice(ctx->device_id) != cudaSuccess) {
        delete ctx;
        return Backend{};
    }
    if (cudaStreamCreateWithFlags(&ctx->stream, cudaStreamNonBlocking)
            != cudaSuccess) {
        delete ctx;
        return Backend{};
    }
    // Pair with a CPU backend for the non-overlap vtable entries.
    // If this fails we still return a GPU backend that will fall back
    // to no-op on classify / vaf, which is acceptable because the
    // pipeline ultimately goes through cpu_classify_batch_with_graph
    // directly from CLI.
    ctx->cpu_delegate = make_cpu_backend();

    return Backend{static_cast<BackendContext>(ctx), &kGpuVTable};
}

}  // namespace branch::backend
