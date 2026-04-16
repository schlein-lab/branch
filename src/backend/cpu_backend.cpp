// BRANCH v0.1 — CPU reference backend implementation.

#include "backend/cpu_backend.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace branch::backend {

struct CpuBackendContext {
    // v0.1: no state. Future: thread-pool handle, scratch buffers, etc.
};

namespace {

void cpu_destroy(BackendContext ctx) {
    delete static_cast<CpuBackendContext*>(ctx);
}

void cpu_compute_overlaps(BackendContext /*ctx*/,
                          const ReadBatch* /*batch*/,
                          std::span<OverlapPair> /*out_pairs*/,
                          std::size_t* out_count) {
    // v0.1: stub. Real minimizer-based overlap lands in v0.2.
    *out_count = 0;
}

void cpu_classify_batch(BackendContext /*ctx*/,
                        std::span<const BubbleCandidate> candidates,
                        std::span<ClassificationResult> out_results) {
    // v0.1: stub. Real cascade in src/classify/cascade.hpp; a future
    // revision wires this to call classify_one() per candidate so the
    // VTable dispatch goes through the public API.
    const std::size_t n = std::min(candidates.size(), out_results.size());
    for (std::size_t i = 0; i < n; ++i) {
        out_results[i].label = BubbleClass::Unknown;
        out_results[i].confidence = 0.0f;
    }
}

void cpu_estimate_vaf_batch(BackendContext /*ctx*/,
                            std::span<const std::uint32_t> branch_ids,
                            std::span<VAFEstimate> out_estimates) {
    const std::size_t n = std::min(branch_ids.size(), out_estimates.size());
    for (std::size_t i = 0; i < n; ++i) {
        out_estimates[i].point = 0.5f;
        out_estimates[i].ci_low = 0.0f;
        out_estimates[i].ci_high = 1.0f;
    }
}

const char* cpu_name(BackendContext /*ctx*/) {
    return kCpuBackendName;
}

constexpr BackendVTable kCpuVTable{
    .destroy = &cpu_destroy,
    .compute_overlaps = &cpu_compute_overlaps,
    .classify_batch = &cpu_classify_batch,
    .estimate_vaf_batch = &cpu_estimate_vaf_batch,
    .name = &cpu_name,
};

}  // namespace

Backend make_cpu_backend() {
    auto* ctx = new CpuBackendContext{};
    return Backend(ctx, &kCpuVTable);
}

}  // namespace branch::backend
