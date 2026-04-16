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
#include <string_view>
#include <vector>

#include <span>

#include "backend/backend_vtable.hpp"
#include "classify/features.hpp"
#include "graph/delta_read.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::backend {

// Per-instance state owned by the CPU backend. Hidden behind a
// forward declaration so the header stays minimal; the full struct
// is opaque to callers outside of the CPU backend .cpp.
struct CpuBackendContext;

// ReadBatch: container of sequences for overlap computation.
// Each read is a (sequence, read_id) pair.
struct ReadBatch {
    struct Read {
        std::string_view seq;
        graph::ReadId id;
    };
    std::vector<Read> reads;
};

// Construct a CPU-backed Backend. Never returns an empty Backend;
// context is always allocated. Throws std::bad_alloc on OOM.
Backend make_cpu_backend();

// Identifier constant for logging + registry.
inline constexpr const char* kCpuBackendName = "cpu-reference";

// v0.2 convenience: classify a batch of BubbleCandidates with the full
// feature-extraction + cascade pipeline. Unlike the vtable entry
// (Backend::classify_batch), this version takes a graph reference so it
// can run extract_features() per candidate before invoking the cascade.
//
// Use this as the real entry point from pipeline code; the vtable
// variant is retained for GPU-backend parity and assumes features are
// already populated on each candidate.
void cpu_classify_batch_with_graph(
    std::span<const branch::classify::BubbleCandidate> candidates,
    const branch::graph::LosslessGraph& graph,
    std::span<ClassificationResult> out_results);

// v0.2 convenience: Wilson-95%-CI VAF estimator driven directly by
// BubbleCandidate.read_support_{branch,alt}. Emits one VAFEstimate per
// candidate. On total == 0 returns the uninformative prior
// {0.5, 0.0, 1.0}; see vaf_stats.hpp::wilson_ci.
void cpu_estimate_vaf_batch_candidates(
    std::span<const branch::classify::BubbleCandidate> candidates,
    std::span<VAFEstimate> out_estimates);

}  // namespace branch::backend
