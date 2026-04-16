#pragma once

// BRANCH v0.1 — Backend VTable (Type-Erased).
//
// Per C++ architect consultation: classical virtual-inheritance interfaces
// per concern (IOverlapComputer / IAligner / IClassifier / IVAFEstimator)
// are collapsed into ONE physical VTable with ~8–10 function pointers
// operating on BATCH APIs. This removes per-element virtual dispatch
// overhead on hot paths and enforces GPU-friendly batching at the
// interface level.
//
// Explicitly NOT a virtual-inheritance class hierarchy. Do not add an
// IBackend with virtual methods — that design was rejected upstream.

#include <cstddef>
#include <cstdint>
#include <span>

// BubbleCandidate lives in classify/features.hpp. Included here (before
// any namespace) because classify_batch takes a span of values.
#include "classify/features.hpp"

namespace branch::backend {

using BackendContext = void*;
using BubbleCandidate = branch::classify::BubbleCandidate;

// Forward-declared because only a pointer is passed.
struct ReadBatch;

// Output of classification.
enum class BubbleClass : std::uint8_t {
    Branch = 0,
    Duplication = 1,
    Mixed = 2,
    NonSeparable = 3,
    Unknown = 255,
};

struct ClassificationResult {
    BubbleClass label;
    float confidence;   // [0.0, 1.0]
    // Cascade stage index that produced the label (0..3), or 255 if
    // not emitted by the cascade (e.g. legacy stub paths). Populated
    // by CPU backend from classify::StageResult.stage_index.
    std::uint8_t stage{255};
};

// VAF estimate with a simple confidence interval.
struct VAFEstimate {
    float point;    // [0.0, 1.0]
    float ci_low;   // [0.0, 1.0]
    float ci_high;  // [0.0, 1.0]
};

// Overlap pair: two reads overlap starting at given offsets.
//
// v0.2: strand distinguishes same-strand from opposite-strand overlaps.
// 0 = reads A and B overlap on the same strand (minimizer hits agreed
//     on canonical orientation).
// 1 = reads overlap on opposite strands (one side is reverse-complement
//     relative to the other).
struct OverlapPair {
    std::uint32_t read_a;
    std::uint32_t read_b;
    std::uint32_t offset_a;
    std::uint32_t offset_b;
    std::uint32_t overlap_len;
    std::int16_t diff_count;
    std::uint8_t strand;  // 0 = same-strand, 1 = opposite-strand
    std::uint8_t _pad;    // explicit padding for 24-byte alignment
};

static_assert(sizeof(OverlapPair) == 24,
              "OverlapPair must pack to 24 bytes for SIMD-friendly batches");

// The VTable. All backend entry points are batch-oriented; single-element
// dispatch is never exposed.
struct BackendVTable {
    // Free the backend context (called by Backend's destructor).
    void (*destroy)(BackendContext ctx);

    // Compute candidate overlap pairs for a read batch.
    void (*compute_overlaps)(BackendContext ctx,
                             const ReadBatch* batch,
                             std::span<OverlapPair> out_pairs,
                             std::size_t* out_count);

    // Classify a batch of bubble candidates.
    void (*classify_batch)(BackendContext ctx,
                           std::span<const BubbleCandidate> candidates,
                           std::span<ClassificationResult> out_results);

    // Estimate VAF for a batch of branches.
    void (*estimate_vaf_batch)(BackendContext ctx,
                               std::span<const std::uint32_t> branch_ids,
                               std::span<VAFEstimate> out_estimates);

    // Backend self-identifier string (static lifetime).
    const char* (*name)(BackendContext ctx);
};

// Value-type handle. Owns ctx_; vtable_ is non-owning and has static
// storage duration (one VTable per backend kind).
class Backend {
public:
    Backend() = default;
    Backend(BackendContext ctx, const BackendVTable* vtable) noexcept
        : ctx_(ctx), vtable_(vtable) {}

    Backend(const Backend&) = delete;
    Backend& operator=(const Backend&) = delete;

    Backend(Backend&& other) noexcept
        : ctx_(other.ctx_), vtable_(other.vtable_) {
        other.ctx_ = nullptr;
        other.vtable_ = nullptr;
    }
    Backend& operator=(Backend&& other) noexcept {
        if (this != &other) {
            release();
            ctx_ = other.ctx_;
            vtable_ = other.vtable_;
            other.ctx_ = nullptr;
            other.vtable_ = nullptr;
        }
        return *this;
    }
    ~Backend() { release(); }

    bool empty() const noexcept { return ctx_ == nullptr || vtable_ == nullptr; }

    void compute_overlaps(const ReadBatch* batch,
                          std::span<OverlapPair> out_pairs,
                          std::size_t* out_count) {
        vtable_->compute_overlaps(ctx_, batch, out_pairs, out_count);
    }

    void classify_batch(std::span<const BubbleCandidate> candidates,
                        std::span<ClassificationResult> out_results) {
        vtable_->classify_batch(ctx_, candidates, out_results);
    }

    void estimate_vaf_batch(std::span<const std::uint32_t> branch_ids,
                            std::span<VAFEstimate> out_estimates) {
        vtable_->estimate_vaf_batch(ctx_, branch_ids, out_estimates);
    }

    const char* name() const { return vtable_->name(ctx_); }

private:
    void release() noexcept {
        if (vtable_ && ctx_) {
            vtable_->destroy(ctx_);
        }
        ctx_ = nullptr;
        vtable_ = nullptr;
    }

    BackendContext ctx_{nullptr};
    const BackendVTable* vtable_{nullptr};
};

}  // namespace branch::backend
