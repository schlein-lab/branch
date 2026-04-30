#pragma once

// BRANCH v0.1 — Minimizer-sketch for overlap computation.
//
// A MinimizerSketch is a compact summary of a read used by the overlap
// kernel (CPU or GPU) to find candidate read-pairs without comparing
// the full sequences. The sketch stores minimizer hashes with their
// in-read positions; overlap candidates are read-pairs that share
// ≥ threshold minimizers in compatible orientation.
//
// Parameters are fixed at v0.1 to keep the interface stable:
//   k = 21 (HiFi-appropriate k-mer size)
//   w = 19 (window size; together with k this yields ~one minimizer
//           per 10 bases on random sequence)
//
// The layout is deliberately the same on CPU and GPU (see
// branch::gpu::MinimizerEntry) so the sketch can be uploaded to the
// device without transformation.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "graph/delta_read.hpp"

namespace branch::graph {

inline constexpr std::size_t kMinimizerK = 21;
inline constexpr std::size_t kMinimizerW = 19;

// ---------------------------------------------------------------------------
// MinimizerProfile — read-tech tunable knobs for sketch + overlap.
//
// `kMinimizerK` / `kMinimizerW` constants above stay valid as the HiFi
// defaults (still consumed by GFA/BED writers as cosmetic ref-window-K).
// The actual overlap path now takes a runtime profile so ONT-tuned
// values can be plugged in via the assemble-CLI `--read-tech ont` flag.
//
// Why these defaults:
//   HiFi  — k=21, w=19 from minimap2/HiCanu HiFi presets. ~0.1% error
//           tolerates a long minimizer; bucket cap unnecessary because
//           HiFi minimizers are highly informative.
//   ONT   — k=15, w=10 from minimap2 `map-ont` preset. R10.4.1 SUP at
//           ~1-2% error rate produces many wrong minimizers at k=21,
//           and the resulting bucket explosion was the documented OOM
//           failure mode (peak RSS 26 GiB on 37 Mbp HG002 IGH chunk).
//           bucket_cap=64 drops the high-multiplicity hashes before
//           they detonate the all-pairs match table; min_qual reserved
//           for a future Q-score mask once ReadBatch carries qualities.
struct MinimizerProfile {
    std::size_t k             = 21;   // HiFi default
    std::size_t w             = 19;   // HiFi default
    std::size_t max_bucket_size = 0;  // 0 = no cap; cap drops over-replicated minimizers
    int         min_qual      = 0;    // Q-score floor; 0 = no filter (reserved)
    const char* name          = "hifi";
};

inline constexpr MinimizerProfile kProfileHiFi{21, 19, 0,  0, "hifi"};
inline constexpr MinimizerProfile kProfileONT {15, 10, 64, 0, "ont"};


// ---------------------------------------------------------------------------
// TechProfile (v0.5) — single source for every tech-specific threshold.
//
// Extends MinimizerProfile with the parameters needed by the v0.5
// whole-genome pipeline (Phase 0..5, see docs/v0.5.md). The minimizer
// fields are duplicated rather than inherited so existing call sites
// using MinimizerProfile keep working without any change. New v0.5
// modules should consume TechProfile directly.
struct TechProfile {
    // Phase 0 — minimizer + k-mer histogram
    std::size_t minimizer_k        = 21;
    std::size_t minimizer_w        = 19;
    std::size_t bucket_cap         = 0;      // 0 = no cap (HiFi); 64 (ONT)
    int         expected_coverage  = 30;     // user-overridable per run
    double      bimodality_z       = 4.0;    // peak-separation z-threshold

    // Phase 1 — master-tiling + curation
    double      min_base_agreement = 0.95;
    int         min_coverage_branch = 2;
    int         min_branch_length_bp = 1;    // 1-bp branches as feature

    // Phase 2 — transitive branch attachment
    double      sketch_jaccard_lo  = 0.20;
    int         max_transitive_depth = 3;

    // Phase 3 — vPCR auto-primer design
    int         primer_min_window_bp = 18;
    int         primer_max_window_bp = 25;
    double      primer_max_entropy = 0.5;

    const char* name = "hifi";
};

inline constexpr TechProfile kTechProfileHiFi{
    21, 19,  0, 30, 4.0,
    0.95, 2, 1,
    0.20, 3,
    18, 25, 0.5,
    "hifi"
};

inline constexpr TechProfile kTechProfileONT{
    15, 10, 64, 35, 2.5,
    0.70, 3, 1,
    0.05, 5,
    25, 40, 0.8,
    "ont"
};

// Convenience: derive a MinimizerProfile (for the existing assemble
// backend) from the v0.5 TechProfile.
inline constexpr MinimizerProfile minimizer_view(const TechProfile& t) noexcept {
    return MinimizerProfile{
        t.minimizer_k, t.minimizer_w, t.bucket_cap, 0, t.name
    };
}



// Minimizer hit — one minimizer from one read.
//
// v0.2 adds a strand bit so the overlap backend can tell same-strand
// from opposite-strand matches. The struct is still 16 bytes: `pos`
// shrinks from 32 to 24 bits (HiFi reads fit comfortably in 16 MB)
// to reclaim one byte for the strand.
//
// Layout is kept identical to branch::gpu::MinimizerEntry on purpose,
// so the sketch can be uploaded to a GPU without transformation.
struct MinimizerHit {
    std::uint64_t hash;       // canonical 2-bit-packed k-mer hash
    std::uint32_t read_id;
    std::uint32_t pos : 24;   // 0-based position in the read (<= 16 MB)
    std::uint32_t strand : 8; // 0 = forward canonical, 1 = reverse canonical
};

static_assert(sizeof(MinimizerHit) == 16,
              "MinimizerHit must pack to 16 bytes; identical layout to "
              "branch::gpu::MinimizerEntry for zero-copy GPU upload");

// A sketch of a single read. Typically populated by a Sketcher (added
// in a future module); here we only define the storage type so the
// overlap backend can iterate over it uniformly on CPU and GPU.
struct ReadSketch {
    ReadId read_id{};
    std::uint32_t read_length{};
    std::vector<MinimizerHit> hits;
};

// Non-owning view over a contiguous block of minimizer hits. Used in
// hot paths to pass around device-resident memory without template
// parameterisation.
struct MinimizerHitsView {
    std::span<const MinimizerHit> hits;
};

}  // namespace branch::graph
