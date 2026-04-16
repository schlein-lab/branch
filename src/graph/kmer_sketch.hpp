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
