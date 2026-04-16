#pragma once

// BRANCH v0.1 — Minimizer sketcher.
//
// Produces a ReadSketch from a DNA sequence. The sketcher iterates a
// rolling k-mer window, computes canonical hashes, and emits the
// minimizer (smallest hash) per w-base window. Results are appended
// to an output vector so callers can sketch many reads in sequence
// without re-allocating.
//
// Parameters are fixed at v0.1 (k=21, w=19 — see kmer_sketch.hpp) so
// the GPU overlap kernel can rely on a stable layout.

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

#include "graph/kmer_sketch.hpp"

namespace branch::graph {

// Encodes one ASCII base to a 2-bit value in {0,1,2,3} for A/C/G/T,
// or 4 for anything else (N, lowercase, etc. — the sketcher treats
// any non-ACGT base as a window break).
[[nodiscard]] inline std::uint8_t encode_base(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// 64-bit splitmix hash. Low-bias, independent of the seed, so the
// same k-mer produces the same hash across runs.
[[nodiscard]] inline std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Reverse-complement a 2-bit-packed k-mer.
[[nodiscard]] inline std::uint64_t rc_kmer(std::uint64_t kmer, std::size_t k) noexcept {
    std::uint64_t rc = 0;
    for (std::size_t i = 0; i < k; ++i) {
        rc = (rc << 2) | ((~kmer) & 3ULL);
        kmer >>= 2;
    }
    return rc;
}

// Canonical k-mer hash: min(hash(kmer), hash(rc(kmer))) so forward
// and reverse-complement strands collide in the sketch.
[[nodiscard]] inline std::uint64_t canonical_hash(std::uint64_t kmer, std::size_t k) noexcept {
    auto hf = splitmix64(kmer);
    auto hr = splitmix64(rc_kmer(kmer, k));
    return hf < hr ? hf : hr;
}

// Append the minimizer sketch of `seq` to `out`. Resets and appends
// one MinimizerHit per window. `read_id` is stamped onto each hit so
// the hits can be merged into a single per-batch table later.
//
// Ambiguous bases (N, lowercase, anything non-ACGT) break the window
// and restart hashing on the next valid base — matching the behaviour
// of most production HiFi overlap tools.
void sketch_read(std::string_view seq,
                 ReadId read_id,
                 std::vector<MinimizerHit>& out) noexcept;

}  // namespace branch::graph
