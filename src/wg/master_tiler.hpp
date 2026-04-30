// BRANCH v0.5 — Phase 1.3: master-read tiling.
//
// For one haplotype's read set, build a sparse master tiling: pick real
// reads as tiles, place them sequentially with a midpoint cutover where
// they overlap. Every output base is a real read base; no synthetic
// consensus.
//
// First-cut algorithm — minimizer-shared-position greedy:
//   1. minimize every read at (k=21/15, w=19/10) per profile
//   2. seed = longest read; emit as tile 0 starting at offset 0
//   3. extend: among unused reads, pick the one sharing the most
//      minimizers with the current frontier whose match-positions
//      suggest "downstream of frontier"; place it with overlap at the
//      midpoint of the shared minimizer range
//   4. repeat until no candidate has min_share minimizers in common;
//      then start a new contig from the longest remaining read

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"

namespace branch::wg {

struct MasterTile {
    std::uint32_t master_read_idx;   // index into the per-hap reads vector
    std::uint32_t start_bp;           // global offset in the haplotype output
    std::uint32_t length_bp;          // bp contributed by this tile
    std::uint32_t source_offset_bp;   // first bp of source read used (≥0; midpoint cutover)
    std::uint32_t q_mean_x10;
};

class MasterTiler {
public:
    explicit MasterTiler(const ::branch::graph::TechProfile& profile);

    /// Build master tiles from a per-haplotype reads vector. `reads_by_idx`
    /// is the (id, seq) list for this haplotype. Tiles are emitted in
    /// `out_tiles`; each tile points to one element of `reads_by_idx`.
    /// Returns total bp covered by tiles.
    std::uint64_t build_tiling(
        const std::vector<std::pair<std::string, std::string>>& reads_by_idx,
        std::vector<MasterTile>& out_tiles) const;

    /// Render the haplotype sequence as the concatenation of tile-source
    /// substrings (start = source_offset_bp for length length_bp).
    [[nodiscard]] std::string render_haplotype_seq(
        const std::vector<MasterTile>& tiles,
        const std::vector<std::pair<std::string, std::string>>& reads_by_idx) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
