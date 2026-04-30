// BRANCH v0.5 — Phase 1.3: master-read tiling.
//
// Given an overlap graph (existing branch::graph::LosslessGraph), find
// a sparsest tiling path with maximum cumulative overlap quality. Each
// position in the haplotype output belongs to exactly one master read;
// adjacent master tiles meet at a midpoint cutover with no synthetic
// consensus.

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"

namespace branch::wg {

struct MasterTile {
    std::uint32_t master_read_id;   // ReadId of this tile's source read
    std::uint32_t start_bp;          // start position in the haplotype output
    std::uint32_t length_bp;
    std::uint32_t q_mean_x10;        // mean Q-score × 10 (for #vph emit)
};

class MasterTiler {
public:
    explicit MasterTiler(const ::branch::graph::TechProfile& profile);

    /// Build the tiling from a per-haplotype overlap-graph (typically
    /// the GFA output of branch assemble run on hap1_reads.fq).
    /// `gfa_path` is the input; tiles are populated in `out_tiles`.
    /// Returns total bp covered.
    std::uint64_t build_tiling(const std::string& gfa_path,
                               std::vector<MasterTile>& out_tiles) const;

    /// Emit master-tile sequence as a FASTA-equivalent string (just the
    /// concatenation of source-read sequences with midpoint cutover).
    [[nodiscard]] std::string render_haplotype_seq(
        const std::vector<MasterTile>& tiles,
        const std::vector<std::pair<std::string, std::string>>& reads_by_id) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
