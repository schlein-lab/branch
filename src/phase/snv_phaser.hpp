#pragma once

// BRANCH — SNV-aware read phasing.
//
// Once every read carries its own deltas (positional differences from
// the unitig consensus, see graph::DeltaOp), positions where multiple
// reads disagree with the consensus become candidate variant sites.
// A position is INFORMATIVE when:
//   - at least `min_alt_reads` reads carry the alt allele, AND
//   - at least `min_ref_reads` reads carry the reference allele.
//
// Reads that span two informative SNVs in a consistent (ref/alt) pair
// link those SNVs into the same PHASE BLOCK. Within a phase block we
// can split reads into two haplotypes by majority-allele clustering;
// this gives us per-haplotype read counts that BubbleAltVAF can use to
// distinguish:
//   - germline het (50/50 split, both haps see the same alt count)
//   - somatic mosaic (unbalanced split, one hap dominates)
//
// This module is intentionally minimal: a greedy single-pass phaser
// (HapCUT2-lite). It does NOT solve the optimal-min-error-correction
// problem. Use it as a coarse haplotype assignment to be paired with
// the bubble detector's per-alt read counts.

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "graph/delta_read.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::graph { class ReadDB; }

namespace branch::phase {

struct InformativeSNV {
    branch::graph::NodeId unitig{};
    std::uint32_t pos{};                     // 0-based offset within unitig consensus
    char ref_base{'.'};
    char alt_base{'.'};
    std::vector<branch::graph::ReadId> ref_reads;
    std::vector<branch::graph::ReadId> alt_reads;
};

struct PhaseBlock {
    // Member SNVs in unitig+position order.
    std::vector<InformativeSNV> snvs;
    // For each read that carries at least one allele in this block:
    // 0 = haplotype 0, 1 = haplotype 1, 255 = ambiguous.
    std::unordered_map<branch::graph::ReadId, std::uint8_t> read_haplotype;
};

struct SnvPhaserConfig {
    // Minimum reads carrying the alt allele for a position to be
    // informative. Setting to 2 collapses HiFi singletons (~1 % error
    // rate gives mostly n=1 false positives).
    std::uint32_t min_alt_reads{2};
    // Minimum reads carrying the ref allele.
    std::uint32_t min_ref_reads{2};
    // Two SNVs are linked into the same phase block if at least
    // `min_link_reads` reads carry consistent allele calls at both.
    std::uint32_t min_link_reads{2};
};

// Identify informative SNVs across the whole ReadDB. Returns one
// entry per informative position, sorted by (unitig, pos).
[[nodiscard]] std::vector<InformativeSNV>
find_informative_snvs(const branch::graph::ReadDB& db,
                      const SnvPhaserConfig& cfg = {});

// Greedy single-pass phasing. Walks SNVs in position order, links
// consecutive SNVs that share `min_link_reads` reads with consistent
// alleles into the same PhaseBlock, and assigns each read to a
// haplotype within the block by majority vote across its SNV calls.
// Reads with conflicting calls get haplotype = 255 (ambiguous) and
// are excluded from per-haplotype counts.
[[nodiscard]] std::vector<PhaseBlock>
build_phase_blocks(const std::vector<InformativeSNV>& snvs,
                   const SnvPhaserConfig& cfg = {});

}  // namespace branch::phase
