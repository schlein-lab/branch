#pragma once
//
// refiner — EM-style read-remap loop.
//
// After the build pass tags each read with its best-guess haplotype
// using the het-pair vote, the refiner re-checks every read against the
// current consensus haplotypes and bubble paths. Reads whose alternative
// home scores significantly higher are moved. Iterates until <0.1% of
// reads change in a pass.
//
// This catches early-iteration mistakes where a read was pinned to
// haplotype 1 because hap2 was sparsely seeded at that point — once
// hap2 grows, the read may align better there.

#include "types.hpp"
#include "fastq_index.hpp"
#include <cstdint>
#include <utility>
#include <vector>

namespace branch::wg::phaser {

struct RefinerOpts {
    double min_score_margin = 0.05;  // remap only if alt > current + margin
    double change_floor     = 0.001; // <0.1% reads moving = converged
    int    max_iters        = 10;
};

struct RefinerProgress {
    std::size_t reads_remapped     = 0;
    std::size_t reads_total        = 0;
    int         iter_n             = 0;
    bool        converged          = false;
    // Per-iteration count of reads that moved. Length == iter_n.
    std::vector<std::size_t> per_iter_changes;
};

class Refiner {
public:
    // Iterates until convergence or max_iters. Updates assignments and
    // unitig::member_reads in place. Does not modify unitig sequences;
    // a separate consensus rebuild step is the caller's responsibility
    // (and only needed when remapped fraction > some threshold).
    RefinerProgress refine(
        std::vector<UnitigNode>& unitigs,
        std::vector<Bubble>& bubbles,
        std::vector<ReadAssignment>& assignments,
        FastqIndex& reads,
        const RefinerOpts& opts = {});
};

}  // namespace branch::wg::phaser
