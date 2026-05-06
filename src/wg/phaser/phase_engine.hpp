#pragma once
//
// phase_engine — assigns ReadTag {H1_ONLY, H2_ONLY, SHARED, BRANCH,
// UNCERTAIN} to each read using het-pair voting per bubble.
//
// Iterative: tightens or loosens the voting threshold based on how much
// signal previous iterations resolved. A bubble that gets a clear vote
// (≥ratio_threshold for one side) becomes "phased"; un-phased bubbles
// stay as UNCERTAIN regions.

#include "types.hpp"
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace branch::wg::phaser {

struct HetPair {
    std::uint64_t  kmer_a;     // allele-A k-mer canonical hash
    std::uint64_t  kmer_b;     // allele-B k-mer canonical hash
    NodeId         anchor_node;
    Position       anchor_pos_bp;
};

struct PhaseEngineOpts {
    int    iter_n = 1;
    int    min_supporting_overlaps = 10;
    double min_ratio = 0.95;

    // When set, force-assign every still-untagged read at the end of this
    // pass with UNCERTAIN tag. Used in the final iteration.
    bool   force_assign_remainder = false;
};

struct PhaseEngineProgress {
    std::size_t reads_newly_tagged = 0;
    std::size_t bubbles_newly_phased = 0;
    std::size_t reads_still_untagged = 0;
};

class PhaseEngine {
public:
    // Called once after overlap_graph; precomputes per-read het-pair
    // signatures so subsequent iterations don't re-extract them.
    void prepare(const std::vector<std::pair<ReadId, std::string>>& reads,
                 const std::vector<HetPair>& het_pairs);

    // Run one iteration of tagging. Updates `assignments` in place
    // (only writes UNTAGGED → tagged; never overwrites). Marks any
    // bubble that gets resolved.
    PhaseEngineProgress run_iteration(
        std::vector<UnitigNode>& unitigs,
        std::vector<Bubble>& bubbles,
        std::vector<ReadAssignment>& assignments,
        const PhaseEngineOpts& opts);

private:
    // Sparse: read_signatures_[i] is a list of (hetpair_idx, allele 0/1).
    std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>>
        read_signatures_;
};

}  // namespace branch::wg::phaser
