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
#include "fastq_index.hpp"
#include "bubble_voter.hpp"
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
    // Streams sequences from the FastqIndex sequentially.
    void prepare(FastqIndex& reads,
                 const std::vector<HetPair>& het_pairs);

    // Run one iteration of tagging. Updates `assignments` in place
    // (only writes UNTAGGED → tagged; never overwrites). Marks any
    // bubble that gets resolved.
    PhaseEngineProgress run_iteration(
        std::vector<UnitigNode>& unitigs,
        std::vector<Bubble>& bubbles,
        std::vector<ReadAssignment>& assignments,
        const PhaseEngineOpts& opts);

    // Bridge: take per-read voter outcomes and write them into
    // assignments. This OVERRIDES home_node-based tags so that every
    // read assigned/bridged/marked-novel by the voter gets the
    // corresponding ReadTag, regardless of where its home_node sat.
    //
    // Priority (per read):
    //   ASSIGNED → tag from primary bubble's alt kind (H1/H2/BRANCH)
    //   BRIDGE   → BRANCH_BRIDGE (CNV-junction)
    //   NOVEL    → BRANCH_NOVEL  (deep sub-branch candidate)
    //   ABSENT   → don't touch (home_node-based tagging took over,
    //              usually SHARED in homozygous regions)
    static void apply_voter_outcomes(
        const std::vector<UnitigNode>& unitigs,
        const std::vector<Bubble>& bubbles,
        const std::vector<VoterReadOutcome>& outcomes,
        std::vector<ReadAssignment>& assignments);

private:
    // Sparse: read_signatures_[i] is a list of (hetpair_idx, allele 0/1).
    std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>>
        read_signatures_;

    // hetpair_anchors_[hp_idx] = HetPair.anchor_node. Used by run_iteration
    // to attribute each (signature, allele) to the bubble whose alt-path
    // contains the anchor. This decouples voting from home_node — at WG
    // scale a read's home_node is one of many unitigs it spans, so
    // home_node-based voting starves bubble alts of evidence (v11: 12.7K
    // votes over 2550 bubbles → 0 phased).
    std::vector<NodeId> hetpair_anchors_;
};

}  // namespace branch::wg::phaser
