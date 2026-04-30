// BRANCH v0.5 — Phase 2: transitive branch attachment.
//
// Inputs from upstream phases:
//   - direct branch candidates (curator-emitted, anchored at hap:pos)
//   - ambiguous reads from the router (no hap1/hap2 majority)
//
// Output:
//   - AttachedBranch: each branch with anchor metadata + transitive depth
//   - OrphanRead: ambig reads that did not anchor anywhere
//
// Algorithm:
//   1. Direct: each curator candidate becomes one AttachedBranch (depth=0)
//   2. Sketch every ambig read with MinHash-of-minimizers
//   3. Build inverted index minimizer → (read_idx, side='ambig'/'branch')
//   4. For each ambig read: count shared sketch minimizers with branches
//      and with haplotypes. If max-shared ≥ jaccard_lo * sketch_size →
//      attach as transitive depth=1.
//   5. Remaining ambig reads cluster among themselves (BFS hops up to
//      max_transitive_depth); if a cluster eventually anchors via any
//      member, all members get the cluster's depth + their own hop.
//   6. Reads still unanchored → OrphanRead list (Phase 5 input).

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"
#include "wg/master_curator.hpp"

namespace branch::wg {

struct AttachedBranch {
    std::string branch_name;
    std::string seq;
    std::uint32_t anchor_hap;
    std::uint32_t anchor_pos;
    std::uint32_t transitive_depth;  // 0 = direct, 1+ = via other branches
    double confidence;
};

struct OrphanRead {
    std::string read_id;
    std::string seq;
};

class BranchAttacher {
public:
    explicit BranchAttacher(const ::branch::graph::TechProfile& profile);

    /// `hap_seqs` are the curated haplotype sequences (Phase 1.5 output);
    /// `candidates` are the curator-emitted direct branch candidates;
    /// `ambig_reads` are router output where neither hap clearly won.
    void attach(const std::vector<std::string>& hap_seqs,
                const std::vector<BranchCandidate>& candidates,
                const std::vector<std::pair<std::string, std::string>>& ambig_reads,
                std::vector<AttachedBranch>& out_branches,
                std::vector<OrphanRead>& out_orphans) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
