// BRANCH v0.5 — Phase 2: transitive branch attachment.

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

    void attach(const std::vector<BranchCandidate>& candidates,
                const std::vector<std::pair<std::string, std::string>>& ambig_reads,
                std::vector<AttachedBranch>& out_branches,
                std::vector<OrphanRead>& out_orphans) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
