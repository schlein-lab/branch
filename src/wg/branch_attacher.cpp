// BRANCH v0.5 — Phase 2 implementation (skeleton + direct-only first cut).

#include "wg/branch_attacher.hpp"
#include <sstream>

namespace branch::wg {

BranchAttacher::BranchAttacher(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void BranchAttacher::attach(
    const std::vector<BranchCandidate>& candidates,
    const std::vector<std::pair<std::string, std::string>>& ambig_reads,
    std::vector<AttachedBranch>& out_branches,
    std::vector<OrphanRead>& out_orphans) const {
    out_branches.clear();
    out_orphans.clear();

    // Direct branches: each curator-emitted candidate becomes one attached
    // branch with depth=0. Confidence proxied by Q-weighted evidence.
    std::uint32_t idx = 0;
    for (const auto& c : candidates) {
        AttachedBranch b{};
        std::ostringstream name;
        name << "branch_" << std::dec << ++idx;
        b.branch_name = name.str();
        b.seq = c.alt_seq;
        b.anchor_hap = c.hap_idx;
        b.anchor_pos = c.pos_in_hap;
        b.transitive_depth = 0;
        b.confidence = c.q_weighted_evidence;
        out_branches.push_back(std::move(b));
    }

    // Transitive attachment + orphan re-clustering: out of scope for
    // this first cut. Ambig reads with no direct match → orphan list.
    for (const auto& [id, seq] : ambig_reads) {
        out_orphans.push_back(OrphanRead{id, seq});
    }
}

}  // namespace branch::wg
