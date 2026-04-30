// BRANCH v0.5 — Phase 1.5 implementation (skeleton).
//
// First-cut: emits zero curation events and zero branches. The full
// implementation requires a per-base pile-up reader for the per-hap
// alignment file (Phase 1.4 output); that file format is not yet
// finalised. This skeleton lets the orchestrator wire the function
// call without the inner work.

#include "wg/master_curator.hpp"

namespace branch::wg {

MasterCurator::MasterCurator(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void MasterCurator::curate(
    std::vector<MasterTile>& /*tiles*/,
    const std::string& /*reads_aligned_to_master_path*/,
    std::vector<CurationEvent>& out_events,
    std::vector<BranchCandidate>& out_branches) const {
    out_events.clear();
    out_branches.clear();
    // TODO: implement pile-up reader + Q-weighted evidence + thresholding.
    // Until then, no curation is performed and no branches are emitted.
}

}  // namespace branch::wg
