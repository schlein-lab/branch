// BRANCH v0.5 — Phase 3 implementation (skeleton, leans on existing src/vpcr/).
//
// First cut: empty design + zero counts. The existing src/vpcr/ module
// already contains the primer-search + counting infrastructure; future
// commit will wire VpcrAutoDesigner to delegate to it instead of
// reimplementing.

#include "wg/vpcr_wg.hpp"

namespace branch::wg {

VpcrAutoDesigner::VpcrAutoDesigner(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void VpcrAutoDesigner::design(
    const std::vector<std::string>& /*hap_seqs*/,
    std::vector<VpcrAmplicon>& out) const {
    out.clear();
    // TODO: window-sliding + entropy + primer-pairing
}

void VpcrAutoDesigner::count(
    std::vector<VpcrAmplicon>& /*amplicons*/,
    const std::vector<std::pair<std::string, std::string>>& /*reads*/,
    int /*expected_single_copy_count*/) const {
    // TODO: bracketed primer scan
}

}  // namespace branch::wg
