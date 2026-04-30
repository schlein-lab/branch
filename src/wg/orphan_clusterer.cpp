// BRANCH v0.5 — Phase 5 implementation (singleton-cluster first cut).

#include "wg/orphan_clusterer.hpp"

namespace branch::wg {

OrphanClusterer::OrphanClusterer(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void OrphanClusterer::cluster(const std::vector<OrphanRead>& orphans,
                              std::vector<OrphanCluster>& out) const {
    // First-cut: each orphan is its own singleton cluster. Future
    // commit: minHash sketch + Jaccard cluster + consensus hint.
    out.clear();
    for (std::size_t i = 0; i < orphans.size(); ++i) {
        OrphanCluster c{};
        c.cluster_id = static_cast<std::uint32_t>(i);
        c.member_read_ids.push_back(orphans[i].read_id);
        c.consensus_hint.clear();
        c.ref_attempt_chrom = "unmapped";
        c.ref_attempt_score = 0.0;
        out.push_back(std::move(c));
    }
}

}  // namespace branch::wg
