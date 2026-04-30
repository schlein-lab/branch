// BRANCH v0.5 — Phase 5: orphan re-clustering.

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"
#include "wg/branch_attacher.hpp"

namespace branch::wg {

struct OrphanCluster {
    std::uint32_t cluster_id;
    std::vector<std::string> member_read_ids;
    std::string consensus_hint;
    std::string ref_attempt_chrom;
    double ref_attempt_score;
};

class OrphanClusterer {
public:
    explicit OrphanClusterer(const ::branch::graph::TechProfile& profile);

    void cluster(const std::vector<OrphanRead>& orphans,
                 std::vector<OrphanCluster>& out) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
