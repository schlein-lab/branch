// BRANCH v0.5 — Phase 1.5: master-read pile-up curation.
//
// For each master tile:
//   1. minimize the master sequence
//   2. for each routed read: find shared minimizer hits → estimate
//      offset of read against the tile (anchor pile-up, no full SW)
//   3. walk the read base-by-base at that offset; for every position
//      where read-base != master-base, accumulate evidence
//   4. when a position's alt-base count crosses min_coverage_branch and
//      Q-weighted evidence exceeds threshold → emit BranchCandidate;
//      when master-base agreement < min_base_agreement → emit
//      CurationEvent (master gets corrected to majority)

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"
#include "wg/master_tiler.hpp"

namespace branch::wg {

struct CurationEvent {
    std::uint32_t hap_idx;        // 0 or 1
    std::uint32_t pos_in_hap;
    char original_base;
    char curated_base;
    std::uint32_t n_supporting_reads;
};

struct BranchCandidate {
    std::uint32_t hap_idx;
    std::uint32_t pos_in_hap;
    std::string alt_seq;          // alt allele (≥ 1 bp)
    std::vector<std::uint32_t> supporting_read_ids;
    double q_weighted_evidence;   // sum of sigmoid(Q-12)/5 across supporters
};

class MasterCurator {
public:
    explicit MasterCurator(const ::branch::graph::TechProfile& profile);

    /// Curate one haplotype against the read pool routed to it.
    void curate(std::uint32_t hap_idx,
                const std::vector<MasterTile>& tiles,
                const std::string& haplotype_seq,
                const std::vector<std::pair<std::string, std::string>>& reads,
                std::vector<CurationEvent>& out_events,
                std::vector<BranchCandidate>& out_branches) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
