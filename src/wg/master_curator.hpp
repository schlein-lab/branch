// BRANCH v0.5 — Phase 1.5: master-read pile-up curation.

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

    void curate(std::vector<MasterTile>& tiles,
                const std::string& reads_aligned_to_master_path,
                std::vector<CurationEvent>& out_events,
                std::vector<BranchCandidate>& out_branches) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
