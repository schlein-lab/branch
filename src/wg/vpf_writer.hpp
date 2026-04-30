// BRANCH v0.5 — extended-GFA writer.
//
// Writes the 3 standard files (HG002.fa.gz, HG002.bed.gz, HG002.gfa.zst)
// + the indices the indexing tools expect.

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "wg/master_tiler.hpp"
#include "wg/master_curator.hpp"
#include "wg/branch_attacher.hpp"
#include "wg/vpcr_wg.hpp"
#include "wg/orphan_clusterer.hpp"
#include "wg/haplotype_router.hpp"

namespace branch::wg {

struct WgRunMeta {
    std::string sample;
    std::string assembler_version;
    std::string read_tech;
    std::uint64_t input_reads;
    std::uint64_t input_bases;
    int detected_coverage;
    int detected_het_coverage;
    double bimodality_z;
};

struct WgOutputs {
    WgRunMeta meta;
    std::vector<std::string> haplotype_seqs;        // hap1, hap2 (sequence bytes)
    std::vector<std::vector<MasterTile>> haplotype_tiles;  // per-hap
    std::vector<CurationEvent> curation_events;
    std::vector<AttachedBranch> branches;
    std::vector<OrphanCluster> orphan_clusters;
    std::vector<VpcrAmplicon> vpcr_amplicons;
    std::vector<std::pair<std::string, std::string>> read_disposition;
                                          // (read_id, "phase=X outcome=Y reason=Z")
    std::vector<std::pair<std::string, double>> phase_log;
                                          // (phase_name, wall_seconds)
};

class VpfWriter {
public:
    explicit VpfWriter(const std::string& out_prefix);

    /// Write all 3 outputs (.fa.gz, .bed.gz, .gfa.zst) + indices.
    /// On error throws std::runtime_error.
    void write(const WgOutputs& out);

private:
    std::string prefix_;
};

}  // namespace branch::wg
