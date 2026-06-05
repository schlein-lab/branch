#pragma once
//
// gfa_writer — emits the assembly graph in GFA1 format with P-lines for
// the h1 and h2 walks.
//
// Output structure:
//
//   H  VN:Z:1.0
//   S  shared_<id>     <seq>     LN:i:<len>  KIND:Z:shared
//   S  bubble_h1_<id>  <seq>     LN:i:<len>  KIND:Z:bubble_h1
//   S  bubble_h2_<id>  <seq>     LN:i:<len>  KIND:Z:bubble_h2
//   S  branch_<id>     <seq>     LN:i:<len>  KIND:Z:branch
//   L  <a> + <b> +  <overlap>M
//   P  hap1            <segment_path>+    <overlaps_csv>
//   P  hap2            <segment_path>+    <overlaps_csv>
//
// Downstream extraction:
//   gfatools gfa2fa -p hap1 out.gfa > out.h1.fa
//   gfatools gfa2fa -p hap2 out.gfa > out.h2.fa

#include "types.hpp"
#include "fastq_index.hpp"
#include <string>
#include <utility>
#include <vector>

namespace branch::wg::phaser {

struct GfaWriterOpts {
    std::string out_gfa_path;
    // Optional FASTA output paths for direct emission without gfatools.
    std::string out_h1_fa_path;
    std::string out_h2_fa_path;
    // Optional FASTA output for BRANCH-tagged unitigs (CNV duplicate
    // copies identified by bubble_classifier, type=BRANCH/BRANCH_IN_BRANCH).
    // Each unitig emits one contig: ">branch_<unitig_id>". v0.10+.
    std::string out_branch_fa_path;
    // v0.10 Variante B: read-level FASTA outputs for the new tag classes.
    // - bridge.fa: reads spanning ≥2 alts of a bubble (CNV-junction reads).
    //   Each emits one record: ">bridge_<read_name>".
    // - novel.fa: reads with sequence outside any known bubble alt (deep
    //   sub-branch candidates, possibly LINEAGE_GAP descendants).
    //   Each emits one record: ">novel_<read_name>".
    std::string out_bridge_fa_path;
    std::string out_novel_fa_path;
    // Optional read-assignment TSV (read_id, tag, confidence, home_node).
    std::string out_assignments_tsv_path;
    // Approximate overlap-cut between consecutive reads in a unitig
    // (kept consistent with overlap_graph collapse logic).
    int approx_overlap_bp = 5 * 23;
};

class GfaWriter {
public:
    // Emits all configured outputs. Skips empty paths.
    // `reads` is used to materialize unitig sequences on-the-fly via
    // FastqIndex random access, one unitig at a time, so peak RSS is
    // bounded.
    void write(
        const std::vector<UnitigNode>& unitigs,
        const std::vector<Bubble>& bubbles,
        const std::vector<ReadAssignment>& assignments,
        const std::vector<std::string>& read_id_strings,
        FastqIndex& reads,
        const GfaWriterOpts& opts);
};

}  // namespace branch::wg::phaser
