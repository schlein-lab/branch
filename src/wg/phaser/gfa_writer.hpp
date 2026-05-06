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
#include <string>
#include <vector>

namespace branch::wg::phaser {

struct GfaWriterOpts {
    std::string out_gfa_path;
    // Optional FASTA output paths for direct emission without gfatools.
    std::string out_h1_fa_path;
    std::string out_h2_fa_path;
    // Optional read-assignment TSV (read_id, tag, confidence, home_node).
    std::string out_assignments_tsv_path;
};

class GfaWriter {
public:
    // Emits all configured outputs. Skips empty paths.
    void write(
        const std::vector<UnitigNode>& unitigs,
        const std::vector<Bubble>& bubbles,
        const std::vector<ReadAssignment>& assignments,
        const std::vector<std::string>& read_id_strings,
        const GfaWriterOpts& opts);
};

}  // namespace branch::wg::phaser
