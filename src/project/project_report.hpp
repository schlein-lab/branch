// BRANCH v0.4 — JSON report writer for `branch project`.

#ifndef BRANCH_PROJECT_PROJECT_REPORT_HPP
#define BRANCH_PROJECT_PROJECT_REPORT_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "linear_mapper.hpp"
#include "pangenome_mapper.hpp"
#include "somatic_delta.hpp"

namespace branch::project {

struct BranchEntry {
    std::string branch_id;
    std::int64_t length_bp = 0;
    double vaf = -1.0;                  // -1 = unknown
    double coverage = -1.0;
    std::vector<LinearMapping> linear_mappings;
    // v0.4.2: HPRC pangenome-graph alignments per branch (one entry per
    // pangenome ref the user passed via `--ref-pangenome`).
    std::vector<PangenomeMapping> pangenome_mappings;
    // v0.4.3: edit distance between this branch and the closest known
    // reference sequence (linear ref window at the mapped coords). One
    // entry per linear ref the branch aligned to. CIGAR-encoded.
    std::vector<SomaticDelta> somatic_deltas;
    bool unannotated = false;           // true if no high-quality mapping anywhere
};

// Write the branch report to JSON. One top-level object with "sample",
// "version" ("0.4.3"), and "branches" array. Returns false on IO error.
bool write_branch_report_json(
    const std::string& out_path,
    const std::string& sample_name,
    const std::vector<BranchEntry>& branches,
    std::string* err_out);

}  // namespace branch::project

#endif
