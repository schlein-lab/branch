// BRANCH v0.4 — JSON report writer for `branch project`.

#ifndef BRANCH_PROJECT_PROJECT_REPORT_HPP
#define BRANCH_PROJECT_PROJECT_REPORT_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "linear_mapper.hpp"

namespace branch::project {

struct BranchEntry {
    std::string branch_id;
    std::int64_t length_bp = 0;
    double vaf = -1.0;               // -1 = unknown
    double coverage = -1.0;
    std::vector<LinearMapping> linear_mappings;
    bool unannotated = false;        // true if no high-quality mapping anywhere
};

// Write the branch report to JSON. One top-level object with "sample",
// "version" ("0.4.1"), and "branches" array. Returns false on IO error.
bool write_branch_report_json(
    const std::string& out_path,
    const std::string& sample_name,
    const std::vector<BranchEntry>& branches,
    std::string* err_out);

}  // namespace branch::project

#endif
