// BRANCH v0.4 — Linear reference mapper for `branch project`.
// Shell-out to minimap2 per branch, parse PAF, return mappings.

#ifndef BRANCH_PROJECT_LINEAR_MAPPER_HPP
#define BRANCH_PROJECT_LINEAR_MAPPER_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace branch::project {

struct LinearMapping {
    std::string branch_id;    // query name (branch ID from FASTA)
    std::string ref_name;     // logical name (CHM13, GRCh38, ...)
    std::string target;       // chrom / contig
    std::int64_t target_start = 0;
    std::int64_t target_end = 0;
    int mapq = 0;
    char strand = '+';
    std::int64_t query_len = 0;
    std::int64_t query_start = 0;
    std::int64_t query_end = 0;
};

struct LinearRef {
    std::string name;    // user-tag (CHM13, GRCh38)
    std::string path;    // .fa or .mmi
};

struct LinearMapOptions {
    std::string minimap2_path = "minimap2";
    std::string preset = "asm20";        // HiFi-vs-ref default
    int threads = 4;
    int min_mapq = 0;                    // emit all; caller filters UNANNOTATED
};

// Map one branch FASTA against all given linear refs. Returns mappings sorted
// by (ref_name, query_start). Shell-out to minimap2 preset (asm20) — one call
// per (branch-FASTA x ref) pair. Caller is responsible for iterating over
// branches; this maps EVERY branch in fasta_path against EVERY ref.
std::vector<LinearMapping> map_branches_linear(
    const std::string& fasta_path,
    const std::vector<LinearRef>& refs,
    const LinearMapOptions& opts,
    std::string* err_out);

}  // namespace branch::project

#endif
