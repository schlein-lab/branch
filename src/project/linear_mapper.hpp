// BRANCH v0.4 — Linear reference mapper for `branch project`.
// Maps branch contigs against linear references via the LLmap classical engine.

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
    std::string preset = "asm20";        // contig-vs-ref divergence preset
    int threads = 4;
    int min_mapq = 0;                    // emit all; caller filters UNANNOTATED
};

// Map one branch FASTA against all given linear refs via the LLmap classical
// engine. Returns mappings sorted by (ref_name, branch_id, query_start). Maps
// EVERY branch in fasta_path against EVERY ref (one minimizer index per ref).
std::vector<LinearMapping> map_branches_linear(
    const std::string& fasta_path,
    const std::vector<LinearRef>& refs,
    const LinearMapOptions& opts,
    std::string* err_out);

}  // namespace branch::project

#endif
