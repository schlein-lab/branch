// BRANCH v0.4 — Pangenome mapper for `branch project`.
// Shell-out to GraphAligner or minigraph, parse GAF, return mappings.

#ifndef BRANCH_PROJECT_PANGENOME_MAPPER_HPP
#define BRANCH_PROJECT_PANGENOME_MAPPER_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace branch::project {

struct PangenomeMapping {
    std::string branch_id;     // query name (branch ID from FASTA)
    std::string ref_name;      // logical name (HPRC-CHM13, HPRC-GRCh38, ...)
    std::string path_names;    // HPRC haplotype path (e.g., ">HG00733#1#chr1...")
    std::int64_t path_start = 0;
    std::int64_t path_end = 0;
    int mapq = 0;
    double identity = 0.0;
    std::int64_t query_len = 0;
};

struct PangenomeRef {
    std::string name;       // user-tag (HPRC-CHM13, HPRC-GRCh38)
    std::string gbz_path;   // .gbz / .gfa
};

struct PangenomeMapOptions {
    std::string tool = "GraphAligner";      // "GraphAligner" or "minigraph"
    std::string tool_path = "GraphAligner"; // binary path
    int threads = 4;
    int min_mapq = 0;                       // caller filters UNANNOTATED
};

// Map one branch FASTA against all given pangenome refs. Returns mappings
// sorted by (ref_name, branch_id, path_start). Shell-out to GraphAligner
// (default) or minigraph.
std::vector<PangenomeMapping> map_branches_pangenome(
    const std::string& fasta_path,
    const std::vector<PangenomeRef>& refs,
    const PangenomeMapOptions& opts,
    std::string* err_out);

}  // namespace branch::project

#endif
