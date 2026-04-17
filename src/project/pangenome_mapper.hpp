// BRANCH v0.4 — Pangenome mapping via minigraph/GraphAligner.
//
// Maps branch consensus sequences to HPRC v1.1 pangenome graphs (GBZ format).
// Identifies the closest haplotype path for each branch.

#pragma once

#include <string>
#include <vector>

namespace branch::project {

/// Pangenome mapping backend selection.
enum class PangenomeMapper {
    Minigraph,      ///< minigraph -c (fast, CIGAR-free)
    GraphAligner,   ///< GraphAligner (high precision for HiFi)
};

/// Represents a pangenome reference configuration.
struct PangenomeRef {
    std::string name;      ///< Reference backbone name (e.g., "CHM13-based", "GRCh38-based")
    std::string gbz_path;  ///< Path to .gbz file
};

/// Result of pangenome mapping for a single branch.
struct PangenomeMapResult {
    std::string branch_id;
    std::string closest_path;   ///< HPRC haplotype path (e.g., "HG00733#1#chr1:...")
    int mapq{0};
    double identity{0.0};
    std::string gaf_line;       ///< Raw GAF record
};

/// Configuration for pangenome mapping.
struct PangenomeMapperConfig {
    std::vector<PangenomeRef> references;
    std::string fasta_path;           ///< Branch consensus FASTA
    std::string output_gaf;           ///< Output GAF path
    PangenomeMapper mapper{PangenomeMapper::Minigraph};
    int threads{1};
    std::string minigraph_path{"minigraph"};
    std::string graphaligner_path{"GraphAligner"};
};

/// Run pangenome mapping on all branches.
///
/// @param config  Mapping configuration
/// @return true on success, false on error
///
/// TODO: Implement in v0.4
bool run_pangenome_mapping(const PangenomeMapperConfig& config);

/// Parse a GAF file and extract per-branch mapping statistics.
///
/// @param gaf_path  Path to GAF file
/// @return Vector of mapping results
///
/// TODO: Implement in v0.4
std::vector<PangenomeMapResult> parse_pangenome_gaf(const std::string& gaf_path);

}  // namespace branch::project
