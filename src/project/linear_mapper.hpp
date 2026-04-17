// BRANCH v0.4 — Linear reference mapping via minimap2 shell-out.
//
// Maps branch consensus sequences to linear references (CHM13, GRCh38)
// using minimap2 asm20 preset. Results are merged into a single PAF
// with rf:Z: tags indicating the reference source.

#pragma once

#include <string>
#include <vector>

namespace branch::project {

/// Represents a linear reference configuration.
struct LinearRef {
    std::string name;      ///< Reference name (e.g., "CHM13", "GRCh38")
    std::string path;      ///< Path to .mmi or .fasta file
};

/// Result of linear mapping for a single branch.
struct LinearMapResult {
    std::string branch_id;
    std::string ref_name;
    int mapq{0};
    double identity{0.0};
    std::string paf_line;  ///< Raw PAF record
};

/// Configuration for linear mapping.
struct LinearMapperConfig {
    std::vector<LinearRef> references;
    std::string fasta_path;        ///< Branch consensus FASTA
    std::string output_paf;        ///< Output PAF path
    int threads{1};
    std::string minimap2_path{"minimap2"};  ///< Path to minimap2 binary
};

/// Run minimap2 on all branches against all linear references.
/// Merges results into a single PAF with rf:Z: tags.
///
/// @param config  Mapping configuration
/// @return true on success, false on error
///
/// TODO: Implement in v0.4
bool run_linear_mapping(const LinearMapperConfig& config);

/// Parse a PAF file and extract per-branch mapping statistics.
///
/// @param paf_path  Path to PAF file
/// @return Vector of mapping results
///
/// TODO: Implement in v0.4
std::vector<LinearMapResult> parse_linear_paf(const std::string& paf_path);

}  // namespace branch::project
