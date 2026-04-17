// BRANCH v0.4 — Project report generation (JSON + VCF).
//
// Generates comprehensive output files from projection results:
// - branch-report.json: Detailed per-branch analysis
// - somatic.vcf: VCF with SNV/INDEL/SV calls

#pragma once

#include <string>
#include <vector>

#include "linear_mapper.hpp"
#include "pangenome_mapper.hpp"
#include "somatic_delta.hpp"

namespace branch::project {

/// Combined projection result for a single branch.
struct BranchProjection {
    std::string branch_id;
    std::uint32_t length_bp{0};
    double coverage{0.0};
    double vaf{0.0};

    // Linear mapping results (per reference)
    std::vector<LinearMapResult> linear_results;

    // Pangenome mapping result
    PangenomeMapResult pangenome_result;

    // Somatic delta result
    SomaticDeltaResult delta_result;

    // Classification
    bool unannotated{false};  ///< True if no confident annotation
};

/// Summary statistics for the projection.
struct ProjectionSummary {
    std::uint32_t total_branches{0};
    std::uint32_t annotated{0};
    std::uint32_t unannotated{0};
    std::uint32_t somatic_snvs{0};
    std::uint32_t somatic_indels{0};
    std::uint32_t somatic_svs{0};
};

/// Configuration for report generation.
struct ReportConfig {
    std::string sample_name;
    std::string out_prefix;
    int min_mapq{20};       ///< Minimum MAPQ for annotation
};

/// Write branch-report.json.
///
/// @param projections  Per-branch projection results
/// @param summary      Summary statistics
/// @param config       Report configuration
/// @return true on success
///
/// TODO: Implement in v0.4
bool write_json_report(
    const std::vector<BranchProjection>& projections,
    const ProjectionSummary& summary,
    const ReportConfig& config);

/// Write somatic.vcf with variant calls.
///
/// @param projections  Per-branch projection results
/// @param config       Report configuration
/// @return true on success
///
/// TODO: Implement in v0.4
bool write_vcf(
    const std::vector<BranchProjection>& projections,
    const ReportConfig& config);

/// Compute summary statistics from projections.
///
/// @param projections  Per-branch projection results
/// @return Summary statistics
///
/// TODO: Implement in v0.4
ProjectionSummary compute_summary(const std::vector<BranchProjection>& projections);

}  // namespace branch::project
