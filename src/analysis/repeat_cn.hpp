// BRANCH — Repeat family copy-number analysis
// Computes CN estimates for repeat families from per-base coverage.

#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace branch::analysis {

/// Entry from a repeat annotation BED file (e.g., RepeatMasker output).
/// BED format: chrom<TAB>start<TAB>end<TAB>class/family
struct RepeatEntry {
    std::string family;   // e.g., "Alu"
    std::string class_;   // e.g., "SINE"
    std::string chrom;
    size_t start = 0;
    size_t end = 0;
};

/// Parse a repeat annotation BED file.
/// Expected BED format: chrom<TAB>start<TAB>end<TAB>class/family
/// The 4th column contains "class/family", which is split on '/'.
std::vector<RepeatEntry> parse_repeat_bed(const std::string& path);

/// Statistics for a repeat family's copy-number.
struct CnStat {
    size_t total_bp = 0;
    double mean_cov = 0.0;
    double cn_estimate = 0.0;
};

/// Compute CN statistics per repeat family.
/// @param entries Parsed repeat BED entries.
/// @param per_base_cov Map from chrom -> vector of per-base coverage depths.
/// @param baseline_cov Haploid baseline coverage (diploid = 2x baseline).
/// @return Map from family name to CnStat.
std::map<std::string, CnStat> compute_family_cn(
    const std::vector<RepeatEntry>& entries,
    const std::map<std::string, std::vector<uint32_t>>& per_base_cov,
    double baseline_cov);

/// Write CN statistics to a TSV file.
/// Columns: family, total_bp, mean_cov, cn_estimate
void write_cn_tsv(const std::map<std::string, CnStat>& stats,
                  const std::string& out_path);

}  // namespace branch::analysis
