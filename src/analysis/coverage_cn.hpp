#pragma once

// BRANCH v0.2 — Coverage-based CN inference.
//
// Given mosdepth per-region coverage and a set of "known single-copy"
// baseline regions, produce per-target CN estimates (mean_coverage /
// haploid_baseline) and paralog-cluster totals. This is NOT a full
// BRANCH call — it operates on linear, collapsed alignments — but it
// is the natural v0.2 first step and the ground-truth sanity check
// against which future graph-based calls will be compared.

#include <cstddef>
#include <string>
#include <vector>

#include "io/mosdepth_reader.hpp"

namespace branch::analysis {

struct CNEstimate {
    std::string name;
    std::string chrom;
    std::uint32_t start{};
    std::uint32_t end{};
    float mean_coverage{};
    float relative_cn{};  // mean_coverage / haploid_baseline
};

// Compute the median mean_coverage of the provided baseline regions.
// Median is robust to outliers (one failed region doesn't skew).
[[nodiscard]] float median_baseline_coverage(
    const std::vector<branch::io::RegionCoverage>& baseline_regions);

// Produce CN estimates for every target region using the given
// haploid baseline. Caller typically uses 0.5 * diploid_median or
// the per-locus single-copy baseline.
[[nodiscard]] std::vector<CNEstimate> estimate_cn(
    const std::vector<branch::io::RegionCoverage>& targets,
    float haploid_baseline);

// Sum the relative CN of target regions matching any of the provided
// name prefixes. Used for paralog-cluster totals like IGHG1+2+3+4.
struct ClusterSummary {
    std::string cluster_name;
    std::vector<std::string> member_names;
    float total_relative_cn{};
    float member_count{};
    float mean_per_member{};
};

[[nodiscard]] ClusterSummary summarise_cluster(
    const std::vector<CNEstimate>& estimates,
    const std::string& cluster_name,
    const std::vector<std::string>& name_prefixes);

}  // namespace branch::analysis
