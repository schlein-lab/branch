#pragma once

// BRANCH v0.1 — IGH locus analysis primitives.
//
// Small analytical building block used by early IGHG-family CNV probes
// while the full BRANCH pipeline is still scaffolding. Operates on
// mosdepth per-region-coverage output and produces paralog-aware
// signals that do NOT rely on the aligner's MAPQ filtering — every
// read-mass in a region counts, multi-mapping reads included.
//
// This is v0.1: rule-based estimators only. When v0.3 LightGBM
// classifier arrives, these signals become features, not decisions.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace branch::analysis {

// A single per-region coverage record, e.g. one row of mosdepth.regions.bed.gz
struct RegionCoverage {
    std::string chrom;
    std::uint32_t start;
    std::uint32_t end;
    std::string name;
    float mean_coverage;
};

// Paralog-cluster CN estimate. Sums mean_coverage across the paralog
// family and divides by the single-copy haploid baseline to get the
// family-wide total copy count. Individual gene CNs are undecidable
// from linear coverage alone — that is BRANCH's job proper — so this
// function only returns the cluster total.
struct ParalogClusterCN {
    std::vector<std::string> gene_names;
    float total_coverage_sum;
    float haploid_baseline;
    float total_copy_count_estimate;  // = total_coverage_sum / haploid_baseline
    float per_gene_mean_if_uniform;   // = total_copy_count_estimate / gene_names.size()
};

// Sum the coverage of a named subset of regions and divide by the
// provided haploid baseline (e.g. median coverage of single-copy
// control loci, or genome-wide mean / 2 for a diploid sample).
[[nodiscard]] ParalogClusterCN estimate_paralog_cluster_cn(
    const std::vector<RegionCoverage>& regions,
    const std::vector<std::string>& gene_filter,
    float haploid_baseline);

// Per-gene relative coverage ratio: mean_coverage / haploid_baseline.
// A ratio near 1 means one haploid copy; near 0 means apparent
// absence (but could be paralog-mapping loss). Low MAPQ upstream
// means this ratio should be interpreted at the *cluster* level
// rather than per gene.
struct PerGeneRelativeCN {
    std::string name;
    float relative_cn;
};

[[nodiscard]] std::vector<PerGeneRelativeCN> per_gene_relative_cn(
    const std::vector<RegionCoverage>& regions,
    float haploid_baseline);

}  // namespace branch::analysis
