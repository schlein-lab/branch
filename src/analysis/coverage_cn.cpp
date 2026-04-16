#include "analysis/coverage_cn.hpp"

#include <algorithm>
#include <cmath>

namespace branch::analysis {

float median_baseline_coverage(
    const std::vector<branch::io::RegionCoverage>& baseline_regions) {
    if (baseline_regions.empty()) return 0.0f;
    std::vector<float> covs;
    covs.reserve(baseline_regions.size());
    for (const auto& r : baseline_regions) covs.push_back(r.mean_coverage);
    std::sort(covs.begin(), covs.end());
    const std::size_t n = covs.size();
    if (n % 2 == 1) return covs[n / 2];
    return 0.5f * (covs[n / 2 - 1] + covs[n / 2]);
}

std::vector<CNEstimate> estimate_cn(
    const std::vector<branch::io::RegionCoverage>& targets,
    float haploid_baseline) {
    std::vector<CNEstimate> out;
    out.reserve(targets.size());
    const float denom = haploid_baseline > 0.0f ? haploid_baseline : 1.0f;
    for (const auto& t : targets) {
        CNEstimate est;
        est.name = t.name;
        est.chrom = t.chrom;
        est.start = t.start;
        est.end = t.end;
        est.mean_coverage = t.mean_coverage;
        est.relative_cn = t.mean_coverage / denom;
        out.push_back(std::move(est));
    }
    return out;
}

ClusterSummary summarise_cluster(
    const std::vector<CNEstimate>& estimates,
    const std::string& cluster_name,
    const std::vector<std::string>& name_prefixes) {
    ClusterSummary s;
    s.cluster_name = cluster_name;
    for (const auto& est : estimates) {
        for (const auto& p : name_prefixes) {
            if (est.name.size() >= p.size() &&
                est.name.compare(0, p.size(), p) == 0) {
                s.member_names.push_back(est.name);
                s.total_relative_cn += est.relative_cn;
                break;
            }
        }
    }
    s.member_count = static_cast<float>(s.member_names.size());
    if (s.member_count > 0.0f) {
        s.mean_per_member = s.total_relative_cn / s.member_count;
    }
    return s;
}

}  // namespace branch::analysis
