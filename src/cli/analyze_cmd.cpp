// BRANCH v0.2 — `branch analyze` subcommand.
//
// Consumes a mosdepth regions.bed file (targets + optional baseline
// regions flagged by their names) and emits per-target CN estimates
// plus paralog-cluster summaries to stdout in TSV.

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "analysis/coverage_cn.hpp"
#include "io/mosdepth_reader.hpp"

namespace branch::cli {

namespace {

void print_analyze_usage(std::ostream& os) {
    os << "branch analyze — CN inference from mosdepth regions\n"
          "\nUsage:\n"
          "  branch analyze --regions <regions.bed>\n"
          "                 [--baseline-prefix <prefix>] [--baseline-cov <float>]\n"
          "                 [--cluster <name:prefix,prefix,...>]\n"
          "\nOne of --baseline-cov (explicit haploid baseline) or\n"
          "--baseline-prefix (median of regions whose name begins with prefix)\n"
          "must be provided.\n"
          "\nExample:\n"
          "  branch analyze --regions mosdepth.regions.bed \\\n"
          "                 --baseline-prefix control_ \\\n"
          "                 --cluster IGHG:IGHG1,IGHG2,IGHG3,IGHG4\n";
}

// Simple flag parser: captures --key value pairs.
struct Args {
    std::string regions;
    std::string baseline_prefix;
    float baseline_cov = 0.0f;
    std::vector<std::pair<std::string, std::vector<std::string>>> clusters;
    bool ok = false;
    std::string err;
};

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string_view k = argv[i];
        auto needs_val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                a.err = std::string("missing value for ") + name;
                return nullptr;
            }
            return argv[++i];
        };
        if (k == "--regions") {
            auto v = needs_val("--regions");
            if (!v) return a;
            a.regions = v;
        } else if (k == "--baseline-prefix") {
            auto v = needs_val("--baseline-prefix");
            if (!v) return a;
            a.baseline_prefix = v;
        } else if (k == "--baseline-cov") {
            auto v = needs_val("--baseline-cov");
            if (!v) return a;
            a.baseline_cov = std::strtof(v, nullptr);
        } else if (k == "--cluster") {
            auto v = needs_val("--cluster");
            if (!v) return a;
            std::string s = v;
            auto colon = s.find(':');
            if (colon == std::string::npos) {
                a.err = "--cluster requires NAME:PREFIX1,PREFIX2,...";
                return a;
            }
            a.clusters.emplace_back(s.substr(0, colon), split_csv(s.substr(colon + 1)));
        } else if (k == "--help" || k == "-h") {
            a.err = "HELP";
            return a;
        } else {
            a.err = std::string("unknown arg: ") + std::string(k);
            return a;
        }
    }
    if (a.regions.empty()) {
        a.err = "--regions is required";
        return a;
    }
    if (a.baseline_cov <= 0.0f && a.baseline_prefix.empty()) {
        a.err = "one of --baseline-cov or --baseline-prefix is required";
        return a;
    }
    a.ok = true;
    return a;
}

}  // namespace

int run_analyze(int argc, char** argv) {
    Args a = parse_args(argc, argv);
    if (!a.ok) {
        if (a.err == "HELP") {
            print_analyze_usage(std::cout);
            return 0;
        }
        std::cerr << "branch analyze: " << a.err << "\n\n";
        print_analyze_usage(std::cerr);
        return 2;
    }

    auto regions = branch::io::read_mosdepth_regions(a.regions);
    if (regions.empty()) {
        std::cerr << "branch analyze: failed to read regions from " << a.regions << "\n";
        return 3;
    }

    // Resolve baseline.
    float baseline = a.baseline_cov;
    if (baseline <= 0.0f) {
        std::vector<branch::io::RegionCoverage> baseline_regions;
        for (const auto& r : regions) {
            if (r.name.size() >= a.baseline_prefix.size() &&
                r.name.compare(0, a.baseline_prefix.size(), a.baseline_prefix) == 0) {
                baseline_regions.push_back(r);
            }
        }
        if (baseline_regions.empty()) {
            std::cerr << "branch analyze: no baseline regions matched prefix '"
                      << a.baseline_prefix << "'\n";
            return 4;
        }
        baseline = branch::analysis::median_baseline_coverage(baseline_regions);
    }

    std::cout << "# branch analyze\n";
    std::cout << "# baseline_haploid_coverage=" << baseline << "\n";
    std::cout << "# regions=" << regions.size() << "\n";
    std::cout << "chrom\tstart\tend\tname\tmean_cov\trelative_cn\n";
    auto est = branch::analysis::estimate_cn(regions, baseline);
    for (const auto& e : est) {
        std::cout << e.chrom << '\t' << e.start << '\t' << e.end << '\t'
                  << e.name << '\t' << e.mean_coverage << '\t' << e.relative_cn << '\n';
    }

    // Cluster summaries.
    for (const auto& [cluster_name, prefixes] : a.clusters) {
        auto s = branch::analysis::summarise_cluster(est, cluster_name, prefixes);
        std::cout << "# cluster " << s.cluster_name
                  << " members=" << s.member_count
                  << " total_rel_cn=" << s.total_relative_cn
                  << " mean_per_member=" << s.mean_per_member << "\n";
    }
    return 0;
}

}  // namespace branch::cli
