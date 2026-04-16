// BRANCH — Repeat family copy-number analysis implementation.

#include "analysis/repeat_cn.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace branch::analysis {

std::vector<RepeatEntry> parse_repeat_bed(const std::string& path) {
    std::vector<RepeatEntry> entries;
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open repeat BED file: " + path);
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chrom, start_str, end_str, class_family;
        if (!(iss >> chrom >> start_str >> end_str >> class_family)) {
            continue;  // Skip malformed lines
        }

        RepeatEntry e;
        e.chrom = chrom;
        e.start = std::stoull(start_str);
        e.end = std::stoull(end_str);

        // Split class/family on '/'
        auto slash = class_family.find('/');
        if (slash != std::string::npos) {
            e.class_ = class_family.substr(0, slash);
            e.family = class_family.substr(slash + 1);
        } else {
            e.class_ = class_family;
            e.family = class_family;
        }

        entries.push_back(std::move(e));
    }

    return entries;
}

std::map<std::string, CnStat> compute_family_cn(
    const std::vector<RepeatEntry>& entries,
    const std::map<std::string, std::vector<uint32_t>>& per_base_cov,
    double baseline_cov) {

    // Accumulate total coverage and bp per family
    std::map<std::string, uint64_t> family_total_cov;
    std::map<std::string, size_t> family_total_bp;

    for (const auto& e : entries) {
        auto it = per_base_cov.find(e.chrom);
        if (it == per_base_cov.end()) continue;

        const auto& cov = it->second;
        size_t effective_end = std::min(e.end, cov.size());
        if (e.start >= effective_end) continue;

        uint64_t sum = 0;
        for (size_t i = e.start; i < effective_end; ++i) {
            sum += cov[i];
        }

        size_t bp = effective_end - e.start;
        family_total_cov[e.family] += sum;
        family_total_bp[e.family] += bp;
    }

    // Compute final stats
    std::map<std::string, CnStat> result;
    for (const auto& [family, total_bp] : family_total_bp) {
        if (total_bp == 0) continue;

        CnStat stat;
        stat.total_bp = total_bp;
        stat.mean_cov = static_cast<double>(family_total_cov[family]) /
                        static_cast<double>(total_bp);
        stat.cn_estimate = (baseline_cov > 0.0)
                               ? (stat.mean_cov / baseline_cov * 2.0)
                               : 0.0;
        result[family] = stat;
    }

    return result;
}

void write_cn_tsv(const std::map<std::string, CnStat>& stats,
                  const std::string& out_path) {
    std::ofstream out(out_path);
    if (!out) {
        throw std::runtime_error("Cannot open output file: " + out_path);
    }

    out << "family\ttotal_bp\tmean_cov\tcn_estimate\n";
    for (const auto& [family, s] : stats) {
        out << family << '\t' << s.total_bp << '\t' << s.mean_cov << '\t'
            << s.cn_estimate << '\n';
    }
}

}  // namespace branch::analysis
