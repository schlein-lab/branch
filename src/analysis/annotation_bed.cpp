#include "analysis/annotation_bed.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <stdexcept>

namespace branch::analysis {

namespace {

std::vector<std::uint32_t> parse_alt_supports(const std::string& tok) {
    std::vector<std::uint32_t> out;
    if (tok.empty() || tok == ".") return out;
    std::string cur;
    for (char c : tok) {
        if (c == ',') {
            if (!cur.empty()) {
                try { out.push_back(static_cast<std::uint32_t>(std::stoul(cur))); }
                catch (const std::exception&) { /* skip malformed entry */ }
                cur.clear();
            }
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        try { out.push_back(static_cast<std::uint32_t>(std::stoul(cur))); }
        catch (const std::exception&) {}
    }
    return out;
}

}  // namespace

std::vector<BedEntry> load_bed(const std::string& path) {
    std::vector<BedEntry> entries;
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("Cannot open BED file: " + path);
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        BedEntry e;
        if (!(iss >> e.chrom >> e.start >> e.end)) continue;
        iss >> e.name; // optional
        // Optional 5th column: per-bubble confidence. "." = unset.
        std::string conf_tok;
        if (iss >> conf_tok && conf_tok != ".") {
            try {
                e.confidence = std::stof(conf_tok);
            } catch (const std::exception&) {
                // Leave unset; malformed confidence is non-fatal.
            }
        }
        // Optional 6th column: comma-separated per-alt read supports.
        // Columns 7 (min_vaf) and 8 (total_support) are derived and skipped
        // on load — they round-trip via recomputation from alt_read_supports.
        std::string alt_tok;
        if (iss >> alt_tok) {
            e.alt_read_supports = parse_alt_supports(alt_tok);
        }
        entries.push_back(std::move(e));
    }
    return entries;
}

void write_bed_entry(std::ostream& out, const BedEntry& entry) {
    out << entry.chrom << '\t' << entry.start << '\t' << entry.end
        << '\t' << (entry.name.empty() ? std::string(".") : entry.name);
    if (entry.confidence.has_value() && !std::isnan(*entry.confidence)) {
        std::ostringstream conf;
        conf << std::fixed << std::setprecision(3) << *entry.confidence;
        out << '\t' << conf.str();
    } else {
        out << "\t.";
    }

    if (!entry.alt_read_supports.empty()) {
        // col 6: comma-separated alt supports
        out << '\t';
        for (std::size_t i = 0; i < entry.alt_read_supports.size(); ++i) {
            if (i) out << ',';
            out << entry.alt_read_supports[i];
        }

        // col 7: min/sum VAF (only defined for >=2 alts with non-zero sum)
        const auto& sup = entry.alt_read_supports;
        std::uint64_t sum = 0;
        for (auto v : sup) sum += v;
        if (sup.size() >= 2 && sum > 0) {
            std::uint32_t mn = *std::min_element(sup.begin(), sup.end());
            double vaf = static_cast<double>(mn) / static_cast<double>(sum);
            std::ostringstream vs;
            vs << std::fixed << std::setprecision(4) << vaf;
            out << '\t' << vs.str();
        } else {
            out << "\t.";
        }

        // col 8: total alt support
        out << '\t' << sum;

        // col 9-10: phased per-alt counts + hap_skew. Emitted only
        // when phasing is available; otherwise the columns are absent
        // (the BED stays at 8 columns, callers that don't expect 9-10
        // see no surprise).
        if (!entry.alt_read_supports_phased.empty() &&
            entry.alt_read_supports_phased.size() == entry.alt_read_supports.size()) {
            out << '\t';
            std::uint64_t total_phased = 0;
            std::uint32_t max_skew_abs = 0;
            for (std::size_t i = 0; i < entry.alt_read_supports_phased.size(); ++i) {
                if (i) out << ',';
                const auto h0 = entry.alt_read_supports_phased[i][0];
                const auto h1 = entry.alt_read_supports_phased[i][1];
                out << h0 << '|' << h1;
                total_phased += h0 + h1;
                const std::uint32_t diff = h0 > h1 ? h0 - h1 : h1 - h0;
                if (diff > max_skew_abs) max_skew_abs = diff;
            }
            if (total_phased >= 2) {
                std::ostringstream ss;
                ss << std::fixed << std::setprecision(4)
                   << (static_cast<double>(max_skew_abs) / static_cast<double>(total_phased));
                out << '\t' << ss.str();
            } else {
                out << "\t.";
            }
        }
    }
    out << '\n';
}

IntervalIndex::IntervalIndex(const std::vector<BedEntry>& entries) {
    for (const auto& e : entries) {
        by_chrom_[e.chrom].push_back(e);
    }
}

std::vector<const BedEntry*> IntervalIndex::query(const std::string& chrom, size_t start, size_t end) const {
    std::vector<const BedEntry*> hits;
    auto it = by_chrom_.find(chrom);
    if (it == by_chrom_.end()) return hits;
    for (const auto& e : it->second) {
        // Overlap: e.start < end && e.end > start
        if (e.start < end && e.end > start) {
            hits.push_back(&e);
        }
    }
    return hits;
}

} // namespace branch::analysis
