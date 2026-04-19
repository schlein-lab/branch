#include "analysis/annotation_bed.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <sstream>
#include <stdexcept>

namespace branch::analysis {

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
