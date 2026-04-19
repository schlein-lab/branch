#pragma once
#include <map>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

namespace branch::analysis {

struct BedEntry {
    std::string chrom;
    std::string name;
    size_t start;
    size_t end;
    // P1.2: per-bubble classifier confidence in [0, 1]. NaN / unset means
    // "no confidence column present"; serialisation falls back to a dot
    // ('.') placeholder, matching UCSC BED conventions.
    std::optional<float> confidence;
};

std::vector<BedEntry> load_bed(const std::string& path);

// P1.2: serialise a BED entry with optional confidence as a 5- or 4-
// column BED line: "<chrom>\t<start>\t<end>\t<name>[\t<confidence>]".
// Confidence is written with 3 decimals; absent -> ".".
void write_bed_entry(std::ostream& out, const BedEntry& entry);

class IntervalIndex {
public:
    explicit IntervalIndex(const std::vector<BedEntry>& entries);
    std::vector<const BedEntry*> query(const std::string& chrom, size_t start, size_t end) const;
    bool empty() const { return by_chrom_.empty(); }

private:
    std::map<std::string, std::vector<BedEntry>> by_chrom_;
};

} // namespace branch::analysis
