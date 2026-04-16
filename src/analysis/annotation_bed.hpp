#pragma once
#include <string>
#include <vector>
#include <map>

namespace branch::analysis {

struct BedEntry {
    std::string chrom;
    std::string name;
    size_t start;
    size_t end;
};

std::vector<BedEntry> load_bed(const std::string& path);

class IntervalIndex {
public:
    explicit IntervalIndex(const std::vector<BedEntry>& entries);
    std::vector<const BedEntry*> query(const std::string& chrom, size_t start, size_t end) const;
    bool empty() const { return by_chrom_.empty(); }

private:
    std::map<std::string, std::vector<BedEntry>> by_chrom_;
};

} // namespace branch::analysis
