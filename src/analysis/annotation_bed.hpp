#pragma once
#include <array>
#include <cstdint>
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

    // P2.4: per-alt cumulative read support (one entry per bubble alt path).
    // Empty vector means "no per-alt support column present".
    // Values are summed edge supports along each alt's reconstructed path
    // (see detect::Bubble::alts[].total_read_support) — a proxy for reads
    // traversing the alt, not a bottleneck. Sufficient to rank alts by
    // relative depth and compute a VAF estimate from min/sum.
    std::vector<std::uint32_t> alt_read_supports;

    // Per-alt haplotype-resolved read counts. Outer index = alt, inner
    // index = haplotype (0 = hap0, 1 = hap1). Populated when analyze
    // was given --reads <gaf> AND phasing produced a phase block that
    // covers the bubble. Empty otherwise.
    std::vector<std::array<std::uint32_t, 2>> alt_read_supports_phased;
};

std::vector<BedEntry> load_bed(const std::string& path);

// Serialise a BED entry. Columns:
//   1 chrom, 2 start, 3 end, 4 name,
//   5 confidence (3 decimals, "." if unset),
//   6 alt_read_supports (comma-sep uint32, "." if empty),
//   7 min_vaf = min(alt_supports)/sum(alt_supports), 4 decimals, "." if <2 alts,
//   8 total_alt_support = sum(alt_supports), "." if empty.
//   9 alt_read_supports_phased (per alt: "hap0|hap1" pairs joined by ',',
//     e.g. "12|3,2|4" — empty/"." when phasing unavailable).
//  10 hap_skew = max(per-alt |hap0 - hap1|) / total_phased_reads, 4
//     decimals, "." when fewer than 2 phased reads in any alt.
// Columns 6-8 only emitted when alt_read_supports is non-empty;
// columns 9-10 only emitted when alt_read_supports_phased is non-empty.
// Older 5-column readers stay compatible.
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
