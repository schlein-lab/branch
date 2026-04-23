#pragma once

// BRANCH v0.2 — Minimal mosdepth output reader.
//
// Reads a mosdepth regions file into memory. Supports both plain
// `.bed` and gzip-compressed `.bed.gz` (auto-detected by suffix).
// Expected layout: 4-column BED (chrom / start / end / mean_coverage)
// or 5-column (chrom / start / end / name / mean_coverage).

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace branch::io {

struct RegionCoverage {
    std::string chrom;
    std::uint32_t start{};
    std::uint32_t end{};
    std::string name;
    float mean_coverage{};
};

// Parse mosdepth regions.bed text lines from the given file path.
// Returns empty vector on read failure. Order is preserved.
[[nodiscard]] std::vector<RegionCoverage>
read_mosdepth_regions(const std::string& path);

}  // namespace branch::io
