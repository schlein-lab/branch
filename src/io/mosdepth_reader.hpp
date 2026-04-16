#pragma once

// BRANCH v0.2 — Minimal mosdepth output reader.
//
// Reads an uncompressed mosdepth *.regions.bed file (4-column BED
// variant: chrom / start / end / name / mean_coverage) into memory.
// For now we take the uncompressed path; gzip support can be added
// without changing the public API.

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
