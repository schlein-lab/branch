#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace branch::consensus {

struct ConsensusConfig {
    // Minimum number of aligned sequences required to emit a base
    // (positions with fewer non-gap symbols are output as 'N').
    std::size_t min_coverage{2};
    // Minimum fraction of aligned non-gap symbols needed for the
    // majority base to be kept; otherwise 'N' is emitted.
    float min_majority_fraction{0.5f};
    // Character used to represent missing data in inputs.
    char gap_char{'-'};
    // Character emitted when coverage or majority is insufficient.
    char ambiguous_char{'N'};
};

struct ConsensusResult {
    std::string sequence;                 // one char per column
    std::vector<std::uint32_t> depths;    // non-gap count per column
    std::vector<float> majority_fraction; // [0,1] per column
};

// Build a per-position majority-base consensus from a set of already-
// aligned sequences (all of equal length). Non-ACGT, non-gap symbols
// are treated as 'N'-equivalents (do not contribute to majority).
[[nodiscard]] ConsensusResult
build_consensus(const std::vector<std::string_view>& aligned_sequences,
                const ConsensusConfig& cfg = {});

}  // namespace branch::consensus
