#include "consensus/majority_voter.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace branch::consensus {

namespace {

// Map an input character to an index into a 4-slot A/C/G/T counter,
// or -1 if it is a gap / ambiguous symbol that must not contribute
// to the majority vote.
constexpr int base_index(char c) noexcept {
    switch (c) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

constexpr char index_base(int i) noexcept {
    constexpr char kBases[4] = {'A', 'C', 'G', 'T'};
    return kBases[i];
}

}  // namespace

ConsensusResult
build_consensus(const std::vector<std::string_view>& aligned_sequences,
                const ConsensusConfig& cfg) {
    ConsensusResult result;

    if (aligned_sequences.empty()) {
        return result;
    }

    const std::size_t ncols = aligned_sequences.front().size();
    for (const auto& s : aligned_sequences) {
        if (s.size() != ncols) {
            throw std::invalid_argument(
                "build_consensus: aligned sequences must all have equal length");
        }
    }

    result.sequence.resize(ncols, cfg.ambiguous_char);
    result.depths.assign(ncols, 0u);
    result.majority_fraction.assign(ncols, 0.0f);

    for (std::size_t col = 0; col < ncols; ++col) {
        std::array<std::uint32_t, 4> counts{0u, 0u, 0u, 0u};
        std::uint32_t non_gap_total = 0;

        for (const auto& s : aligned_sequences) {
            const char c = s[col];
            if (c == cfg.gap_char) {
                continue;  // gap: does not count toward coverage
            }
            const int idx = base_index(c);
            if (idx < 0) {
                // Non-ACGT, non-gap symbol (e.g., 'N', 'R'): treat as
                // ambiguous — skip per contract ("do not contribute to
                // majority").
                continue;
            }
            ++counts[static_cast<std::size_t>(idx)];
            ++non_gap_total;
        }

        result.depths[col] = non_gap_total;

        if (non_gap_total == 0) {
            result.majority_fraction[col] = 0.0f;
            result.sequence[col] = cfg.ambiguous_char;
            continue;
        }

        int best_idx = 0;
        std::uint32_t best_count = counts[0];
        for (int i = 1; i < 4; ++i) {
            if (counts[static_cast<std::size_t>(i)] > best_count) {
                best_count = counts[static_cast<std::size_t>(i)];
                best_idx = i;
            }
        }

        const float frac = static_cast<float>(best_count) /
                           static_cast<float>(non_gap_total);
        result.majority_fraction[col] = frac;

        if (non_gap_total < cfg.min_coverage) {
            result.sequence[col] = cfg.ambiguous_char;
            continue;
        }

        if (frac < cfg.min_majority_fraction) {
            result.sequence[col] = cfg.ambiguous_char;
            continue;
        }

        // Detect tie at the top: if another base equals best_count,
        // the column is ambiguous regardless of the fraction threshold.
        bool tie = false;
        for (int i = 0; i < 4; ++i) {
            if (i == best_idx) {
                continue;
            }
            if (counts[static_cast<std::size_t>(i)] == best_count) {
                tie = true;
                break;
            }
        }
        if (tie) {
            result.sequence[col] = cfg.ambiguous_char;
            continue;
        }

        result.sequence[col] = index_base(best_idx);
    }

    return result;
}

}  // namespace branch::consensus
