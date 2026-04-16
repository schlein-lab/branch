#pragma once

// BRANCH v0.1 — Delta-Graph Read representation.
//
// A Read in BRANCH is NOT stored as an independent sequence. After graph
// integration, every read becomes a path through the string-graph plus a
// small list of positional differences (deltas) from the path consensus.
// This is the single most important design decision separating BRANCH
// from classical string-graph assemblers like hifiasm.

#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <vector>

namespace branch::graph {

using NodeId = std::uint32_t;
using ReadId = std::uint32_t;

// A positional delta of a single read from its graph path consensus.
// Deltas are read-local: offset is measured from the start of the read's
// concatenated path, not from genome coordinates.
struct PositionalDelta {
    std::uint32_t offset;           // read-local position (0..read_length)
    std::uint8_t length;            // length of the replacement (1..255)
    std::uint8_t alt_base;          // 2-bit packed base for length==1
    std::uint16_t alt_bases_index;  // index into external alt-bases pool
                                    // for multi-base insertions; 0 == inline
};

static_assert(sizeof(PositionalDelta) == 8,
              "PositionalDelta must pack to 8 bytes for cache efficiency");

// A ReadPath describes how a single HiFi read traverses the graph.
// The deltas encode deviations from the concatenated consensus of the
// traversed nodes: typically sequencing errors (<1% on HiFi) and rare
// true variants (branch points where the read diverges).
struct ReadPath {
    ReadId read_id{};
    std::uint32_t read_length{};
    std::vector<NodeId> path;
    std::vector<PositionalDelta> deltas;
};

// Compact view used in hot paths. Backed by contiguous arrays in
// a Graph-owned storage pool; no allocations per read access.
struct ReadPathView {
    ReadId read_id{};
    std::uint32_t read_length{};
    std::span<const NodeId> path;
    std::span<const PositionalDelta> deltas;
};

}  // namespace branch::graph
