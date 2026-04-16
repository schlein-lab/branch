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

// Forward-declare the backend OverlapPair so this header does not
// pull in the backend target (avoids a branch_graph <-> branch_backend
// CMake cycle: branch_backend depends on branch_graph).
namespace branch::backend { struct OverlapPair; }

namespace branch::graph {

class LosslessGraph;

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

// Populate read paths from overlap-based graph construction.
//
// v0.2 best-effort mapping: for each read, emit a ReadPath containing
// its own node and the nodes of all reads it overlaps with, ordered by
// the overlap's offset along this read. `deltas` is left empty because
// Node currently stores only `length_bp` (no consensus sequence). This
// is sufficient for:
//   - Debugging: which reads landed on which branch of a bubble.
//   - Downstream VAF: count reads per bubble branch.
//
// TODO(v0.3): compute_delta(ref, read_sub) once Node has consensus.
//
// Parameters:
//   n_reads      Number of input reads. Output vector has n_reads entries
//                in read_id order; reads with no overlaps get a ReadPath
//                containing only their own node (still useful downstream).
//   read_to_node Mapping from read_id -> NodeId as produced by build_graph.
//   read_lengths Read lengths in bp, indexed by read_id. May be empty; if
//                empty the ReadPath.read_length is set to 0.
//   overlaps     The same overlap pairs passed to build_graph.
//   out          Output read paths, one per read_id in [0, n_reads).
void populate_read_paths(
    std::size_t n_reads,
    std::span<const NodeId> read_to_node,
    std::span<const std::uint32_t> read_lengths,
    std::span<const branch::backend::OverlapPair> overlaps,
    std::vector<ReadPath>& out);

}  // namespace branch::graph
