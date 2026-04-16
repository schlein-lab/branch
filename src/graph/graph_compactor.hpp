#pragma once

// BRANCH v0.2 — Unitig compaction.
//
// Collapses linear chains of nodes (in-degree ≤ 1, out-degree ≤ 1,
// single successor has in-degree 1) into single unitig nodes. The
// result is a graph with the same branching structure but vastly
// fewer nodes — essential for downstream bubble detection and
// classifier work on real-world data.
//
// v0.2 ships the simple/safe variant:
//   - Respects existing branch/merge points.
//   - Preserves per-edge read_support by summing along the chain.
//   - Node length_bp is summed along the chain.
//   - Node copy_count is taken from the chain start (CN-aware
//     inference happens later, in a separate pass).
//
// The compactor returns a new LosslessGraph; the input is not
// modified. A mapping from old NodeId to new NodeId is returned so
// external structures (e.g. ReverseIndex) can be rewritten.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "graph/lossless_graph.hpp"

namespace branch::graph {

struct CompactionResult {
    LosslessGraph compacted;
    // old_to_new[old_id] is the id of the unitig that contains old_id
    // in the new graph.
    std::vector<NodeId> old_to_new;
};

[[nodiscard]] CompactionResult compact_unitigs(const LosslessGraph& input);

}  // namespace branch::graph
