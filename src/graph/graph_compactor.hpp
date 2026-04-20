#pragma once

// BRANCH v0.2 — Unitig compaction with consensus building.
//
// Collapses linear chains of nodes (in-degree ≤ 1, out-degree ≤ 1,
// single successor has in-degree 1) into single unitig nodes, then
// runs two structural merge passes that handle cases linear-chain
// compaction cannot reach:
//
//   1. Closed-twin merge. Two singleton nodes with identical closed
//      undirected neighborhood AND identical length_bp are merged
//      into a single unitig. In a read-overlap graph these are
//      reads from the SAME allele covering the SAME locus — they
//      form a tournament/clique where every member has in>1 and
//      out>1 (so the linear-chain pass cannot touch them).
//
//   2. Isolated-node drop. Nodes untouched by any edge are residue
//      from graph_filter's containment drop; they carry no topology
//      and are not emitted as unitigs.
//
// The result is a graph with the same branching structure but vastly
// fewer nodes — essential for downstream bubble detection and
// classifier work on real-world data.
//
// v0.2 ships the simple/safe variant:
//   - Respects existing branch/merge points.
//   - Preserves per-edge read_support on inter-unitig edges.
//   - Node length_bp is summed along the chain; closed-twin clusters
//     take the representative (shared) length.
//   - Node copy_count is taken from the chain start (CN-aware
//     inference happens later, in a separate pass).
//   - Node consensus is built from member sequences via MSA.
//
// The compactor returns a new LosslessGraph; the input is not
// modified. A mapping from old NodeId to new NodeId is returned so
// external structures (e.g. ReverseIndex) can be rewritten. Dropped
// isolated nodes map to `std::numeric_limits<NodeId>::max()` —
// consumers must treat that value as "dropped from the compacted
// graph".

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "graph/lossless_graph.hpp"

namespace branch::graph {

struct CompactionResult {
    LosslessGraph compacted;
    // old_to_new[old_id] is the id of the unitig that contains old_id
    // in the new graph.
    std::vector<NodeId> old_to_new;
};

// Compact unitigs, building consensus from Node::consensus fields.
[[nodiscard]] CompactionResult compact_unitigs(const LosslessGraph& input);

// Compact unitigs with explicit sequences for consensus building.
// node_sequences[i] is the sequence for node i (may be empty).
[[nodiscard]] CompactionResult compact_unitigs_with_sequences(
    const LosslessGraph& input,
    const std::vector<std::string>& node_sequences);

}  // namespace branch::graph
