#pragma once

// BRANCH v0.1 — Graph builder from read overlaps.
//
// Consumes minimizer-sketched reads and a set of computed overlap
// pairs, emits a LosslessGraph. The builder is intentionally minimal
// in v0.1: it allocates one node per read (no unitig compaction yet)
// and one directed edge per overlap pair. Unitig compaction, bubble
// detection, and CN-aware merging arrive in subsequent modules.
//
// The builder is a plain function rather than a stateful class so it
// can be invoked from any backend (CPU or GPU) after the backend has
// produced its overlap candidates.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "graph/kmer_sketch.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::graph {

struct ReadMeta {
    ReadId read_id;
    std::uint32_t length_bp;
};

// Build a LosslessGraph from a list of read-metas and overlap pairs.
// Allocates one node per read (in the order they appear in `reads`,
// so node_id == reads[i].read_id only if IDs are 0..N-1 dense — no
// assumption is made; returned map gives the translation).
struct BuildResult {
    LosslessGraph graph;
    std::vector<NodeId> read_to_node;  // indexed by read_id
};

BuildResult build_graph(std::span<const ReadMeta> reads,
                        std::span<const branch::backend::OverlapPair> overlaps);

}  // namespace branch::graph
