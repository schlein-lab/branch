#pragma once

// BRANCH v0.2 — Bubble detection in a LosslessGraph.
//
// A "bubble" is a subgraph where exactly one entry node fans out into
// multiple parallel paths that reconverge at a single exit node. In
// BRANCH terms every bubble is a classification candidate: its
// parallel alternatives are what downstream cascade has to label as
// branch, duplication, mixed, or non-separable.
//
// v0.2 ships a structural detector only — it enumerates bubbles by
// local graph topology. It does not yet extract path sequences or
// coverage per alt-path; those arrive with the CN-aware pass.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "graph/lossless_graph.hpp"

namespace branch::detect {

struct AltPath {
    std::vector<branch::graph::NodeId> nodes;  // intermediate nodes between entry and exit
    std::uint32_t total_read_support{};
};

struct Bubble {
    branch::graph::NodeId entry;
    branch::graph::NodeId exit;
    std::vector<AltPath> alts;
    std::uint32_t total_read_support{};
};

struct BubbleDetectorConfig {
    // Maximum number of intermediate nodes per alt path. Longer paths
    // are more likely structural variation than true bubbles; v0.2
    // hardcodes a conservative cap so detection stays O(V) per entry.
    std::uint32_t max_alt_path_length{8};
};

// Detect bubbles in the graph. Returns a vector of Bubble records,
// each describing one entry+exit pair and the parallel alt paths
// between them.
[[nodiscard]] std::vector<Bubble> detect_bubbles(
    const branch::graph::LosslessGraph& graph,
    const BubbleDetectorConfig& cfg = {});

}  // namespace branch::detect
