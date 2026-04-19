#pragma once

// BRANCH v0.3 — Conservation-state container.
//
// Internal working state for the coverage-conservation solver. Kept in
// its own header so tests and the solver share a single definition and
// the public API in `coverage_conservation.hpp` stays compact.
//
// The solver treats each node as having an "inflow" (Σ read_support on
// incoming edges) and an "outflow" (Σ read_support on outgoing edges).
// For a read-flow-conserving assembly these must match, modulo source /
// sink nodes and scaling by the node's expected copy number: a 2x
// duplication legitimately doubles outflow.
//
// NOTE: state is scratch memory. The public ConservationReport is the
// stable artifact; this struct only leaks out for unit tests that want
// to inspect the intermediate flow field.

#include <cstdint>
#include <vector>

#include "graph/lossless_graph.hpp"

namespace branch::analysis::conservation {

// Per-node flow accounting. Values are floats because the iterative
// reallocation pass produces non-integer intermediate estimates.
struct NodeFlow {
    branch::graph::NodeId id{};
    float inflow{0.0f};         // Σ incoming read_support
    float outflow{0.0f};        // Σ outgoing read_support
    float expected_cn{1.0f};    // expected copy number (1 = haploid)
    float residual{0.0f};       // outflow − inflow · expected_cn
    std::uint32_t in_degree{0};
    std::uint32_t out_degree{0};
};

// Mutable flow state for a single solver iteration. `edge_flow` mirrors
// the graph's edge vector 1:1; reallocation updates these floats rather
// than touching the graph's integer read_support field.
struct SolverState {
    std::vector<NodeFlow> nodes;
    std::vector<float> edge_flow;   // indexed by edge position
    float global_residual{0.0f};    // Σ |residual| across non-boundary nodes
    std::uint32_t iteration{0};
};

}  // namespace branch::analysis::conservation
