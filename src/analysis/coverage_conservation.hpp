#pragma once

// BRANCH v0.3 — Coverage-Conservation-Axiom (P1.1).
//
// Read-flow conservation on the lossless assembly graph. At every
// branch point we require:
//
//   Σ(incoming read_support)  ≈  Σ(outgoing read_support / expected_cn)
//
// Violations of this identity are the primary signal for low-frequency
// CNV: if a 2x-duplicated segment carries only 1x the expected flow,
// the duplication is incomplete (mosaic) or the branch is spurious.
//
// Two passes are run in sequence:
//   1. Local  — per-node and per-bubble check, flags deficit / surplus.
//   2. Global — iterative reallocation until Σ|residual| converges or
//               the iteration cap is hit (default 10). This is the
//               fixed-point stage; it will not mutate the input graph.
//
// The solver is stateless with respect to the graph: it reads
// `Node.read_support`, `Edge.read_support`, and an optional
// per-node CN hint (from `analysis::coverage_cn::CNEstimate`), and
// emits a `ConservationReport` summarising violations and solver
// status. Callers integrate the report downstream (classifier,
// VAF-band, repeat_cn).

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "analysis/coverage_cn.hpp"
#include "analysis/conservation_state.hpp"
#include "detect/bubble_detector.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::analysis {

// Tunables for the conservation solver. Every threshold that gates a
// decision in the hot path lives here; no magic numbers downstream.
struct ConservationConfig {
    // Max fixed-point iterations. Safety cap to keep worst case O(V·k).
    std::uint32_t max_iterations{10};

    // Convergence threshold on Σ|residual| (absolute). Solver stops
    // once the global residual falls below this value.
    float convergence_epsilon{1e-3f};

    // Minimum |residual| to flag a node locally as a violation. Below
    // this we treat the node as conserved and drop it from the report.
    float local_violation_tolerance{0.5f};

    // Minimum |bubble_residual| to flag a whole bubble. Bubbles sum
    // alt-path supports, so tolerance is looser than per-node.
    float bubble_violation_tolerance{1.0f};

    // Relative damping factor for global reallocation. In [0, 1]; 0
    // freezes flows, 1 applies the full correction each step. Smaller
    // values stabilise oscillations on near-cyclic subgraphs.
    float damping{0.5f};
};

// One flagged non-conservation at a single node. Sign of `residual`
// indicates direction:
//   residual > 0  → outflow exceeds expected inflow*CN (surplus)
//   residual < 0  → outflow is short of expected inflow*CN (deficit)
struct NodeViolation {
    branch::graph::NodeId node;
    float inflow{0.0f};
    float outflow{0.0f};
    float expected_cn{1.0f};
    float residual{0.0f};
};

// Conservation accounting for one bubble. `alt_flow_sum` aggregates
// all alt-path supports; entry/exit inflow come from the parent graph.
struct BubbleViolation {
    branch::graph::NodeId entry;
    branch::graph::NodeId exit;
    float entry_outflow{0.0f};
    float exit_inflow{0.0f};
    float alt_flow_sum{0.0f};
    float residual{0.0f};           // exit_inflow − entry_outflow (scaled)
};

struct ConservationReport {
    std::vector<NodeViolation> per_node_violations;
    std::vector<BubbleViolation> per_bubble_violations;
    std::uint32_t iterations_used{0};
    bool converged{false};
    float final_global_residual{0.0f};
};

// Public entry point. Runs local + global passes on the graph and
// returns the report. `bubbles` may be empty (skip bubble-level pass);
// `cn_hints` is indexed by NodeId — entries outside its size are
// treated as expected_cn = 1.0.
[[nodiscard]] ConservationReport run_conservation(
    const branch::graph::LosslessGraph& graph,
    const std::vector<branch::detect::Bubble>& bubbles,
    const std::vector<CNEstimate>& cn_hints,
    const ConservationConfig& cfg = {});

// Exposed for tests: builds the initial solver state without running
// any iterations. Callers that want to inspect pre-reallocation flows
// call this directly.
[[nodiscard]] conservation::SolverState
build_solver_state(
    const branch::graph::LosslessGraph& graph,
    const std::vector<CNEstimate>& cn_hints);

// Exposed for tests: run a single Gauss-Seidel-style reallocation
// sweep on an already-initialised state. Returns the updated global
// residual. No-op on graphs with no violations.
float reallocate_once(
    conservation::SolverState& state,
    const branch::graph::LosslessGraph& graph,
    const ConservationConfig& cfg);

}  // namespace branch::analysis
