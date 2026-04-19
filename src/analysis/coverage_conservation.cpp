// BRANCH v0.3 — Coverage-Conservation-Axiom implementation.
//
// Algorithm sketch (fixed-point flow balancing):
//
//   1. Build NodeFlow records: sum incoming/outgoing edge read_support,
//      assign expected_cn from cn_hints (fall back to Node.copy_count,
//      then to 1.0). Record per-node residual.
//
//   2. Local pass — snapshot every node / bubble whose |residual| is
//      above the respective tolerance. This is the report that ships
//      without mutating any flow.
//
//   3. Global pass — fixed-point reallocation. On each sweep, distribute
//      a damped fraction of each node's residual back onto the
//      edge_flow[] mirror, proportional to current edge flow. Nodes
//      with in_degree == 0 or out_degree == 0 are sources/sinks and
//      are excluded from residual calculation (their imbalance is
//      expected — it's where reads enter/leave the graph).
//
//   4. Stop when Σ|residual| drops below convergence_epsilon or the
//      iteration cap is reached.
//
// Determinism: iteration order is edge-index order (matches the graph's
// `edges()` vector) and node-id order. No hash maps in the hot path.

#include "analysis/coverage_conservation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace branch::analysis {

using conservation::NodeFlow;
using conservation::SolverState;

namespace {

// Expected-CN lookup: cn_hints is keyed by name / chrom range, not by
// NodeId. We do not have a node-to-region map here, so the v0.3 policy
// is: if a hint's relative_cn is set and the hint list is the same
// length as the node vector, use it positionally; otherwise fall back
// to Node.copy_count. This mirrors the convention used by repeat_cn
// for per-region CN assignment and keeps the solver decoupled from
// chrom coordinates until the mapper lands.
float expected_cn_for(
    const branch::graph::Node& node,
    std::size_t idx,
    const std::vector<CNEstimate>& cn_hints) {
    if (idx < cn_hints.size() && cn_hints[idx].relative_cn > 0.0f) {
        return cn_hints[idx].relative_cn;
    }
    return static_cast<float>(node.copy_count > 0u ? node.copy_count : 1u);
}

// Σ|residual| over non-boundary nodes. Sources (in_degree == 0) and
// sinks (out_degree == 0) are intentionally excluded — their mismatch
// is the point at which reads enter / leave the graph and is not a
// violation of conservation.
float compute_global_residual(const SolverState& state) {
    float sum = 0.0f;
    for (const auto& nf : state.nodes) {
        if (nf.in_degree == 0 || nf.out_degree == 0) continue;
        sum += std::fabs(nf.residual);
    }
    return sum;
}

// Recompute every NodeFlow.inflow/outflow/residual from edge_flow[].
// Used between reallocation sweeps.
void refresh_node_flows(
    SolverState& state,
    const branch::graph::LosslessGraph& graph) {
    for (auto& nf : state.nodes) {
        nf.inflow = 0.0f;
        nf.outflow = 0.0f;
    }
    const auto& edges = graph.edges();
    for (std::size_t ei = 0; ei < edges.size(); ++ei) {
        const auto& e = edges[ei];
        const float f = state.edge_flow[ei];
        if (e.to < state.nodes.size())   state.nodes[e.to].inflow   += f;
        if (e.from < state.nodes.size()) state.nodes[e.from].outflow += f;
    }
    for (auto& nf : state.nodes) {
        // residual = outflow − inflow · expected_cn. A node that
        // legitimately duplicates (expected_cn = 2) must produce 2x
        // its inflow on its outgoing edges to be conservation-balanced.
        nf.residual = nf.outflow - nf.inflow * nf.expected_cn;
    }
    state.global_residual = compute_global_residual(state);
}

}  // namespace

conservation::SolverState build_solver_state(
    const branch::graph::LosslessGraph& graph,
    const std::vector<CNEstimate>& cn_hints) {
    SolverState st;
    st.nodes.resize(graph.node_count());
    st.edge_flow.resize(graph.edge_count(), 0.0f);

    // Initialise per-node CN + id + degrees.
    for (std::size_t i = 0; i < graph.nodes().size(); ++i) {
        const auto& n = graph.nodes()[i];
        auto& nf = st.nodes[i];
        nf.id          = n.id;
        nf.expected_cn = expected_cn_for(n, i, cn_hints);
    }
    for (std::size_t ei = 0; ei < graph.edges().size(); ++ei) {
        const auto& e = graph.edges()[ei];
        st.edge_flow[ei] = static_cast<float>(e.read_support);
        if (e.from < st.nodes.size()) st.nodes[e.from].out_degree++;
        if (e.to   < st.nodes.size()) st.nodes[e.to].in_degree++;
    }

    refresh_node_flows(st, graph);
    return st;
}

float reallocate_once(
    SolverState& state,
    const branch::graph::LosslessGraph& graph,
    const ConservationConfig& cfg) {
    // Strategy: for every interior node with non-trivial residual,
    // push the correction uniformly across its outgoing edges (or pull
    // from incoming edges if residual is negative). Uniform split is
    // the conservative default — weighted-by-current-flow variants
    // biased toward spurious high-count edges in early experiments.
    const auto& edges = graph.edges();
    std::vector<std::vector<std::size_t>> out_edges(state.nodes.size());
    std::vector<std::vector<std::size_t>> in_edges(state.nodes.size());
    for (std::size_t ei = 0; ei < edges.size(); ++ei) {
        const auto& e = edges[ei];
        if (e.from < out_edges.size()) out_edges[e.from].push_back(ei);
        if (e.to   < in_edges.size())  in_edges[e.to].push_back(ei);
    }

    const float damp = std::clamp(cfg.damping, 0.0f, 1.0f);

    for (auto& nf : state.nodes) {
        if (nf.in_degree == 0 || nf.out_degree == 0) continue;
        if (std::fabs(nf.residual) < cfg.local_violation_tolerance) continue;

        const float correction = damp * nf.residual;
        // Positive residual → outflow too high. Shave correction off the
        // outgoing edges uniformly.
        if (correction > 0.0f) {
            const auto& oes = out_edges[nf.id];
            if (oes.empty()) continue;
            const float per_edge = correction / static_cast<float>(oes.size());
            for (std::size_t ei : oes) {
                state.edge_flow[ei] = std::max(0.0f, state.edge_flow[ei] - per_edge);
            }
        } else {
            // Negative residual → outflow too low. Add to outgoing edges
            // (reads that ought to be there but aren't, proportional
            // distribution — uniform in v0.3).
            const auto& oes = out_edges[nf.id];
            if (oes.empty()) continue;
            const float per_edge = -correction / static_cast<float>(oes.size());
            for (std::size_t ei : oes) {
                state.edge_flow[ei] += per_edge;
            }
        }
    }

    refresh_node_flows(state, graph);
    state.iteration++;
    return state.global_residual;
}

ConservationReport run_conservation(
    const branch::graph::LosslessGraph& graph,
    const std::vector<branch::detect::Bubble>& bubbles,
    const std::vector<CNEstimate>& cn_hints,
    const ConservationConfig& cfg) {
    ConservationReport report;

    SolverState state = build_solver_state(graph, cn_hints);

    // --- Local pass: snapshot violations from the pre-solver state. ---
    for (const auto& nf : state.nodes) {
        if (nf.in_degree == 0 || nf.out_degree == 0) continue;
        if (std::fabs(nf.residual) < cfg.local_violation_tolerance) continue;
        report.per_node_violations.push_back(NodeViolation{
            .node        = nf.id,
            .inflow      = nf.inflow,
            .outflow     = nf.outflow,
            .expected_cn = nf.expected_cn,
            .residual    = nf.residual,
        });
    }

    // Bubble-level check: at each (entry, exit) the sum of alt
    // read-supports must match entry.outflow ≈ exit.inflow scaled by
    // the entry's expected_cn. Deviations flag mosaic / ghost paths.
    for (const auto& b : bubbles) {
        if (b.entry >= state.nodes.size() || b.exit >= state.nodes.size()) continue;
        const auto& entry = state.nodes[b.entry];
        const auto& exit_n = state.nodes[b.exit];
        float alt_sum = 0.0f;
        for (const auto& a : b.alts) alt_sum += static_cast<float>(a.total_read_support);
        // Only residual vs entry_outflow scaled by CN is the reliable
        // bubble-level signal — alt_sum double-counts when alt paths
        // traverse further interior nodes, and bubble_detector can
        // emit overlapping (entry, deep-exit) pairs by construction.
        const float residual = exit_n.inflow - entry.outflow * entry.expected_cn;
        if (std::fabs(residual) < cfg.bubble_violation_tolerance) {
            continue;
        }
        report.per_bubble_violations.push_back(BubbleViolation{
            .entry         = b.entry,
            .exit          = b.exit,
            .entry_outflow = entry.outflow,
            .exit_inflow   = exit_n.inflow,
            .alt_flow_sum  = alt_sum,
            .residual      = residual,
        });
    }

    // --- Global pass: iterate until convergence or cap. ---
    float prev = state.global_residual;
    for (std::uint32_t it = 0; it < cfg.max_iterations; ++it) {
        if (prev < cfg.convergence_epsilon) {
            report.converged = true;
            break;
        }
        const float now = reallocate_once(state, graph, cfg);
        if (!(now < prev)) {
            // No progress — treat as converged to avoid wasted sweeps
            // on oscillating subgraphs. `damping` < 1 usually prevents
            // this but we guard anyway.
            report.converged = now < cfg.convergence_epsilon;
            break;
        }
        prev = now;
    }
    // Final check if we fell out of the loop by exhausting iterations.
    if (state.global_residual < cfg.convergence_epsilon) {
        report.converged = true;
    }

    report.iterations_used       = state.iteration;
    report.final_global_residual = state.global_residual;

    // Deterministic ordering of output vectors.
    std::sort(report.per_node_violations.begin(), report.per_node_violations.end(),
              [](const NodeViolation& a, const NodeViolation& b) {
                  return a.node < b.node;
              });
    std::sort(report.per_bubble_violations.begin(), report.per_bubble_violations.end(),
              [](const BubbleViolation& a, const BubbleViolation& b) {
                  if (a.entry != b.entry) return a.entry < b.entry;
                  return a.exit < b.exit;
              });

    return report;
}

}  // namespace branch::analysis
