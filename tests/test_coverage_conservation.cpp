#include <gtest/gtest.h>

#include <vector>

#include "analysis/coverage_conservation.hpp"
#include "detect/bubble_detector.hpp"
#include "graph/lossless_graph.hpp"

using branch::analysis::ConservationConfig;
using branch::analysis::ConservationReport;
using branch::analysis::build_solver_state;
using branch::analysis::reallocate_once;
using branch::analysis::run_conservation;
using branch::analysis::CNEstimate;
using branch::detect::Bubble;
using branch::detect::detect_bubbles;
using branch::graph::LosslessGraph;
using branch::graph::NodeId;

namespace {

// Source → mid → sink with matched read_support. No violations expected.
LosslessGraph make_balanced_linear() {
    LosslessGraph g;
    NodeId a = g.add_node(100);
    NodeId b = g.add_node(100);
    NodeId c = g.add_node(100);
    g.add_edge(a, b, 10);
    g.add_edge(b, c, 10);
    return g;
}

// Source → entry → (alt1, alt2) → exit → sink. Both alts carry 10; entry
// sees 20 in and 20 out. All interior nodes conserve when expected_cn = 1.
LosslessGraph make_balanced_bubble() {
    LosslessGraph g;
    NodeId src   = g.add_node(100);
    NodeId entry = g.add_node(100);
    NodeId alt1  = g.add_node(100);
    NodeId alt2  = g.add_node(100);
    NodeId exit  = g.add_node(100);
    NodeId sink  = g.add_node(100);
    g.add_edge(src, entry, 20);
    g.add_edge(entry, alt1, 10);
    g.add_edge(entry, alt2, 10);
    g.add_edge(alt1, exit, 10);
    g.add_edge(alt2, exit, 10);
    g.add_edge(exit, sink, 20);
    return g;
}

// Same shape as balanced_bubble but with a 2x-duplicated middle: node
// `mid` has expected_cn = 2 and must therefore emit 2x its inflow. If
// our cn_hints align, conservation is preserved; if they don't, the
// solver flags a surplus on the duplicate.
LosslessGraph make_duplication_surplus() {
    LosslessGraph g;
    NodeId src   = g.add_node(100);
    NodeId entry = g.add_node(100);  // inflow 10, outflow 20 → needs CN=2
    NodeId mid1  = g.add_node(100);
    NodeId mid2  = g.add_node(100);
    NodeId exit  = g.add_node(100);
    NodeId sink  = g.add_node(100);
    g.add_edge(src, entry, 10);
    g.add_edge(entry, mid1, 10);
    g.add_edge(entry, mid2, 10);
    g.add_edge(mid1, exit, 10);
    g.add_edge(mid2, exit, 10);
    g.add_edge(exit, sink, 20);
    return g;
}

}  // namespace

TEST(CoverageConservation, NoViolationOnBalancedGraph) {
    auto g = make_balanced_linear();
    auto report = run_conservation(g, {}, {});
    EXPECT_TRUE(report.per_node_violations.empty());
    EXPECT_TRUE(report.per_bubble_violations.empty());
    EXPECT_TRUE(report.converged);
    EXPECT_LT(report.final_global_residual, 1e-3f);
}

TEST(CoverageConservation, BalancedBubbleConserves) {
    auto g = make_balanced_bubble();
    auto bubbles = detect_bubbles(g);
    auto report = run_conservation(g, bubbles, {});
    EXPECT_TRUE(report.per_node_violations.empty());
    EXPECT_TRUE(report.per_bubble_violations.empty());
    EXPECT_TRUE(report.converged);
}

TEST(CoverageConservation, DuplicationSurplusWithoutCNHintIsFlagged) {
    // The entry node has in=10, out=20; without a CN=2 hint this looks
    // like a surplus of +10 on the entry, which is exactly what we
    // want the local pass to flag.
    auto g = make_duplication_surplus();
    auto report = run_conservation(g, {}, {});
    ASSERT_FALSE(report.per_node_violations.empty());
    bool saw_entry_surplus = false;
    for (const auto& v : report.per_node_violations) {
        if (v.residual > 0.0f) saw_entry_surplus = true;
    }
    EXPECT_TRUE(saw_entry_surplus);
}

TEST(CoverageConservation, DuplicationConservesWithCNHint) {
    auto g = make_duplication_surplus();
    // cn_hints is indexed by NodeId; node ids are 0..5 in insertion
    // order from make_duplication_surplus. Entry is at index 1 and
    // legitimately carries CN=2.
    std::vector<CNEstimate> hints(g.node_count());
    for (auto& h : hints) h.relative_cn = 1.0f;
    hints[1].relative_cn = 2.0f;
    auto report = run_conservation(g, {}, hints);
    EXPECT_TRUE(report.per_node_violations.empty())
        << "node_violations=" << report.per_node_violations.size()
        << " residual=" << report.final_global_residual;
}

TEST(CoverageConservation, IterationBoundRespected) {
    auto g = make_duplication_surplus();
    ConservationConfig cfg;
    cfg.max_iterations = 3;
    cfg.convergence_epsilon = 1e-12f;  // force hitting the cap
    cfg.damping = 0.1f;
    auto report = run_conservation(g, {}, {});
    // Solver must never report more iterations than the cap.
    EXPECT_LE(report.iterations_used, 10u);

    ConservationReport r2 = run_conservation(g, {}, {}, cfg);
    EXPECT_LE(r2.iterations_used, cfg.max_iterations);
}

TEST(CoverageConservation, SolverReducesGlobalResidualMonotonically) {
    // With damping < 1 on a finite graph, each reallocation sweep must
    // reduce (or leave unchanged) the global |residual|.
    auto g = make_duplication_surplus();
    auto state = build_solver_state(g, {});
    const float initial = state.global_residual;
    ConservationConfig cfg;
    const float after_one = reallocate_once(state, g, cfg);
    EXPECT_LE(after_one, initial);
}

TEST(CoverageConservation, EmptyGraphIsTriviallyConserved) {
    LosslessGraph g;
    auto report = run_conservation(g, {}, {});
    EXPECT_TRUE(report.converged);
    EXPECT_EQ(report.iterations_used, 0u);
    EXPECT_FLOAT_EQ(report.final_global_residual, 0.0f);
}
