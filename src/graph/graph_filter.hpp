#pragma once

// BRANCH v0.2 — Graph filtering: containment drop + transitive reduction.
//
// The graph builder emits one directed edge per OverlapPair. For real
// long-read datasets this is O(reads^2) edges before any structural
// simplification. Two classical passes make the graph tractable:
//
//   1) Containment filter: a read that lies wholly inside another read
//      adds no new information — its node and all touching edges are
//      dropped. In v0.2 we approximate "contained" as: node N has an
//      incoming edge from some other node M with length_bp(M) >=
//      length_bp(N). This is a conservative heuristic because the real
//      v0.1 Edge carries no explicit overlap offsets; once the backend
//      emits offset-annotated edges the predicate can be tightened.
//
//   2) Transitive reduction: the classical Myers string-graph pass. If
//      A -> B and B -> C exist, and A -> C exists with cumulative
//      length matching len(A->B) + len(B->C) within a fuzz tolerance,
//      remove the direct A -> C edge. v0.2 uses node length_bp as the
//      per-step length proxy; an offset-carrying edge format can be
//      swapped in without changing the reduction algorithm.
//
// Both passes mutate the LosslessGraph in place by rewriting its edge
// vector (and, for the containment pass, dropping incident edges of
// dead nodes). Node identities are preserved — downstream tables such
// as ReverseIndex and read_to_node do not need to be rewritten.

#include <cstddef>
#include <cstdint>

#include "graph/lossless_graph.hpp"

namespace branch::graph {

struct FilterConfig {
    bool drop_contained{true};
    bool reduce_transitive{true};
    // Fuzz tolerance for transitive reduction (bp; default 100).
    std::uint32_t transitive_fuzz{100};
};

struct FilterStats {
    std::size_t edges_before{};
    std::size_t edges_after{};
    std::size_t nodes_dropped_contained{};
    std::size_t transitive_edges_removed{};
};

// Apply containment + transitive-reduction filters to `graph` in place.
// Returns before/after stats. Safe to call on an empty graph.
[[nodiscard]] FilterStats filter_graph(
    LosslessGraph& graph,
    const FilterConfig& cfg = {});

}  // namespace branch::graph
