// BRANCH v0.2 — Graph filter implementation.
//
// Two passes:
//   (1) Containment drop. Conservative heuristic: a node N is
//       "contained" if it has at least one incoming edge from a node M
//       where length_bp(M) >= length_bp(N) AND M also has an edge to
//       at least one other successor of N (i.e. M already witnesses
//       everything N witnesses). Without per-edge overlap offsets this
//       is the tightest predicate we can safely apply; it still
//       removes the obvious case of a short read wholly engulfed by a
//       longer one. Dropped nodes are not erased from the node vector
//       (that would renumber NodeIds); instead, every edge touching a
//       dropped node is removed.
//
//   (2) Transitive reduction (Myers 2005, "The fragment assembly
//       string graph"). For each node A, look at its successors B.
//       If there exists a two-hop path A -> B -> C whose cumulative
//       length matches a direct edge A -> C within `transitive_fuzz`
//       bp, remove the direct edge. The length proxy in v0.2 is
//       node length_bp; a later revision that attaches explicit
//       overlap offsets to Edge can swap the metric in one place.
//
// Both passes operate on the same `LosslessGraph::replace_edges` bulk
// mutator so that the node vector — and therefore all external
// NodeId-indexed structures — stays stable.

#include "graph/graph_filter.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

namespace branch::graph {

namespace {

// Build per-node outgoing-edge index lists into `edges`. Returns a
// CSR-style pair: `offsets` of size N+1, `indices` of size edges.size()
// where edges[indices[offsets[v] .. offsets[v+1]]] are v's out-edges.
struct OutAdjacency {
    std::vector<std::size_t> offsets;
    std::vector<std::size_t> indices;
};

OutAdjacency build_out_adjacency(std::size_t n,
                                 const std::vector<Edge>& edges) {
    OutAdjacency adj;
    adj.offsets.assign(n + 1, 0);
    for (const auto& e : edges) {
        if (e.from < n) ++adj.offsets[e.from + 1];
    }
    for (std::size_t i = 1; i <= n; ++i) {
        adj.offsets[i] += adj.offsets[i - 1];
    }
    adj.indices.assign(edges.size(), 0);
    std::vector<std::size_t> cursor(adj.offsets.begin(), adj.offsets.end() - 1);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        if (e.from < n) {
            adj.indices[cursor[e.from]++] = i;
        }
    }
    return adj;
}

// Pass 1: mark contained nodes.
//
// v0.2 (without per-edge overlap offsets) predicate: a node N is
// contained iff there exists a predecessor M such that
//   (i)  length_bp(M) > length_bp(N)                 (strictly longer)
//   (ii) every direct successor of N is also a direct successor of M
//        (M already witnesses all of N's forward context; N is
//        redundant).
// A leaf N (out_degree == 0) trivially satisfies (ii).
//
// Condition (i) is STRICT so two equally-long sibling successors of a
// branching parent are never mutually contained. When edges gain
// explicit overlap offsets in a later revision, condition (ii) can be
// tightened to also check offset agreement; the shape of the
// predicate remains.
std::vector<std::uint8_t> mark_contained(const LosslessGraph& g,
                                         const std::vector<Edge>& edges) {
    const std::size_t n = g.node_count();
    std::vector<std::uint8_t> dropped(n, 0);
    if (n == 0) return dropped;

    std::vector<std::vector<std::size_t>> in_edges(n);
    std::vector<std::vector<std::size_t>> out_edges(n);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        if (e.from < n && e.to < n) {
            in_edges[e.to].push_back(i);
            out_edges[e.from].push_back(i);
        }
    }

    auto has_out_edge = [&](NodeId pred, NodeId succ) -> bool {
        for (std::size_t ei : out_edges[pred]) {
            if (edges[ei].to == succ) return true;
        }
        return false;
    };

    for (NodeId v = 0; v < static_cast<NodeId>(n); ++v) {
        if (in_edges[v].empty()) continue;
        const std::uint32_t vlen = g.node(v).length_bp;

        for (std::size_t eidx : in_edges[v]) {
            NodeId pred = edges[eidx].from;
            if (pred == v) continue;                     // self-loop
            if (dropped[pred]) continue;                 // stale predecessor
            const std::uint32_t plen = g.node(pred).length_bp;
            if (plen <= vlen) continue;                  // condition (i) strict

            // Condition (ii): every successor of v is also a successor of pred.
            bool covers_all = true;
            for (std::size_t vo : out_edges[v]) {
                NodeId w = edges[vo].to;
                if (w == pred) continue;                 // pred<->v cycle edge
                if (!has_out_edge(pred, w)) {
                    covers_all = false;
                    break;
                }
            }
            if (covers_all) {
                dropped[v] = 1;
                break;
            }
        }
    }
    return dropped;
}

// Pass 2: classical Myers transitive reduction with a length-sum fuzz.
// Works on the already-filtered `edges` vector. Returns a new edge
// vector with transitive edges removed, and writes the count of
// removed edges to *out_removed.
std::vector<Edge> reduce_transitive_edges(const LosslessGraph& g,
                                          const std::vector<Edge>& edges,
                                          std::uint32_t fuzz,
                                          std::size_t* out_removed) {
    *out_removed = 0;
    const std::size_t n = g.node_count();
    if (n == 0 || edges.empty()) return edges;

    OutAdjacency adj = build_out_adjacency(n, edges);

    // Mark, per edge, whether it is a transitive (redundant) edge.
    std::vector<std::uint8_t> is_transitive(edges.size(), 0);

    // For quick "direct A->C length" lookup for a fixed source A, we
    // build a flat "direct_to_len" table keyed by destination NodeId.
    // We only need it to contain A's direct successors at any time.
    // Using sentinel max() to mean "no direct edge".
    const std::uint64_t kNoEdge = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::uint64_t> direct_to_len(n, kNoEdge);
    std::vector<std::size_t>   direct_to_edge(n, 0);

    for (NodeId a = 0; a < static_cast<NodeId>(n); ++a) {
        const std::size_t o_begin = adj.offsets[a];
        const std::size_t o_end   = adj.offsets[a + 1];
        if (o_end - o_begin < 2) continue;  // need >=2 successors to reduce

        // Populate direct_to_len for every successor of A.
        for (std::size_t k = o_begin; k < o_end; ++k) {
            const Edge& e = edges[adj.indices[k]];
            // Per-step length proxy: length_bp of the target node.
            // (With offset-aware edges this would be the overlap
            // "advance" distance instead.)
            direct_to_len[e.to] = g.node(e.to).length_bp;
            direct_to_edge[e.to] = adj.indices[k];
        }

        // For each pair (A->B, B->C), test if A has a direct A->C
        // whose length is within [ len(A->B) + len(B->C) - fuzz,
        //                         len(A->B) + len(B->C) + fuzz ].
        // If so, mark the direct A->C as transitive.
        for (std::size_t kb = o_begin; kb < o_end; ++kb) {
            const Edge& eab = edges[adj.indices[kb]];
            NodeId b = eab.to;
            const std::uint64_t len_ab = g.node(b).length_bp;

            const std::size_t bo_begin = adj.offsets[b];
            const std::size_t bo_end   = adj.offsets[b + 1];
            for (std::size_t kc = bo_begin; kc < bo_end; ++kc) {
                const Edge& ebc = edges[adj.indices[kc]];
                NodeId c = ebc.to;
                if (c == a) continue;       // no self-reduction
                if (c == b) continue;       // defensive

                if (direct_to_len[c] == kNoEdge) continue;  // no direct A->C
                const std::uint64_t len_bc   = g.node(c).length_bp;
                const std::uint64_t sum      = len_ab + len_bc;
                const std::uint64_t direct   = direct_to_len[c];
                const std::uint64_t lo = sum > fuzz ? sum - fuzz : 0;
                const std::uint64_t hi = sum + fuzz;
                if (direct >= lo && direct <= hi) {
                    std::size_t ei = direct_to_edge[c];
                    if (!is_transitive[ei]) {
                        is_transitive[ei] = 1;
                    }
                }
            }
        }

        // Clear direct_to_len entries we populated (keep the vector
        // allocated between iterations for cache reuse).
        for (std::size_t k = o_begin; k < o_end; ++k) {
            const Edge& e = edges[adj.indices[k]];
            direct_to_len[e.to] = kNoEdge;
        }
    }

    std::vector<Edge> out;
    out.reserve(edges.size());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        if (is_transitive[i]) {
            ++(*out_removed);
            continue;
        }
        out.push_back(edges[i]);
    }
    return out;
}

}  // namespace

FilterStats filter_graph(LosslessGraph& graph, const FilterConfig& cfg) {
    FilterStats stats{};
    stats.edges_before = graph.edge_count();

    if (graph.node_count() == 0 || graph.edge_count() == 0) {
        stats.edges_after = graph.edge_count();
        return stats;
    }

    // Snapshot current edges into a working buffer.
    std::vector<Edge> working(graph.edges().begin(), graph.edges().end());

    // -------- Pass 1: containment drop --------
    if (cfg.drop_contained) {
        auto dropped = mark_contained(graph, working);

        std::size_t drop_count = 0;
        for (auto flag : dropped) if (flag) ++drop_count;
        stats.nodes_dropped_contained = drop_count;

        if (drop_count > 0) {
            std::vector<Edge> kept;
            kept.reserve(working.size());
            for (const auto& e : working) {
                const bool from_dead = (e.from < dropped.size() && dropped[e.from]);
                const bool to_dead   = (e.to   < dropped.size() && dropped[e.to]);
                if (from_dead || to_dead) continue;
                kept.push_back(e);
            }
            working = std::move(kept);
        }
    }

    // -------- Pass 2: transitive reduction --------
    if (cfg.reduce_transitive) {
        std::size_t removed = 0;
        working = reduce_transitive_edges(graph, working, cfg.transitive_fuzz,
                                          &removed);
        stats.transitive_edges_removed = removed;
    }

    graph.replace_edges(std::move(working));
    stats.edges_after = graph.edge_count();
    return stats;
}

}  // namespace branch::graph
