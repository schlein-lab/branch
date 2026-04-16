// BRANCH v0.2 — Unitig compaction implementation.
//
// Algorithm (two passes):
//   1) Build in-degree / out-degree arrays and per-node adjacency indices
//      into the input edge vector.
//   2) Identify chain starts: a node N starts a unitig if either
//        - in_degree(N) != 1, or
//        - the unique predecessor P has out_degree(P) != 1.
//      From each start, walk forward while the current node has out_degree
//      exactly 1 AND the successor has in_degree exactly 1. Collect all
//      visited nodes into the same unitig.
//   3) Any node that is not yet assigned to a unitig must be on a pure
//      cycle (every node in the cycle has in=out=1). Break such a cycle
//      arbitrarily at an unvisited node and walk it as its own unitig.
//   4) Build the output graph: one new node per unitig with
//        length_bp = sum of member lengths
//        copy_count = chain start's copy_count
//      Copy only inter-unitig edges (from_unitig != to_unitig), keeping
//      the original read_support of the boundary-crossing edge.

#include "graph/graph_compactor.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace branch::graph {

namespace {

constexpr NodeId kUnassigned = std::numeric_limits<NodeId>::max();

}  // namespace

CompactionResult compact_unitigs(const LosslessGraph& input) {
    CompactionResult result;
    const std::size_t n = input.node_count();
    result.old_to_new.assign(n, kUnassigned);

    if (n == 0) {
        return result;
    }

    // ---- Pass 1: degree counts and adjacency lookup ----
    std::vector<std::uint32_t> in_degree(n, 0);
    std::vector<std::uint32_t> out_degree(n, 0);

    const auto& edges = input.edges();
    for (const auto& e : edges) {
        if (e.from < n) ++out_degree[e.from];
        if (e.to < n)   ++in_degree[e.to];
    }

    // For each node that has out_degree == 1, remember its unique successor
    // edge index so we can walk chains in O(1) per step.
    std::vector<std::size_t> unique_out_edge(n, std::numeric_limits<std::size_t>::max());
    std::vector<std::size_t> unique_in_edge(n, std::numeric_limits<std::size_t>::max());
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        if (e.from < n && out_degree[e.from] == 1) {
            unique_out_edge[e.from] = i;
        }
        if (e.to < n && in_degree[e.to] == 1) {
            unique_in_edge[e.to] = i;
        }
    }

    // Helper: is node v an "internal" chain node? (in=out=1 AND the
    // unique predecessor has out=1 — but that predecessor condition is
    // symmetric to the start-detection below. For deciding "should I
    // extend past v to v's successor" we only need: out_degree[v]==1 and
    // in_degree[succ]==1.)

    // ---- Pass 2: walk unitigs starting from every chain start ----
    std::vector<std::vector<NodeId>> unitigs;
    unitigs.reserve(n);

    auto is_chain_start = [&](NodeId v) -> bool {
        if (in_degree[v] != 1) return true;
        // unique predecessor:
        std::size_t ein = unique_in_edge[v];
        const Edge& pe = edges[ein];
        NodeId pred = pe.from;
        return out_degree[pred] != 1;
    };

    for (NodeId v = 0; v < n; ++v) {
        if (result.old_to_new[v] != kUnassigned) continue;
        if (!is_chain_start(v)) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            members.push_back(cur);
            result.old_to_new[cur] = unitig_id;

            if (out_degree[cur] != 1) break;
            std::size_t eout = unique_out_edge[cur];
            NodeId next = edges[eout].to;
            if (next >= n) break;  // defensive
            if (in_degree[next] != 1) break;
            if (result.old_to_new[next] != kUnassigned) break;  // cycle-safety
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // ---- Pass 2b: pure cycles (every node in=out=1, no start found) ----
    // Any node still unassigned lies on a pure cycle. Break arbitrarily.
    for (NodeId v = 0; v < n; ++v) {
        if (result.old_to_new[v] != kUnassigned) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            if (result.old_to_new[cur] != kUnassigned) break;
            members.push_back(cur);
            result.old_to_new[cur] = unitig_id;
            if (out_degree[cur] != 1) break;  // should not happen on a pure cycle
            std::size_t eout = unique_out_edge[cur];
            NodeId next = edges[eout].to;
            if (next >= n) break;
            if (in_degree[next] != 1) break;
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // ---- Pass 3: build the compacted graph ----
    LosslessGraph& out = result.compacted;
    for (const auto& members : unitigs) {
        std::uint64_t total_length = 0;
        for (NodeId m : members) {
            total_length += input.node(m).length_bp;
        }
        const auto& start_node = input.node(members.front());
        // add_node takes uint32_t; clamp if (unrealistically) overflowed.
        std::uint32_t len32 = (total_length > std::numeric_limits<std::uint32_t>::max())
            ? std::numeric_limits<std::uint32_t>::max()
            : static_cast<std::uint32_t>(total_length);
        out.add_node(len32, start_node.copy_count);
    }

    // Emit only inter-unitig edges. Intra-unitig edges (i.e. the edges
    // absorbed into a chain) are dropped. Inter-unitig edges keep their
    // original read_support verbatim.
    for (const auto& e : edges) {
        if (e.from >= n || e.to >= n) continue;
        NodeId uf = result.old_to_new[e.from];
        NodeId ut = result.old_to_new[e.to];
        if (uf == ut) {
            // Could be an intra-chain edge (drop) or a self-loop on a
            // pure-cycle unitig (also drop, since the cycle was linearised).
            continue;
        }
        out.add_edge(uf, ut, e.read_support);
    }

    return result;
}

}  // namespace branch::graph
