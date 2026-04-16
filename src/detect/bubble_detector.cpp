// BRANCH v0.2 — Bubble detection implementation.
//
// Algorithm (v0.2, simple-bubble-only):
//  For every node with out-degree >= 2:
//    - Enumerate outgoing successors.
//    - For each pair of successors, follow forward edges up to
//      max_alt_path_length steps until they reach a common node.
//    - If a common exit exists, record the bubble.
//
// Limitations: only detects "binary" bubbles between two alternatives
// per entry. Multi-allelic (n>2) and nested bubbles are v0.3 work.
// No transitive-reduction here; duplicate emissions are filtered by
// (entry, exit) uniqueness at the end.

#include "detect/bubble_detector.hpp"

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace branch::detect {

namespace {

using NodeId = branch::graph::NodeId;

struct Succ {
    NodeId to;
    std::uint32_t read_support;
};

// Build out-adjacency list from the graph's edge list.
std::vector<std::vector<Succ>> build_adjacency(
    const branch::graph::LosslessGraph& graph) {
    std::vector<std::vector<Succ>> adj(graph.node_count());
    for (const auto& e : graph.edges()) {
        if (e.from < adj.size()) {
            adj[e.from].push_back(Succ{.to = e.to, .read_support = e.read_support});
        }
    }
    return adj;
}

// Walk forward from `start` up to max_steps, recording every visited
// node and the path to reach it. Returns a map node→(path, support).
struct PathRecord {
    std::vector<NodeId> nodes;
    std::uint32_t support{};
};

std::unordered_map<NodeId, PathRecord> reachable_within(
    NodeId start,
    const std::vector<std::vector<Succ>>& adj,
    std::uint32_t max_steps) {
    std::unordered_map<NodeId, PathRecord> out;
    struct Frame { NodeId node; std::vector<NodeId> path; std::uint32_t support; std::uint32_t steps; };
    std::vector<Frame> stack;
    stack.push_back({start, {}, 0, 0});

    while (!stack.empty()) {
        auto f = std::move(stack.back());
        stack.pop_back();
        if (f.steps > max_steps) continue;

        // Record reaching f.node (overwrite is fine; we just need any path).
        if (f.node != start) {
            auto [it, inserted] = out.try_emplace(f.node,
                PathRecord{.nodes = f.path, .support = f.support});
            (void)it; (void)inserted;
        }

        if (f.node >= adj.size()) continue;
        for (const auto& s : adj[f.node]) {
            auto new_path = f.path;
            new_path.push_back(s.to);
            stack.push_back({s.to, std::move(new_path),
                             f.support + s.read_support,
                             f.steps + 1});
        }
    }
    return out;
}

}  // namespace

std::vector<Bubble> detect_bubbles(
    const branch::graph::LosslessGraph& graph,
    const BubbleDetectorConfig& cfg) {
    std::vector<Bubble> out;
    auto adj = build_adjacency(graph);

    // Deduplicate (entry,exit) so repeat emissions merge.
    struct EntryExitKey {
        NodeId entry;
        NodeId exit;
        bool operator==(const EntryExitKey& o) const {
            return entry == o.entry && exit == o.exit;
        }
    };
    struct EntryExitHash {
        std::size_t operator()(const EntryExitKey& k) const noexcept {
            return (static_cast<std::uint64_t>(k.entry) << 32) ^ k.exit;
        }
    };
    std::unordered_map<EntryExitKey, Bubble, EntryExitHash> bubbles;

    for (NodeId entry = 0; entry < adj.size(); ++entry) {
        if (adj[entry].size() < 2) continue;

        // For each pair of out-successors, look for a common reachable exit.
        for (std::size_t i = 0; i < adj[entry].size(); ++i) {
            for (std::size_t j = i + 1; j < adj[entry].size(); ++j) {
                auto ri = reachable_within(adj[entry][i].to, adj,
                                           cfg.max_alt_path_length);
                auto rj = reachable_within(adj[entry][j].to, adj,
                                           cfg.max_alt_path_length);

                // Iterate the smaller map, look up in the larger.
                auto* small = &ri;
                auto* large = &rj;
                if (ri.size() > rj.size()) std::swap(small, large);

                for (const auto& [exit_node, path_i] : *small) {
                    auto it = large->find(exit_node);
                    if (it == large->end()) continue;
                    if (exit_node == entry) continue;

                    EntryExitKey key{entry, exit_node};
                    auto [bit, inserted] = bubbles.try_emplace(key);
                    if (inserted) {
                        bit->second.entry = entry;
                        bit->second.exit = exit_node;
                    }

                    // Alt from i
                    AltPath alt_i;
                    alt_i.nodes.push_back(adj[entry][i].to);
                    alt_i.nodes.insert(alt_i.nodes.end(),
                                       path_i.nodes.begin(), path_i.nodes.end());
                    alt_i.total_read_support = adj[entry][i].read_support + path_i.support;
                    bit->second.alts.push_back(std::move(alt_i));

                    AltPath alt_j;
                    alt_j.nodes.push_back(adj[entry][j].to);
                    const auto& pj = it->second;
                    alt_j.nodes.insert(alt_j.nodes.end(),
                                       pj.nodes.begin(), pj.nodes.end());
                    alt_j.total_read_support = adj[entry][j].read_support + pj.support;
                    bit->second.alts.push_back(std::move(alt_j));
                    bit->second.total_read_support =
                        adj[entry][i].read_support + adj[entry][j].read_support;
                }
            }
        }
    }

    out.reserve(bubbles.size());
    for (auto& [k, b] : bubbles) {
        // De-duplicate alts (both (i,j) directions may have added same path).
        std::sort(b.alts.begin(), b.alts.end(),
                  [](const AltPath& a, const AltPath& c) { return a.nodes < c.nodes; });
        b.alts.erase(
            std::unique(b.alts.begin(), b.alts.end(),
                        [](const AltPath& a, const AltPath& c) {
                            return a.nodes == c.nodes;
                        }),
            b.alts.end());
        out.push_back(std::move(b));
    }
    return out;
}

}  // namespace branch::detect
