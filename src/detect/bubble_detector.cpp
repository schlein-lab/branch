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
#include <atomic>
#include <unordered_map>
#include <unordered_set>

#include <omp.h>

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

// Walk forward from `start` up to max_steps, recording every reachable
// node with the *shortest* path found and its cumulative read support.
// Returns a map node→(path, support).
//
// Implementation: level-BFS with a visited set. Each node is visited at
// most once, regardless of cycles or graph density, yielding O(V + E)
// per start — bounded by max_steps layers. Previous DFS-with-path-copy
// had no visited set and blew up exponentially on dense/cyclic graphs
// (branching factor^max_steps path copies per entry).
struct PathRecord {
    std::vector<NodeId> nodes;
    std::uint32_t support{};
};

std::unordered_map<NodeId, PathRecord> reachable_within(
    NodeId start,
    const std::vector<std::vector<Succ>>& adj,
    std::uint32_t max_steps) {
    std::unordered_map<NodeId, PathRecord> out;
    if (start >= adj.size()) return out;

    // Use vector-indexed arrays rather than unordered_map: NodeIds are
    // dense (0..adj.size()-1) so O(1) direct access is faster than
    // hash-table inserts/lookups by ~10-50x on dense graphs.
    const std::size_t N = adj.size();
    std::vector<std::uint8_t> visited(N, 0);
    std::vector<NodeId> parent(N, 0);
    std::vector<std::uint32_t> edge_support(N, 0);

    // Track nodes touched so we only iterate visited ones at reconstruction
    // and only clear those slots on drop (not strictly needed here since
    // `visited` lives for one call, but keeps reconstruction linear).
    std::vector<NodeId> touched;
    touched.reserve(128);

    visited[start] = 1;
    std::vector<NodeId> frontier;
    frontier.push_back(start);
    for (std::uint32_t step = 0; step < max_steps && !frontier.empty(); ++step) {
        std::vector<NodeId> next_frontier;
        next_frontier.reserve(frontier.size() * 2);
        for (NodeId u : frontier) {
            for (const auto& s : adj[u]) {
                if (s.to >= N || visited[s.to]) continue;
                visited[s.to] = 1;
                parent[s.to] = u;
                edge_support[s.to] = s.read_support;
                touched.push_back(s.to);
                next_frontier.push_back(s.to);
            }
        }
        frontier = std::move(next_frontier);
    }

    // Reconstruct path [first_succ_of_start, ..., target] + cumulative support.
    out.reserve(touched.size());
    for (NodeId node : touched) {
        std::vector<NodeId> rev_path;
        std::uint32_t total_support = 0;
        NodeId cur = node;
        for (std::uint32_t hops = 0; cur != start && hops <= max_steps + 1; ++hops) {
            rev_path.push_back(cur);
            total_support += edge_support[cur];
            cur = parent[cur];
        }
        std::reverse(rev_path.begin(), rev_path.end());
        out.emplace(node, PathRecord{.nodes = std::move(rev_path), .support = total_support});
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

    // Each entry's work is completely independent of every other entry's
    // work — per-entry we compute reachability from each successor,
    // enumerate pairs, and propose (entry, exit) bubbles. The only
    // shared mutable state is the final `bubbles` map, which we merge
    // sequentially after the parallel region so per-thread insertions
    // don't fight for the same bucket under contention. Determinism is
    // preserved because the alt dedup + (entry, exit) key sort at the
    // end of this function already collapses any order-of-insertion
    // difference into a canonical output.
    using BubbleMap = std::unordered_map<EntryExitKey, Bubble, EntryExitHash>;
    std::atomic<std::size_t> skipped_high_outdeg{0};
    const int n_threads = std::max(1, omp_get_max_threads());
    std::vector<BubbleMap> thread_bubbles(n_threads);

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        BubbleMap& local = thread_bubbles[tid];

        #pragma omp for schedule(dynamic, 16)
        for (NodeId entry = 0; entry < static_cast<NodeId>(adj.size()); ++entry) {
            if (adj[entry].size() < 2) continue;
            if (cfg.max_entry_out_degree > 0 &&
                adj[entry].size() > cfg.max_entry_out_degree) {
                // Repeat-hub entries would explode the O(outdeg^2) pair loop.
                // Skipping them is the correct behaviour — a true biological
                // bubble has a handful of alternatives, not hundreds.
                skipped_high_outdeg.fetch_add(1, std::memory_order_relaxed);
                continue;
            }

            // Precompute reachability per successor *once* per entry so the
            // pair loop below is O(outdeg^2) map lookups rather than
            // O(outdeg^2) BFS calls. On the old code each BFS was recomputed
            // (outdeg-1) times per successor, which combined with the
            // unbounded DFS previously in reachable_within produced runtimes
            // that were effectively quadratic-in-outdeg on top of the
            // exponential path-copy blow-up.
            std::vector<std::unordered_map<NodeId, PathRecord>> succ_reach;
            succ_reach.reserve(adj[entry].size());
            for (const auto& s : adj[entry]) {
                succ_reach.emplace_back(reachable_within(s.to, adj,
                                                          cfg.max_alt_path_length));
            }

            // For each pair of out-successors, look for a common reachable exit.
            for (std::size_t i = 0; i < adj[entry].size(); ++i) {
                for (std::size_t j = i + 1; j < adj[entry].size(); ++j) {
                    const auto& ri = succ_reach[i];
                    const auto& rj = succ_reach[j];

                    // Iterate the smaller map, look up in the larger.
                    auto* small = &ri;
                    auto* large = &rj;
                    if (ri.size() > rj.size()) std::swap(small, large);

                    // Deterministic iteration: extract keys, sort, then iterate.
                    // std::unordered_map iteration order is impl-defined and
                    // allocation-order dependent — iterating directly causes
                    // non-deterministic memory-allocation patterns downstream
                    // (peak-RAM varies across runs on identical input).
                    std::vector<NodeId> exits;
                    exits.reserve(small->size());
                    for (const auto& kv : *small) exits.push_back(kv.first);
                    std::sort(exits.begin(), exits.end());

                    for (NodeId exit_node : exits) {
                        const auto& path_i = small->at(exit_node);
                        auto it = large->find(exit_node);
                        if (it == large->end()) continue;
                        if (exit_node == entry) continue;

                        EntryExitKey key{entry, exit_node};
                        auto [bit, inserted] = local.try_emplace(key);
                        if (inserted) {
                            bit->second.entry = entry;
                            bit->second.exit = exit_node;
                        }

                        // Alt from i
                        AltPath alt_i;
                        alt_i.nodes.push_back(adj[entry][i].to);
                        alt_i.nodes.insert(alt_i.nodes.end(),
                                           path_i.nodes.begin(), path_i.nodes.end());
                        alt_i.total_read_support =
                            adj[entry][i].read_support + path_i.support;
                        bit->second.alts.push_back(std::move(alt_i));

                        AltPath alt_j;
                        alt_j.nodes.push_back(adj[entry][j].to);
                        const auto& pj = it->second;
                        alt_j.nodes.insert(alt_j.nodes.end(),
                                           pj.nodes.begin(), pj.nodes.end());
                        alt_j.total_read_support =
                            adj[entry][j].read_support + pj.support;
                        bit->second.alts.push_back(std::move(alt_j));
                        bit->second.total_read_support =
                            adj[entry][i].read_support + adj[entry][j].read_support;
                    }
                }
            }
        }
    }

    // Sequential merge: fold per-thread maps into the global map. We
    // iterate threads in index order, and for each shared (entry, exit)
    // key we concatenate the alts; the final dedup-sort-unique pass
    // below canonicalises regardless of source-thread order.
    for (auto& tb : thread_bubbles) {
        for (auto& kv : tb) {
            auto [bit, inserted] = bubbles.try_emplace(kv.first);
            if (inserted) {
                bit->second = std::move(kv.second);
            } else {
                // Merge alts; keep the larger total_read_support as a
                // conservative aggregate. Alts will be deduped below.
                for (auto& alt : kv.second.alts) {
                    bit->second.alts.push_back(std::move(alt));
                }
                bit->second.total_read_support = std::max(
                    bit->second.total_read_support,
                    kv.second.total_read_support);
            }
        }
        tb.clear();
    }

    // Deterministic emission: sort (entry, exit) keys before final pass
    // so output vector has a stable order independent of hash bucketing.
    std::vector<EntryExitKey> keys;
    keys.reserve(bubbles.size());
    for (const auto& kv : bubbles) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end(),
              [](const EntryExitKey& a, const EntryExitKey& b) {
                  if (a.entry != b.entry) return a.entry < b.entry;
                  return a.exit < b.exit;
              });

    out.reserve(bubbles.size());
    for (const auto& k : keys) {
        auto& b = bubbles[k];
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
