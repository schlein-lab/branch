#include "bubble_finder.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace branch::wg::phaser {

namespace {

// BFS from `source` for at most `max_depth` unitig hops, collecting all
// nodes reachable. Returns the set of reachable nodes and any node
// where multiple paths re-converge — those are bubble sinks.
struct BfsResult {
    std::unordered_set<NodeId> reachable;
    std::unordered_set<NodeId> sinks;  // re-convergence points
};

BfsResult bfs_local(const std::vector<UnitigNode>& unitigs,
                    NodeId source, int max_depth) {
    BfsResult r;
    std::queue<std::pair<NodeId, int>> q;
    q.emplace(source, 0);
    r.reachable.insert(source);
    std::unordered_map<NodeId, int> parents_count;
    while (!q.empty()) {
        auto [cur, d] = q.front(); q.pop();
        if (d >= max_depth) continue;
        if (cur >= unitigs.size()) continue;
        for (NodeId n : unitigs[cur].next_nodes) {
            if (n >= unitigs.size()) continue;
            ++parents_count[n];
            if (r.reachable.insert(n).second) {
                q.emplace(n, d + 1);
            }
            if (parents_count[n] >= 2) r.sinks.insert(n);
        }
    }
    return r;
}

// Trace one alt-path from `source` toward `sink` choosing edge `which_out`
// at the source (0 = first outgoing, 1 = second, ...).
std::vector<NodeId> trace_alt_path(const std::vector<UnitigNode>& unitigs,
                                   NodeId source, NodeId sink,
                                   std::size_t which_out, int max_steps) {
    std::vector<NodeId> path;
    if (source >= unitigs.size()) return path;
    if (which_out >= unitigs[source].next_nodes.size()) return path;
    NodeId cur = unitigs[source].next_nodes[which_out];
    int steps = 0;
    while (cur != sink && cur < unitigs.size() && steps++ < max_steps) {
        path.push_back(cur);
        if (unitigs[cur].next_nodes.empty()) break;
        cur = unitigs[cur].next_nodes[0];  // greedy
    }
    if (cur == sink) {
        return path;
    }
    return {};
}

}  // namespace

void BubbleFinder::find(std::vector<UnitigNode>& unitigs,
                        std::vector<Bubble>& out_bubbles,
                        const BubbleFinderOpts& opts)
{
    out_bubbles.clear();
    BubbleId next_id = 0;

    // Pre-filter forks. Iterating all 1.95M unitigs at WG scale was 45 min
    // due to cache-cold UnitigNode access; ~99.99% of nodes have <2
    // outgoing edges and would early-continue anyway. Building the fork
    // list once cuts bubble_finder from minutes to milliseconds.
    auto t_start = std::chrono::steady_clock::now();
    std::vector<NodeId> forks;
    forks.reserve(unitigs.size() / 100);
    for (NodeId i = 0; i < unitigs.size(); ++i) {
        if (unitigs[i].next_nodes.size() >= 2) forks.push_back(i);
    }
    std::cerr << "[phaser/bf] " << forks.size() << " fork unitigs of "
              << unitigs.size() << "\n";
    std::fflush(stderr);

    for (std::size_t fi = 0; fi < forks.size(); ++fi) {
        const NodeId src = forks[fi];
        if (fi % 50'000 == 0 && fi > 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            std::fprintf(stderr,
                "[phaser/bf] scan %zu/%zu (%.1f%%) bubbles=%u elapsed=%.0fs\n",
                fi, forks.size(),
                100.0 * static_cast<double>(fi) /
                static_cast<double>(forks.size()),
                next_id, el);
            std::fflush(stderr);
        }
        const auto& s = unitigs[src];
        if (s.next_nodes.size() < 2) continue;  // double-check (defensive)

        BfsResult r = bfs_local(unitigs, src, opts.max_alt_path_unitigs);
        // Pick the closest re-convergence sink.
        NodeId sink = ~0u;
        for (NodeId cand : r.sinks) {
            if (unitigs[cand].prev_nodes.size() < 2) continue;
            sink = cand;
            break;
        }
        if (sink == ~0u) continue;

        // Trace each fork-out path independently.
        std::vector<std::vector<NodeId>> alt_paths;
        for (std::size_t k = 0; k < s.next_nodes.size(); ++k) {
            auto p = trace_alt_path(unitigs, src, sink, k,
                                    opts.max_alt_path_unitigs);
            if (!p.empty()) alt_paths.push_back(std::move(p));
        }
        if (alt_paths.size() < 2) continue;

        // Span check.
        Position span = 0;
        for (NodeId nid : alt_paths.front()) {
            if (nid < unitigs.size())
                span += unitigs[nid].length_bp;
        }
        if (span > opts.max_bubble_span_bp) continue;

        Bubble b;
        b.id = next_id++;
        b.source_anchor = src;
        b.sink_anchor   = sink;
        for (auto& p : alt_paths) {
            for (NodeId nid : p) {
                if (nid < unitigs.size()) {
                    unitigs[nid].bubble = b.id;
                }
            }
        }
        // Save the first node of each alt path as the "label" position
        // for phase_engine to vote on. Anchors stay UNCLASSIFIED here;
        // they'll be tagged SHARED later by phase_engine when it finds
        // them outside any bubble.
        std::vector<NodeId> alt_heads;
        alt_heads.reserve(alt_paths.size());
        for (auto& p : alt_paths) alt_heads.push_back(p.front());
        b.alt_paths = std::move(alt_heads);
        out_bubbles.push_back(std::move(b));
    }

    std::cerr << "[phaser/bf] found " << out_bubbles.size() << " bubbles\n";
}

}  // namespace branch::wg::phaser
