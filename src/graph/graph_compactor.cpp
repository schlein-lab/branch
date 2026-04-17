// BRANCH v0.2 — Unitig compaction implementation with consensus building.
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
//        consensus = MSA-derived consensus from member sequences
//      Copy only inter-unitig edges (from_unitig != to_unitig), keeping
//      the original read_support of the boundary-crossing edge.

#include "graph/graph_compactor.hpp"
#include "consensus/majority_voter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

namespace branch::graph {

namespace {

constexpr NodeId kUnassigned = std::numeric_limits<NodeId>::max();

// Position-wise majority voting for unaligned sequences.
// Sequences of different lengths are padded with conceptual gaps
// (positions past the end contribute nothing to the base count).
// Adequate for well-overlapping unitig members where indels are rare;
// a POA-based aligner (abPOA or equivalent) would improve quality on
// pathological regions but is not vendored yet — revisit when the
// assembler starts hitting high-indel targets.
std::string simple_majority_consensus(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {};
    }
    if (sequences.size() == 1) {
        return sequences[0];
    }
    
    // Find max length
    std::size_t max_len = 0;
    for (const auto& seq : sequences) {
        max_len = std::max(max_len, seq.size());
    }
    
    std::string result;
    result.reserve(max_len);
    
    for (std::size_t pos = 0; pos < max_len; ++pos) {
        // Count bases at this position
        std::array<int, 5> counts{}; // A=0, C=1, G=2, T=3, other=4
        
        for (const auto& seq : sequences) {
            if (pos >= seq.size()) continue;
            char c = seq[pos];
            switch (c) {
                case 'A': case 'a': ++counts[0]; break;
                case 'C': case 'c': ++counts[1]; break;
                case 'G': case 'g': ++counts[2]; break;
                case 'T': case 't': ++counts[3]; break;
                default: ++counts[4]; break;
            }
        }
        
        // Find majority base
        int max_count = 0;
        int max_idx = 4;
        for (int i = 0; i < 4; ++i) {
            if (counts[i] > max_count) {
                max_count = counts[i];
                max_idx = i;
            }
        }
        
        constexpr char bases[] = "ACGTN";
        result += bases[max_idx];
    }
    
    return result;
}

// Build consensus from member node sequences.
// Strategy:
//   - 0 sequences with content -> empty consensus
//   - 1 sequence -> use directly
//   - 2+ sequences -> simple_majority_consensus (position-wise majority
//     with length padding). A POA-based aligner would improve indel-
//     heavy regions but is not yet vendored.
std::string build_unitig_consensus(
    const std::vector<NodeId>& members,
    const LosslessGraph& input) {
    
    // Collect non-empty consensus sequences from member nodes
    std::vector<std::string> sequences;
    sequences.reserve(members.size());
    
    for (NodeId m : members) {
        const auto& node = input.node(m);
        if (!node.consensus.empty()) {
            sequences.push_back(node.consensus);
        }
    }
    
    if (sequences.empty()) {
        return {};
    }
    
    if (sequences.size() == 1) {
        return sequences[0];
    }
    
    return simple_majority_consensus(sequences);
}

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

    // ---- Pass 2: walk unitigs starting from every chain start ----
    std::vector<std::vector<NodeId>> unitigs;
    unitigs.reserve(n);

    auto is_chain_start = [&](NodeId v) -> bool {
        if (in_degree[v] != 1) return true;  // source or branching -> chain start
        // Safe here: in_degree[v] == 1 guarantees unique_in_edge[v] was set above.
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
    for (NodeId v = 0; v < n; ++v) {
        if (result.old_to_new[v] != kUnassigned) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            if (result.old_to_new[cur] != kUnassigned) break;
            members.push_back(cur);
            result.old_to_new[cur] = unitig_id;
            if (out_degree[cur] != 1) break;
            std::size_t eout = unique_out_edge[cur];
            NodeId next = edges[eout].to;
            if (next >= n) break;
            if (in_degree[next] != 1) break;
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // ---- Pass 3: build the compacted graph with consensus ----
    LosslessGraph& out = result.compacted;
    for (const auto& members : unitigs) {
        std::uint64_t total_length = 0;
        for (NodeId m : members) {
            total_length += input.node(m).length_bp;
        }
        const auto& start_node = input.node(members.front());
        std::uint32_t len32 = (total_length > std::numeric_limits<std::uint32_t>::max())
            ? std::numeric_limits<std::uint32_t>::max()
            : static_cast<std::uint32_t>(total_length);
        
        NodeId new_id = out.add_node(len32, start_node.copy_count);
        
        // Build and set consensus for the new unitig node
        std::string consensus = build_unitig_consensus(members, input);
        if (!consensus.empty()) {
            out.node(new_id).consensus = std::move(consensus);
        }
    }

    // Emit only inter-unitig edges
    for (const auto& e : edges) {
        if (e.from >= n || e.to >= n) continue;
        NodeId uf = result.old_to_new[e.from];
        NodeId ut = result.old_to_new[e.to];
        if (uf == ut) {
            continue;
        }
        out.add_edge(uf, ut, e.read_support);
    }

    return result;
}

// New overload: compact with explicit sequences for consensus building
CompactionResult compact_unitigs_with_sequences(
    const LosslessGraph& input,
    const std::vector<std::string>& node_sequences) {
    
    // First, copy sequences into input graph's consensus fields if provided
    // (This is a workaround since input is const - we build consensus from
    // the separate sequences vector)
    
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

    // ---- Pass 2: walk unitigs ----
    std::vector<std::vector<NodeId>> unitigs;
    unitigs.reserve(n);

    auto is_chain_start = [&](NodeId v) -> bool {
        if (in_degree[v] != 1) return true;  // source or branching -> chain start
        // Safe here: in_degree[v] == 1 guarantees unique_in_edge[v] was set above.
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
            if (next >= n) break;
            if (in_degree[next] != 1) break;
            if (result.old_to_new[next] != kUnassigned) break;
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // Pass 2b: pure cycles
    for (NodeId v = 0; v < n; ++v) {
        if (result.old_to_new[v] != kUnassigned) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            if (result.old_to_new[cur] != kUnassigned) break;
            members.push_back(cur);
            result.old_to_new[cur] = unitig_id;
            if (out_degree[cur] != 1) break;
            std::size_t eout = unique_out_edge[cur];
            NodeId next = edges[eout].to;
            if (next >= n) break;
            if (in_degree[next] != 1) break;
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // ---- Pass 3: build compacted graph with consensus from sequences ----
    // Split into three sub-passes so the expensive per-unitig consensus
    // computation can run in parallel:
    //   (3a) Sequential: allocate a new_id per unitig and gather its member
    //        sequences. add_node mutates `out` and must stay single-threaded.
    //   (3b) Parallel: compute simple_majority_consensus(seqs[i]) for every
    //        unitig. Pure per-unitig, no shared state.
    //   (3c) Sequential: assign each consensus into out.node(new_ids[i]).
    LosslessGraph& out = result.compacted;
    const std::size_t u_count = unitigs.size();

    std::vector<NodeId> new_ids(u_count);
    std::vector<std::vector<std::string>> seqs_per_unitig(u_count);

    // (3a)
    for (std::size_t i = 0; i < u_count; ++i) {
        const auto& members = unitigs[i];
        std::uint64_t total_length = 0;
        for (NodeId m : members) {
            total_length += input.node(m).length_bp;
        }
        const auto& start_node = input.node(members.front());
        std::uint32_t len32 = (total_length > std::numeric_limits<std::uint32_t>::max())
            ? std::numeric_limits<std::uint32_t>::max()
            : static_cast<std::uint32_t>(total_length);

        new_ids[i] = out.add_node(len32, start_node.copy_count);

        auto& seqs = seqs_per_unitig[i];
        seqs.reserve(members.size());
        for (NodeId m : members) {
            if (m < node_sequences.size() && !node_sequences[m].empty()) {
                seqs.push_back(node_sequences[m]);
            }
        }
    }

    // (3b) — the expensive bit, now parallel
    std::vector<std::string> consensuses(u_count);
    #pragma omp parallel for schedule(dynamic, 16)
    for (std::size_t i = 0; i < u_count; ++i) {
        if (!seqs_per_unitig[i].empty()) {
            consensuses[i] = simple_majority_consensus(seqs_per_unitig[i]);
        }
    }

    // (3c)
    std::size_t empty_consensus_unitigs = 0;
    for (std::size_t i = 0; i < u_count; ++i) {
        if (!consensuses[i].empty()) {
            out.node(new_ids[i]).consensus = std::move(consensuses[i]);
        } else {
            // Finding 4 cascade: empty consensus here means the downstream
            // write_bed_with_refs path cannot FASTA-dump this node, so
            // minimap2 never sees it and the BED row falls back to chrom=NA.
            // Upstream fix is tracked under Finding 1 (unitig collapse);
            // see feat/compactor-debug. Count-and-warn until then.
            ++empty_consensus_unitigs;
        }
    }
    if (empty_consensus_unitigs > 0) {
        std::cerr << "BRANCH: graph_compactor produced "
                  << empty_consensus_unitigs
                  << " unitigs with empty consensus (no member sequences "
                     "supplied). Downstream BED rows will have chrom=NA for "
                     "these. Root cause tracked in Finding 1 (unitig "
                     "collapse, feat/compactor-debug).\n";
    }

    // Emit inter-unitig edges
    for (const auto& e : edges) {
        if (e.from >= n || e.to >= n) continue;
        NodeId uf = result.old_to_new[e.from];
        NodeId ut = result.old_to_new[e.to];
        if (uf == ut) continue;
        out.add_edge(uf, ut, e.read_support);
    }

    return result;
}

}  // namespace branch::graph
