#include "branch/graph/graph_compactor.hpp"
#include <unordered_set>
#include <unordered_map>
#include <algorithm>

namespace branch {
namespace graph {

std::vector<std::vector<size_t>> GraphCompactor::find_unitig_chains(
    const OverlapGraph& graph) {
    
    std::vector<std::vector<size_t>> chains;
    std::unordered_set<size_t> visited;
    
    const auto& nodes = graph.nodes();
    const auto& edges = graph.edges();
    
    // Build adjacency lists
    std::unordered_map<size_t, std::vector<size_t>> outgoing;
    std::unordered_map<size_t, std::vector<size_t>> incoming;
    
    for (const auto& edge : edges) {
        outgoing[edge.from].push_back(edge.to);
        incoming[edge.to].push_back(edge.from);
    }
    
    // Find chain starts (nodes with in-degree != 1 or out-degree != 1)
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (visited.count(i)) continue;
        
        size_t in_deg = incoming[i].size();
        size_t out_deg = outgoing[i].size();
        
        // Start a chain from branch points or sources
        if (in_deg != 1 || out_deg != 1) {
            // Try to extend forward
            if (out_deg > 0) {
                std::vector<size_t> chain;
                chain.push_back(i);
                visited.insert(i);
                
                size_t current = i;
                while (outgoing[current].size() == 1) {
                    size_t next = outgoing[current][0];
                    if (visited.count(next)) break;
                    if (incoming[next].size() != 1) break;
                    
                    chain.push_back(next);
                    visited.insert(next);
                    current = next;
                }
                
                if (chain.size() > 1) {
                    chains.push_back(chain);
                }
            }
        }
    }
    
    // Handle isolated linear chains (cycles or unvisited linear segments)
    for (size_t i = 0; i < nodes.size(); ++i) {
        if (visited.count(i)) continue;
        
        std::vector<size_t> chain;
        chain.push_back(i);
        visited.insert(i);
        
        size_t current = i;
        while (outgoing[current].size() == 1) {
            size_t next = outgoing[current][0];
            if (visited.count(next)) break;
            chain.push_back(next);
            visited.insert(next);
            current = next;
        }
        
        if (chain.size() > 1) {
            chains.push_back(chain);
        }
    }
    
    return chains;
}

OverlapGraph GraphCompactor::compact_unitigs_with_sequences(
    const OverlapGraph& graph,
    const std::vector<std::string>& /* sequences */) {
    
    auto chains = find_unitig_chains(graph);
    
    if (chains.empty()) {
        return graph;  // Nothing to compact
    }
    
    const auto& old_nodes = graph.nodes();
    const auto& old_edges = graph.edges();
    
    // Map old node IDs to new node IDs (or -1 if compacted away)
    std::unordered_map<size_t, size_t> node_mapping;
    std::unordered_set<size_t> compacted_nodes;
    
    OverlapGraph new_graph;
    
    // First, mark all nodes that are part of chains (except chain heads)
    for (const auto& chain : chains) {
        for (size_t i = 1; i < chain.size(); ++i) {
            compacted_nodes.insert(chain[i]);
        }
    }
    
    // Add non-compacted nodes first
    size_t new_id = 0;
    for (size_t i = 0; i < old_nodes.size(); ++i) {
        if (!compacted_nodes.count(i)) {
            node_mapping[i] = new_id++;
            new_graph.add_node(old_nodes[i].length_bp, old_nodes[i].read_support);
        }
    }
    
    // Process chains: update the chain head node with summed values
    for (const auto& chain : chains) {
        size_t head = chain[0];
        size_t head_new_id = node_mapping[head];
        
        // Sum length and read_support along the chain
        size_t total_length = 0;
        size_t total_read_support = 0;
        for (size_t node_id : chain) {
            total_length += old_nodes[node_id].length_bp;
            total_read_support += old_nodes[node_id].read_support;
        }
        
        // Update the head node with summed values
        // We need to access the node directly - use a workaround
        // by rebuilding: remove and re-add is not supported, so we
        // track the sums and apply them
        auto& nodes_ref = const_cast<std::vector<Node>&>(new_graph.nodes());
        nodes_ref[head_new_id].length_bp = total_length;
        nodes_ref[head_new_id].read_support = total_read_support;
        
        // Map all chain nodes to the head's new ID
        for (size_t node_id : chain) {
            node_mapping[node_id] = head_new_id;
        }
    }
    
    // Add edges with updated node IDs
    std::unordered_set<size_t> added_edges;  // To avoid duplicates
    for (const auto& edge : old_edges) {
        size_t new_from = node_mapping[edge.from];
        size_t new_to = node_mapping[edge.to];
        
        // Skip self-loops created by compaction
        if (new_from == new_to) continue;
        
        // Skip edges internal to chains
        size_t edge_key = new_from * old_nodes.size() + new_to;
        if (added_edges.count(edge_key)) continue;
        added_edges.insert(edge_key);
        
        new_graph.add_edge(new_from, new_to, edge.read_support);
    }
    
    return new_graph;
}

} // namespace graph
} // namespace branch
