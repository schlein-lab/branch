#include "branch/graph/graph_builder.hpp"
#include "branch/graph/overlap_graph.hpp"
#include "branch/overlapper/overlap_pair.hpp"
#include <unordered_map>
#include <utility>

namespace branch {
namespace graph {

// Hash for std::pair<NodeId, NodeId>
struct NodePairHash {
    std::size_t operator()(const std::pair<size_t, size_t>& p) const {
        return std::hash<size_t>()(p.first) ^ (std::hash<size_t>()(p.second) << 1);
    }
};

OverlapGraph GraphBuilder::build_from_overlaps(
    const std::vector<overlapper::OverlapPair>& overlaps,
    size_t num_reads) {
    
    OverlapGraph graph;
    
    // Add nodes for each read with read_support = 1
    for (size_t i = 0; i < num_reads; ++i) {
        graph.add_node(0, 1);  // length=0 (will be set later), read_support=1
    }
    
    // Aggregate edges: count how many times each (from, to) pair appears
    std::unordered_map<std::pair<size_t, size_t>, size_t, NodePairHash> edge_counts;
    
    for (const auto& overlap : overlaps) {
        auto key = std::make_pair(overlap.read_a_id, overlap.read_b_id);
        edge_counts[key]++;
    }
    
    // Add edges with aggregated read_support
    for (auto it = edge_counts.begin(); it != edge_counts.end(); ++it) {
        graph.add_edge(it->first.first, it->first.second, it->second);
    }
    
    return graph;
}

} // namespace graph
} // namespace branch
