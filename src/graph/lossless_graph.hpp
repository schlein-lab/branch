#pragma once

// BRANCH v0.1 - Lossless String Graph.
//
// Central data model: a string-graph that preserves ALL overlaps and does
// not collapse Tandem-Repeats or Segmental-Duplicates. Nodes carry
// CN-awareness (copy-count is a per-node attribute rather than node
// duplication). Edges carry per-edge VAF (variant allele frequency)
// and read-support counts.
//
// This is a plain-data structure for v0.1; an actual builder (from
// overlaps) lives in a separate .cpp module added later.

#include <cstdint>
#include <cstddef>
#include <vector>

#include "graph/delta_read.hpp"

namespace branch::graph {

struct Node {
    NodeId id{};
    std::uint32_t length_bp{};   // consensus sequence length
    std::uint32_t copy_count{1}; // CN-aware: per-node attribute, not duplication
    float copy_count_confidence{1.0f};   // [0, 1]
    std::uint32_t read_support{0};       // number of reads traversing this node
};

struct Edge {
    NodeId from{};
    NodeId to{};
    std::uint32_t read_support{0};       // reads supporting this edge
    float vaf{1.0f};                     // variant allele frequency in [0, 1]
    float vaf_confidence{0.0f};          // [0, 1]
    std::uint16_t flags{0};              // bit-packed: is_branch, is_dup, etc.
    std::uint16_t _pad{0};
};

static_assert(sizeof(Edge) == 24,
              "Edge must pack to 24 bytes for SIMD-friendly iteration");

// Lossless graph container. Nodes and edges are kept in contiguous
// vectors keyed by NodeId. The adjacency lists live in edges_from_
// and edges_to_ flat arrays for cache-friendly traversal.
//
// Intentionally minimal in v0.1: a builder API + traversal helpers
// arrive in subsequent headers.
class LosslessGraph {
public:
    LosslessGraph() = default;

    std::size_t node_count() const noexcept { return nodes_.size(); }
    std::size_t edge_count() const noexcept { return edges_.size(); }

    NodeId add_node(std::uint32_t length_bp, std::uint32_t copy_count = 1) {
        NodeId id = static_cast<NodeId>(nodes_.size());
        nodes_.push_back(Node{
            .id = id,
            .length_bp = length_bp,
            .copy_count = copy_count,
        });
        return id;
    }

    void add_edge(NodeId from, NodeId to, std::uint32_t read_support = 0) {
        edges_.push_back(Edge{.from = from, .to = to, .read_support = read_support});
    }

    const Node& node(NodeId id) const { return nodes_[id]; }
    Node& node(NodeId id) { return nodes_[id]; }

    const std::vector<Node>& nodes() const noexcept { return nodes_; }
    const std::vector<Edge>& edges() const noexcept { return edges_; }

    // Bulk-update helpers used by the classifier and CN-inference passes.
    void set_copy_count(NodeId id, std::uint32_t cc, float confidence) {
        auto& n = nodes_[id];
        n.copy_count = cc;
        n.copy_count_confidence = confidence;
    }

    void set_edge_vaf(std::size_t edge_index, float vaf, float confidence) {
        edges_[edge_index].vaf = vaf;
        edges_[edge_index].vaf_confidence = confidence;
    }

private:
    std::vector<Node> nodes_;
    std::vector<Edge> edges_;
};

}  // namespace branch::graph
