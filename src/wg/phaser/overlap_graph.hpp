#pragma once
//
// overlap_graph — minimizer-based all-vs-all read overlap, then
// transitive reduction → string graph of unitigs.
//
// Hifiasm-style approach: each read becomes a node, edges are computed
// from minimizer-sketch Jaccard similarity. Linear chains in the read
// graph are collapsed into unitigs.

#include "types.hpp"
#include <cstdint>
#include <string>
#include <vector>
#include <utility>

namespace branch::wg::phaser {

struct OverlapGraphOpts {
    int minimizer_k    = 21;
    int minimizer_w    = 11;
    int min_overlap_bp = 1000;
    double min_jaccard = 0.35;
    int n_threads      = 1;
};

struct OverlapEdge {
    ReadId  a;
    ReadId  b;
    int     overlap_bp;
    float   jaccard;
    bool    a_to_b;        // direction in canonical orientation
};

class OverlapGraph {
public:
    OverlapGraph() = default;

    // Build the read overlap graph from raw reads. Sequences are accessed
    // through indices into `reads`; the function does not retain pointers.
    void build(const std::vector<std::pair<ReadId, std::string>>& reads,
               const OverlapGraphOpts& opts);

    // After build(): collapse linear chains into unitigs.
    // The output `unitigs` array becomes the working node set used by
    // bubble_finder, phase_engine, etc.
    void collapse_to_unitigs(
        const std::vector<std::pair<ReadId, std::string>>& reads,
        std::vector<UnitigNode>& out_unitigs) const;

    // Direct access to the read-level overlap edges (pre-collapse), for
    // diagnostics / tests.
    const std::vector<OverlapEdge>& edges() const { return edges_; }

    // n × ~50 bytes graph-state estimate, reported up-front.
    static std::size_t estimate_ram_bytes(std::size_t n_reads);

private:
    std::vector<OverlapEdge> edges_;
    OverlapGraphOpts         opts_{};
};

}  // namespace branch::wg::phaser
