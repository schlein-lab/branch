#pragma once
//
// overlap_graph — minimizer-based all-vs-all read overlap, then
// transitive reduction → string graph of unitigs.
//
// Hifiasm-style approach: each read becomes a node, edges are computed
// from minimizer-sketch Jaccard similarity. Linear chains in the read
// graph are collapsed into unitigs.

#include "types.hpp"
#include "fastq_index.hpp"
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
    // After dedup, keep only the `max_sketch_size` smallest hashes.
    // 0 = no cap. With 256, each sketch is a bottom-k MinHash signature
    // good enough for Jaccard estimation at variance ~1/k. Caps RSS at
    // WG scale (4.3M × 600 mins ≈ 21 GB → 4.3M × 256 = 8.8 GB).
    std::size_t max_sketch_size = 256;
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
    // sequentially during sketching only; the caller may free the vector
    // immediately after build() returns.
    void build(const std::vector<std::pair<ReadId, std::string>>& reads,
               const OverlapGraphOpts& opts);

    // Transitive reduction (Myers, 2005). Removes shortcut edges A→B if
    // a stronger detour A→C→B exists. After this pass the dense raw
    // overlap graph (avg deg ≈ 24) collapses to a near-string graph
    // (avg deg ≈ 2-3), which is what the unitig walker actually needs.
    //
    // Without this step, mutual-best top-1 + degree-2 unitig extraction
    // produces ~2 reads per unitig at WG scale (the v10[a-h] symptom of
    // 1.95M unitigs from 4.3M reads → only 217 bubbles found). With it,
    // we expect ~50K-200K unitigs and ~50K-150K bubbles.
    void transitive_reduce();

    // After build() (and ideally transitive_reduce()): collapse linear
    // chains into unitigs. Only read lengths are needed; FastqIndex
    // provides them at O(1).
    void collapse_to_unitigs(
        FastqIndex& reads,
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
