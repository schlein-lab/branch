// BRANCH v0.1 — Graph builder implementation.

#include "graph/graph_builder.hpp"

#include <algorithm>
#include <limits>

namespace branch::graph {

BuildResult build_graph(std::span<const ReadMeta> reads,
                        std::span<const branch::backend::OverlapPair> overlaps) {
    BuildResult r{};

    // Determine mapping capacity: the largest read_id we'll see.
    ReadId max_id = 0;
    for (const auto& rm : reads) {
        if (rm.read_id > max_id) max_id = rm.read_id;
    }
    r.read_to_node.assign(static_cast<std::size_t>(max_id) + 1,
                          std::numeric_limits<NodeId>::max());

    // Allocate one node per read. Overlap-derived node attributes are
    // filled in later passes (unitig compaction, CN inference).
    //
    // read_support starts at 1: every raw-read node is, by definition,
    // supported by exactly one read (itself). The compactor sums these
    // across collapsed members so the post-compaction unitig's RC:i
    // reflects the real number of reads traversing it.
    for (const auto& rm : reads) {
        NodeId nid = r.graph.add_node(rm.length_bp);
        r.graph.node(nid).read_support = 1;
        r.read_to_node[rm.read_id] = nid;
    }

    // Translate overlap pairs to edges. Count a per-edge read_support
    // of 1 for each pair; later passes aggregate. Overlap direction
    // is implicit in the pair's offset_a/offset_b semantics — v0.1
    // treats every pair as a single undirected edge emitted as an
    // a→b directed edge, matching the LosslessGraph model.
    //
    // v0.2: propagate the strand bit from OverlapPair onto the Edge
    // via flag bit 4 (0x10). graph_io.cpp consumes that bit when
    // emitting GFA L-line orientation.
    for (const auto& op : overlaps) {
        if (op.read_a >= r.read_to_node.size() ||
            op.read_b >= r.read_to_node.size()) {
            continue;
        }
        NodeId na = r.read_to_node[op.read_a];
        NodeId nb = r.read_to_node[op.read_b];
        if (na == std::numeric_limits<NodeId>::max() ||
            nb == std::numeric_limits<NodeId>::max()) {
            continue;
        }
        r.graph.add_edge(na, nb, 1u);
        if (op.strand != 0) {
            auto& e = const_cast<Edge&>(r.graph.edges().back());
            e.flags = static_cast<std::uint16_t>(e.flags | 0x10);
        }
    }

    return r;
}

}  // namespace branch::graph
