// BRANCH v0.2 — ReadPath population from overlap pairs.
//
// Design note (v0.2 best-effort):
// The graph builder currently allocates one Node per read and one Edge
// per OverlapPair. Node does NOT yet carry a consensus sequence, so a
// true base-level read-to-node trace is not possible. Instead, we emit
// a ReadPath per read that lists:
//   1. The read's own node (always present).
//   2. The nodes of all partner reads this read overlapped with,
//      ordered by the overlap offset along the read.
//
// This is useful for downstream VAF estimation and bubble-branch
// debugging. Deltas are left empty until Node gains a consensus
// sequence (v0.3).

#include "graph/delta_read.hpp"

#include <algorithm>
#include <limits>
#include <utility>

#include "backend/backend_vtable.hpp"

namespace branch::graph {

namespace {

// A single (offset_along_this_read, partner_node_id) entry used while
// accumulating per-read path information.
struct OffsetNode {
    std::uint32_t offset;
    NodeId node;
};

}  // namespace

void populate_read_paths(
    std::size_t n_reads,
    std::span<const NodeId> read_to_node,
    std::span<const std::uint32_t> read_lengths,
    std::span<const branch::backend::OverlapPair> overlaps,
    std::vector<ReadPath>& out) {
    out.clear();
    out.resize(n_reads);

    constexpr NodeId kNoNode = std::numeric_limits<NodeId>::max();

    // Initialise each path with the read's own node (if mapped) and its
    // recorded length.
    for (std::size_t rid = 0; rid < n_reads; ++rid) {
        auto& rp = out[rid];
        rp.read_id = static_cast<ReadId>(rid);
        rp.read_length = (rid < read_lengths.size()) ? read_lengths[rid] : 0u;
        rp.path.clear();
        rp.deltas.clear();

        if (rid < read_to_node.size()) {
            NodeId self = read_to_node[rid];
            if (self != kNoNode) {
                rp.path.push_back(self);
            }
        }
    }

    // Per-read partner offsets. Built as flat vectors keyed by read_id
    // to avoid per-read allocations of transient maps.
    std::vector<std::vector<OffsetNode>> per_read_partners(n_reads);

    for (const auto& op : overlaps) {
        const auto a = op.read_a;
        const auto b = op.read_b;

        if (a >= n_reads || b >= n_reads) continue;
        if (a >= read_to_node.size() || b >= read_to_node.size()) continue;

        const NodeId na = read_to_node[a];
        const NodeId nb = read_to_node[b];
        if (na == kNoNode || nb == kNoNode) continue;

        // From read A's perspective, partner is node B, anchored at offset_a.
        per_read_partners[a].push_back(OffsetNode{op.offset_a, nb});
        // And symmetrically for read B.
        per_read_partners[b].push_back(OffsetNode{op.offset_b, na});
    }

    // Order partner nodes along each read by their offset, then append
    // to the read's path. Duplicate partner node IDs are kept (they
    // mean multiple overlaps with the same read; downstream passes can
    // fold them).
    for (std::size_t rid = 0; rid < n_reads; ++rid) {
        auto& partners = per_read_partners[rid];
        std::sort(partners.begin(), partners.end(),
                  [](const OffsetNode& x, const OffsetNode& y) {
                      if (x.offset != y.offset) return x.offset < y.offset;
                      return x.node < y.node;
                  });
        auto& rp = out[rid];
        rp.path.reserve(rp.path.size() + partners.size());
        for (const auto& pn : partners) {
            rp.path.push_back(pn.node);
        }
        // TODO(v0.3): compute_delta(ref, read_sub) once Node has consensus.
    }
}

}  // namespace branch::graph
