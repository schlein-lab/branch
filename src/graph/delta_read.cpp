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
#include <cstdlib>
#include <limits>
#include <string>
#include <utility>
#include <vector>

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
        // deltas stay empty here — populate_read_paths has no access to
        // the read or node-consensus sequences. Callers that want real
        // deltas must run encode_delta(read_seq, concat-of-node-consensus)
        // separately; see encode_delta / reconstruct_read below.
    }
}

// ---------------------------------------------------------------------
// Delta encoding for read reconstruction
// ---------------------------------------------------------------------
//
// encode_delta: returns ops that describe how read_seq differs from
// path_seq, with a hard cap of kMaxOps to avoid pathological cases
// (very divergent pairs or O(n*m) DP memory blow-up).
//
// Strategy:
//  1. Same-length fast path: walk positions, emit 'S' for each mismatch.
//  2. Different-length fallback: banded DP with band = 64 around the
//     diagonal. Traceback emits S/I/D ops. HiFi reads are ~99% identity,
//     so a narrow band is sufficient.
//  3. If the fast path or banded DP yields more than kMaxOps deltas or
//     cannot cover the strings, return an empty vector — the caller must
//     treat this as "delta encoding unavailable, store raw".

namespace {
constexpr std::size_t kMaxOps = 1000;
constexpr int kBand = 64;  // banded DP radius
}  // namespace

std::vector<DeltaOp> encode_delta(std::string_view read_seq,
                                  std::string_view path_seq) {
    std::vector<DeltaOp> ops;

    // Fast path: same length -> only substitutions.
    if (read_seq.size() == path_seq.size()) {
        for (std::size_t i = 0; i < read_seq.size(); ++i) {
            if (read_seq[i] != path_seq[i]) {
                ops.push_back(DeltaOp{
                    .pos = static_cast<std::uint32_t>(i),
                    .op = 'S',
                    .base = read_seq[i],
                });
                if (ops.size() > kMaxOps) return {};
            }
        }
        return ops;
    }

    const int n = static_cast<int>(path_seq.size());
    const int m = static_cast<int>(read_seq.size());
    if (std::abs(n - m) > kBand) return {};  // too divergent for banded DP

    // Banded DP over (i=0..n, j=0..m) with |i-j| <= kBand.
    // dp[i][j] = min edit ops to transform path[0..i] to read[0..j].
    // Stored as (kBand*2+1)-wide row per i.
    const int width = 2 * kBand + 1;
    const int kInf = std::numeric_limits<int>::max() / 2;
    std::vector<int> dp((n + 1) * width, kInf);

    auto at = [&](int i, int j) -> int& {
        const int offset = j - i + kBand;
        if (offset < 0 || offset >= width) {
            static int sink = kInf;
            sink = kInf;
            return sink;
        }
        return dp[i * width + offset];
    };

    at(0, 0) = 0;
    for (int j = 1; j <= std::min(m, kBand); ++j) at(0, j) = j;
    for (int i = 1; i <= std::min(n, kBand); ++i) at(i, 0) = i;

    for (int i = 1; i <= n; ++i) {
        const int jlo = std::max(1, i - kBand);
        const int jhi = std::min(m, i + kBand);
        for (int j = jlo; j <= jhi; ++j) {
            const int cost = (path_seq[i - 1] == read_seq[j - 1]) ? 0 : 1;
            const int sub = at(i - 1, j - 1) + cost;
            const int del = (j - (i - 1) + kBand >= 0 &&
                             j - (i - 1) + kBand < width)
                              ? at(i - 1, j) + 1
                              : kInf;
            const int ins = (j - 1 - i + kBand >= 0 &&
                             j - 1 - i + kBand < width)
                              ? at(i, j - 1) + 1
                              : kInf;
            at(i, j) = std::min({sub, del, ins});
        }
    }

    if (at(n, m) >= kInf) return {};

    // Traceback from (n, m) to (0, 0).
    int i = n;
    int j = m;
    std::vector<DeltaOp> reversed;
    while (i > 0 || j > 0) {
        if (reversed.size() > kMaxOps) return {};
        const int here = at(i, j);
        if (i > 0 && j > 0) {
            const int cost = (path_seq[i - 1] == read_seq[j - 1]) ? 0 : 1;
            if (here == at(i - 1, j - 1) + cost) {
                if (cost == 1) {
                    reversed.push_back(DeltaOp{
                        .pos = static_cast<std::uint32_t>(i - 1),
                        .op = 'S',
                        .base = read_seq[j - 1],
                    });
                }
                --i;
                --j;
                continue;
            }
        }
        if (i > 0 && here == at(i - 1, j) + 1) {
            reversed.push_back(DeltaOp{
                .pos = static_cast<std::uint32_t>(i - 1),
                .op = 'D',
                .base = path_seq[i - 1],
            });
            --i;
            continue;
        }
        if (j > 0 && here == at(i, j - 1) + 1) {
            reversed.push_back(DeltaOp{
                .pos = static_cast<std::uint32_t>(i),
                .op = 'I',
                .base = read_seq[j - 1],
            });
            --j;
            continue;
        }
        // Unreachable in a consistent DP; guard against bugs.
        return {};
    }
    std::reverse(reversed.begin(), reversed.end());
    return reversed;
}

std::string reconstruct_read(std::string_view path_seq,
                             std::span<const DeltaOp> delta) {
    // Walk path_seq left-to-right, applying deltas in order. Deltas are
    // assumed sorted by pos (encode_delta emits in traceback order,
    // which is ascending).
    std::string out;
    out.reserve(path_seq.size() + delta.size());
    std::size_t p = 0;
    for (const auto& d : delta) {
        while (p < d.pos && p < path_seq.size()) {
            out.push_back(path_seq[p++]);
        }
        switch (d.op) {
            case 'S':
                if (p < path_seq.size()) {
                    out.push_back(d.base);
                    ++p;
                }
                break;
            case 'I':
                out.push_back(d.base);
                break;
            case 'D':
                if (p < path_seq.size()) ++p;
                break;
            default:
                break;
        }
    }
    while (p < path_seq.size()) {
        out.push_back(path_seq[p++]);
    }
    return out;
}

}  // namespace branch::graph
