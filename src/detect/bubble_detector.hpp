#pragma once

// BRANCH v0.2 — Bubble detection in a LosslessGraph.
//
// A "bubble" is a subgraph where exactly one entry node fans out into
// multiple parallel paths that reconverge at a single exit node. In
// BRANCH terms every bubble is a classification candidate: its
// parallel alternatives are what downstream cascade has to label as
// branch, duplication, mixed, or non-separable.
//
// v0.2 ships a structural detector only — it enumerates bubbles by
// local graph topology. It does not yet extract path sequences or
// coverage per alt-path; those arrive with the CN-aware pass.

#include <cstddef>
#include <cstdint>
#include <vector>

#include "graph/lossless_graph.hpp"

namespace branch::graph { class ReadDB; }

namespace branch::detect {

struct AltPath {
    std::vector<branch::graph::NodeId> nodes;  // intermediate nodes between entry and exit
    // Topological support: cumulative edge `read_support` along this alt
    // path (the legacy v0.2 metric — saturates at very low values
    // because edge supports are overlap-pair counts, not read counts).
    // Kept for backward-compat callers; use `reads_traversing` when a
    // ReadDB is available to count reads exactly.
    std::uint32_t total_read_support{};
    // Per-alt read count: number of reads in the supplied ReadDB whose
    // path contains the bubble entry, this alt's intermediate nodes
    // (in order), and the bubble exit. Empty / zero when detect_bubbles
    // was called without a ReadDB. This is the quantity downstream VAF
    // analysis should use: it counts reads, not topology events.
    std::vector<branch::graph::ReadId> reads_traversing;
};

struct Bubble {
    branch::graph::NodeId entry;
    branch::graph::NodeId exit;
    std::vector<AltPath> alts;
    std::uint32_t total_read_support{};
    // Total reads (union over all alts, deduplicated) that traverse
    // this bubble entry → exit. Filled when detect_bubbles is given a
    // ReadDB; left at 0 otherwise. Use this as the denominator when
    // computing per-alt VAF from `alts[i].reads_traversing.size()`.
    std::uint32_t total_reads_traversing{};
};

struct BubbleDetectorConfig {
    // Maximum number of intermediate nodes per alt path. Longer paths
    // are more likely structural variation than true bubbles; v0.2
    // hardcodes a conservative cap so detection stays O(V) per entry.
    std::uint32_t max_alt_path_length{8};

    // Cap on out-degree of entry nodes. A node with out-degree O has
    // O*(O-1)/2 successor pairs — above ~20 this is almost certainly an
    // assembly-artefact repeat hub, not a true bubble entry, and the
    // quadratic pair enumeration becomes the dominant runtime cost on
    // dense graphs (observed: 200+ outdegree nodes at IGH repeats
    // produce ~20k pairs each and hours of runtime). Entries exceeding
    // this cap are skipped. 0 disables the cap.
    std::uint32_t max_entry_out_degree{20};
};

// Detect bubbles in the graph. Returns a vector of Bubble records,
// each describing one entry+exit pair and the parallel alt paths
// between them.
[[nodiscard]] std::vector<Bubble> detect_bubbles(
    const branch::graph::LosslessGraph& graph,
    const BubbleDetectorConfig& cfg = {});

// Same, but additionally counts reads traversing each alt by
// intersecting with the supplied ReadDB. After this call every
// AltPath has `reads_traversing` populated with the read IDs whose
// graph path contains the entry, the alt's intermediate nodes (in
// order), and the exit. Use AltPath::reads_traversing.size() as the
// per-alt support for downstream VAF (denominator =
// Bubble::total_reads_traversing).
[[nodiscard]] std::vector<Bubble> detect_bubbles_with_reads(
    const branch::graph::LosslessGraph& graph,
    const branch::graph::ReadDB& reads,
    const BubbleDetectorConfig& cfg = {});

}  // namespace branch::detect
