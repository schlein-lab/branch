#pragma once
//
// BRANCH v0.6 phaser — shared types.
//
// All sub-modules under src/wg/phaser/ exchange information through these
// typed structures. No sub-module reaches into another's internals.

#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace branch::wg::phaser {

using ReadId   = std::uint32_t;
using NodeId   = std::uint32_t;     // index into a UnitigNode array
using BubbleId = std::uint32_t;
using Position = std::int64_t;      // coordinate in bp; signed for ε arithmetic

// Per-read tag emitted by phase_engine. Each read receives EXACTLY ONE tag.
enum class ReadTag : std::uint8_t {
    UNTAGGED  = 0,    // initial state before phasing
    H1_ONLY   = 1,    // het-pair vote consistent with haplotype 1
    H2_ONLY   = 2,    // het-pair vote consistent with haplotype 2
    SHARED    = 3,    // homozygous region — no phase signal possible
    BRANCH    = 4,    // SV-bearing or third-allele read; emitted as branch_N
    UNCERTAIN = 5,    // signal conflict or below-threshold; flagged in output
};

inline const char* read_tag_str(ReadTag t) {
    switch (t) {
        case ReadTag::UNTAGGED:  return "untagged";
        case ReadTag::H1_ONLY:   return "h1_only";
        case ReadTag::H2_ONLY:   return "h2_only";
        case ReadTag::SHARED:    return "shared";
        case ReadTag::BRANCH:    return "branch";
        case ReadTag::UNCERTAIN: return "uncertain";
    }
    return "?";
}

// Every read carries a confidence in [0,1] alongside its tag.
// 1.0 = unanimous consensus among ≥10 supporting overlaps; 0.0 = forced.
struct ReadAssignment {
    ReadId   read_id;
    ReadTag  tag       = ReadTag::UNTAGGED;
    float    confidence = 0.0f;
    NodeId   home_node  = 0;       // unitig where this read sits in the graph
    std::uint32_t branch_idx = 0;  // only meaningful when tag == BRANCH
};

// One node in the assembly string graph. Built by overlap_graph,
// classified by bubble_finder, walked by gfa_writer.
enum class NodeKind : std::uint8_t {
    UNCLASSIFIED = 0,
    SHARED       = 1,    // emitted in both h1 and h2 walks
    BUBBLE_H1    = 2,    // emitted only in h1 walk
    BUBBLE_H2    = 3,    // emitted only in h2 walk
    BRANCH       = 4,    // SV unitig; not on either main walk
};

struct UnitigNode {
    NodeId        id;
    NodeKind      kind   = NodeKind::UNCLASSIFIED;
    std::string   seq;                 // unitig consensus sequence
    std::vector<ReadId> member_reads;  // reads contributing to this unitig
    std::vector<NodeId> next_nodes;    // outgoing edges (forward strand)
    std::vector<NodeId> prev_nodes;    // incoming edges
    // The bubble this node belongs to, if any. UINT32_MAX = not in a bubble.
    BubbleId      bubble = ~0u;
};

// A bubble is a pair (or rarely a fan) of alternative paths between two
// shared anchor nodes. Detected topologically by bubble_finder.
struct Bubble {
    BubbleId  id;
    NodeId    source_anchor;   // last shared node before divergence
    NodeId    sink_anchor;     // first shared node after divergence
    std::vector<NodeId> alt_paths;  // 2 = standard diploid; >2 = multi-allelic
    bool      phased = false;       // set true once phase_engine resolves it
};

// Coverage / size imbalance flag emitted by balance_qc.
enum class FlagSeverity : std::uint8_t {
    NONE    = 0,
    NOTICE  = 1,    // <5% deviation
    WARNING = 2,    // 5-20% deviation
    ERROR   = 3,    // >20% deviation; region marked unphased in output
};

struct ImbalanceFlag {
    std::string  region_label;     // e.g. "chr14:106Mb-107Mb" or unitig set
    Position     start_bp = 0;
    Position     end_bp   = 0;
    double       h1_bp_share = 0.0;
    double       h2_bp_share = 0.0;
    double       sum_coverage_z = 0.0;
    FlagSeverity severity = FlagSeverity::NONE;
    std::string  reason;
};

struct PhaserStats {
    std::size_t n_reads_total       = 0;
    std::size_t n_reads_h1          = 0;
    std::size_t n_reads_h2          = 0;
    std::size_t n_reads_shared      = 0;
    std::size_t n_reads_branch      = 0;
    std::size_t n_reads_uncertain   = 0;
    std::size_t n_unitigs_shared    = 0;
    std::size_t n_unitigs_h1        = 0;
    std::size_t n_unitigs_h2        = 0;
    std::size_t n_unitigs_branch    = 0;
    std::size_t n_bubbles           = 0;
    std::size_t n_bubbles_phased    = 0;
    int         build_iterations    = 0;
    int         refine_iterations   = 0;
    double      h1_bp = 0.0;
    double      h2_bp = 0.0;
    double      shared_bp = 0.0;
    double      branch_bp = 0.0;
};

}  // namespace branch::wg::phaser
