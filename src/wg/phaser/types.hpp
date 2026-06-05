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
//
// Important design principle (v0.10+): reads are NOT discarded as
// "unmappable" when current haplotype models can't place them. Reads
// that span multiple bubble-alts (CNV bridge), or carry sequence
// outside any known alt (novel deep sub-branch evidence), are tagged
// explicitly. The only tag that means "we don't know" is UNCERTAIN,
// and it should be empirically <1% of reads.
enum class ReadTag : std::uint8_t {
    UNTAGGED      = 0,   // initial state before phasing
    H1_ONLY       = 1,   // voter assigned to alt 0 of a HET / parent of BRANCH bubble
    H2_ONLY       = 2,   // voter assigned to alt 1 of a HET bubble
    SHARED        = 3,   // homozygous region — read carries no bubble signal
    BRANCH        = 4,   // assigned to a non-parent child alt of a BRANCH/BIB bubble
    UNCERTAIN     = 5,   // genuine failure: voter ambiguous AND no novel signal
    BRANCH_BRIDGE = 6,   // spans ≥2 alts of same bubble — CNV-junction evidence
    BRANCH_NOVEL  = 7,   // carries sequence not in any known alt — deep sub-branch candidate
    // Reserved for B-3 (molecular-clock validation): read descends from
    // a lineage whose intermediate sub-clones are extinct. VAF-derived
    // depth disagrees with mutations-along-branch-derived depth — the
    // "isotope decay chain" analogue: measured parent + daughter, but
    // the unstable intermediate isotopes already decayed away.
    BRANCH_EXTINCT_PARENT = 8,
};

inline const char* read_tag_str(ReadTag t) {
    switch (t) {
        case ReadTag::UNTAGGED:       return "untagged";
        case ReadTag::H1_ONLY:        return "h1_only";
        case ReadTag::H2_ONLY:        return "h2_only";
        case ReadTag::SHARED:         return "shared";
        case ReadTag::BRANCH:         return "branch";
        case ReadTag::UNCERTAIN:      return "uncertain";
        case ReadTag::BRANCH_BRIDGE:         return "branch_bridge";
        case ReadTag::BRANCH_NOVEL:          return "branch_novel";
        case ReadTag::BRANCH_EXTINCT_PARENT: return "branch_extinct_parent";
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

    // Path through the BranchTree of whichever bubble this read sits on.
    // Empty for reads outside bubbles (SHARED). For BRANCH/H1/H2 reads,
    // each entry is an index into Bubble.tree.nodes; the last entry is
    // the leaf the read most strongly supports. In B-1 paths are length
    // 1 (root → child); B-2 produces longer paths.
    std::vector<std::uint32_t> branch_path;
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
    Position      length_bp = 0;       // authoritative length; populated at
                                       // collapse time (cheap)
    std::string   seq;                 // OPTIONAL — empty after collapse;
                                       // gfa_writer materializes from
                                       // member_reads + reads on-the-fly,
                                       // one unitig at a time, to keep RSS
                                       // bounded at WG scale.
    std::vector<ReadId> member_reads;  // reads contributing to this unitig
    std::vector<NodeId> next_nodes;    // outgoing edges (forward strand)
    std::vector<NodeId> prev_nodes;    // incoming edges
    // The bubble this node belongs to, if any. UINT32_MAX = not in a bubble.
    BubbleId      bubble = ~0u;
};

// Topological classification of a bubble's allele structure based on the
// VAF distribution across its alt-paths. Drives downstream phasing logic:
// HET uses 2-vote MaxCut, BRANCH uses parent-attach via VAF arithmetic,
// BRANCH_IN_BRANCH recurses, etc.
//
// Patterns (with ±vaf_tol slack):
//   HET              [0.5, 0.5]
//   HOMOZ            [1.0]            single dominant allele
//   BRANCH           [0.5, 0.25, 0.25]   CNV duplication on one haplotype
//   BRANCH_IN_BRANCH [0.5, 0.25, 0.125, 0.125]   recursive duplication
//   TRIALLELIC       [0.33, 0.33, 0.33]  three distinct alleles, no CNV
//   NOISY            doesn't match any signature within tol
enum class BubbleType : std::uint8_t {
    UNKNOWN          = 0,
    HET              = 1,
    HOMOZ            = 2,
    BRANCH           = 3,
    BRANCH_IN_BRANCH = 4,
    TRIALLELIC       = 5,
    NOISY            = 6,
    // Reserved for B-3 (mol-clock): VAF cascade doesn't close because
    // an intermediate sub-clone is extinct — the isotope-decay-chain
    // analogue. We measure parent + final daughter, but the unstable
    // intermediates have all decayed. Classical example: two SHM
    // mutations co-occurring at 12.5% VAF with NO single-mutation
    // intermediate visible (Landau Cell 2013, Sottoriva Nature 2015).
    // The two mutations being linked at low VAF proves the intermediate
    // existed; the intermediate-only reads are absent because that
    // sub-clone died out.
    EXTINCT_INTERMEDIATE = 7,
};

inline const char* bubble_type_str(BubbleType t) {
    switch (t) {
        case BubbleType::UNKNOWN:          return "unknown";
        case BubbleType::HET:              return "het";
        case BubbleType::HOMOZ:            return "homoz";
        case BubbleType::BRANCH:           return "branch";
        case BubbleType::BRANCH_IN_BRANCH: return "branch_in_branch";
        case BubbleType::TRIALLELIC:       return "triallelic";
        case BubbleType::NOISY:            return "noisy";
        case BubbleType::EXTINCT_INTERMEDIATE:      return "extinct_intermediate";
    }
    return "?";
}

// Recursive sub-clonal structure within a single bubble's alt-path.
//
// We've historically classified bubbles by ONE top-level VAF pattern
// (e.g. [50, 25, 25] = HET founder + 2 sister sub-clones at depth-1).
// But a bubble's alt-path spans L bp of sequence, and over that length
// further sub-mutations accumulate at deeper depths. The right model is
// a **tree of branches**, where each branch may itself branch.
//
// Each `BranchNode` is one node in that tree:
//   - root (depth=0): the bubble's founder allele
//   - depth-K: a sub-clone K mitotic generations from founder, VAF ~ 2^-(K+1)
//
// In v0.10-α (Stufe B-1), bubble_voter only populates the root and its
// immediate children (one per alt_path). The recursive depth structure
// (B-2) is filled in by a later pass that walks each alt-path
// position-by-position and clusters variant positions by VAF + read
// co-occurrence. Hooks for B-2 are present already so callers can
// traverse `children` without re-checking the version.
struct BranchNode {
    float        vaf            = 0.0f;       // population frequency
    std::uint8_t depth          = 0;          // 0 = founder, 1 = first branch
    std::uint32_t alt_idx       = ~0u;        // which Bubble.alt_paths this is
    std::vector<Position> defining_positions; // variants on this branch
    std::vector<std::uint32_t> children;      // child indices into BranchTree.nodes
    std::vector<ReadId> supporting_reads;

    // Filled by mol-clock pass (B-3). Empty/zero in B-1.
    std::uint32_t n_substitutions          = 0;
    double        est_mitotic_generations  = 0.0;
};

struct BranchTree {
    std::vector<BranchNode> nodes;          // 0 = root by convention
    std::uint32_t           root_node_idx = 0;

    // Sequence-position-keyed VAF profile across this bubble's alt-paths.
    // Co-indexed: per_position_vaf[i] is the population frequency of the
    // minor allele at variant_positions_bp[i]. Empty in B-1; populated
    // by the per-position scan in B-2.
    std::vector<Position> variant_positions_bp;
    std::vector<float>    per_position_vaf;
};

// A bubble is a fan of alternative paths between two shared anchor nodes.
// Detected topologically by bubble_finder, classified by bubble_classifier
// (top-level VAF pattern), structured by bubble_voter (BranchTree).
struct Bubble {
    BubbleId  id;
    NodeId    source_anchor;   // last shared node before divergence
    NodeId    sink_anchor;     // first shared node after divergence
    std::vector<NodeId> alt_paths;  // 2 = standard diploid; >2 = multi-allelic
    bool      phased = false;       // set true once phase_engine resolves it

    // Filled in by bubble_classifier from per-path read counts. Parallel to
    // alt_paths. Empty until classification has run.
    std::vector<std::uint32_t> read_counts;
    std::vector<float>         vafs;        // sum to ~1
    BubbleType type = BubbleType::UNKNOWN;

    // Top-level branch depth: 0 = HET founder, 1 = first sub-clone,
    // 2 = sub-sub-clone, ... matches type to integer K where VAF of the
    // smallest sister-sub-clone ≈ 2^-(K+1).
    std::uint8_t branch_depth = 0;

    // Parent index into alt_paths for BRANCH-typed bubbles: the alt that
    // carries the unduplicated founder allele (highest VAF). UINT32_MAX
    // if not applicable.
    std::uint32_t parent_alt = ~0u;

    // Full clonal structure within this bubble. B-1 populates root + one
    // child per alt_path; B-2 fills in deeper sub-branches.
    BranchTree tree;
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
    std::size_t n_reads_bridge      = 0;  // BRANCH_BRIDGE (multi-alt)
    std::size_t n_reads_novel       = 0;  // BRANCH_NOVEL (off-alt sequence)
    std::size_t n_unitigs_shared    = 0;
    std::size_t n_unitigs_h1        = 0;
    std::size_t n_unitigs_h2        = 0;
    std::size_t n_unitigs_branch    = 0;
    std::size_t n_bubbles           = 0;
    std::size_t n_bubbles_phased    = 0;
    // Per-type bubble counts from bubble_classifier (v0.10).
    std::size_t n_bubbles_het              = 0;
    std::size_t n_bubbles_homoz            = 0;
    std::size_t n_bubbles_branch           = 0;
    std::size_t n_bubbles_branch_in_branch = 0;
    std::size_t n_bubbles_triallelic       = 0;
    std::size_t n_bubbles_noisy            = 0;
    int         build_iterations    = 0;
    int         refine_iterations   = 0;
    double      h1_bp = 0.0;
    double      h2_bp = 0.0;
    double      shared_bp = 0.0;
    double      branch_bp = 0.0;
};

}  // namespace branch::wg::phaser
