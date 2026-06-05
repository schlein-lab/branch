#pragma once
//
// bubble_voter — sequence-similarity read→alt assignment (v0.10 Variante B).
//
// Replaces home_node-based voting in phase_engine. Where the legacy path
// asked "is this read's home_node sitting on a bubble alt?" (filtering
// 99.7% of evidence at WG scale), bubble_voter asks the right question:
// "does this read's sequence k-mer set match this bubble alt's sequence
// k-mer set?" — and answers via Jaccard similarity.
//
// Algorithm:
//   1. Per bubble, per alt_path: materialize sequence by walking the
//      alt unitig chain and stitching member-read sequences (recycles
//      gfa_writer's build_unitig_seq helper).
//   2. Extract canonical k-mers from each alt sequence → kmer_set[bubble][alt].
//   3. Build an inverted index: kmer → list of (bubble_id, alt_idx) it
//      occurs in. Most bubble-specific k-mers occur in exactly one alt;
//      shared k-mers (alts that overlap) occur in multiple.
//   4. For each read, scan its canonical k-mers, accumulate per-(bubble,alt)
//      hits. Read's "support score" for alt i = #(kmer ∈ read ∩ alt[i]).
//   5. Assign read to argmax_i support[i] if best beats second-best by
//      a Jaccard-margin; else "ambiguous" and uncounted.
//   6. Emit per-bubble per-alt read counts → VAFs.
//
// This decouples voting from the per-read home_node and lets every read
// with any k-mer in any alt contribute. The full v16 measurement was 12.7K
// reads voting; the voter should pull in 1-2M.
//
// B-2 hooks: bubble_voter's `BranchTree` output stays flat (root + one
// child per alt) in this version. The per-position VAF cascade that
// resolves sub-branches inside an alt is a follow-up pass.

#include "types.hpp"
#include "fastq_index.hpp"
#include <cstdint>
#include <vector>

namespace branch::wg::phaser {

struct BubbleVoterOpts {
    // k-mer size for similarity matching. Should match the upstream
    // het-pair k-mer choice for consistency; 21 is the v0.10 default.
    int k = 21;
    // Approximate overlap between consecutive reads when stitching an
    // alt sequence (must match overlap_graph collapse rule).
    int approx_overlap_bp = 5 * 23;
    // A read is assigned to alt i only if its support count beats the
    // second-best alt by at least (1 + margin) × second_best. Reads
    // failing this go to bridge_reads, NOT dropped.
    float assign_margin = 0.10f;
    // Reads with fewer than this many alt-k-mer hits across ALL alts
    // are flagged as BRANCH_NOVEL if they had >= min_alt_hits_for_novel
    // hits anywhere (i.e. they overlap a bubble region but carry mostly
    // novel sequence). Otherwise truly off-target → no voter output.
    std::uint32_t min_alt_hits = 3;
    // Fraction of a read's k-mers that match SOME alt's k-mers. Below
    // this fraction, AND with at least min_alt_hits, the read is
    // BRANCH_NOVEL — overlaps a bubble region with deep-sub-clone
    // sequence not yet captured by any known alt. v19 chr14 showed 0
    // NOVEL with threshold=0.10 because most reads have some kmer-hits
    // (avg 92%). Raised to 0.30 so any read with <30% kmer-match to any
    // alt gets BRANCH_NOVEL — i.e. read overlaps bubble region but
    // carries mostly-foreign sequence, possibly an extinct-intermediate
    // descendant.
    float novel_match_frac_threshold = 0.30f;
    // Heartbeat period for stderr logging.
    std::size_t log_every = 100000;
};

// Per-bubble vote tallies and resulting VAF vector. Parallel to
// Bubble.alt_paths. Populated by `vote()` below.
struct BubbleVote {
    std::vector<std::uint32_t> read_counts;     // per alt
    std::vector<float>         vafs;             // per alt, sums to ~1
    // Reads attributed to each alt. Each Read appears in exactly one
    // alt's list (its argmax-match). Empty for unassigned reads.
    std::vector<std::vector<ReadId>> reads_per_alt;
    // Reads spanning ≥2 alts of THIS bubble with comparable hit counts.
    // Biologically: CNV-junction reads / SV-bridge reads.
    std::vector<ReadId>        bridge_reads;
};

// Per-read voter result spanning all bubbles. Each read appears in
// exactly one of: assigned (to some bubble alt), bridge (to some bubble),
// novel (off-alt evidence), or absent (truly outside all bubbles → SHARED).
struct VoterReadOutcome {
    ReadId   read_id;
    enum Kind : std::uint8_t {
        ASSIGNED = 0,   // 1+ bubble assignments → tag from alt's kind
        BRIDGE   = 1,   // 1+ bubble bridge memberships → BRANCH_BRIDGE
        NOVEL    = 2,   // overlaps bubble region, mostly novel sequence
        ABSENT   = 3,   // no hits anywhere → out-of-bubble (SHARED)
    } kind = ABSENT;
    // Primary bubble + alt for ASSIGNED reads. UINT32_MAX otherwise.
    std::uint32_t primary_bubble = ~0u;
    std::uint8_t  primary_alt    = 0xFF;
};

class BubbleVoter {
public:
    // Materialize alt sequences + run voting + emit per-bubble VAFs.
    //
    // `reads` is needed to materialize alt sequences from member_reads
    // AND to scan each read's own k-mer set during voting. The FastqIndex
    // is opened twice in this call — once per pass — to keep peak RSS
    // bounded.
    //
    // The voter does NOT drop "ambiguous" or "noisy" reads. Every read
    // gets an outcome in `out_read_outcomes` (parallel to FastqIndex
    // read order): ASSIGNED, BRIDGE, NOVEL, or ABSENT. Phase_engine
    // bridges these to ReadTag.
    void vote(
        const std::vector<UnitigNode>& unitigs,
        std::vector<Bubble>& bubbles,
        FastqIndex& reads,
        std::vector<BubbleVote>& out_votes,
        std::vector<VoterReadOutcome>& out_read_outcomes,
        const BubbleVoterOpts& opts = {});

    // After vote(), each Bubble.tree gets a root + one child per alt,
    // populated from read_counts. Use this to populate `tree` for
    // downstream phase_engine consumption.
    static void populate_branch_trees(
        std::vector<Bubble>& bubbles,
        const std::vector<BubbleVote>& votes);
};

}  // namespace branch::wg::phaser
