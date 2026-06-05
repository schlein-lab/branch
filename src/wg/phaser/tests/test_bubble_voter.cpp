// Unit tests for BubbleVoter (v0.10 Variante B).
//
// Tests the clonally branched haplotype lineage voter:
//   - Pass 1: alt-seq materialization + kmer index build
//   - Pass 2: read→alt assignment via Jaccard
//   - Pass 3: VAF computation
//   - populate_branch_trees: flat tree (root + children) construction

#include "wg/phaser/bubble_voter.hpp"
#include "wg/phaser/fastq_index.hpp"
#include "wg/phaser/types.hpp"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace branch::wg::phaser;

namespace {

// Populate `idx` from a list of read sequences in-place. Each read gets
// ReadId = its index in the input vector.
void build_fastq(FastqIndex& idx, const std::vector<std::string>& seqs) {
    std::vector<std::pair<ReadId, std::string>> reads_vec;
    reads_vec.reserve(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        reads_vec.emplace_back(static_cast<ReadId>(i), seqs[i]);
    }
    idx.build_from_memory(std::move(reads_vec), {});
}

// Mini-graph with one bubble: 0 (source) → {1, 2} → 3 (sink).
// Unitig 1's sequence comes from reads[0]; unitig 2's from reads[1].
// Reads 2..N+1 are the "test reads" we feed into the voter — their
// match against alt 1 vs alt 2 determines assignment.
//
// We construct a HET-shaped scenario: alt 1 sequence ≠ alt 2 sequence,
// half the test reads carry alt 1's k-mers, half carry alt 2's.
void test_het_50_50() {
    // Two divergent sequences sharing prefix/suffix but with a unique
    // 21-bp k-mer in the middle that distinguishes them.
    const std::string base_left  = "ACGTACGTACGTACGTACGTACGTAC";  // 26 bp
    const std::string base_right = "GGGGGGGGGGGGGGGGGGGGGGGGGG";  // 26 bp
    const std::string alt1_mid   = "AAAAAATTTTTTAAAAAATTTTTT";    // 24 bp, unique
    const std::string alt2_mid   = "CCCCCCGGGGGGCCCCCCGGGGGG";    // 24 bp, unique
    const std::string alt1_seq = base_left + alt1_mid + base_right;
    const std::string alt2_seq = base_left + alt2_mid + base_right;

    // Reads:
    //   0 = alt1 seq (member of unitig 1)
    //   1 = alt2 seq (member of unitig 2)
    //   2..5 = 4 test reads similar to alt1
    //   6..9 = 4 test reads similar to alt2
    std::vector<std::string> seqs = {
        alt1_seq, alt2_seq,
        alt1_seq, alt1_seq, alt1_seq, alt1_seq,
        alt2_seq, alt2_seq, alt2_seq, alt2_seq,
    };
    FastqIndex idx;
    build_fastq(idx, seqs);

    // Build unitigs: 4 nodes (source, alt1, alt2, sink)
    std::vector<UnitigNode> unitigs(4);
    for (NodeId i = 0; i < 4; ++i) {
        unitigs[i].id = i;
        unitigs[i].length_bp = 100;
    }
    unitigs[0].next_nodes = {1, 2};
    unitigs[1].next_nodes = {3};
    unitigs[2].next_nodes = {3};
    unitigs[1].prev_nodes = {0};
    unitigs[2].prev_nodes = {0};
    // alt unitigs' sequences come from their first member_read
    unitigs[1].member_reads = {0};
    unitigs[2].member_reads = {1};

    Bubble b;
    b.id = 0;
    b.source_anchor = 0;
    b.sink_anchor = 3;
    b.alt_paths = {1, 2};
    std::vector<Bubble> bubbles = {b};

    BubbleVoter v;
    BubbleVoterOpts opts;
    opts.k = 21;
    opts.approx_overlap_bp = 0;
    opts.min_alt_hits = 1;
    opts.assign_margin = 0.10f;
    opts.log_every = 0;

    std::vector<BubbleVote> votes;
    std::vector<VoterReadOutcome> outcomes;
    v.vote(unitigs, bubbles, idx, votes, outcomes, opts);

    assert(votes.size() == 1);
    assert(votes[0].read_counts.size() == 2);
    // Expect alt 1 to absorb reads 0,2,3,4,5 = 5 reads; alt 2 the other 5.
    // (Reads 0 and 1 are the "founder" alt sequences themselves but also
    // get scanned and assigned to their own alt.)
    std::cerr << "[test_het] alt1_count=" << votes[0].read_counts[0]
              << " alt2_count=" << votes[0].read_counts[1]
              << " vafs=" << votes[0].vafs[0] << "," << votes[0].vafs[1] << "\n";
    assert(votes[0].read_counts[0] == 5);
    assert(votes[0].read_counts[1] == 5);
    // VAFs ≈ 0.5 / 0.5
    assert(votes[0].vafs[0] > 0.45f && votes[0].vafs[0] < 0.55f);
    assert(votes[0].vafs[1] > 0.45f && votes[0].vafs[1] < 0.55f);
}

void test_branch_50_25_25() {
    // Three alts, one founder + two sister sub-clones.
    const std::string base_left  = "ACGTACGTACGTACGTACGTACGTAC";
    const std::string base_right = "GGGGGGGGGGGGGGGGGGGGGGGGGG";
    const std::string p_mid = "AAAAAATTTTTTAAAAAATTTTTT";   // parent (founder)
    const std::string c1_mid = "CCCCCCGGGGGGCCCCCCGGGGGG";  // sub-clone 1
    const std::string c2_mid = "TGTGTGAATCATGAATGTGTGCAT";  // sub-clone 2

    const std::string parent_seq   = base_left + p_mid  + base_right;
    const std::string child1_seq   = base_left + c1_mid + base_right;
    const std::string child2_seq   = base_left + c2_mid + base_right;

    // 8 reads total, target VAFs 0.5/0.25/0.25:
    //   reads 0..3 → parent (4 reads = 50%)
    //   reads 4,5 → child1 (2 reads = 25%)
    //   reads 6,7 → child2 (2 reads = 25%)
    std::vector<std::string> seqs = {
        parent_seq, parent_seq, parent_seq, parent_seq,
        child1_seq, child1_seq,
        child2_seq, child2_seq,
    };
    FastqIndex idx;
    build_fastq(idx, seqs);

    // 5 unitigs: source(0) → {parent(1), child1(2), child2(3)} → sink(4)
    std::vector<UnitigNode> unitigs(5);
    for (NodeId i = 0; i < 5; ++i) {
        unitigs[i].id = i; unitigs[i].length_bp = 100;
    }
    unitigs[0].next_nodes = {1, 2, 3};
    unitigs[1].next_nodes = {4};
    unitigs[2].next_nodes = {4};
    unitigs[3].next_nodes = {4};
    unitigs[1].prev_nodes = {0};
    unitigs[2].prev_nodes = {0};
    unitigs[3].prev_nodes = {0};
    unitigs[1].member_reads = {0};
    unitigs[2].member_reads = {4};
    unitigs[3].member_reads = {6};

    Bubble b;
    b.id = 0;
    b.source_anchor = 0;
    b.sink_anchor = 4;
    b.alt_paths = {1, 2, 3};
    std::vector<Bubble> bubbles = {b};

    BubbleVoter v;
    BubbleVoterOpts opts;
    opts.k = 21;
    opts.approx_overlap_bp = 0;
    opts.min_alt_hits = 1;
    opts.assign_margin = 0.10f;
    opts.log_every = 0;

    std::vector<BubbleVote> votes;
    std::vector<VoterReadOutcome> outcomes;
    v.vote(unitigs, bubbles, idx, votes, outcomes, opts);

    std::cerr << "[test_branch] counts="
              << votes[0].read_counts[0] << "/"
              << votes[0].read_counts[1] << "/"
              << votes[0].read_counts[2]
              << " vafs="
              << votes[0].vafs[0] << "/"
              << votes[0].vafs[1] << "/"
              << votes[0].vafs[2] << "\n";
    assert(votes[0].read_counts[0] == 4);  // parent
    assert(votes[0].read_counts[1] == 2);  // child1
    assert(votes[0].read_counts[2] == 2);  // child2
    // VAFs ≈ 0.5 / 0.25 / 0.25
    assert(votes[0].vafs[0] > 0.45f && votes[0].vafs[0] < 0.55f);
    assert(votes[0].vafs[1] > 0.20f && votes[0].vafs[1] < 0.30f);
    assert(votes[0].vafs[2] > 0.20f && votes[0].vafs[2] < 0.30f);

    // populate_branch_trees: root + 3 children
    BubbleVoter::populate_branch_trees(bubbles, votes);
    assert(bubbles[0].tree.nodes.size() == 4);
    assert(bubbles[0].tree.nodes[0].depth == 0);
    assert(bubbles[0].tree.nodes[0].children.size() == 3);
    assert(bubbles[0].tree.nodes[1].depth == 1);
    assert(bubbles[0].tree.nodes[1].alt_idx == 0);
    assert(bubbles[0].tree.nodes[2].alt_idx == 1);
    assert(bubbles[0].tree.nodes[3].alt_idx == 2);
}

}  // namespace

int main() {
    test_het_50_50();
    test_branch_50_25_25();
    std::cout << "test_bubble_voter: all tests passed\n";
    return 0;
}
