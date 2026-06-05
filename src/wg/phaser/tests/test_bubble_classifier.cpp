// Unit tests for BubbleClassifier (v0.10).
//
// Pattern-matching of VAF distributions to bubble types:
//   HET              [0.5, 0.5]
//   HOMOZ            [1.0]
//   BRANCH           [0.5, 0.25, 0.25]   CNV duplication on one hap
//   BRANCH_IN_BRANCH [0.5, 0.25, 0.125, 0.125]
//   TRIALLELIC       [0.33, 0.33, 0.33]
//   NOISY            doesn't match anything within tol

#include "wg/phaser/bubble_classifier.hpp"
#include "wg/phaser/types.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using namespace branch::wg::phaser;

namespace {

void test_classify_vaf_het() {
    auto t = classify_vaf_vector({0.5f, 0.5f}, 0.07f);
    assert(t == BubbleType::HET);
}

void test_classify_vaf_branch() {
    auto t = classify_vaf_vector({0.5f, 0.25f, 0.25f}, 0.07f);
    assert(t == BubbleType::BRANCH);
}

void test_classify_vaf_branch_slack() {
    // Within tol — real read counts won't be exactly 0.5/0.25/0.25.
    auto t = classify_vaf_vector({0.52f, 0.26f, 0.22f}, 0.07f);
    assert(t == BubbleType::BRANCH);
}

void test_classify_vaf_branch_in_branch() {
    auto t = classify_vaf_vector({0.5f, 0.25f, 0.125f, 0.125f}, 0.07f);
    assert(t == BubbleType::BRANCH_IN_BRANCH);
}

void test_classify_vaf_triallelic() {
    auto t = classify_vaf_vector({1.0f / 3, 1.0f / 3, 1.0f / 3}, 0.07f);
    assert(t == BubbleType::TRIALLELIC);
}

void test_classify_vaf_homoz() {
    auto t = classify_vaf_vector({1.0f}, 0.07f);
    assert(t == BubbleType::HOMOZ);
}

void test_classify_vaf_noisy() {
    auto t = classify_vaf_vector({0.7f, 0.3f}, 0.07f);
    assert(t == BubbleType::NOISY);
}

void test_classify_vaf_noisy_4way_uniform() {
    // 4-way uniform [0.25,0.25,0.25,0.25] doesn't match any canonical
    // pattern (no 0.5 anchor, no triallele) → NOISY.
    auto t = classify_vaf_vector({0.25f, 0.25f, 0.25f, 0.25f}, 0.07f);
    assert(t == BubbleType::NOISY);
}

// End-to-end: build minimal unitig graph with one bubble, populate
// assignments such that 4 reads land on alt 0 and 2+2 reads on alts 1 & 2.
// Expect BRANCH classification with parent_alt = 0.
void test_end_to_end_branch_pattern() {
    // Graph: 0 (source anchor) → {1, 2, 3} → 4 (sink)
    std::vector<UnitigNode> u(5);
    for (NodeId i = 0; i < 5; ++i) { u[i].id = i; u[i].length_bp = 100; }
    u[0].next_nodes = {1, 2, 3};
    u[1].next_nodes = {4};  u[2].next_nodes = {4};  u[3].next_nodes = {4};
    u[1].prev_nodes = {0};  u[2].prev_nodes = {0};  u[3].prev_nodes = {0};

    Bubble b;
    b.id = 0;
    b.source_anchor = 0;
    b.sink_anchor = 4;
    b.alt_paths = {1, 2, 3};
    std::vector<Bubble> bubbles = {b};

    // 8 reads total: 4 on alt 0 (node 1), 2 on alt 1 (node 2), 2 on alt 2 (node 3)
    // → VAFs [0.5, 0.25, 0.25] = BRANCH
    std::vector<ReadAssignment> as;
    auto mk = [](ReadId rid, NodeId home) {
        ReadAssignment a; a.read_id = rid; a.home_node = home;
        a.tag = ReadTag::UNTAGGED; a.confidence = 0; return a;
    };
    for (int i = 0; i < 4; ++i) as.push_back(mk(i, 1));      // alt 0
    for (int i = 4; i < 6; ++i) as.push_back(mk(i, 2));      // alt 1
    for (int i = 6; i < 8; ++i) as.push_back(mk(i, 3));      // alt 2

    BubbleClassifier bc;
    BubbleClassifierOpts opts;
    opts.min_total_reads = 4;  // 8 total, easily passes
    opts.log_every = 0;
    bc.classify(u, as, bubbles, opts);

    assert(bubbles[0].type == BubbleType::BRANCH);
    assert(bubbles[0].read_counts.size() == 3);
    assert(bubbles[0].read_counts[0] == 4);
    assert(bubbles[0].read_counts[1] == 2);
    assert(bubbles[0].read_counts[2] == 2);
    assert(bubbles[0].parent_alt == 0);
}

void test_end_to_end_het_pattern() {
    // 0 → {1, 2} → 3, 4 reads each → [0.5, 0.5] = HET
    std::vector<UnitigNode> u(4);
    for (NodeId i = 0; i < 4; ++i) { u[i].id = i; u[i].length_bp = 100; }
    u[0].next_nodes = {1, 2};
    u[1].next_nodes = {3};  u[2].next_nodes = {3};
    u[1].prev_nodes = {0};  u[2].prev_nodes = {0};

    Bubble b;
    b.id = 0;
    b.source_anchor = 0;
    b.sink_anchor = 3;
    b.alt_paths = {1, 2};
    std::vector<Bubble> bubbles = {b};

    std::vector<ReadAssignment> as;
    auto mk = [](ReadId rid, NodeId home) {
        ReadAssignment a; a.read_id = rid; a.home_node = home;
        a.tag = ReadTag::UNTAGGED; a.confidence = 0; return a;
    };
    for (int i = 0; i < 4; ++i) as.push_back(mk(i, 1));
    for (int i = 4; i < 8; ++i) as.push_back(mk(i, 2));

    BubbleClassifier bc;
    BubbleClassifierOpts opts;
    opts.min_total_reads = 4;
    opts.log_every = 0;
    bc.classify(u, as, bubbles, opts);

    assert(bubbles[0].type == BubbleType::HET);
    assert(bubbles[0].parent_alt == ~0u);  // not set for HET
}

void test_end_to_end_below_threshold_is_unknown() {
    std::vector<UnitigNode> u(4);
    for (NodeId i = 0; i < 4; ++i) { u[i].id = i; u[i].length_bp = 100; }
    u[0].next_nodes = {1, 2};
    u[1].next_nodes = {3};  u[2].next_nodes = {3};
    u[1].prev_nodes = {0};  u[2].prev_nodes = {0};

    Bubble b;
    b.id = 0; b.source_anchor = 0; b.sink_anchor = 3;
    b.alt_paths = {1, 2};
    std::vector<Bubble> bubbles = {b};

    // Only 3 reads total; below default min_total_reads=10
    std::vector<ReadAssignment> as;
    auto mk = [](ReadId rid, NodeId home) {
        ReadAssignment a; a.read_id = rid; a.home_node = home;
        a.tag = ReadTag::UNTAGGED; a.confidence = 0; return a;
    };
    as.push_back(mk(0, 1));
    as.push_back(mk(1, 1));
    as.push_back(mk(2, 2));

    BubbleClassifier bc;
    BubbleClassifierOpts opts;
    opts.log_every = 0;
    bc.classify(u, as, bubbles, opts);
    assert(bubbles[0].type == BubbleType::UNKNOWN);
}

}  // namespace

int main() {
    test_classify_vaf_het();
    test_classify_vaf_branch();
    test_classify_vaf_branch_slack();
    test_classify_vaf_branch_in_branch();
    test_classify_vaf_triallelic();
    test_classify_vaf_homoz();
    test_classify_vaf_noisy();
    test_classify_vaf_noisy_4way_uniform();
    test_end_to_end_branch_pattern();
    test_end_to_end_het_pattern();
    test_end_to_end_below_threshold_is_unknown();
    std::cout << "test_bubble_classifier: all tests passed\n";
    return 0;
}
