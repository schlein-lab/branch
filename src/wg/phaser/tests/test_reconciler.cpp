// Unit test for Reconciler.
//
// Tests the v0.8 dual-phaser merge logic:
//   1. Both agree → AGREED, full confidence
//   2. Disagree, one weak → resolve to stronger
//   3. Disagree, both weak → FLAGGED
//   4. Disagree, both strong → resolve heuristically (until neural lands)

#include "wg/phaser/reconcile/reconciler.hpp"
#include "wg/phaser/types.hpp"

#include <cassert>
#include <iostream>

using namespace branch::wg::phaser;

namespace {

ReadAssignment mk(ReadId rid, ReadTag t, float conf, NodeId home = 0) {
    ReadAssignment a;
    a.read_id = rid;
    a.tag = t;
    a.confidence = conf;
    a.home_node = home;
    return a;
}

void test_agreement_keeps_consensus() {
    std::vector<ReadAssignment> a = {
        mk(0, ReadTag::H1_ONLY, 0.9f),
        mk(1, ReadTag::H2_ONLY, 0.8f),
        mk(2, ReadTag::SHARED,  0.7f),
    };
    std::vector<ReadAssignment> b = a;  // identical

    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts ropts;
    ropts.verbose = false;
    r.reconcile(a, b, out, ropts);
    assert(out.size() == 3);
    for (const auto& x : out) {
        assert(x.source == PhaseSource::AGREED);
        assert(x.tag_final == x.tag_a);
    }
    assert(r.stats().n_agreed == 3);
    assert(r.stats().n_flagged == 0);
}

void test_weak_disagreement_flags() {
    // Both phasers are uncertain → flag.
    std::vector<ReadAssignment> a = { mk(0, ReadTag::H1_ONLY, 0.3f) };
    std::vector<ReadAssignment> b = { mk(0, ReadTag::H2_ONLY, 0.4f) };
    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts ropts;
    ropts.verbose = false;
    r.reconcile(a, b, out, ropts);
    assert(out.size() == 1);
    assert(out[0].source == PhaseSource::FLAGGED);
    assert(out[0].tag_final == ReadTag::UNCERTAIN);
    assert(r.stats().n_flagged == 1);
}

void test_strong_disagreement_picks_stronger() {
    // One side is confident, other is weak → pick confident.
    std::vector<ReadAssignment> a = { mk(0, ReadTag::H1_ONLY, 0.95f) };
    std::vector<ReadAssignment> b = { mk(0, ReadTag::H2_ONLY, 0.30f) };
    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts ropts;
    ropts.verbose = false;
    r.reconcile(a, b, out, ropts);
    assert(out.size() == 1);
    assert(out[0].source == PhaseSource::NEURAL_A);
    assert(out[0].tag_final == ReadTag::H1_ONLY);
    assert(r.stats().n_neural_a == 1);
    assert(r.stats().n_flagged == 0);
}

void test_strong_disagreement_picks_b_when_higher() {
    std::vector<ReadAssignment> a = { mk(0, ReadTag::H1_ONLY, 0.55f) };
    std::vector<ReadAssignment> b = { mk(0, ReadTag::H2_ONLY, 0.95f) };
    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts ropts;
    ropts.verbose = false;
    r.reconcile(a, b, out, ropts);
    assert(out[0].source == PhaseSource::NEURAL_B);
    assert(out[0].tag_final == ReadTag::H2_ONLY);
}

void test_confusion_matrix_populated() {
    std::vector<ReadAssignment> a = {
        mk(0, ReadTag::H1_ONLY, 0.9f),
        mk(1, ReadTag::H1_ONLY, 0.9f),
        mk(2, ReadTag::SHARED,  0.9f),
    };
    std::vector<ReadAssignment> b = {
        mk(0, ReadTag::H1_ONLY, 0.9f),
        mk(1, ReadTag::H2_ONLY, 0.9f),
        mk(2, ReadTag::H1_ONLY, 0.9f),
    };
    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts ropts;
    ropts.verbose = false;
    r.reconcile(a, b, out, ropts);
    // Confusion: (1,1)=1, (1,2)=1, (3,1)=1
    assert(r.stats().confusion[1][1] == 1);
    assert(r.stats().confusion[1][2] == 1);
    assert(r.stats().confusion[3][1] == 1);
}

}  // namespace

int main() {
    test_agreement_keeps_consensus();
    std::cout << "agreement ok\n";
    test_weak_disagreement_flags();
    std::cout << "weak-disagreement → flag ok\n";
    test_strong_disagreement_picks_stronger();
    std::cout << "strong-disagreement-A ok\n";
    test_strong_disagreement_picks_b_when_higher();
    std::cout << "strong-disagreement-B ok\n";
    test_confusion_matrix_populated();
    std::cout << "confusion matrix ok\n";
    std::cout << "test_reconciler: OK\n";
    return 0;
}
