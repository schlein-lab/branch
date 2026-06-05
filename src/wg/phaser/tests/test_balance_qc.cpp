// Unit test for balance_qc — pure logic, no I/O.
//
// Validates that balanced and imbalanced unitig sets produce the
// expected severity flags.

#include "wg/phaser/balance_qc.hpp"

#include <cassert>
#include <iostream>
#include <string>

using namespace branch::wg::phaser;

namespace {

UnitigNode mk(NodeId id, NodeKind k, std::size_t len) {
    UnitigNode u;
    u.id = id;
    u.kind = k;
    u.length_bp = static_cast<Position>(len);
    return u;
}

void test_balanced_global_is_clean() {
    std::vector<UnitigNode> u = {
        mk(0, NodeKind::SHARED,    1'000'000),
        mk(1, NodeKind::BUBBLE_H1,   500'000),
        mk(2, NodeKind::BUBBLE_H2,   500'000),
    };
    std::vector<Bubble> bubbles;
    std::vector<ReadAssignment> a;
    std::vector<ImbalanceFlag> flags;
    BalanceQc qc;
    qc.check(u, bubbles, a, flags);
    bool any_global = false;
    for (const auto& f : flags) {
        if (f.region_label == "GLOBAL") any_global = true;
    }
    assert(!any_global && "balanced data should not produce a GLOBAL flag");
}

void test_imbalance_warning() {
    std::vector<UnitigNode> u = {
        mk(0, NodeKind::BUBBLE_H1, 1'000'000),
        mk(1, NodeKind::BUBBLE_H2,   850'000),  // 15% smaller
    };
    std::vector<Bubble> bubbles;
    std::vector<ReadAssignment> a;
    std::vector<ImbalanceFlag> flags;
    BalanceQc qc;
    qc.check(u, bubbles, a, flags);
    bool found = false;
    for (const auto& f : flags) {
        if (f.region_label == "GLOBAL" && f.severity == FlagSeverity::WARNING) {
            found = true;
        }
    }
    assert(found && "15% imbalance should produce a WARNING");
}

void test_imbalance_error() {
    std::vector<UnitigNode> u = {
        mk(0, NodeKind::BUBBLE_H1, 1'000'000),
        mk(1, NodeKind::BUBBLE_H2,   500'000),  // 50% smaller
    };
    std::vector<Bubble> bubbles;
    std::vector<ReadAssignment> a;
    std::vector<ImbalanceFlag> flags;
    BalanceQc qc;
    qc.check(u, bubbles, a, flags);
    bool found = false;
    for (const auto& f : flags) {
        if (f.region_label == "GLOBAL" && f.severity == FlagSeverity::ERROR) {
            found = true;
        }
    }
    assert(found && "50% imbalance should produce an ERROR");
}

void test_summarize_aggregates_correctly() {
    std::vector<UnitigNode> u = {
        mk(0, NodeKind::SHARED,    100),
        mk(1, NodeKind::BUBBLE_H1,  50),
        mk(2, NodeKind::BUBBLE_H2,  60),
        mk(3, NodeKind::BRANCH,     20),
    };
    double h1, h2, sh, br;
    BalanceQc::summarize(u, h1, h2, sh, br);
    assert(h1 == 50.0);
    assert(h2 == 60.0);
    assert(sh == 100.0);
    assert(br == 20.0);
}

}  // namespace

int main() {
    test_balanced_global_is_clean();
    test_imbalance_warning();
    test_imbalance_error();
    test_summarize_aggregates_correctly();
    std::cout << "test_balance_qc: OK\n";
    return 0;
}
