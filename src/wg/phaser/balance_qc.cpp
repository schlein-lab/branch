#include "balance_qc.hpp"

#include <cmath>
#include <sstream>
#include <unordered_map>

namespace branch::wg::phaser {

void BalanceQc::summarize(
    const std::vector<UnitigNode>& unitigs,
    double& out_h1_bp, double& out_h2_bp,
    double& out_shared_bp, double& out_branch_bp)
{
    out_h1_bp = out_h2_bp = out_shared_bp = out_branch_bp = 0.0;
    for (const auto& u : unitigs) {
        const double L = static_cast<double>(u.seq.size());
        switch (u.kind) {
            case NodeKind::SHARED:    out_shared_bp += L; break;
            case NodeKind::BUBBLE_H1: out_h1_bp     += L; break;
            case NodeKind::BUBBLE_H2: out_h2_bp     += L; break;
            case NodeKind::BRANCH:    out_branch_bp += L; break;
            case NodeKind::UNCLASSIFIED: break;
        }
    }
}

namespace {

FlagSeverity classify_imbalance(double h1_bp, double h2_bp,
                                double warn, double err) {
    const double mx = std::max(h1_bp, h2_bp);
    if (mx <= 0.0) return FlagSeverity::NONE;
    const double dev = std::abs(h1_bp - h2_bp) / mx;
    if (dev > err)  return FlagSeverity::ERROR;
    if (dev > warn) return FlagSeverity::WARNING;
    if (dev > warn / 2.0) return FlagSeverity::NOTICE;
    return FlagSeverity::NONE;
}

}  // namespace

void BalanceQc::check(
    const std::vector<UnitigNode>& unitigs,
    const std::vector<Bubble>& bubbles,
    const std::vector<ReadAssignment>& assignments,
    std::vector<ImbalanceFlag>& out_flags,
    const BalanceQcOpts& opts)
{
    out_flags.clear();
    (void)bubbles;  // currently summary-only; per-bubble flags added later
    (void)assignments;

    // Global summary flag.
    double g_h1, g_h2, g_shared, g_branch;
    summarize(unitigs, g_h1, g_h2, g_shared, g_branch);
    {
        ImbalanceFlag f;
        f.region_label = "GLOBAL";
        f.start_bp = 0;
        f.end_bp   = 0;
        const double total = g_h1 + g_h2 + 1e-9;
        f.h1_bp_share = g_h1 / total;
        f.h2_bp_share = g_h2 / total;
        f.severity = classify_imbalance(g_h1, g_h2,
                                        opts.imbalance_warn,
                                        opts.imbalance_err);
        if (f.severity != FlagSeverity::NONE) {
            std::ostringstream r;
            r << "h1=" << static_cast<long long>(g_h1)
              << " h2=" << static_cast<long long>(g_h2)
              << " shared=" << static_cast<long long>(g_shared)
              << " branch=" << static_cast<long long>(g_branch);
            f.reason = r.str();
            out_flags.push_back(std::move(f));
        }
    }

    // Per-region sliding flags. We use the unitig array order as a proxy
    // for genomic order; a real chromosome-aware grouping needs the
    // graph walker, which is gfa_writer's domain. For v1 this is "good
    // enough" — it catches gross local imbalances inside the unitig
    // stream.
    Position cursor = 0;
    Position window_start = 0;
    double  win_h1 = 0.0, win_h2 = 0.0, win_shared = 0.0;
    auto    flush_window = [&](Position end_bp) {
        if (end_bp <= window_start) return;
        ImbalanceFlag f;
        std::ostringstream lab;
        lab << "win_" << window_start << "_" << end_bp;
        f.region_label = lab.str();
        f.start_bp = window_start;
        f.end_bp   = end_bp;
        const double total = win_h1 + win_h2 + 1e-9;
        f.h1_bp_share = win_h1 / total;
        f.h2_bp_share = win_h2 / total;
        f.severity = classify_imbalance(win_h1, win_h2,
                                        opts.imbalance_warn,
                                        opts.imbalance_err);
        if (f.severity >= FlagSeverity::WARNING) {
            std::ostringstream r;
            r << "h1=" << static_cast<long long>(win_h1)
              << " h2=" << static_cast<long long>(win_h2)
              << " shared=" << static_cast<long long>(win_shared);
            f.reason = r.str();
            out_flags.push_back(std::move(f));
        }
    };

    for (const auto& u : unitigs) {
        const Position L = static_cast<Position>(u.seq.size());
        const double Ld = static_cast<double>(L);
        switch (u.kind) {
            case NodeKind::SHARED:    win_shared += Ld; break;
            case NodeKind::BUBBLE_H1: win_h1     += Ld; break;
            case NodeKind::BUBBLE_H2: win_h2     += Ld; break;
            default: break;
        }
        cursor += L;
        if (cursor - window_start >= opts.region_window_bp) {
            flush_window(cursor);
            window_start = cursor;
            win_h1 = win_h2 = win_shared = 0.0;
        }
    }
    flush_window(cursor);
}

}  // namespace branch::wg::phaser
