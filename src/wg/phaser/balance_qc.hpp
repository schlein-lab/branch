#pragma once
//
// balance_qc — diagnose, never fix.
//
// We do NOT move reads to balance haplotype size. Instead we report
// regions where balance is off, with severity flags, so downstream
// users can investigate (artificial CNV? phasing error? real
// hemizygosity?). Honest output beats silent stuffing.

#include "types.hpp"
#include <vector>

namespace branch::wg::phaser {

struct BalanceQcOpts {
    double imbalance_warn = 0.05;    // 5%
    double imbalance_err  = 0.20;    // 20%
    double cov_uniform_warn_z = 2.0;
    double cov_uniform_err_z  = 4.0;

    // Window size for per-region checks. Smaller = finer flags but more
    // noise in low-coverage areas.
    Position region_window_bp = 1'000'000;  // 1 Mb
};

class BalanceQc {
public:
    // Read-only over unitigs/assignments. Emits per-region flags.
    void check(
        const std::vector<UnitigNode>& unitigs,
        const std::vector<Bubble>& bubbles,
        const std::vector<ReadAssignment>& assignments,
        std::vector<ImbalanceFlag>& out_flags,
        const BalanceQcOpts& opts = {});

    // Aggregate per-haplotype span in bp. Useful for the global summary
    // line that goes into stderr at the end of run().
    static void summarize(
        const std::vector<UnitigNode>& unitigs,
        double& out_h1_bp, double& out_h2_bp,
        double& out_shared_bp, double& out_branch_bp);
};

}  // namespace branch::wg::phaser
