#pragma once
//
// Reconciler — merge two independent phaser outputs.
//
// Takes assignments from de-novo phaser (track A) and pangenome-anchored
// phaser (track B). For each read:
//   - if both agree: AGREED, full confidence
//   - if they disagree and a neural voter resolves it: NEURAL_*, model confidence
//   - if disagreement remains unresolved: FLAGGED (UNCERTAIN with provenance)
//
// The neural voter is optional. Without it, all disagreements become
// FLAGGED — still useful information, since FLAGGED regions are exactly
// where novel structural variants / mosaic / cancer-CNV live.

#include "../types.hpp"
#include "../types_dual.hpp"
#include "mlp.hpp"
#include <array>
#include <string>
#include <vector>

namespace branch::wg::phaser {

struct ReconcilerOpts {
    // Optional neural voter model. Empty path = no neural voting; all
    // disagreements stay FLAGGED.
    std::string neural_model_path;

    // Disagreement statistics to stderr.
    bool verbose = true;
};

struct ReconcilerStats {
    std::size_t n_total = 0;
    std::size_t n_agreed = 0;
    std::size_t n_neural_a = 0;      // model sided with de-novo
    std::size_t n_neural_b = 0;      // model sided with anchored
    std::size_t n_heuristic_a = 0;   // confidence-fallback sided with de-novo
    std::size_t n_heuristic_b = 0;   // confidence-fallback sided with anchored
    std::size_t n_flagged = 0;
    bool model_used = false;         // true if a neural model adjudicated
    // Confusion matrix: tag_a × tag_b
    std::array<std::array<std::size_t, 6>, 6> confusion{};
};

class Reconciler {
public:
    void reconcile(
        const std::vector<ReadAssignment>& assignments_a,
        const std::vector<ReadAssignment>& assignments_b,
        std::vector<ReconciledAssignment>& out,
        const ReconcilerOpts& opts = {});

    const ReconcilerStats& stats() const { return stats_; }

private:
    ReconcilerStats stats_;

    // Neural voter (v0.9). When ReconcilerOpts::neural_model_path points at a
    // valid BRANCH_MLP weights file, the model adjudicates each disagreement;
    // otherwise the confidence-fallback below is used and labelled HEURISTIC_*.
    MLP  voter_;
    bool neural_loaded_ = false;

    // Model-based adjudication: features → softmax over {A, B, flag}.
    void neural_resolve(
        const ReadAssignment& a,
        const ReadAssignment& b,
        ReconciledAssignment& out);

    // Confidence fallback: side with the higher-confidence track, FLAG only
    // when both are below threshold. Labelled HEURISTIC_*, never NEURAL_*.
    void heuristic_resolve(
        const ReadAssignment& a,
        const ReadAssignment& b,
        ReconciledAssignment& out);
};

}  // namespace branch::wg::phaser
