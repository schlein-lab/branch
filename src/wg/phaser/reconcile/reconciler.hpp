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
    std::size_t n_neural_a = 0;
    std::size_t n_neural_b = 0;
    std::size_t n_flagged = 0;
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

    // Neural model placeholder. v0.9 implementation will load a TorchScript
    // or ONNX model here; v0.8 ships with a heuristic fallback that
    // sides with whichever phaser had higher per-read confidence and
    // FLAGS only when both confidences are below a threshold.
    bool neural_loaded_ = false;
    void heuristic_resolve(
        const ReadAssignment& a,
        const ReadAssignment& b,
        ReconciledAssignment& out);
};

}  // namespace branch::wg::phaser
