#include "reconciler.hpp"

#include <algorithm>
#include <iostream>

namespace branch::wg::phaser {

namespace {

constexpr float kFlagConfidenceThreshold = 0.5f;

// Index a ReadTag into 0..5 for confusion matrix counting.
inline std::size_t tag_idx(ReadTag t) {
    switch (t) {
        case ReadTag::UNTAGGED:              return 0;
        case ReadTag::H1_ONLY:               return 1;
        case ReadTag::H2_ONLY:               return 2;
        case ReadTag::SHARED:                return 3;
        case ReadTag::BRANCH:                return 4;
        case ReadTag::UNCERTAIN:             return 5;
        // v0.10 read-tag classes: lumped to BRANCH idx for confusion mtx.
        case ReadTag::BRANCH_BRIDGE:         return 4;
        case ReadTag::BRANCH_NOVEL:          return 4;
        case ReadTag::BRANCH_EXTINCT_PARENT: return 4;
    }
    return 0;
}

}  // namespace

void Reconciler::heuristic_resolve(
    const ReadAssignment& a,
    const ReadAssignment& b,
    ReconciledAssignment& out)
{
    // No neural model loaded — fall back to confidence-weighted picking.
    // FLAG when both confidences are weak (genuine disagreement we can't
    // adjudicate). This is the v0.8 default until v0.9's neural voter
    // ships.
    const bool a_weak = a.confidence < kFlagConfidenceThreshold;
    const bool b_weak = b.confidence < kFlagConfidenceThreshold;

    if (a_weak && b_weak) {
        out.source = PhaseSource::FLAGGED;
        out.tag_final = ReadTag::UNCERTAIN;
        out.confidence = 0.0f;
        return;
    }
    if (a.confidence > b.confidence) {
        out.source = PhaseSource::NEURAL_A;  // heuristic stand-in
        out.tag_final = a.tag;
        out.confidence = a.confidence;
    } else {
        out.source = PhaseSource::NEURAL_B;
        out.tag_final = b.tag;
        out.confidence = b.confidence;
    }
}

void Reconciler::reconcile(
    const std::vector<ReadAssignment>& assignments_a,
    const std::vector<ReadAssignment>& assignments_b,
    std::vector<ReconciledAssignment>& out,
    const ReconcilerOpts& opts)
{
    stats_ = {};

    if (assignments_a.size() != assignments_b.size()) {
        std::cerr << "[reconcile] WARN: size mismatch: A="
                  << assignments_a.size() << " B=" << assignments_b.size()
                  << " — using min, others marked FLAGGED\n";
    }

    if (!opts.neural_model_path.empty()) {
        // v0.9: load TorchScript / ONNX neural voter.
        // For v0.8 we ship the heuristic fallback only.
        std::cerr << "[reconcile] neural model path provided ("
                  << opts.neural_model_path
                  << ") but neural inference not yet wired in v0.8; "
                     "using heuristic resolver.\n";
    }
    neural_loaded_ = false;

    const std::size_t n = std::min(assignments_a.size(), assignments_b.size());
    out.clear();
    out.reserve(n);
    stats_.n_total = n;

    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = assignments_a[i];
        const auto& b = assignments_b[i];
        ReconciledAssignment r;
        r.read_id     = a.read_id;
        r.tag_a       = a.tag;
        r.tag_b       = b.tag;
        r.home_node_a = a.home_node;
        r.home_node_b = b.home_node;

        ++stats_.confusion[tag_idx(a.tag)][tag_idx(b.tag)];

        if (a.tag == b.tag) {
            r.source     = PhaseSource::AGREED;
            r.tag_final  = a.tag;
            r.confidence = std::max(a.confidence, b.confidence);
            ++stats_.n_agreed;
        } else {
            heuristic_resolve(a, b, r);
            switch (r.source) {
                case PhaseSource::NEURAL_A: ++stats_.n_neural_a; break;
                case PhaseSource::NEURAL_B: ++stats_.n_neural_b; break;
                case PhaseSource::FLAGGED:  ++stats_.n_flagged;  break;
                default: break;
            }
        }

        out.push_back(std::move(r));
    }

    if (opts.verbose) {
        const double pct = stats_.n_total > 0
            ? 100.0 * static_cast<double>(stats_.n_agreed) /
                      static_cast<double>(stats_.n_total)
            : 0.0;
        std::cerr << "[reconcile] " << stats_.n_total << " reads | "
                  << "agreed=" << stats_.n_agreed << " (" << pct << "%) | "
                  << "resolved_to_A=" << stats_.n_neural_a
                  << " resolved_to_B=" << stats_.n_neural_b
                  << " flagged=" << stats_.n_flagged << "\n";
    }
}

}  // namespace branch::wg::phaser
