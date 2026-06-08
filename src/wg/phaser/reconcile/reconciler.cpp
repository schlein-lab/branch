#include "reconciler.hpp"
#include "neural_features.hpp"

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

void Reconciler::neural_resolve(
    const ReadAssignment& a,
    const ReadAssignment& b,
    ReconciledAssignment& out)
{
    const std::vector<float> probs = voter_.forward(make_feature_vector(a, b));
    if (probs.size() != kVoterClasses) {
        // Defensive: a malformed model output falls back to the heuristic
        // rather than emitting a garbage call.
        heuristic_resolve(a, b, out);
        return;
    }

    // argmax over {PickA, PickB, Flag}.
    std::uint32_t arg = 0;
    for (std::uint32_t i = 1; i < kVoterClasses; ++i) {
        if (probs[i] > probs[arg]) arg = i;
    }
    const float conf = probs[arg];

    switch (static_cast<VoterClass>(arg)) {
        case VoterClass::PickA:
            out.source = PhaseSource::NEURAL_A;
            out.tag_final = a.tag;
            out.confidence = conf;
            break;
        case VoterClass::PickB:
            out.source = PhaseSource::NEURAL_B;
            out.tag_final = b.tag;
            out.confidence = conf;
            break;
        case VoterClass::Flag:
            out.source = PhaseSource::FLAGGED;
            out.tag_final = ReadTag::UNCERTAIN;
            out.confidence = conf;
            break;
    }
}

void Reconciler::heuristic_resolve(
    const ReadAssignment& a,
    const ReadAssignment& b,
    ReconciledAssignment& out)
{
    // No neural model loaded — fall back to confidence-weighted picking.
    // FLAG when both confidences are weak (genuine disagreement we can't
    // adjudicate). Labelled HEURISTIC_* so it is never confused with a real
    // model decision.
    const bool a_weak = a.confidence < kFlagConfidenceThreshold;
    const bool b_weak = b.confidence < kFlagConfidenceThreshold;

    if (a_weak && b_weak) {
        out.source = PhaseSource::FLAGGED;
        out.tag_final = ReadTag::UNCERTAIN;
        out.confidence = 0.0f;
        return;
    }
    if (a.confidence > b.confidence) {
        out.source = PhaseSource::HEURISTIC_A;
        out.tag_final = a.tag;
        out.confidence = a.confidence;
    } else {
        out.source = PhaseSource::HEURISTIC_B;
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

    neural_loaded_ = false;
    if (!opts.neural_model_path.empty()) {
        std::string err;
        if (voter_.load(opts.neural_model_path, &err) &&
            voter_.input_dim() == kFeatureDim &&
            voter_.output_dim() == kVoterClasses) {
            neural_loaded_ = true;
            std::cerr << "[reconcile] neural voter loaded ("
                      << opts.neural_model_path << ", in=" << voter_.input_dim()
                      << " out=" << voter_.output_dim() << ")\n";
        } else {
            std::cerr << "[reconcile] WARN: could not use neural model '"
                      << opts.neural_model_path << "' (";
            if (!err.empty()) std::cerr << err;
            else std::cerr << "expected in=" << kFeatureDim << " out="
                           << kVoterClasses << ", got in=" << voter_.input_dim()
                           << " out=" << voter_.output_dim();
            std::cerr << "); falling back to heuristic resolver.\n";
        }
    }
    stats_.model_used = neural_loaded_;

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
        r.conf_a      = a.confidence;
        r.conf_b      = b.confidence;
        r.home_node_a = a.home_node;
        r.home_node_b = b.home_node;

        ++stats_.confusion[tag_idx(a.tag)][tag_idx(b.tag)];

        if (a.tag == b.tag) {
            r.source     = PhaseSource::AGREED;
            r.tag_final  = a.tag;
            r.confidence = std::max(a.confidence, b.confidence);
            ++stats_.n_agreed;
        } else {
            if (neural_loaded_) neural_resolve(a, b, r);
            else                heuristic_resolve(a, b, r);
            switch (r.source) {
                case PhaseSource::NEURAL_A:    ++stats_.n_neural_a;    break;
                case PhaseSource::NEURAL_B:    ++stats_.n_neural_b;    break;
                case PhaseSource::HEURISTIC_A: ++stats_.n_heuristic_a; break;
                case PhaseSource::HEURISTIC_B: ++stats_.n_heuristic_b; break;
                case PhaseSource::FLAGGED:     ++stats_.n_flagged;     break;
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
                  << (stats_.model_used ? "neural" : "heuristic")
                  << " resolved_to_A="
                  << (stats_.model_used ? stats_.n_neural_a : stats_.n_heuristic_a)
                  << " resolved_to_B="
                  << (stats_.model_used ? stats_.n_neural_b : stats_.n_heuristic_b)
                  << " flagged=" << stats_.n_flagged << "\n";
    }
}

}  // namespace branch::wg::phaser
