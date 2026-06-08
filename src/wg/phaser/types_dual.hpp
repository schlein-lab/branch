#pragma once
//
// Types for v0.8/v0.9 dual-phaser architecture.
//
// Two independent phasers (de-novo + pangenome-anchored) produce parallel
// assignments. A reconciliation pass merges them, exposing disagreements
// as either neural-resolved confident calls or biologically-meaningful
// UNCERTAIN flags.

#include "types.hpp"
#include <string>
#include <vector>

namespace branch::wg::phaser {

// Which model produced this read assignment.
enum class PhaseSource : std::uint8_t {
    DENOVO      = 0,   // overlap-graph + bubble + het-pair vote
    ANCHORED    = 1,   // pangenome-anchored phaser
    AGREED      = 2,   // both models agreed
    NEURAL_A    = 3,   // disagreement; neural voter sided with DENOVO
    NEURAL_B    = 4,   // disagreement; neural voter sided with ANCHORED
    FLAGGED     = 5,   // disagreement; voter unsure → biologically-meaningful uncertain
    HEURISTIC_A = 6,   // disagreement; confidence-fallback sided with DENOVO (no model)
    HEURISTIC_B = 7,   // disagreement; confidence-fallback sided with ANCHORED (no model)
};

inline const char* phase_source_str(PhaseSource s) {
    switch (s) {
        case PhaseSource::DENOVO:      return "denovo";
        case PhaseSource::ANCHORED:    return "anchored";
        case PhaseSource::AGREED:      return "agreed";
        case PhaseSource::NEURAL_A:    return "neural_denovo";
        case PhaseSource::NEURAL_B:    return "neural_anchored";
        case PhaseSource::FLAGGED:     return "flagged";
        case PhaseSource::HEURISTIC_A: return "heuristic_denovo";
        case PhaseSource::HEURISTIC_B: return "heuristic_anchored";
    }
    return "?";
}

// Per-read result after reconciliation. Carries provenance so downstream
// callers (analysis, visualization, paper figures) can distinguish
// confident calls from neural-resolved from genuinely-uncertain regions.
struct ReconciledAssignment {
    ReadId   read_id;
    ReadTag  tag_a;      // de-novo phaser call
    ReadTag  tag_b;      // anchored phaser call
    ReadTag  tag_final;  // consensus / neural / flagged
    PhaseSource source;
    float    confidence;     // [0,1]; 1.0 = both agree; <0.5 = uncertain
    float    conf_a = 0.0f;  // de-novo track per-read confidence
    float    conf_b = 0.0f;  // anchored track per-read confidence
    NodeId   home_node_a;
    NodeId   home_node_b;
};

}  // namespace branch::wg::phaser
