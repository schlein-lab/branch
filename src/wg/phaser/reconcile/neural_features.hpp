#pragma once
//
// neural_features — the canonical feature layout for the dual-phaser neural
// voter. THIS IS A CONTRACT: scripts/train_neural_voter.py must build feature
// vectors with exactly this layout and class ordering, or the trained weights
// will be meaningless at inference time. Keep the two in lockstep.
//
// A disagreement between the de-novo track (A) and the anchored track (B) is
// encoded as a fixed-width float vector; the model emits a 3-way distribution
// over {pick A, pick B, flag}.

#include "../types.hpp"

#include <cstdint>
#include <vector>

namespace branch::wg::phaser {

// Number of ReadTag classes (UNTAGGED..BRANCH_EXTINCT_PARENT). Update together
// with the ReadTag enum in types.hpp.
inline constexpr std::uint32_t kNumReadTags = 9;

// Feature vector: one-hot(tag_a) | one-hot(tag_b) | conf_a | conf_b |
//                 (conf_a - conf_b) | (home_node_a == home_node_b)
inline constexpr std::uint32_t kFeatureDim = 2 * kNumReadTags + 4;

// Output classes (argmax of the model's softmax).
enum class VoterClass : std::uint8_t {
    PickA = 0,  // trust the de-novo track
    PickB = 1,  // trust the anchored track
    Flag  = 2,  // genuine disagreement → FLAGGED / UNCERTAIN
};
inline constexpr std::uint32_t kVoterClasses = 3;

inline std::uint32_t read_tag_index(ReadTag t) {
    const auto v = static_cast<std::uint32_t>(t);
    return v < kNumReadTags ? v : 0u;
}

// Build the feature vector for one disagreeing read pair. The layout MUST
// match scripts/train_neural_voter.py exactly.
inline std::vector<float> make_feature_vector(const ReadAssignment& a,
                                              const ReadAssignment& b) {
    std::vector<float> f(kFeatureDim, 0.0f);
    f[read_tag_index(a.tag)] = 1.0f;                 // one-hot tag_a [0..8]
    f[kNumReadTags + read_tag_index(b.tag)] = 1.0f;  // one-hot tag_b [9..17]
    f[2 * kNumReadTags + 0] = a.confidence;          // [18]
    f[2 * kNumReadTags + 1] = b.confidence;          // [19]
    f[2 * kNumReadTags + 2] = a.confidence - b.confidence;            // [20]
    f[2 * kNumReadTags + 3] = (a.home_node == b.home_node) ? 1.0f : 0.0f;  // [21]
    return f;
}

}  // namespace branch::wg::phaser
