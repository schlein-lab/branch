#pragma once

// BRANCH v0.1 - Classification features.
//
// Per ML expert consultation: 8-12 handcrafted features feed the
// LightGBM-Cascade at classify time. Features are extracted per
// bubble-candidate from graph + coverage + reference alignment state.
// Order below is stable and used as model input feature index.

#include <cstdint>
#include <array>

namespace branch::classify {

enum class Feature : std::uint8_t {
    FlankJaccardK31 = 0,       // k=31 Jaccard of left+right flanks
    DepthRatioDiploid = 1,     // observed_depth / expected_diploid_depth
    ReadSpanCoverageRatio = 2, // fraction of bubble-spanning reads
    ReadSpanCoverageIQR = 3,   // IQR of per-position span coverage
    RefAlignScoreLeft = 4,     // minimap2 score of left alt vs. ref
    RefAlignScoreRight = 5,    // minimap2 score of right alt vs. ref
    PathLengthDifference = 6,  // |len(alt1) - len(alt2)| / min(len)
    GcContentDivergence = 7,   // |gc(alt1) - gc(alt2)|
    SegdupAnnotationFlag = 8,  // 1 if bubble in UCSC-SegDup, else 0
    RepeatAnnotationFlag = 9,  // 1 if bubble in RepeatMasker, else 0
    AlleleFrequencyEstimate = 10,  // depth-ratio-derived AF
    BubbleLengthBp = 11,       // len(longest alt) in bp
};

inline constexpr std::size_t kNumFeatures = 12;

using FeatureVector = std::array<float, kNumFeatures>;

struct BubbleCandidate {
    std::uint32_t bubble_id;
    std::uint32_t chrom_id;
    std::uint32_t start;
    std::uint32_t end;
    FeatureVector features;

    // Graph linkage to the structural bubble this candidate is derived
    // from. NodeIds (branch::graph::NodeId == std::uint32_t) of the
    // detect::Bubble's entry and exit. Zero means "unknown / not set".
    std::uint32_t entry_node{0};
    std::uint32_t exit_node{0};

    // Support on the two alt paths of the bubble (branch vs. alt).
    // Populated by the feature extractor from detect::Bubble.alts.
    // Used by the VAF estimator: vaf = branch / (branch + alt).
    std::uint32_t read_support_branch{0};
    std::uint32_t read_support_alt{0};

    // Cached length of the longest alt path in bp — passed into
    // Feature::BubbleLengthBp by the extractor, kept here as the raw
    // source value for downstream passes (e.g. VAF / phasing).
    std::uint32_t bubble_length_bp{0};

    // Annotation flags (segdup / repeat) propagated from the reference-
    // track intersector. Zero is the "unannotated" default and means the
    // feature extractor will emit 0.0 for the corresponding Feature slot.
    std::uint8_t segdup_flag{0};
    std::uint8_t repeat_flag{0};
};

static_assert(sizeof(FeatureVector) == kNumFeatures * sizeof(float),
              "FeatureVector must be densely packed");

}  // namespace branch::classify
