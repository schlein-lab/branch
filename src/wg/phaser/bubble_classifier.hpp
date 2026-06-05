#pragma once
//
// bubble_classifier — labels each Bubble as HET / BRANCH / TRIALLELIC / etc
// based purely on the VAF distribution across its alt-paths.
//
// The classifier is *independent* of het-pair signal: it only counts how
// many reads have home_node sitting on each alt-path of each bubble, then
// pattern-matches the resulting fraction vector against canonical
// signatures (see BubbleType in types.hpp).
//
// This produces information whether or not phase_engine successfully
// resolved the bubble — i.e. it works at WG scale even when phasing fails.
// CNV duplications (the [0.5, 0.25, 0.25] pattern) are the primary biology
// target for IGHG; classifier output is what feeds downstream branch-tag
// assignment.

#include "types.hpp"
#include <vector>

namespace branch::wg::phaser {

struct BubbleClassifierOpts {
    // Slack around canonical VAF patterns. 0.12 means [0.5,0.5] matches
    // up to [0.38,0.62] — wide enough to absorb real-world HiFi coverage
    // jitter where bubbles rarely hit exactly 50/50 but cluster around it.
    // v18 measured 2492/2550 bubbles as NOISY with the old 0.07 tol;
    // most of those are real HET/BRANCH bubbles with skewed coverage.
    float  vaf_tol = 0.12f;
    // Bubbles with fewer than this many supporting reads (summed over all
    // alts) are left as UNKNOWN — too noisy to classify reliably.
    // 4 is the floor for distinguishing [0.5,0.5] from [0.75,0.25] noise;
    // anything sparser is statistically meaningless.
    std::uint32_t  min_total_reads = 4;
    // Heartbeat period for stderr logging during classification.
    std::size_t  log_every = 50000;
};

class BubbleClassifier {
public:
    // Legacy path: counts reads-per-alt from home_node, then classifies.
    // Replaced at WG scale by bubble_voter + classify_from_vafs().
    void classify(const std::vector<UnitigNode>& unitigs,
                  const std::vector<ReadAssignment>& assignments,
                  std::vector<Bubble>& bubbles,
                  const BubbleClassifierOpts& opts = {});

    // Run only the pattern-match step. Caller has already populated
    // each Bubble.vafs (e.g. from bubble_voter); this function reads
    // those, stamps Bubble.type and Bubble.branch_depth, and logs a
    // tally. read_counts is preserved as-is.
    void classify_from_vafs(std::vector<Bubble>& bubbles,
                            const BubbleClassifierOpts& opts = {});
};

// Match a sorted-descending VAF vector against canonical signatures.
// Exposed for unit tests. `vafs` must already be sorted high→low.
BubbleType classify_vaf_vector(const std::vector<float>& vafs, float tol);

}  // namespace branch::wg::phaser
