#pragma once
#include <vector>
#include <cstddef>

namespace branch::classify {

struct SubBubble {
    std::vector<size_t> read_ids;
    double estimated_vaf;
};

// Forward declaration - BubbleCandidate is defined elsewhere
struct BubbleCandidate;

/// Decompose a Mixed bubble into sub-bubbles via SNP-vector clustering.
/// Returns empty if no valid decomposition found (< 2 clusters or min_cluster_size not met).
std::vector<SubBubble> decompose_mixed(
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    int max_depth = 3,
    size_t min_cluster_size = 3
);

} // namespace branch::classify
