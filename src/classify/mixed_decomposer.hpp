#pragma once
#include <cstddef>
#include <functional>
#include <vector>

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

// P1.2: Post-decompose re-classification hook.
//
// Signature: called once per sub-bubble that came out of decompose_mixed().
// Consumers (e.g. the hierarchical disambiguator) can plug into this hook
// to re-classify sub-paths without the decomposer itself having to know
// about FeatureVectors or BubbleClass.
//
// `sub_index` is the position in the returned vector, `sub` the sub-bubble.
using SubBubbleReclassifyHook =
    std::function<void(std::size_t sub_index, const SubBubble& sub)>;

/// Same as decompose_mixed() but invokes `hook` (if non-null) once per
/// sub-bubble produced. Mutates nothing; purely a re-classification
/// callout. Returns the same vector decompose_mixed() would return.
std::vector<SubBubble> decompose_mixed_with_hook(
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    const SubBubbleReclassifyHook& hook,
    int max_depth = 3,
    std::size_t min_cluster_size = 3
);

} // namespace branch::classify
