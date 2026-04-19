#include "classify/mixed_decomposer.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace branch::classify {

namespace {

// Hamming distance between two SNP vectors
int hamming_distance(const std::vector<char>& a, const std::vector<char>& b) {
    if (a.size() != b.size()) return std::numeric_limits<int>::max();
    int dist = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) ++dist;
    }
    return dist;
}

// Single-linkage clustering: find closest pair of clusters and merge
struct Cluster {
    std::vector<size_t> members;  // indices into original read set
};

// Find minimum distance between any two reads in different clusters
std::pair<size_t, size_t> find_closest_clusters(
    const std::vector<Cluster>& clusters,
    const std::vector<std::vector<char>>& snp_vecs,
    int& min_dist
) {
    min_dist = std::numeric_limits<int>::max();
    size_t ci = 0, cj = 0;
    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size(); ++j) {
            for (size_t ri : clusters[i].members) {
                for (size_t rj : clusters[j].members) {
                    int d = hamming_distance(snp_vecs[ri], snp_vecs[rj]);
                    if (d < min_dist) {
                        min_dist = d;
                        ci = i;
                        cj = j;
                    }
                }
            }
        }
    }
    return {ci, cj};
}

} // anonymous namespace

std::vector<SubBubble> decompose_mixed(
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    [[maybe_unused]] int max_depth,
    size_t min_cluster_size
) {
    if (per_read_snp_vectors.size() < 2 * min_cluster_size) {
        return {};  // Not enough reads for 2 clusters
    }

    // Initialize: each read is its own cluster
    std::vector<Cluster> clusters;
    for (size_t i = 0; i < per_read_snp_vectors.size(); ++i) {
        clusters.push_back({{i}});
    }

    // Agglomerative clustering until we have 2 clusters
    while (clusters.size() > 2) {
        int min_dist;
        auto [ci, cj] = find_closest_clusters(clusters, per_read_snp_vectors, min_dist);
        
        // Merge cj into ci
        for (size_t m : clusters[cj].members) {
            clusters[ci].members.push_back(m);
        }
        clusters.erase(clusters.begin() + cj);
    }

    // Check if both clusters meet minimum size
    if (clusters.size() < 2) return {};
    if (clusters[0].members.size() < min_cluster_size ||
        clusters[1].members.size() < min_cluster_size) {
        return {};
    }

    // Build SubBubbles
    std::vector<SubBubble> result;
    size_t total_reads = per_read_snp_vectors.size();
    for (const auto& c : clusters) {
        SubBubble sb;
        sb.read_ids = c.members;
        sb.estimated_vaf = static_cast<double>(c.members.size()) / static_cast<double>(total_reads);
        result.push_back(std::move(sb));
    }

    return result;
}

std::vector<SubBubble> decompose_mixed_with_hook(
    const std::vector<std::vector<char>>& per_read_snp_vectors,
    const SubBubbleReclassifyHook& hook,
    int max_depth,
    std::size_t min_cluster_size
) {
    auto subs = decompose_mixed(per_read_snp_vectors, max_depth, min_cluster_size);
    if (hook) {
        for (std::size_t i = 0; i < subs.size(); ++i) {
            hook(i, subs[i]);
        }
    }
    return subs;
}

} // namespace branch::classify
