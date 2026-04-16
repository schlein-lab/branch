// BRANCH v0.3 — Feature extraction implementation.
//
// Populates 8 of the 12 FeatureVector slots from BubbleCandidate and
// LosslessGraph. Features 3 (ReadSpanCoverageIQR) and 7 (GcContentDivergence)
// are now implemented.

#include "classify/feature_extractor.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace branch::classify {

namespace {

constexpr std::size_t kFlankK = 31;
constexpr std::size_t kFlankLength = 5 * kFlankK;  // 155bp

std::unordered_set<std::string_view> extract_kmers(
    std::string_view seq, std::size_t k) {
    std::unordered_set<std::string_view> kmers;
    if (seq.size() < k) return kmers;
    kmers.reserve(seq.size() - k + 1);
    for (std::size_t i = 0; i + k <= seq.size(); ++i) {
        kmers.insert(seq.substr(i, k));
    }
    return kmers;
}

float jaccard_similarity(
    const std::unordered_set<std::string_view>& a,
    const std::unordered_set<std::string_view>& b) {
    if (a.empty() && b.empty()) return 0.0f;
    std::size_t intersection = 0;
    for (const auto& kmer : a) {
        if (b.count(kmer) > 0) ++intersection;
    }
    const std::size_t union_size = a.size() + b.size() - intersection;
    if (union_size == 0) return 0.0f;
    return static_cast<float>(intersection) / static_cast<float>(union_size);
}

float compute_flank_jaccard_k31(
    const branch::graph::LosslessGraph& graph,
    std::uint32_t entry_node,
    std::uint32_t exit_node) {
    if (entry_node >= graph.node_count() || exit_node >= graph.node_count()) {
        return 0.0f;
    }
    const auto& entry_consensus = graph.node(entry_node).consensus;
    const auto& exit_consensus = graph.node(exit_node).consensus;
    if (entry_consensus.empty() || exit_consensus.empty()) return 0.0f;

    std::string_view left_flank;
    if (entry_consensus.size() >= kFlankLength) {
        left_flank = std::string_view(entry_consensus).substr(
            entry_consensus.size() - kFlankLength, kFlankLength);
    } else {
        left_flank = entry_consensus;
    }

    std::string_view right_flank;
    if (exit_consensus.size() >= kFlankLength) {
        right_flank = std::string_view(exit_consensus).substr(0, kFlankLength);
    } else {
        right_flank = exit_consensus;
    }

    auto left_kmers = extract_kmers(left_flank, kFlankK);
    auto right_kmers = extract_kmers(right_flank, kFlankK);
    return jaccard_similarity(left_kmers, right_kmers);
}

float diploid_baseline(const branch::graph::LosslessGraph& graph) {
    const auto& edges = graph.edges();
    if (edges.empty()) return 1.0f;
    std::uint64_t sum = 0;
    for (const auto& e : edges) sum += e.read_support;
    const float mean = static_cast<float>(sum) / static_cast<float>(edges.size());
    return std::max(mean, 1.0f);
}

// Compute GC content: (G+C) / (A+C+G+T). Returns 0.0f if empty.
float compute_gc_content(std::string_view seq) {
    if (seq.empty()) return 0.0f;
    std::size_t gc = 0, total = 0;
    for (char c : seq) {
        switch (c) {
            case 'G': case 'g': case 'C': case 'c': ++gc; ++total; break;
            case 'A': case 'a': case 'T': case 't': ++total; break;
            default: break;
        }
    }
    if (total == 0) return 0.0f;
    return static_cast<float>(gc) / static_cast<float>(total);
}

// Compute IQR (Q3 - Q1) from read spans. Returns 0.0f if <4 elements.
float compute_read_span_iqr(std::vector<std::uint32_t> spans) {
    if (spans.size() < 4) return 0.0f;
    std::sort(spans.begin(), spans.end());
    const std::size_t n = spans.size();
    const float q1_idx = 0.25f * static_cast<float>(n - 1);
    const float q3_idx = 0.75f * static_cast<float>(n - 1);
    auto percentile = [&](float idx) -> float {
        const std::size_t lo = static_cast<std::size_t>(idx);
        const std::size_t hi = std::min(lo + 1, n - 1);
        const float frac = idx - static_cast<float>(lo);
        return static_cast<float>(spans[lo]) * (1.0f - frac) +
               static_cast<float>(spans[hi]) * frac;
    };
    return percentile(q3_idx) - percentile(q1_idx);
}

}  // namespace

FeatureVector extract_features(
    const BubbleCandidate& candidate,
    const branch::graph::LosslessGraph& graph) {

    FeatureVector f{};

    // FlankJaccardK31 (index 0)
    f[static_cast<std::size_t>(Feature::FlankJaccardK31)] =
        compute_flank_jaccard_k31(graph, candidate.entry_node, candidate.exit_node);

    // DepthRatioDiploid (index 1)
    const std::uint32_t total_support =
        candidate.read_support_branch + candidate.read_support_alt;
    const float baseline = diploid_baseline(graph);
    f[static_cast<std::size_t>(Feature::DepthRatioDiploid)] =
        static_cast<float>(total_support) / baseline;

    // ReadSpanCoverageRatio (index 2)
    float span_ratio = 0.0f;
    if (candidate.entry_node < graph.node_count()) {
        const std::uint32_t entry_reads =
            graph.node(candidate.entry_node).read_support;
        if (entry_reads > 0) {
            span_ratio = static_cast<float>(total_support) /
                         static_cast<float>(entry_reads);
            if (span_ratio > 1.0f) span_ratio = 1.0f;
        }
    }
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageRatio)] = span_ratio;

    // ReadSpanCoverageIQR (index 3) — NEW
    f[static_cast<std::size_t>(Feature::ReadSpanCoverageIQR)] =
        compute_read_span_iqr(candidate.read_spans);

    // GcContentDivergence (index 7) — NEW
    float gc_div = 0.0f;
    if (candidate.entry_node < graph.node_count() &&
        candidate.exit_node < graph.node_count()) {
        const auto& entry_cons = graph.node(candidate.entry_node).consensus;
        const auto& exit_cons = graph.node(candidate.exit_node).consensus;
        float gc_entry = compute_gc_content(entry_cons);
        float gc_exit = compute_gc_content(exit_cons);
        gc_div = std::abs(gc_entry - gc_exit);
    }
    f[static_cast<std::size_t>(Feature::GcContentDivergence)] = gc_div;

    // Annotation flags (indices 8, 9)
    f[static_cast<std::size_t>(Feature::SegdupAnnotationFlag)] =
        candidate.segdup_flag ? 1.0f : 0.0f;
    f[static_cast<std::size_t>(Feature::RepeatAnnotationFlag)] =
        candidate.repeat_flag ? 1.0f : 0.0f;

    // BubbleLengthBp (index 11)
    f[static_cast<std::size_t>(Feature::BubbleLengthBp)] =
        static_cast<float>(candidate.bubble_length_bp);

    return f;
}

}  // namespace branch::classify
