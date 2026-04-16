// BRANCH v0.1 — CPU reference backend implementation.

#include "backend/cpu_backend.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "graph/minimizer_sketcher.hpp"

namespace branch::backend {

// Definition of the forward-declared context struct
struct CpuBackendContext {
    // v0.1: no state. Future: thread-pool handle, scratch buffers, etc.
};

// Configuration for overlap detection.
struct OverlapConfig {
    std::size_t min_matches{5};       // Minimum shared minimizers
    std::int32_t offset_tolerance{100}; // Max offset deviation in bp
};

namespace {

// Hash bucket entry: (read_id, position in read)
struct BucketEntry {
    std::uint32_t read_id;
    std::uint32_t pos;
};

// Candidate pair key for deduplication
struct PairKey {
    std::uint32_t a;
    std::uint32_t b;

    bool operator==(const PairKey& o) const { return a == o.a && b == o.b; }
};

struct PairKeyHash {
    std::size_t operator()(const PairKey& k) const {
        return std::hash<std::uint64_t>{}(
            (static_cast<std::uint64_t>(k.a) << 32) | k.b);
    }
};

// Match info for a candidate pair
struct MatchInfo {
    std::int32_t offset;  // posA - posB
};

// Compute overlaps using minimizer bucketing.
std::size_t compute_overlaps_impl(
    const ReadBatch& batch,
    std::span<OverlapPair> out_pairs,
    const OverlapConfig& cfg) {

    if (batch.reads.empty() || out_pairs.empty()) {
        return 0;
    }

    // Step 1: Sketch all reads
    std::vector<graph::MinimizerHit> all_hits;
    all_hits.reserve(batch.reads.size() * 100);

    for (const auto& read : batch.reads) {
        graph::sketch_read(read.seq, read.id, all_hits);
    }

    if (all_hits.empty()) {
        return 0;
    }

    // Step 2: Build hash buckets
    std::unordered_map<std::uint64_t, std::vector<BucketEntry>> buckets;
    buckets.reserve(all_hits.size());

    for (const auto& hit : all_hits) {
        buckets[hit.hash].push_back(BucketEntry{hit.read_id, hit.pos});
    }

    // Step 3+4: Collect matches per candidate pair
    std::unordered_map<PairKey, std::vector<MatchInfo>, PairKeyHash> pair_matches;

    for (const auto& [hash, entries] : buckets) {
        if (entries.size() < 2) continue;

        for (std::size_t i = 0; i < entries.size(); ++i) {
            for (std::size_t j = i + 1; j < entries.size(); ++j) {
                const auto& ea = entries[i];
                const auto& eb = entries[j];

                if (ea.read_id == eb.read_id) continue;

                PairKey key;
                std::int32_t offset;
                if (ea.read_id < eb.read_id) {
                    key = {ea.read_id, eb.read_id};
                    offset = static_cast<std::int32_t>(ea.pos) -
                             static_cast<std::int32_t>(eb.pos);
                } else {
                    key = {eb.read_id, ea.read_id};
                    offset = static_cast<std::int32_t>(eb.pos) -
                             static_cast<std::int32_t>(ea.pos);
                }

                pair_matches[key].push_back(MatchInfo{offset});
            }
        }
    }

    // Step 5: Filter pairs by match count and offset consistency
    std::size_t out_idx = 0;

    for (const auto& [pair, matches] : pair_matches) {
        if (matches.size() < cfg.min_matches) continue;

        std::int32_t best_offset = matches[0].offset;
        std::size_t best_count = 0;

        for (const auto& m : matches) {
            std::size_t count = 0;
            for (const auto& other : matches) {
                if (std::abs(m.offset - other.offset) <= cfg.offset_tolerance) {
                    ++count;
                }
            }
            if (count > best_count) {
                best_count = count;
                best_offset = m.offset;
            }
        }

        if (best_count < cfg.min_matches) continue;

        if (out_idx >= out_pairs.size()) break;

        std::uint32_t offset_a = static_cast<std::uint32_t>(std::max(0, -best_offset));
        std::uint32_t offset_b = static_cast<std::uint32_t>(std::max(0, best_offset));

        out_pairs[out_idx] = OverlapPair{
            .read_a = pair.a,
            .read_b = pair.b,
            .offset_a = offset_a,
            .offset_b = offset_b,
            .overlap_len = static_cast<std::uint32_t>(best_count * 10),
            .diff_count = 0,
            ._pad = 0,
        };
        ++out_idx;
    }

    return out_idx;
}

void cpu_destroy(BackendContext ctx) {
    delete static_cast<CpuBackendContext*>(ctx);
}

void cpu_compute_overlaps(BackendContext /*ctx*/,
                          const ReadBatch* batch,
                          std::span<OverlapPair> out_pairs,
                          std::size_t* out_count) {
    if (!batch || out_pairs.empty()) {
        *out_count = 0;
        return;
    }

    OverlapConfig cfg;
    *out_count = compute_overlaps_impl(*batch, out_pairs, cfg);
}

void cpu_classify_batch(BackendContext /*ctx*/,
                        std::span<const BubbleCandidate> candidates,
                        std::span<ClassificationResult> out_results) {
    const std::size_t n = std::min(candidates.size(), out_results.size());
    for (std::size_t i = 0; i < n; ++i) {
        out_results[i].label = BubbleClass::Unknown;
        out_results[i].confidence = 0.0f;
    }
}

void cpu_estimate_vaf_batch(BackendContext /*ctx*/,
                            std::span<const std::uint32_t> branch_ids,
                            std::span<VAFEstimate> out_estimates) {
    const std::size_t n = std::min(branch_ids.size(), out_estimates.size());
    for (std::size_t i = 0; i < n; ++i) {
        out_estimates[i].point = 0.5f;
        out_estimates[i].ci_low = 0.0f;
        out_estimates[i].ci_high = 1.0f;
    }
}

const char* cpu_name(BackendContext /*ctx*/) {
    return kCpuBackendName;
}

constexpr BackendVTable kCpuVTable{
    .destroy = &cpu_destroy,
    .compute_overlaps = &cpu_compute_overlaps,
    .classify_batch = &cpu_classify_batch,
    .estimate_vaf_batch = &cpu_estimate_vaf_batch,
    .name = &cpu_name,
};

}  // namespace (anonymous)

Backend make_cpu_backend() {
    auto* ctx = new CpuBackendContext{};
    return Backend(ctx, &kCpuVTable);
}

}  // namespace branch::backend
