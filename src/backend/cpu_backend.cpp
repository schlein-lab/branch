// BRANCH v0.1 — CPU reference backend implementation.
// v0.2: classify_batch now drives the real cascade (via
// classify::classify_one + feature_extractor), and a candidate-based
// VAF helper is exposed alongside the legacy ID-based vtable entry.

#include "backend/cpu_backend.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <future>
#include <limits>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <vector>

#include "backend/vaf_stats.hpp"
#include "classify/cascade.hpp"
#include "classify/feature_extractor.hpp"
#include "graph/minimizer_sketcher.hpp"

namespace branch::backend {

// Definition of the forward-declared context struct
struct CpuBackendContext {
    // Stateless — worker threads are spawned per compute_overlaps call
    // via std::async and joined before returning, so no per-instance
    // pool handle is required.
};

// Configuration for overlap detection.
struct OverlapConfig {
    std::size_t min_matches{5};         // Minimum shared minimizers
    std::int32_t offset_tolerance{100}; // Max offset deviation in bp
    // v0.2 strand-majority filter: within the offset-consistent
    // cluster, ≥ this fraction of matches must agree on XOR strand.
    // Mixed clusters (spurious cross-strand collisions) are dropped.
    float strand_majority{0.80f};
    // Worker threads used to sketch the read batch in parallel.
    // 0 = single-threaded (legacy behaviour, deterministic ordering of
    // hits). 1..N = std::thread fan-out over read chunks.
    unsigned int threads{0};
};

namespace {

// Hash bucket entry: (read_id, position in read, strand bit).
// The strand bit comes from the sketcher's canonical-hash picker
// (0 = forward k-mer was canonical, 1 = reverse-complement was).
// Two hits from the same strand XOR to 0; opposite strands XOR to 1 —
// that XOR is the pair-level strand emitted on the OverlapPair.
struct BucketEntry {
    std::uint32_t read_id;
    std::uint32_t pos;
    std::uint8_t strand;
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

// Match info for a candidate pair.
//
// Positions are preserved through the clustering step so the emission
// phase can derive `overlap_len` from the actual span of consistent
// minimizer matches rather than a proxy like `match_count * const`.
struct MatchInfo {
    std::int32_t offset;  // posA - posB (A/B are in canonical pair order)
    std::uint32_t pos_a;  // position in the read with the smaller read_id
    std::uint32_t pos_b;  // position in the read with the larger read_id
    std::uint8_t pair_strand;  // XOR of hit-strands: 0 = same, 1 = opposite
};

// Compute overlaps using minimizer bucketing.
std::size_t compute_overlaps_impl(
    const ReadBatch& batch,
    std::span<OverlapPair> out_pairs,
    const OverlapConfig& cfg) {

    if (batch.reads.empty() || out_pairs.empty()) {
        return 0;
    }

    // Step 1: Sketch all reads. With cfg.threads > 1 the batch is
    // partitioned into chunks sketched in parallel; the per-chunk
    // MinimizerHit vectors are concatenated afterwards, preserving the
    // per-read hit order (read_id grows monotonically across chunks).
    std::vector<graph::MinimizerHit> all_hits;
    all_hits.reserve(batch.reads.size() * 100);

    const unsigned int req_threads =
        (cfg.threads == 0) ? 1u : cfg.threads;
    const unsigned int n_threads = std::min<unsigned int>(
        req_threads, static_cast<unsigned int>(batch.reads.size()));

    if (n_threads <= 1) {
        for (const auto& read : batch.reads) {
            graph::sketch_read(read.seq, read.id, all_hits);
        }
    } else {
        const std::size_t n_reads = batch.reads.size();
        const std::size_t chunk =
            (n_reads + n_threads - 1) / n_threads;

        std::vector<std::future<std::vector<graph::MinimizerHit>>> futs;
        futs.reserve(n_threads);

        for (unsigned int t = 0; t < n_threads; ++t) {
            const std::size_t lo = static_cast<std::size_t>(t) * chunk;
            const std::size_t hi = std::min(lo + chunk, n_reads);
            if (lo >= hi) break;
            futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
                std::vector<graph::MinimizerHit> local;
                local.reserve((hi - lo) * 100);
                for (std::size_t i = lo; i < hi; ++i) {
                    const auto& r = batch.reads[i];
                    graph::sketch_read(r.seq, r.id, local);
                }
                return local;
            }));
        }

        for (auto& f : futs) {
            auto chunk_hits = f.get();
            all_hits.insert(all_hits.end(),
                            std::make_move_iterator(chunk_hits.begin()),
                            std::make_move_iterator(chunk_hits.end()));
        }
    }

    if (all_hits.empty()) {
        return 0;
    }

    // Step 2: Build hash buckets. Sized for the number of *unique*
    // hashes rather than total hits to avoid a wasteful over-reserve
    // (all_hits >> unique hashes when reads share many k-mers). An
    // approximate factor of 4 is a conservative lower-bound seen on
    // HiFi data; the map will rehash once or twice in the worst case
    // instead of never, which is a small, deterministic RAM hit.
    std::unordered_map<std::uint64_t, std::vector<BucketEntry>> buckets;
    buckets.reserve(all_hits.size() / 4 + 16);

    for (const auto& hit : all_hits) {
        buckets[hit.hash].push_back(BucketEntry{
            hit.read_id,
            static_cast<std::uint32_t>(hit.pos),
            static_cast<std::uint8_t>(hit.strand),
        });
    }

    // Snapshot bucket keys in sorted order so Step 3+4 iterates the
    // buckets deterministically — without this, two runs on the same
    // input can iterate buckets in different orders (hash seed, load
    // factor) and emit OverlapPairs in shuffled positions, making
    // output files irreproducible.
    std::vector<std::uint64_t> bucket_keys;
    bucket_keys.reserve(buckets.size());
    for (const auto& kv : buckets) bucket_keys.push_back(kv.first);
    std::sort(bucket_keys.begin(), bucket_keys.end());

    // Step 3+4: Collect matches per candidate pair. Reserve to the
    // number of distinct pairs a well-behaved HiFi batch generates
    // (~same scale as buckets.size()) to keep rehash spikes bounded.
    //
    // Parallel path: shard bucket_keys across cfg.threads workers.
    // Each worker builds a thread-local pair_matches map over its
    // bucket shard, then results merge into the global map. Because
    // each bucket key is disjoint (no two workers touch the same
    // bucket), this is lock-free and deterministic: the final output
    // is the union of all shards, identical regardless of the shard
    // boundaries or the merge order (sorting at Step 5 normalises).
    std::unordered_map<PairKey, std::vector<MatchInfo>, PairKeyHash> pair_matches;
    pair_matches.reserve(buckets.size());

    auto process_bucket_range =
        [&](std::size_t lo, std::size_t hi,
            std::unordered_map<PairKey, std::vector<MatchInfo>, PairKeyHash>& out) {
        for (std::size_t bi = lo; bi < hi; ++bi) {
            const auto& entries = buckets[bucket_keys[bi]];
            if (entries.size() < 2) continue;
            for (std::size_t i = 0; i < entries.size(); ++i) {
                for (std::size_t j = i + 1; j < entries.size(); ++j) {
                    const auto& ea = entries[i];
                    const auto& eb = entries[j];
                    if (ea.read_id == eb.read_id) continue;

                    PairKey key;
                    std::int32_t offset;
                    std::uint32_t pos_a;
                    std::uint32_t pos_b;
                    if (ea.read_id < eb.read_id) {
                        key = {ea.read_id, eb.read_id};
                        pos_a = ea.pos;
                        pos_b = eb.pos;
                        offset = static_cast<std::int32_t>(ea.pos) -
                                 static_cast<std::int32_t>(eb.pos);
                    } else {
                        key = {eb.read_id, ea.read_id};
                        pos_a = eb.pos;
                        pos_b = ea.pos;
                        offset = static_cast<std::int32_t>(eb.pos) -
                                 static_cast<std::int32_t>(ea.pos);
                    }
                    const std::uint8_t pair_strand =
                        static_cast<std::uint8_t>((ea.strand ^ eb.strand) & 1);
                    out[key].push_back(
                        MatchInfo{offset, pos_a, pos_b, pair_strand});
                }
            }
        }
    };

    const unsigned int match_threads = (n_threads <= 1) ? 1u : n_threads;
    if (match_threads <= 1) {
        process_bucket_range(0, bucket_keys.size(), pair_matches);
    } else {
        const std::size_t total = bucket_keys.size();
        const std::size_t chunk = (total + match_threads - 1) / match_threads;

        std::vector<std::future<std::unordered_map<PairKey, std::vector<MatchInfo>, PairKeyHash>>> futs;
        futs.reserve(match_threads);

        for (unsigned int t = 0; t < match_threads; ++t) {
            const std::size_t lo = static_cast<std::size_t>(t) * chunk;
            const std::size_t hi = std::min(lo + chunk, total);
            if (lo >= hi) break;
            futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
                std::unordered_map<PairKey, std::vector<MatchInfo>, PairKeyHash> local;
                local.reserve((hi - lo) * 2);
                process_bucket_range(lo, hi, local);
                return local;
            }));
        }

        for (auto& f : futs) {
            auto shard = f.get();
            for (auto& kv : shard) {
                auto& sink = pair_matches[kv.first];
                sink.insert(sink.end(),
                            std::make_move_iterator(kv.second.begin()),
                            std::make_move_iterator(kv.second.end()));
            }
        }
    }

    // Step 5: Filter pairs by match count and offset consistency.
    // Iterate pair_matches in sorted (read_a, read_b) order so the
    // emitted OverlapPairs appear in the same order on every run.
    std::vector<PairKey> pair_keys;
    pair_keys.reserve(pair_matches.size());
    for (const auto& kv : pair_matches) pair_keys.push_back(kv.first);
    std::sort(pair_keys.begin(), pair_keys.end(),
              [](const PairKey& x, const PairKey& y) {
                  if (x.a != y.a) return x.a < y.a;
                  return x.b < y.b;
              });

    std::size_t out_idx = 0;

    for (const auto& pair : pair_keys) {
        const auto& matches = pair_matches[pair];
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

        // Strand-majority filter: inside the offset-consistent cluster,
        // the majority of matches must agree on the XOR-strand bit. A
        // cluster that splits roughly 50/50 between same- and opposite-
        // strand matches is an artefact of hash-collision noise rather
        // than a real read overlap.
        std::size_t cluster_total = 0;
        std::size_t same_strand = 0;
        std::size_t opp_strand = 0;
        for (const auto& m : matches) {
            if (std::abs(m.offset - best_offset) > cfg.offset_tolerance) continue;
            ++cluster_total;
            if (m.pair_strand == 0) ++same_strand; else ++opp_strand;
        }
        if (cluster_total == 0) continue;
        const bool same_wins = same_strand >= opp_strand;
        const std::size_t winner = same_wins ? same_strand : opp_strand;
        const float frac =
            static_cast<float>(winner) / static_cast<float>(cluster_total);
        if (frac < cfg.strand_majority) continue;
        const std::uint8_t cluster_strand = same_wins ? 0 : 1;

        if (out_idx >= out_pairs.size()) break;

        // Convert best_offset (= posA - posB of matching minimizers) to the
        // BRANCH dovetail convention:
        //
        //   offset_x = 0-based start of overlap within read x
        //
        // If best_offset >= 0 the minimizer sits later in A than in B, which
        // means B starts first and A joins the overlap at position
        // best_offset inside itself (B is fully consumed from its 0 bp onward).
        // If best_offset <  0, symmetrically, A starts first and B joins at
        // -best_offset inside itself.
        std::uint32_t offset_a = static_cast<std::uint32_t>(std::max(0, best_offset));
        std::uint32_t offset_b = static_cast<std::uint32_t>(std::max(0, -best_offset));

        // Estimate overlap_len as the actual span of the overlapping region
        // in bp, derived from the positions of the matches that agree with
        // best_offset (within offset_tolerance). We take the max span across
        // the two reads and add k to account for the final minimizer's kmer.
        std::uint32_t min_pos_a = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t max_pos_a = 0;
        std::uint32_t min_pos_b = std::numeric_limits<std::uint32_t>::max();
        std::uint32_t max_pos_b = 0;
        for (const auto& m : matches) {
            if (std::abs(m.offset - best_offset) > cfg.offset_tolerance) continue;
            min_pos_a = std::min(min_pos_a, m.pos_a);
            max_pos_a = std::max(max_pos_a, m.pos_a);
            min_pos_b = std::min(min_pos_b, m.pos_b);
            max_pos_b = std::max(max_pos_b, m.pos_b);
        }
        const std::uint32_t span_a = max_pos_a - min_pos_a;
        const std::uint32_t span_b = max_pos_b - min_pos_b;
        const std::uint32_t overlap_len =
            std::max(span_a, span_b) +
            static_cast<std::uint32_t>(graph::kMinimizerK);

        out_pairs[out_idx] = OverlapPair{
            .read_a = pair.a,
            .read_b = pair.b,
            .offset_a = offset_a,
            .offset_b = offset_b,
            .overlap_len = overlap_len,
            .diff_count = 0,
            .strand = cluster_strand,
            ._pad = 0,
        };
        ++out_idx;
    }

    return out_idx;
}

void cpu_destroy(BackendContext ctx) {
    delete static_cast<CpuBackendContext*>(ctx);
}

// Process-wide tuning knob set by CLI / test harness through
// set_cpu_overlap_threads(). 0 means single-threaded (deterministic).
std::atomic<unsigned int> g_cpu_overlap_threads{0};

void cpu_compute_overlaps(BackendContext /*ctx*/,
                          const ReadBatch* batch,
                          std::span<OverlapPair> out_pairs,
                          std::size_t* out_count) {
    if (!batch || out_pairs.empty()) {
        *out_count = 0;
        return;
    }

    OverlapConfig cfg;
    cfg.threads = g_cpu_overlap_threads.load(std::memory_order_relaxed);
    *out_count = compute_overlaps_impl(*batch, out_pairs, cfg);
}

void cpu_classify_batch(BackendContext /*ctx*/,
                        std::span<const BubbleCandidate> candidates,
                        std::span<ClassificationResult> out_results) {
    // Pre-calibrated thresholds (see docs/CODE_REVIEW_classify.md):
    // DepthRatio 2.0x is biologisch korrekt for diploid->duplication.
    // v0.3 will swap the static rule-based stages for LightGBM.
    static const branch::classify::CascadeConfig cascade_config{};

    const std::size_t n = std::min(candidates.size(), out_results.size());
    for (std::size_t i = 0; i < n; ++i) {
        // Contract: candidates[i].features is expected to be populated
        // upstream via classify::extract_features(candidate, graph).
        // The vtable cannot carry a LosslessGraph reference, so the
        // extraction step lives outside this batch entry point.
        const auto r = branch::classify::classify_one(
            candidates[i].features, cascade_config);
        out_results[i] = ClassificationResult{
            .label = r.label,
            .confidence = r.confidence,
            .stage = r.stage_index,
        };
    }
}

void cpu_estimate_vaf_batch(BackendContext /*ctx*/,
                            std::span<const std::uint32_t> branch_ids,
                            std::span<VAFEstimate> out_estimates) {
    // Legacy ID-based entry: v0.2 has no branch-id -> read-support
    // lookup table yet, so this returns the uninformative prior per
    // branch. The real VAF path goes through the candidate-based
    // overload cpu_estimate_vaf_batch_candidates(...) below, which
    // reads counts directly off the BubbleCandidate.
    const std::size_t n = std::min(branch_ids.size(), out_estimates.size());
    for (std::size_t i = 0; i < n; ++i) {
        out_estimates[i] = wilson_ci(0, 0);  // {0.5, 0.0, 1.0}
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

void set_cpu_overlap_threads(unsigned int threads) noexcept {
    g_cpu_overlap_threads.store(threads, std::memory_order_relaxed);
}

unsigned int get_cpu_overlap_threads() noexcept {
    return g_cpu_overlap_threads.load(std::memory_order_relaxed);
}

void cpu_classify_batch_with_graph(
    std::span<const branch::classify::BubbleCandidate> candidates,
    const branch::graph::LosslessGraph& graph,
    std::span<ClassificationResult> out_results) {
    static const branch::classify::CascadeConfig cascade_config{};
    const std::size_t n = std::min(candidates.size(), out_results.size());
    for (std::size_t i = 0; i < n; ++i) {
        const auto features =
            branch::classify::extract_features(candidates[i], graph);
        const auto r =
            branch::classify::classify_one(features, cascade_config);
        out_results[i] = ClassificationResult{
            .label = r.label,
            .confidence = r.confidence,
            .stage = r.stage_index,
        };
    }
}

void cpu_estimate_vaf_batch_candidates(
    std::span<const branch::classify::BubbleCandidate> candidates,
    std::span<VAFEstimate> out_estimates) {
    const std::size_t n = std::min(candidates.size(), out_estimates.size());
    for (std::size_t i = 0; i < n; ++i) {
        const std::uint32_t branch_reads = candidates[i].read_support_branch;
        const std::uint32_t alt_reads = candidates[i].read_support_alt;
        const std::uint32_t total = branch_reads + alt_reads;
        out_estimates[i] = wilson_ci(branch_reads, total);
    }
}

}  // namespace branch::backend
