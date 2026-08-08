// BRANCH — CSR overlap engine implementation. See header for design
// rationale. The per-pair clustering in emit_pairs_for_query() is a
// line-for-line port of compute_overlaps_impl() steps 5a-5d
// (cpu_backend.cpp) so both engines emit identical OverlapPairs.

#include "backend/csr_overlap_engine.hpp"

#include <algorithm>
#include <cstdlib>
#include <future>
#include <limits>

#include "graph/minimizer_sketcher.hpp"

namespace branch::backend {

namespace {

// Match candidate collected for one query read. `target` is always the
// smaller read id (the query is the larger one by construction).
//
// `hash` and `entry_idx` reproduce the legacy engine's match order for
// bit-identical tie-breaking in best_offset selection: legacy iterates
// buckets in ascending hash order and enumerates bucket entries in
// insertion (= read) order. Sorting CandHits by (target, hash,
// entry_idx) with a stable sort — collection order is query-hit order —
// yields exactly that sequence per candidate pair.
struct CandHit {
    std::uint64_t hash;
    std::uint32_t target;
    std::uint32_t entry_idx;  // rank of the target occurrence in its bucket
    std::int32_t offset;      // pos_target - pos_query  (== pos_a - pos_b)
    std::uint32_t pos_a;      // position in target (smaller id)
    std::uint32_t pos_b;      // position in query (larger id)
    std::uint8_t pair_strand;
};

// Sketch a chunk of reads. Shared helper for both index build passes
// and the query pass.
std::vector<graph::MinimizerHit> sketch_chunk(
    const ReadBatch& batch, std::size_t lo, std::size_t hi,
    const graph::MinimizerProfile& profile) {
    std::vector<graph::MinimizerHit> local;
    local.reserve((hi - lo) * 100);
    for (std::size_t i = lo; i < hi; ++i) {
        const auto& r = batch.reads[i];
        graph::sketch_read(r.seq, r.id, local, profile);
    }
    return local;
}

}  // namespace

std::size_t CsrMinimizerIndex::find_slot(std::uint64_t key) const noexcept {
    std::size_t i = static_cast<std::size_t>(key) & mask_;
    while (true) {
        const Slot& s = slots_[i];
        if (s.key == key) return i;
        if (s.key == kEmptyKey) return slots_.size();  // not found
        i = (i + 1) & mask_;
    }
}

std::size_t CsrMinimizerIndex::find_or_insert_slot(std::uint64_t key) noexcept {
    std::size_t i = static_cast<std::size_t>(key) & mask_;
    while (true) {
        Slot& s = slots_[i];
        if (s.key == key) return i;
        if (s.key == kEmptyKey) {
            s.key = key;
            s.count = 0;
            s.offset = 0;
            ++n_keys_;
            return i;
        }
        i = (i + 1) & mask_;
    }
}

void CsrMinimizerIndex::build(const ReadBatch& batch,
                              const graph::MinimizerProfile& profile,
                              std::size_t seed_cap,
                              unsigned int threads,
                              CsrIndexStats* stats) {
    slots_.clear();
    entries_.clear();
    n_keys_ = 0;

    if (batch.reads.empty()) {
        mask_ = 0;
        if (stats) *stats = CsrIndexStats{};
        return;
    }

    const unsigned int n_threads =
        std::max(1u, std::min<unsigned int>(
                         threads, static_cast<unsigned int>(batch.reads.size())));
    const std::size_t n_reads = batch.reads.size();
    const std::size_t chunk = (n_reads + n_threads - 1) / n_threads;

    // ---- Pass A: sketch (parallel) + count (sequential consumer) ----
    // Total hit count is needed to size the table before inserting, so
    // chunks are sketched first and kept; peak extra memory is one
    // MinimizerHit (16 B) per hit for the duration of the build. This
    // is the same transient the legacy engine keeps permanently.
    std::vector<std::future<std::vector<graph::MinimizerHit>>> futs;
    futs.reserve(n_threads);
    for (unsigned int t = 0; t < n_threads; ++t) {
        const std::size_t lo = static_cast<std::size_t>(t) * chunk;
        const std::size_t hi = std::min(lo + chunk, n_reads);
        if (lo >= hi) break;
        futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
            return sketch_chunk(batch, lo, hi, profile);
        }));
    }
    std::vector<std::vector<graph::MinimizerHit>> chunks;
    chunks.reserve(futs.size());
    std::uint64_t total_hits = 0;
    for (auto& f : futs) {
        chunks.push_back(f.get());
        total_hits += chunks.back().size();
    }

    // Power-of-two table, sized so even the worst case (every hit a
    // distinct hash) stays at <= 50% load factor. Real HiFi batches
    // have ~10-30% distinct, so the effective load factor is lower.
    std::size_t cap = 16;
    while (cap < total_hits * 2) cap <<= 1;
    mask_ = cap - 1;
    slots_.assign(cap, Slot{kEmptyKey, 0, 0});

    for (const auto& ch : chunks) {
        for (const auto& hit : ch) {
            Slot& s = slots_[find_or_insert_slot(hit.hash)];
            ++s.count;
        }
    }

    // ---- Cap + prefix sum (slot order = deterministic) ----
    std::uint64_t n_capped = 0, dropped = 0, kept = 0;
    std::vector<CsrIndexStats::CappedHash> capped_report;
    for (auto& s : slots_) {
        if (s.key == kEmptyKey) continue;
        if (seed_cap > 0 && s.count > seed_cap) {
            ++n_capped;
            dropped += s.count;
            capped_report.push_back({s.key, s.count});
            s.offset = kCapped;
        }
    }
    std::sort(capped_report.begin(), capped_report.end(),
              [](const auto& x, const auto& y) {
                  if (x.count != y.count) return x.count > y.count;
                  return x.hash < y.hash;
              });
    if (capped_report.size() > CsrIndexStats::kMaxReportedCapped) {
        capped_report.resize(CsrIndexStats::kMaxReportedCapped);
    }

    std::uint64_t run = 0;
    for (auto& s : slots_) {
        if (s.key == kEmptyKey || s.offset == kCapped) continue;
        s.offset = static_cast<std::uint32_t>(run);
        run += s.count;
        kept += s.count;
    }
    entries_.resize(run);

    // ---- Pass B: scatter (sequential, reuses sketched chunks) ----
    std::vector<std::uint32_t> fill(slots_.size(), 0);
    for (const auto& ch : chunks) {
        for (const auto& hit : ch) {
            const std::size_t si = find_slot(hit.hash);
            Slot& s = slots_[si];
            if (s.offset == kCapped) continue;
            entries_[s.offset + fill[si]++] = CsrEntry::make(
                hit.read_id, static_cast<std::uint32_t>(hit.pos),
                static_cast<std::uint8_t>(hit.strand));
        }
    }

    if (stats) {
        stats->n_reads = n_reads;
        stats->n_hits_total = total_hits;
        stats->n_distinct_hashes = n_keys_;
        stats->n_entries_indexed = kept;
        stats->n_hashes_capped = n_capped;
        stats->n_entries_dropped = dropped;
        stats->capped_hashes = std::move(capped_report);
    }
}

std::span<const CsrEntry> CsrMinimizerIndex::lookup(
    std::uint64_t hash) const noexcept {
    if (slots_.empty()) return {};
    const std::size_t si = find_slot(hash);
    if (si == slots_.size()) return {};
    const Slot& s = slots_[si];
    if (s.offset == kCapped) return {};
    return {entries_.data() + s.offset, s.count};
}

namespace {

// Port of compute_overlaps_impl steps 5a-5d for one candidate pair.
// `matches` spans the CandHits of a single (query, target) pair.
bool cluster_and_emit(std::span<const CandHit> matches,
                      std::uint32_t target, std::uint32_t query,
                      const CsrOverlapConfig& cfg, OverlapPair& out) {
    if (matches.size() < cfg.min_matches) return false;

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
    if (best_count < cfg.min_matches) return false;

    std::size_t cluster_total = 0, same_strand = 0, opp_strand = 0;
    for (const auto& m : matches) {
        if (std::abs(m.offset - best_offset) > cfg.offset_tolerance) continue;
        ++cluster_total;
        if (m.pair_strand == 0) ++same_strand; else ++opp_strand;
    }
    if (cluster_total == 0) return false;
    const bool same_wins = same_strand >= opp_strand;
    const std::size_t winner = same_wins ? same_strand : opp_strand;
    const float frac =
        static_cast<float>(winner) / static_cast<float>(cluster_total);
    if (frac < cfg.strand_majority) return false;
    const std::uint8_t cluster_strand = same_wins ? 0 : 1;

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
        static_cast<std::uint32_t>(graph::current_profile().k);

    out = OverlapPair{
        .read_a = target,
        .read_b = query,
        .offset_a = static_cast<std::uint32_t>(std::max(0, best_offset)),
        .offset_b = static_cast<std::uint32_t>(std::max(0, -best_offset)),
        .overlap_len = overlap_len,
        .diff_count = 0,
        .strand = cluster_strand,
        ._pad = 0,
    };
    return true;
}

// Process queries [lo, hi): sketch each, gather candidates against the
// index, cluster per target, emit pairs into `out`.
void process_query_range(const ReadBatch& batch, std::size_t lo, std::size_t hi,
                         const CsrMinimizerIndex& index,
                         const graph::MinimizerProfile& profile,
                         const CsrOverlapConfig& cfg,
                         std::vector<OverlapPair>& out) {
    std::vector<graph::MinimizerHit> hits;
    std::vector<CandHit> scratch;
    for (std::size_t qi = lo; qi < hi; ++qi) {
        const auto& q = batch.reads[qi];
        hits.clear();
        scratch.clear();
        graph::sketch_read(q.seq, q.id, hits, profile);

        for (const auto& h : hits) {
            const auto entries = index.lookup(h.hash);
            for (std::size_t ei = 0; ei < entries.size(); ++ei) {
                const CsrEntry& e = entries[ei];
                const std::uint32_t t = e.read_id();
                if (t >= q.id) continue;  // pair handled when other read queries
                scratch.push_back(CandHit{
                    .hash = h.hash,
                    .target = t,
                    .entry_idx = static_cast<std::uint32_t>(ei),
                    .offset = static_cast<std::int32_t>(e.pos()) -
                              static_cast<std::int32_t>(h.pos),
                    .pos_a = e.pos(),
                    .pos_b = static_cast<std::uint32_t>(h.pos),
                    .pair_strand = static_cast<std::uint8_t>(
                        (e.strand() ^ h.strand) & 1),
                });
            }
        }
        if (scratch.empty()) continue;

        // Stable sort keeps query-hit collection order as the final
        // tie-breaker, completing the legacy-identical match order
        // (see CandHit comment).
        std::stable_sort(scratch.begin(), scratch.end(),
                         [](const CandHit& x, const CandHit& y) {
                             if (x.target != y.target) return x.target < y.target;
                             if (x.hash != y.hash) return x.hash < y.hash;
                             return x.entry_idx < y.entry_idx;
                         });

        std::size_t g0 = 0;
        for (std::size_t i = 1; i <= scratch.size(); ++i) {
            if (i == scratch.size() || scratch[i].target != scratch[g0].target) {
                OverlapPair pair;
                if (cluster_and_emit({scratch.data() + g0, i - g0},
                                     scratch[g0].target, q.id, cfg, pair)) {
                    out.push_back(pair);
                }
                g0 = i;
            }
        }
    }
}

}  // namespace

std::size_t compute_overlaps_csr(const ReadBatch& batch,
                                 const CsrOverlapConfig& cfg,
                                 std::vector<OverlapPair>& out,
                                 CsrIndexStats* stats) {
    if (batch.reads.empty()) {
        if (stats) *stats = CsrIndexStats{};
        return 0;
    }

    const auto& profile = graph::current_profile();
    const unsigned int n_threads = std::max(1u, cfg.threads);

    CsrMinimizerIndex index;
    index.build(batch, profile, cfg.seed_cap, n_threads, stats);

    const std::size_t n_reads = batch.reads.size();
    const std::size_t before = out.size();

    if (n_threads <= 1 || n_reads < 2 * n_threads) {
        process_query_range(batch, 0, n_reads, index, profile, cfg, out);
    } else {
        const std::size_t chunk = (n_reads + n_threads - 1) / n_threads;
        std::vector<std::future<std::vector<OverlapPair>>> futs;
        futs.reserve(n_threads);
        for (unsigned int t = 0; t < n_threads; ++t) {
            const std::size_t lo = static_cast<std::size_t>(t) * chunk;
            const std::size_t hi = std::min(lo + chunk, n_reads);
            if (lo >= hi) break;
            futs.push_back(std::async(std::launch::async, [&, lo, hi]() {
                std::vector<OverlapPair> local;
                process_query_range(batch, lo, hi, index, profile, cfg, local);
                return local;
            }));
        }
        for (auto& f : futs) {
            auto local = f.get();
            out.insert(out.end(), local.begin(), local.end());
        }
    }

    std::sort(out.begin() + static_cast<std::ptrdiff_t>(before), out.end(),
              [](const OverlapPair& x, const OverlapPair& y) {
                  if (x.read_a != y.read_a) return x.read_a < y.read_a;
                  if (x.read_b != y.read_b) return x.read_b < y.read_b;
                  return x.offset_a < y.offset_a;
              });
    return out.size() - before;
}

}  // namespace branch::backend
