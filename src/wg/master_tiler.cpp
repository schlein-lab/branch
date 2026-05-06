// BRANCH v0.5 — Phase 1.3 implementation.
//
// Greedy real-read tiling driven by shared minimizers. The algorithm
// is intentionally lightweight: full sparsest-tiling DP (~2-5 % length
// improvement on real data) is a v0.6 follow-up.
//
// Memory: O(total minimizer count) for the inverted index. With w=10
// (ONT preset) that's ~one minimizer per 10 bp, so a 100k-read input
// at ~10 kbp/read yields ~100 M minimizer entries → ~1.6 GB. Within
// the SLURM 768 GB std-partition envelope by a wide margin.

#include "wg/master_tiler.hpp"
#include "wg/kmer_hist.hpp"  // canonical_kmer_bits

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <unordered_map>
#include <vector>

namespace branch::wg {

namespace {

inline std::uint8_t encode_base(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// Per-read minimizer table: (canonical_bits, position_in_read).
struct ReadMinimizer { std::uint64_t bits; std::uint32_t pos; };

void compute_minimizers(std::string_view seq, std::size_t k, std::size_t w,
                        std::vector<ReadMinimizer>& out) {
    out.clear();
    if (seq.size() < k) return;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);

    // window-of-w sliding minimum of canonical hashes; emit on minimum-change.
    std::vector<std::pair<std::uint64_t, std::uint32_t>> wbuf;  // (bits, pos)
    wbuf.reserve(w);
    std::uint64_t kmer = 0;
    std::size_t filled = 0;
    std::uint64_t prev_min_bits = ~0ULL;

    for (std::size_t i = 0; i < seq.size(); ++i) {
        std::uint8_t b = encode_base(seq[i]);
        if (b > 3) {
            filled = 0; kmer = 0; wbuf.clear(); prev_min_bits = ~0ULL; continue;
        }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;
        std::uint64_t cb = canonical_kmer_bits(kmer, k);
        std::uint32_t pos = static_cast<std::uint32_t>(i - k + 1);
        wbuf.emplace_back(cb, pos);
        if (wbuf.size() > w) wbuf.erase(wbuf.begin());

        // Find min in window.
        if (wbuf.size() == w) {
            auto mn = wbuf[0];
            for (std::size_t j = 1; j < wbuf.size(); ++j) {
                if (wbuf[j].first < mn.first) mn = wbuf[j];
            }
            if (mn.first != prev_min_bits) {
                out.push_back(ReadMinimizer{mn.first, mn.second});
                prev_min_bits = mn.first;
            }
        }
    }
}

}  // namespace

MasterTiler::MasterTiler(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

std::uint64_t MasterTiler::build_tiling(
    const std::vector<std::pair<std::string, std::string>>& reads_by_idx,
    std::vector<MasterTile>& out_tiles,
    const std::string& partial_dump_path) const {
    out_tiles.clear();
    if (reads_by_idx.empty()) return 0;

    const auto t_start = std::chrono::steady_clock::now();
    const std::size_t k = profile_.minimizer_k;
    const std::size_t w = std::max<std::size_t>(profile_.minimizer_w, 5);

    // Open partial-dump TSV (append mode + line-buffered) if requested.
    // Each emitted tile is written as one row and flushed immediately so
    // a walltime-killed run leaves a usable partial output on disk.
    std::FILE* partial = nullptr;
    if (!partial_dump_path.empty()) {
        partial = std::fopen(partial_dump_path.c_str(), "w");
        if (partial) {
            std::setvbuf(partial, nullptr, _IOLBF, 0);
            std::fprintf(partial,
                "tile_idx\tmaster_read_idx\tstart_bp\tlength_bp\t"
                "source_offset_bp\tread_id\n");
        }
    }
    auto close_partial = [&] {
        if (partial) { std::fclose(partial); partial = nullptr; }
    };

    // 1. Compute per-read minimizer lists.
    std::vector<std::vector<ReadMinimizer>> mins(reads_by_idx.size());
    for (std::size_t i = 0; i < reads_by_idx.size(); ++i) {
        compute_minimizers(reads_by_idx[i].second, k, w, mins[i]);
        if ((i + 1) % 100'000 == 0 || i + 1 == reads_by_idx.size()) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            std::fprintf(stderr,
                "[wg/tiler] minimize %zu/%zu reads, elapsed=%.0fs\n",
                i + 1, reads_by_idx.size(), el);
        }
    }

    // 2. Inverted index: minimizer → list of (read_idx, pos).
    struct Hit { std::uint32_t read_idx; std::uint32_t pos; };
    std::unordered_map<std::uint64_t, std::vector<Hit>> idx;
    idx.reserve(reads_by_idx.size() * 16);
    for (std::uint32_t i = 0; i < reads_by_idx.size(); ++i) {
        for (const auto& m : mins[i]) idx[m.bits].push_back(Hit{i, m.pos});
    }

    // 3. Sort reads by length descending; the longest read seeds the first contig.
    std::vector<std::uint32_t> order(reads_by_idx.size());
    for (std::uint32_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) {
        return reads_by_idx[a].second.size() > reads_by_idx[b].second.size();
    });

    std::vector<bool> used(reads_by_idx.size(), false);
    std::uint32_t cursor_bp = 0;
    std::size_t processed = 0;
    auto t_loop_start = std::chrono::steady_clock::now();

    // Simple greedy: walk `order`, for each read either start a new contig
    // or extend the current one if it shares enough minimizers with the
    // previous tile.
    for (std::uint32_t r : order) {
        ++processed;
        if (processed % 100'000 == 0 || processed == order.size()) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_loop_start).count();
            const double frac = static_cast<double>(processed)
                              / static_cast<double>(order.size());
            const double eta = el / std::max(frac, 1e-6) - el;
            std::fprintf(stderr,
                "[wg/tiler] tile %zu/%zu reads (%.1f%%) tiles=%zu cursor_bp=%u "
                "elapsed=%.0fs ETA=%.0fs\n",
                processed, order.size(), frac * 100.0,
                out_tiles.size(), cursor_bp, el, eta);
        }
        if (used[r]) continue;
        const std::string& seq = reads_by_idx[r].second;
        if (seq.empty()) continue;

        // Decide on overlap with the most recent tile (if any).
        std::uint32_t source_offset = 0;
        std::uint32_t length = static_cast<std::uint32_t>(seq.size());

        if (!out_tiles.empty()) {
            const auto& prev = out_tiles.back();
            const auto& prev_seq = reads_by_idx[prev.master_read_idx].second;

            // Count shared minimizers between r and prev's source range.
            std::size_t shared = 0;
            std::uint64_t sum_prev_pos = 0, sum_cur_pos = 0;
            std::size_t sum_n = 0;
            for (const auto& m : mins[r]) {
                auto it = idx.find(m.bits);
                if (it == idx.end()) continue;
                for (const auto& h : it->second) {
                    if (h.read_idx == prev.master_read_idx) {
                        ++shared;
                        // Accumulate position centroids for midpoint cutover.
                        sum_prev_pos += h.pos;
                        sum_cur_pos  += m.pos;
                        ++sum_n;
                    }
                }
            }
            // Need at least 5 shared minimizers AND the average match
            // position in the previous tile is in its second half (i.e.
            // we are extending forward, not backward).
            if (shared >= 5 && sum_n > 0) {
                std::uint32_t avg_prev = static_cast<std::uint32_t>(sum_prev_pos / sum_n);
                std::uint32_t avg_cur  = static_cast<std::uint32_t>(sum_cur_pos / sum_n);
                if (avg_prev >= prev_seq.size() / 2) {
                    // Midpoint cutover: previous tile ends at avg_prev,
                    // current tile starts at avg_cur.
                    std::uint32_t prev_orig_len = prev.length_bp;
                    std::uint32_t prev_orig_off = prev.source_offset_bp;
                    // Recompute prev.length_bp so it ends at avg_prev:
                    if (avg_prev > prev_orig_off) {
                        std::uint32_t new_len = avg_prev - prev_orig_off;
                        if (new_len < prev_orig_len) {
                            // Shrink prev tile to the cutover point.
                            cursor_bp -= (prev_orig_len - new_len);
                            out_tiles.back().length_bp = new_len;
                        }
                    }
                    source_offset = avg_cur;
                    length = static_cast<std::uint32_t>(seq.size()) - source_offset;
                }
            }
        }

        if (length == 0) continue;
        MasterTile t{};
        t.master_read_idx = r;
        t.start_bp = cursor_bp;
        t.source_offset_bp = source_offset;
        t.length_bp = length;
        t.q_mean_x10 = static_cast<std::uint32_t>(
            (profile_.name && profile_.name[0] == 'h') ? 350 : 230);
        out_tiles.push_back(t);
        if (partial) {
            const auto& rid = reads_by_idx[r].first;
            std::fprintf(partial, "%zu\t%u\t%u\t%u\t%u\t%s\n",
                out_tiles.size() - 1,
                t.master_read_idx, t.start_bp, t.length_bp,
                t.source_offset_bp, rid.c_str());
        }
        cursor_bp += length;
        used[r] = true;
    }

    close_partial();
    return cursor_bp;
}

std::string MasterTiler::render_haplotype_seq(
    const std::vector<MasterTile>& tiles,
    const std::vector<std::pair<std::string, std::string>>& reads_by_idx) const {
    std::string out;
    std::size_t total = 0;
    for (const auto& t : tiles) total += t.length_bp;
    out.reserve(total);
    for (const auto& t : tiles) {
        if (t.master_read_idx >= reads_by_idx.size()) continue;
        const auto& seq = reads_by_idx[t.master_read_idx].second;
        if (t.source_offset_bp >= seq.size()) continue;
        std::size_t end = std::min<std::size_t>(seq.size(),
            static_cast<std::size_t>(t.source_offset_bp) + t.length_bp);
        out.append(seq, t.source_offset_bp, end - t.source_offset_bp);
    }
    return out;
}

}  // namespace branch::wg
