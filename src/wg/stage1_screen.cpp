// BRANCH v0.5 — Phase 3 Stage 1 implementation.
//
// For each requested window size, slide across the haplotype with 50 %
// overlap and compute four signals per window. Standardise against the
// haplotype's median + MAD. OR the flags. Hand the windows back so the
// caller can either persist them as a coverage track or merge them into
// CandidateRegions.
//
// All four signals use only the Phase 0 k-mer counter + the haplotype
// sequence — no read-level walk, no alignment, no annotation. Cheap.
//
// CPU implementation (OpenMP-parallel over windows). The companion GPU
// kernel `wg_stage1_kernel.cu` evaluates the same algorithm on H100;
// the host-side `Stage1Screener` orchestrates whichever backend is
// available (selected at run time by the orchestrator).

#include "wg/stage1_screen.hpp"
#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

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

double dimer_entropy_of(std::string_view w) {
    if (w.size() < 2) return 0.0;
    std::array<std::uint32_t, 16> c{};
    std::uint32_t total = 0;
    for (std::size_t i = 0; i + 1 < w.size(); ++i) {
        std::uint8_t a = encode_base(w[i]);
        std::uint8_t b = encode_base(w[i + 1]);
        if (a > 3 || b > 3) continue;
        ++c[a * 4 + b];
        ++total;
    }
    if (total == 0) return 0.0;
    double H = 0.0;
    for (std::uint32_t v : c) {
        if (v == 0) continue;
        double p = static_cast<double>(v) / total;
        H -= p * std::log2(p);
    }
    return H / 4.0;  // normalise to [0,1]
}

double median_of(std::vector<float> v) {
    if (v.empty()) return 0.0;
    auto m = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + m, v.end());
    return v[m];
}

double mad_of(const std::vector<float>& v, double med) {
    if (v.empty()) return 1.0;
    std::vector<float> dev;
    dev.reserve(v.size());
    for (float x : v) dev.push_back(static_cast<float>(std::abs(x - med)));
    double m = median_of(std::move(dev));
    return m > 0 ? 1.4826 * m : 1.0;  // 1.4826 → MAD-to-σ for normal data
}

}  // namespace

Stage1Screener::Stage1Screener(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void Stage1Screener::screen_haplotype(
    std::uint32_t hap_idx,
    const std::string& hap_seq,
    const std::unordered_map<std::uint64_t, std::uint32_t>& counter,
    std::vector<Stage1Window>& out_windows,
    const std::vector<PinnedLocus>& pinned) const {
    out_windows.clear();
    if (hap_seq.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::uint64_t kmask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);

    // ---- Pre-compute per-base canonical kmer hits + minimizer flags ----
    // base_kmer_count[i] = counter[canonical_bits(seq[i-k+1..i])] for valid i.
    // base_kmer_canon[i] = the canonical bits at position i (for distinct-count).
    // base_is_min[i]     = 1 if the canonical kmer at i is a window-minimum.
    std::vector<std::uint32_t> base_kmer_count(hap_seq.size(), 0);
    std::vector<std::uint64_t> base_kmer_canon(hap_seq.size(), 0);
    std::vector<std::uint8_t>  base_is_min(hap_seq.size(), 0);

    {
        const std::size_t w_min = std::max<std::size_t>(profile_.minimizer_w, 5);
        std::vector<std::pair<std::uint64_t, std::size_t>> wbuf;
        wbuf.reserve(w_min);
        std::uint64_t kmer = 0;
        std::size_t filled = 0;
        std::uint64_t prev_min = ~0ULL;
        for (std::size_t i = 0; i < hap_seq.size(); ++i) {
            std::uint8_t b = encode_base(hap_seq[i]);
            if (b > 3) {
                filled = 0; kmer = 0; wbuf.clear(); prev_min = ~0ULL;
                continue;
            }
            kmer = ((kmer << 2) | b) & kmask;
            ++filled;
            if (filled < k) continue;
            std::uint64_t cb = canonical_kmer_bits(kmer, k);
            base_kmer_canon[i] = cb;
            auto it = counter.find(cb);
            if (it != counter.end()) base_kmer_count[i] = it->second;

            wbuf.emplace_back(cb, i);
            if (wbuf.size() > w_min) wbuf.erase(wbuf.begin());
            if (wbuf.size() == w_min) {
                auto mn = wbuf[0];
                for (std::size_t j = 1; j < wbuf.size(); ++j) {
                    if (wbuf[j].first < mn.first) mn = wbuf[j];
                }
                if (mn.first != prev_min) {
                    base_is_min[mn.second] = 1;
                    prev_min = mn.first;
                }
            }
        }
    }

    // ---- For each window-size, build windows ----
    struct WindowSlot {
        std::uint32_t start_bp;
        std::uint32_t length_bp;
        std::size_t   ws_index;  // index into kStage1WindowSizes (window-size class)
    };

    std::vector<WindowSlot> slots;
    for (std::size_t ws = 0; ws < std::size(kStage1WindowSizes); ++ws) {
        std::uint32_t W = kStage1WindowSizes[ws];
        std::uint32_t step = static_cast<std::uint32_t>(W * (1.0 - kStage1OverlapFraction));
        if (step == 0) step = 1;
        for (std::uint32_t i = 0; i + W <= hap_seq.size(); i += step) {
            slots.push_back(WindowSlot{i, W, ws});
        }
    }

    out_windows.assign(slots.size(), Stage1Window{});

    // ---- Compute per-window signals (parallel) ----
    #pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t s = 0; s < slots.size(); ++s) {
        const auto& sl = slots[s];
        Stage1Window w{};
        w.hap_idx = hap_idx;
        w.start_bp = sl.start_bp;
        w.length_bp = sl.length_bp;

        // mean_kmer_count: average over valid (i ≥ k-1) bases in the window
        std::uint64_t sum = 0;
        std::uint32_t n = 0;
        std::unordered_set<std::uint64_t> distinct;
        std::uint32_t n_min = 0;
        std::size_t lo = sl.start_bp;
        std::size_t hi = sl.start_bp + sl.length_bp;
        for (std::size_t i = lo; i < hi; ++i) {
            if (base_kmer_canon[i] != 0 || base_kmer_count[i] != 0) {
                sum += base_kmer_count[i];
                ++n;
                distinct.insert(base_kmer_canon[i]);
            }
            if (base_is_min[i]) ++n_min;
        }
        w.mean_kmer_count = n > 0
            ? static_cast<float>(static_cast<double>(sum) / static_cast<double>(n))
            : 0.0f;
        w.distinct_kmer_count = static_cast<std::uint32_t>(distinct.size());
        w.dimer_entropy = static_cast<float>(
            dimer_entropy_of(std::string_view(hap_seq.data() + lo, sl.length_bp)));
        w.minimizer_density = sl.length_bp > 0
            ? static_cast<float>(n_min) / static_cast<float>(sl.length_bp)
            : 0.0f;
        out_windows[s] = w;
    }

    // ---- Per-window-size standardisation + flagging ----
    for (std::size_t ws = 0; ws < std::size(kStage1WindowSizes); ++ws) {
        std::vector<float> means, distincts, mindens;
        for (std::size_t s = 0; s < slots.size(); ++s) {
            if (slots[s].ws_index != ws) continue;
            means.push_back(out_windows[s].mean_kmer_count);
            distincts.push_back(static_cast<float>(out_windows[s].distinct_kmer_count));
            mindens.push_back(out_windows[s].minimizer_density);
        }
        if (means.empty()) continue;
        double m_med = median_of(means);
        double m_mad = mad_of(means, m_med);
        double d_med = median_of(distincts);
        double d_mad = mad_of(distincts, d_med);
        double n_med = median_of(mindens);
        double n_mad = mad_of(mindens, n_med);

        for (std::size_t s = 0; s < slots.size(); ++s) {
            if (slots[s].ws_index != ws) continue;
            auto& w = out_windows[s];
            w.z_mean = m_mad > 0
                ? static_cast<float>((w.mean_kmer_count - m_med) / m_mad)
                : 0.0f;
            w.z_distinct = d_mad > 0
                ? static_cast<float>((static_cast<double>(w.distinct_kmer_count) - d_med) / d_mad)
                : 0.0f;
            float z_mindens = n_mad > 0
                ? static_cast<float>((w.minimizer_density - n_med) / n_mad)
                : 0.0f;

            std::uint8_t f = 0;
            if (w.z_mean >=  z_threshold) f |= kStage1Elevated;
            if (w.z_mean <= -z_threshold) f |= kStage1Depleted;
            if (w.dimer_entropy <= rep_entropy) f |= kStage1Repeat;
            if (std::abs(w.z_distinct) >= z_threshold) f |= kStage1ComplexOut;
            if (z_mindens <= -kStage1MinDensityZ) f |= kStage1MinDrop;
            w.flags = f;
        }
    }

    // ---- Pinned paralog loci ----
    for (const auto& p : pinned) {
        for (auto& w : out_windows) {
            if (w.start_bp + w.length_bp <= p.start_bp) continue;
            if (w.start_bp >= p.end_bp) continue;
            w.flags |= kStage1Paralog;
        }
    }
}

void Stage1Screener::merge_into_regions(
    const std::vector<Stage1Window>& windows,
    std::vector<CandidateRegion>& out_regions) const {
    out_regions.clear();
    if (windows.empty()) return;

    // Keep only flagged windows, sort by (hap, start).
    struct W { std::uint32_t hap; std::uint32_t s; std::uint32_t e; std::uint8_t fl; float az; };
    std::vector<W> flagged;
    flagged.reserve(windows.size());
    for (const auto& w : windows) {
        if (w.flags == 0) continue;
        flagged.push_back(W{w.hap_idx, w.start_bp, w.start_bp + w.length_bp,
                            w.flags, std::abs(w.z_mean)});
    }
    if (flagged.empty()) return;
    std::sort(flagged.begin(), flagged.end(), [](const W& a, const W& b) {
        return std::tie(a.hap, a.s) < std::tie(b.hap, b.s);
    });

    CandidateRegion cur{};
    cur.hap_idx = flagged[0].hap;
    cur.start_bp = flagged[0].s;
    cur.end_bp = flagged[0].e;
    cur.flags = flagged[0].fl;
    cur.max_abs_z = flagged[0].az;
    cur.n_windows = 1;
    for (std::size_t i = 1; i < flagged.size(); ++i) {
        const auto& f = flagged[i];
        bool same_hap = f.hap == cur.hap_idx;
        bool gap_ok = static_cast<std::int64_t>(f.s) - static_cast<std::int64_t>(cur.end_bp)
                      <= static_cast<std::int64_t>(kStage1MaxRegionGap);
        if (same_hap && gap_ok) {
            cur.end_bp = std::max(cur.end_bp, f.e);
            cur.flags |= f.fl;
            cur.max_abs_z = std::max(cur.max_abs_z, f.az);
            ++cur.n_windows;
        } else {
            out_regions.push_back(cur);
            cur = CandidateRegion{};
            cur.hap_idx = f.hap;
            cur.start_bp = f.s;
            cur.end_bp = f.e;
            cur.flags = f.fl;
            cur.max_abs_z = f.az;
            cur.n_windows = 1;
        }
    }
    out_regions.push_back(cur);
}

}  // namespace branch::wg
