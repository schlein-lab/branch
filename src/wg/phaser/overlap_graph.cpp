#include "overlap_graph.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace branch::wg::phaser {

namespace {

inline std::uint64_t hash64(std::uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33; return x;
}

inline std::uint8_t base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;  // ambiguous; breaks the k-mer
    }
}

// Extract canonical (k,w)-minimizer set from a sequence.
// Returns a vector of (hash, position) pairs, deduplicated by hash.
void extract_minimizers(const std::string& seq, int k, int w,
                        std::vector<std::uint64_t>& out_hashes) {
    out_hashes.clear();
    if (static_cast<int>(seq.size()) < k) return;
    const int len = static_cast<int>(seq.size());
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t fwd = 0, rev = 0;
    int valid = 0;
    std::vector<std::pair<std::uint64_t, int>> window;
    window.reserve(static_cast<std::size_t>(w + 1));

    auto push_window = [&](std::uint64_t h, int pos) {
        while (!window.empty() && window.back().first > h) window.pop_back();
        window.emplace_back(h, pos);
        while (!window.empty() && window.front().second <= pos - w)
            window.erase(window.begin());
    };

    // Tight reserve: minimizers occur at most every w bases AFTER the
    // dedup step. Over-reserving here cost the WG run ~50 GB (4.3M reads
    // × 13 KB/sketch), so we reserve only the post-dedup budget and
    // shrink_to_fit at the end. ~600 mins per 18kb HiFi read at w=11.
    out_hashes.reserve(static_cast<std::size_t>(len / std::max(w, 1) / 2) + 16);
    std::uint64_t last_emitted = ~0ULL;
    for (int i = 0; i < len; ++i) {
        const std::uint8_t b = base_code(seq[i]);
        if (b == 4) { fwd = 0; rev = 0; valid = 0; window.clear(); continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << (2 * (k - 1)));
        ++valid;
        if (valid >= k) {
            const std::uint64_t canon = std::min(fwd, rev);
            const std::uint64_t h = hash64(canon);
            push_window(h, i);
            if (valid >= k + w - 1) {
                const std::uint64_t cur = window.front().first;
                if (cur != last_emitted) {
                    out_hashes.push_back(cur);
                    last_emitted = cur;
                }
            }
        }
    }
    std::sort(out_hashes.begin(), out_hashes.end());
    out_hashes.erase(std::unique(out_hashes.begin(), out_hashes.end()),
                     out_hashes.end());
    // Release any over-reserved capacity. Critical at WG scale: 4.3M
    // reads × ~10 KB excess capacity = ~40 GB wasted RSS otherwise.
    out_hashes.shrink_to_fit();
}

// Jaccard on two sorted minimizer hash arrays.
float sorted_jaccard(const std::vector<std::uint64_t>& a,
                     const std::vector<std::uint64_t>& b) {
    if (a.empty() || b.empty()) return 0.0f;
    std::size_t i = 0, j = 0, inter = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) { ++inter; ++i; ++j; }
        else if (a[i] < b[j]) ++i; else ++j;
    }
    const std::size_t uni = a.size() + b.size() - inter;
    return uni == 0 ? 0.0f : static_cast<float>(inter) / static_cast<float>(uni);
}

}  // namespace

std::size_t OverlapGraph::estimate_ram_bytes(std::size_t n_reads) {
    // Sketches: ~100 minimizer hashes per ~10kbp read at k=21,w=11.
    // Each hash 8B → 800 B/read. Plus edge list: ~30 edges/read × 24 B.
    return n_reads * (800 + 30 * 24);
}

void OverlapGraph::build(
    const std::vector<std::pair<ReadId, std::string>>& reads,
    const OverlapGraphOpts& opts)
{
    opts_ = opts;
    edges_.clear();
    const std::size_t N = reads.size();

    // Sketch each read. Heartbeat every 100k reads.
    std::vector<std::vector<std::uint64_t>> sketches(N);
    std::atomic<std::size_t> sketched{0};
    auto t_sk_start = std::chrono::steady_clock::now();
    std::cerr << "[phaser/og] sketching " << N << " reads (k="
              << opts.minimizer_k << " w=" << opts.minimizer_w << ")\n";
    std::fflush(stderr);
    #ifdef _OPENMP
    #pragma omp parallel for num_threads(opts.n_threads) schedule(dynamic, 64)
    #endif
    for (std::size_t i = 0; i < N; ++i) {
        extract_minimizers(reads[i].second,
                           opts.minimizer_k, opts.minimizer_w,
                           sketches[i]);
        const std::size_t done = sketched.fetch_add(1) + 1;
        if (done % 100'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_sk_start).count();
            const double eta = (done > 0)
                ? el * (static_cast<double>(N) - static_cast<double>(done)) /
                       static_cast<double>(done)
                : 0.0;
            std::fprintf(stderr,
                "[phaser/og] sketch %zu/%zu (%.1f%%) elapsed=%.0fs ETA=%.0fs\n",
                done, N, 100.0 * static_cast<double>(done) /
                static_cast<double>(N), el, eta);
            std::fflush(stderr);
        }
    }
    {
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_sk_start).count();
        std::cerr << "[phaser/og] sketches built in " << el << "s\n";
        std::fflush(stderr);
    }

    // Inverted index: minimizer hash → list of read indices carrying it.
    auto t_idx_start = std::chrono::steady_clock::now();
    std::cerr << "[phaser/og] building inverted minimizer index\n";
    std::fflush(stderr);
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> idx;
    idx.reserve(N * 16);
    for (std::size_t i = 0; i < N; ++i) {
        for (auto h : sketches[i]) idx[h].push_back(static_cast<std::uint32_t>(i));
        if ((i + 1) % 500'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_idx_start).count();
            std::fprintf(stderr,
                "[phaser/og] idx %zu/%zu (%.1f%%) entries=%zu elapsed=%.0fs\n",
                i + 1, N, 100.0 * static_cast<double>(i + 1) /
                static_cast<double>(N), idx.size(), el);
            std::fflush(stderr);
        }
    }
    {
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_idx_start).count();
        std::cerr << "[phaser/og] index built in " << el << "s ("
                  << idx.size() << " unique minimizers)\n";
        std::fflush(stderr);
    }

    // Skip ultra-frequent minimizers (likely repeats).
    const std::size_t freq_cap = std::max<std::size_t>(64, N / 100);

    // For each read: candidate set = reads sharing ≥ K_THRESH minimizers.
    const int min_shared = std::max(2, opts.min_overlap_bp / 200);

    // Output buffer per thread, then merged.
    #ifdef _OPENMP
    int nthreads_actual = std::max(1, opts.n_threads);
    #else
    int nthreads_actual = 1;
    #endif
    std::vector<std::vector<OverlapEdge>> per_thread(nthreads_actual);

    auto t_ed_start = std::chrono::steady_clock::now();
    std::atomic<std::size_t> processed{0};
    std::cerr << "[phaser/og] computing pairwise overlap edges\n";
    std::fflush(stderr);
    #ifdef _OPENMP
    #pragma omp parallel num_threads(opts.n_threads)
    #endif
    {
        #ifdef _OPENMP
        const int tid = omp_get_thread_num();
        #else
        const int tid = 0;
        #endif
        std::unordered_map<std::uint32_t, std::uint32_t> shared;
        shared.reserve(64);

        #ifdef _OPENMP
        #pragma omp for schedule(dynamic, 32)
        #endif
        for (std::size_t i = 0; i < N; ++i) {
            shared.clear();
            for (auto h : sketches[i]) {
                auto it = idx.find(h);
                if (it == idx.end()) continue;
                if (it->second.size() > freq_cap) continue;
                for (auto j : it->second) {
                    if (j <= i) continue;  // canonical ordering, no self
                    ++shared[j];
                }
            }
            for (auto& [j, cnt] : shared) {
                if (static_cast<int>(cnt) < min_shared) continue;
                const float jc = sorted_jaccard(sketches[i], sketches[j]);
                if (jc < opts.min_jaccard) continue;
                OverlapEdge e;
                e.a = reads[i].first;
                e.b = reads[j].first;
                e.overlap_bp = static_cast<int>(cnt) * opts.minimizer_w;
                e.jaccard = jc;
                e.a_to_b = true;
                per_thread[tid].push_back(e);
            }
            const std::size_t done = processed.fetch_add(1) + 1;
            if (done % 100'000 == 0) {
                const double el = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t_ed_start).count();
                const double eta = (done > 0)
                    ? el * (static_cast<double>(N) - static_cast<double>(done)) /
                           static_cast<double>(done)
                    : 0.0;
                std::fprintf(stderr,
                    "[phaser/og] edges %zu/%zu (%.1f%%) elapsed=%.0fs ETA=%.0fs\n",
                    done, N, 100.0 * static_cast<double>(done) /
                    static_cast<double>(N), el, eta);
                std::fflush(stderr);
            }
        }
    }
    {
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_ed_start).count();
        std::cerr << "[phaser/og] edges computed in " << el << "s\n";
        std::fflush(stderr);
    }

    for (auto& v : per_thread) {
        edges_.insert(edges_.end(), v.begin(), v.end());
    }
    std::cerr << "[phaser/og] N=" << N << " edges=" << edges_.size() << "\n";
}

void OverlapGraph::collapse_to_unitigs(
    const std::vector<std::pair<ReadId, std::string>>& reads,
    std::vector<UnitigNode>& out_unitigs) const
{
    out_unitigs.clear();

    // Map ReadId → index in reads[]
    std::unordered_map<ReadId, std::size_t> rd_to_idx;
    rd_to_idx.reserve(reads.size() * 2);
    for (std::size_t i = 0; i < reads.size(); ++i) {
        rd_to_idx[reads[i].first] = i;
    }

    // Build adjacency: read_idx → list of (neighbor read_idx, ovl_bp, jc)
    std::vector<std::vector<std::pair<std::size_t, float>>>
        adj(reads.size());
    for (const auto& e : edges_) {
        auto a = rd_to_idx.find(e.a);
        auto b = rd_to_idx.find(e.b);
        if (a == rd_to_idx.end() || b == rd_to_idx.end()) continue;
        adj[a->second].emplace_back(b->second, e.jaccard);
        adj[b->second].emplace_back(a->second, e.jaccard);
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end(),
                  [](auto& x, auto& y){ return x.second > y.second; });
    }

    // Greedy unitig walk: each unvisited read becomes the seed of a new
    // unitig; we extend forward/backward via the highest-jaccard
    // neighbor that itself sees us as its top-1 (mutual-best) and has
    // degree-2 in the overlap graph (interior of a chain).
    std::vector<bool> seen(reads.size(), false);
    auto top_neighbor = [&](std::size_t i) -> std::size_t {
        if (adj[i].empty()) return SIZE_MAX;
        return adj[i][0].first;
    };
    auto degree_2 = [&](std::size_t i) {
        return adj[i].size() == 2;
    };

    auto t_co_start = std::chrono::steady_clock::now();
    std::cerr << "[phaser/og] collapsing into unitigs\n";
    std::fflush(stderr);
    NodeId next_node_id = 0;
    for (std::size_t seed = 0; seed < reads.size(); ++seed) {
        if (seed % 250'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_co_start).count();
            std::fprintf(stderr,
                "[phaser/og] collapse %zu/%zu (%.1f%%) unitigs=%u elapsed=%.0fs\n",
                seed, reads.size(), 100.0 * static_cast<double>(seed) /
                static_cast<double>(reads.size()),
                next_node_id, el);
            std::fflush(stderr);
        }
        if (seen[seed]) continue;

        std::vector<std::size_t> chain;
        chain.push_back(seed);
        seen[seed] = true;

        // Extend forward (top neighbor of current that is mutual + d2)
        std::size_t cur = seed;
        while (true) {
            std::size_t nxt = top_neighbor(cur);
            if (nxt == SIZE_MAX || seen[nxt]) break;
            if (!degree_2(cur) && cur != seed) break;
            if (top_neighbor(nxt) != cur && !chain.empty()) {
                // mutual-best required for clean chain
                if (chain.size() == 1) break;
                break;
            }
            chain.push_back(nxt);
            seen[nxt] = true;
            cur = nxt;
        }
        // Backward extend
        cur = seed;
        while (true) {
            // Find the second-best neighbor that hasn't been used.
            std::size_t back = SIZE_MAX;
            for (auto& [n, jc] : adj[cur]) {
                if (!seen[n]) { back = n; break; }
            }
            if (back == SIZE_MAX) break;
            if (!degree_2(cur) && cur != seed) break;
            chain.insert(chain.begin(), back);
            seen[back] = true;
            cur = back;
        }

        UnitigNode u;
        u.id = next_node_id++;
        u.kind = NodeKind::UNCLASSIFIED;
        u.member_reads.reserve(chain.size());
        // Concatenate sequences with simple overlap-cut at minimizer_w*5 bp
        // (a v1-stub: real consensus uses POA; this is enough to get
        // reasonable lengths and unit testability).
        const int approx_overlap = opts_.minimizer_w * 5;
        std::string concat;
        concat.reserve(chain.size() * 18000);
        for (std::size_t k = 0; k < chain.size(); ++k) {
            const auto& seq = reads[chain[k]].second;
            u.member_reads.push_back(reads[chain[k]].first);
            if (k == 0) concat += seq;
            else if (static_cast<int>(seq.size()) > approx_overlap) {
                concat.append(seq.begin() + approx_overlap, seq.end());
            }
        }
        u.seq = std::move(concat);
        out_unitigs.push_back(std::move(u));
    }

    // Build node-level next/prev edges by scanning original adjacency
    // for cross-unitig links.
    std::unordered_map<std::size_t, NodeId> read_idx_to_node(reads.size() * 2);
    for (std::size_t ni = 0; ni < out_unitigs.size(); ++ni) {
        for (ReadId rid : out_unitigs[ni].member_reads) {
            auto it = rd_to_idx.find(rid);
            if (it != rd_to_idx.end())
                read_idx_to_node[it->second] = static_cast<NodeId>(ni);
        }
    }
    for (std::size_t i = 0; i < reads.size(); ++i) {
        auto it_self = read_idx_to_node.find(i);
        if (it_self == read_idx_to_node.end()) continue;
        for (auto& [j, jc] : adj[i]) {
            auto it_other = read_idx_to_node.find(j);
            if (it_other == read_idx_to_node.end()) continue;
            if (it_self->second == it_other->second) continue;
            auto& self = out_unitigs[it_self->second];
            auto& othr = out_unitigs[it_other->second];
            if (std::find(self.next_nodes.begin(), self.next_nodes.end(),
                          othr.id) == self.next_nodes.end()) {
                self.next_nodes.push_back(othr.id);
            }
            if (std::find(othr.prev_nodes.begin(), othr.prev_nodes.end(),
                          self.id) == othr.prev_nodes.end()) {
                othr.prev_nodes.push_back(self.id);
            }
        }
    }

    std::cerr << "[phaser/og] collapsed " << reads.size()
              << " reads into " << out_unitigs.size() << " unitigs\n";
}

}  // namespace branch::wg::phaser
