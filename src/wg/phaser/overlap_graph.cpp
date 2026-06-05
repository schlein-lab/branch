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
// If max_keep > 0, the result is truncated to its `max_keep` smallest
// hashes (bottom-k MinHash signature; reduces RSS at WG scale).
void extract_minimizers(const std::string& seq, int k, int w,
                        std::vector<std::uint64_t>& out_hashes,
                        std::size_t max_keep = 0) {
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
    // Bottom-k MinHash truncation: keep only the `max_keep` smallest
    // hashes. The list is already sorted ascending so this is the
    // bottom-k subset by hash value, which is the canonical MinHash
    // signature.
    if (max_keep > 0 && out_hashes.size() > max_keep) {
        out_hashes.resize(max_keep);
    }
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
                           sketches[i], opts.max_sketch_size);
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

    // Inverted index as a flat sorted array of (hash, read_id) postings.
    // Replaces the prior unordered_map<u64, vector<u32>> which had
    // ~56 bytes/entry overhead and crashed v10c at 23 % build with
    // 233 M entries × 56 B ≈ 13 GB just for index nodes (plus per-vector
    // heap overhead). Flat sorted postings cost 12 B/entry with no
    // pointer chasing and one contiguous heap allocation.
    auto t_idx_start = std::chrono::steady_clock::now();
    std::size_t total_postings = 0;
    for (const auto& s : sketches) total_postings += s.size();
    std::cerr << "[phaser/og] building flat posting array ("
              << total_postings << " postings)\n";
    std::fflush(stderr);
    struct Posting { std::uint64_t hash; std::uint32_t read_id; };
    std::vector<Posting> postings;
    postings.reserve(total_postings);
    for (std::size_t i = 0; i < N; ++i) {
        for (auto h : sketches[i]) {
            postings.push_back({h, static_cast<std::uint32_t>(i)});
        }
        if ((i + 1) % 500'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_idx_start).count();
            std::fprintf(stderr,
                "[phaser/og] postings %zu/%zu (%.1f%%) total=%zu elapsed=%.0fs\n",
                i + 1, N, 100.0 * static_cast<double>(i + 1) /
                static_cast<double>(N), postings.size(), el);
            std::fflush(stderr);
        }
    }
    std::cerr << "[phaser/og] sorting postings\n";
    std::fflush(stderr);
    std::sort(postings.begin(), postings.end(),
              [](const Posting& a, const Posting& b) {
                  return a.hash != b.hash ? a.hash < b.hash
                                          : a.read_id < b.read_id;
              });
    {
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_idx_start).count();
        std::cerr << "[phaser/og] postings ready in " << el << "s ("
                  << postings.size() << " postings)\n";
        std::fflush(stderr);
    }

    // Skip ultra-frequent minimizers (likely repeats). At expected
    // diploid coverage ~24×, a real-overlap minimizer appears in
    // ~24-50 reads. A block in the hundreds is already a repeat region;
    // the previous N/100 cap (43 k for HG002) let blocks of ~10 k reads
    // through and exploded the candidate-pair set into hundreds of
    // millions, OOM-killing v10d at the aggregation step. 200 is a
    // hard cap that keeps real overlaps and drops repeats.
    const std::size_t freq_cap = 200;

    // For each read: candidate set = reads sharing ≥ K_THRESH minimizers.
    const int min_shared = std::max(2, opts.min_overlap_bp / 200);

    // Output buffer per thread, then merged.
    #ifdef _OPENMP
    int nthreads_actual = std::max(1, opts.n_threads);
    #else
    int nthreads_actual = 1;
    #endif
    std::vector<std::vector<OverlapEdge>> per_thread(nthreads_actual);

    // Streaming edge computation: never materialize all pair candidates
    // simultaneously. For each read i, binary-search each of its
    // sketches in the sorted postings, accumulate a per-i shared-count
    // map (small, ~hundreds of entries), filter by min_shared + jaccard,
    // emit edges directly. v10e OOMed at 4 G materialized pair keys —
    // this version keeps only a per-thread tiny map plus the edge list.
    auto t_ed_start = std::chrono::steady_clock::now();
    std::cerr << "[phaser/og] streaming edge computation from sorted postings\n";
    std::fflush(stderr);
    std::atomic<std::size_t> processed{0};

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
        shared.reserve(256);

        #ifdef _OPENMP
        #pragma omp for schedule(dynamic, 64)
        #endif
        for (std::size_t i = 0; i < N; ++i) {
            shared.clear();
            const std::uint32_t i32 = static_cast<std::uint32_t>(i);
            for (auto h : sketches[i]) {
                // Find the contiguous range of postings with this hash.
                auto cmp_lo = [](const Posting& p, std::uint64_t v) {
                    return p.hash < v;
                };
                auto begin = std::lower_bound(postings.begin(),
                                              postings.end(), h, cmp_lo);
                auto end_it = begin;
                while (end_it != postings.end() && end_it->hash == h) ++end_it;
                const std::size_t block =
                    static_cast<std::size_t>(end_it - begin);
                if (block > freq_cap || block < 2) continue;
                for (auto it = begin; it != end_it; ++it) {
                    if (it->read_id <= i32) continue;  // canonical, no self
                    ++shared[it->read_id];
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
    // Free postings (no longer needed; sketches stay until edges merged).
    std::vector<Posting>().swap(postings);
    {
        const double el = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_ed_start).count();
        std::cerr << "[phaser/og] streaming edges done in " << el << "s\n";
        std::fflush(stderr);
    }

    for (auto& v : per_thread) {
        edges_.insert(edges_.end(), v.begin(), v.end());
    }
    std::cerr << "[phaser/og] N=" << N << " edges=" << edges_.size() << "\n";
}

void OverlapGraph::transitive_reduce() {
    if (edges_.empty()) return;

    auto t_start = std::chrono::steady_clock::now();
    const std::size_t n_in = edges_.size();
    std::cerr << "[phaser/og] transitive reduction on " << n_in << " edges\n";
    std::fflush(stderr);

    // Find max ReadId so we can size the adjacency table.
    ReadId max_id = 0;
    for (const auto& e : edges_) {
        if (e.a > max_id) max_id = e.a;
        if (e.b > max_id) max_id = e.b;
    }
    const std::size_t N = static_cast<std::size_t>(max_id) + 1;

    // Adjacency list: per node, vector of (neighbor, jaccard, edge_idx),
    // sorted by jaccard descending. Same edge appears in both endpoints'
    // lists; edge_idx points back to the edges_ entry.
    struct Adj { ReadId b; float jaccard; std::uint32_t edge_idx; };
    std::vector<std::vector<Adj>> adj(N);
    for (std::uint32_t i = 0; i < edges_.size(); ++i) {
        const auto& e = edges_[i];
        adj[e.a].push_back({e.b, e.jaccard, i});
        adj[e.b].push_back({e.a, e.jaccard, i});
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end(),
                  [](const Adj& x, const Adj& y) {
                      return x.jaccard > y.jaccard;
                  });
    }

    // Mark redundant edges. For each direct edge (a → b) with weight
    // j_ab, look for any detour neighbor c of a with j_ac > j_ab such
    // that c also reaches b with j_cb >= j_ab. If found, (a,b) is
    // dominated by the detour and can be deleted.
    //
    // Two passes write to the same `drop` flags from both endpoints; we
    // use a vector<uint8_t> instead of <bool> for safe parallel writes.
    std::vector<std::uint8_t> drop(edges_.size(), 0);

    std::atomic<std::size_t> processed{0};
    auto t_loop = std::chrono::steady_clock::now();

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 1024)
    #endif
    for (std::size_t a = 0; a < N; ++a) {
        const auto& a_neighbors = adj[a];
        for (std::size_t i = 0; i < a_neighbors.size(); ++i) {
            const ReadId b = a_neighbors[i].b;
            const float j_ab = a_neighbors[i].jaccard;
            const std::uint32_t edge_idx = a_neighbors[i].edge_idx;
            if (drop[edge_idx]) continue;

            // Detour candidates are stronger neighbors of a (positions 0..i-1).
            for (std::size_t k = 0; k < i; ++k) {
                const ReadId c = a_neighbors[k].b;
                if (c == b) continue;
                // Linear membership scan in adj[c] for b. deg(c) is
                // small (avg 24) so this stays cache-friendly.
                bool detour_strong = false;
                for (const auto& cn : adj[c]) {
                    if (cn.b == b) {
                        if (cn.jaccard >= j_ab) detour_strong = true;
                        break;
                    }
                }
                if (detour_strong) {
                    drop[edge_idx] = 1;
                    break;
                }
            }
        }
        const std::size_t done = processed.fetch_add(1) + 1;
        if (done % 200'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_loop).count();
            std::fprintf(stderr,
                "[phaser/og] TR scan %zu/%zu (%.1f%%) elapsed=%.0fs\n",
                done, N,
                100.0 * static_cast<double>(done) /
                static_cast<double>(N), el);
            std::fflush(stderr);
        }
    }

    // Compact edges_.
    std::vector<OverlapEdge> kept;
    kept.reserve(edges_.size());
    for (std::size_t i = 0; i < edges_.size(); ++i) {
        if (!drop[i]) kept.push_back(edges_[i]);
    }
    const std::size_t n_out = kept.size();
    edges_ = std::move(kept);

    const double el = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count();
    std::cerr << "[phaser/og] TR done in " << el << "s — kept "
              << n_out << " / " << n_in << " edges ("
              << (100.0 * static_cast<double>(n_out) /
                  static_cast<double>(n_in))
              << "% retained, "
              << (100.0 * static_cast<double>(n_in - n_out) /
                  static_cast<double>(n_in))
              << "% reduced)\n";
    std::fflush(stderr);
}

void OverlapGraph::collapse_to_unitigs(
    FastqIndex& reads,
    std::vector<UnitigNode>& out_unitigs) const
{
    out_unitigs.clear();

    // ReadId == sequential index into FastqIndex, so the previous
    // ReadId→idx map is identity. We just use ReadId as the array index.
    const std::size_t N = reads.n_reads();

    // Build adjacency: read_idx → list of (neighbor read_idx, ovl_bp, jc)
    std::vector<std::vector<std::pair<std::size_t, float>>> adj(N);
    for (const auto& e : edges_) {
        if (e.a >= N || e.b >= N) continue;
        adj[e.a].emplace_back(static_cast<std::size_t>(e.b), e.jaccard);
        adj[e.b].emplace_back(static_cast<std::size_t>(e.a), e.jaccard);
    }
    for (auto& v : adj) {
        std::sort(v.begin(), v.end(),
                  [](auto& x, auto& y){ return x.second > y.second; });
    }

    // Greedy unitig walk: each unvisited read becomes the seed of a new
    // unitig; we extend forward/backward via the highest-jaccard
    // neighbor that itself sees us as its top-1 (mutual-best) and has
    // degree-2 in the overlap graph (interior of a chain).
    std::vector<bool> seen(N, false);
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
    for (std::size_t seed = 0; seed < N; ++seed) {
        if (seed % 250'000 == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_co_start).count();
            std::fprintf(stderr,
                "[phaser/og] collapse %zu/%zu (%.1f%%) unitigs=%u elapsed=%.0fs\n",
                seed, N, 100.0 * static_cast<double>(seed) /
                static_cast<double>(N),
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
        // Compute length only — DO NOT materialize sequence here.
        // gfa_writer reconstructs sequences on-the-fly per unitig at
        // output time, reading them lazily from the FastqIndex.
        const Position approx_overlap =
            static_cast<Position>(opts_.minimizer_w) * 5;
        Position L = 0;
        for (std::size_t k = 0; k < chain.size(); ++k) {
            const ReadId rid = static_cast<ReadId>(chain[k]);
            const Position rlen = reads.seq_length(rid);
            u.member_reads.push_back(rid);
            if (k == 0) {
                L += rlen;
            } else if (rlen > approx_overlap) {
                L += rlen - approx_overlap;
            }
        }
        u.length_bp = L;
        out_unitigs.push_back(std::move(u));
    }

    // Build node-level next/prev edges by scanning original adjacency
    // for cross-unitig links. ReadId == array index throughout.
    auto t_ne_start = std::chrono::steady_clock::now();
    std::cerr << "[phaser/og] building next/prev edges\n";
    std::fflush(stderr);
    std::vector<NodeId> read_to_node(N, ~0u);
    for (std::size_t ni = 0; ni < out_unitigs.size(); ++ni) {
        for (ReadId rid : out_unitigs[ni].member_reads) {
            if (rid < N) read_to_node[rid] = static_cast<NodeId>(ni);
        }
    }
    for (std::size_t i = 0; i < N; ++i) {
        if (i % 500'000 == 0 && i > 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_ne_start).count();
            std::fprintf(stderr,
                "[phaser/og] next/prev %zu/%zu (%.1f%%) elapsed=%.0fs\n",
                i, N,
                100.0 * static_cast<double>(i) /
                static_cast<double>(N), el);
            std::fflush(stderr);
        }
        const NodeId self_node = read_to_node[i];
        if (self_node == ~0u) continue;
        for (auto& [j, jc] : adj[i]) {
            const NodeId other_node = read_to_node[j];
            if (other_node == ~0u) continue;
            if (self_node == other_node) continue;
            auto& self = out_unitigs[self_node];
            auto& othr = out_unitigs[other_node];
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

    std::cerr << "[phaser/og] collapsed " << N
              << " reads into " << out_unitigs.size() << " unitigs\n";
}

}  // namespace branch::wg::phaser
