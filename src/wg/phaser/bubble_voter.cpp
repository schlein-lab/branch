#include "bubble_voter.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace branch::wg::phaser {

namespace {

inline std::uint8_t base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 4;
    }
}

// Canonical k-mers as raw 2-bit-packed bits (NOT hashed) — must stay
// byte-identical with router's canonical_kmer_bits and with
// phase_engine::scan_canonical_kmers. Any extra hashing here breaks the
// inverted index lookup.
void scan_canonical_kmers(const std::string& seq, int k,
                          std::vector<std::uint64_t>& out) {
    out.clear();
    if (static_cast<int>(seq.size()) < k) return;
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t fwd = 0, rev = 0;
    int valid = 0;
    out.reserve(seq.size());
    for (char c : seq) {
        const std::uint8_t b = base_code(c);
        if (b == 4) { fwd = 0; rev = 0; valid = 0; continue; }
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << (2 * (k - 1)));
        ++valid;
        if (valid >= k) {
            out.push_back(std::min(fwd, rev));
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

// Walk forward from start until reaching sink_anchor, collecting the
// alt-internal unitigs (not including sink itself).
std::vector<NodeId> walk_alt_path(
    const std::vector<UnitigNode>& unitigs,
    NodeId start, NodeId sink)
{
    std::vector<NodeId> out;
    NodeId cur = start;
    int guard = 0;
    while (cur != sink && cur < unitigs.size() && guard++ < 64) {
        out.push_back(cur);
        if (unitigs[cur].next_nodes.empty()) break;
        cur = unitigs[cur].next_nodes[0];
    }
    return out;
}

}  // namespace

void BubbleVoter::vote(
    const std::vector<UnitigNode>& unitigs,
    std::vector<Bubble>& bubbles,
    FastqIndex& reads,
    std::vector<BubbleVote>& out_votes,
    std::vector<VoterReadOutcome>& out_read_outcomes,
    const BubbleVoterOpts& opts)
{
    auto t0 = std::chrono::steady_clock::now();
    out_votes.assign(bubbles.size(), {});

#ifdef _OPENMP
    const int n_threads = omp_get_max_threads();
#else
    const int n_threads = 1;
#endif
    std::cerr << "[phaser/voter] running on " << n_threads << " threads\n";
    std::fflush(stderr);

    // ---- Pass 1: materialize alt-sequences + extract canonical k-mers.
    //
    // Parallelized: each thread processes a chunk of bubbles, accumulates
    // its k-mers locally, then we merge into the global inverted index.
    // The index maps each canonical k-mer (raw 2-bit bits) to the set of
    // (bubble_id, alt_idx) tuples whose alt-sequence contains it. K-mers
    // shared across alts are kept (they contribute weak evidence to each).
    std::unordered_map<std::uint64_t,
        std::vector<std::pair<std::uint32_t, std::uint8_t>>> kmer_index;
    kmer_index.reserve(2'000'000);

    // Resize per-bubble output slots up-front (so parallel pass 2 can
    // write without resizing).
    for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
        const auto& b = bubbles[bi];
        out_votes[bi].read_counts.assign(b.alt_paths.size(), 0);
        out_votes[bi].vafs.assign(b.alt_paths.size(), 0.0f);
        out_votes[bi].reads_per_alt.assign(b.alt_paths.size(), {});
    }

    std::atomic<std::size_t> total_alt_bp{0};
    std::atomic<std::size_t> total_alt_kmers{0};
    std::vector<std::vector<std::tuple<std::uint64_t, std::uint32_t,
                                       std::uint8_t>>> tl_kmer_lists(n_threads);

    #pragma omp parallel
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        auto thread_reader = reads.open_thread_reader();
        auto& my_kmers_out = tl_kmer_lists[tid];
        my_kmers_out.reserve(64 * 1024);
        std::string alt_seq;
        std::string tmp_seq;
        std::vector<std::uint64_t> kmers;

        #pragma omp for schedule(dynamic, 32)
        for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
            const auto& b = bubbles[bi];
            if (b.alt_paths.empty()) continue;
            for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
                auto chain = walk_alt_path(unitigs, b.alt_paths[a],
                                           b.sink_anchor);
                if (chain.empty()) continue;
                // Inline materialize with thread-safe reader.
                alt_seq.clear();
                bool first_read_of_alt = true;
                for (NodeId nid : chain) {
                    const auto& u = unitigs[nid];
                    for (std::size_t k = 0; k < u.member_reads.size(); ++k) {
                        reads.read_seq_via(thread_reader,
                                           u.member_reads[k], tmp_seq);
                        if (tmp_seq.empty()) continue;
                        if (first_read_of_alt) {
                            alt_seq += tmp_seq;
                            first_read_of_alt = false;
                        } else if (static_cast<int>(tmp_seq.size()) >
                                   opts.approx_overlap_bp) {
                            alt_seq.append(
                                tmp_seq.begin() + opts.approx_overlap_bp,
                                tmp_seq.end());
                        }
                    }
                }
                if (alt_seq.empty()) continue;
                total_alt_bp.fetch_add(alt_seq.size(),
                                       std::memory_order_relaxed);
                scan_canonical_kmers(alt_seq, opts.k, kmers);
                total_alt_kmers.fetch_add(kmers.size(),
                                          std::memory_order_relaxed);
                for (auto kmer : kmers) {
                    my_kmers_out.emplace_back(kmer,
                        static_cast<std::uint32_t>(bi),
                        static_cast<std::uint8_t>(a));
                }
            }
        }
    }
    // Merge per-thread kmer lists into the global inverted index.
    for (const auto& tl : tl_kmer_lists) {
        for (const auto& [kmer, bid, alt] : tl) {
            kmer_index[kmer].emplace_back(bid, alt);
        }
    }
    std::vector<std::vector<std::tuple<std::uint64_t, std::uint32_t,
                                       std::uint8_t>>>().swap(tl_kmer_lists);
    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[phaser/voter] pass 1 done — "
              << total_alt_kmers << " alt-kmers indexed (" << kmer_index.size()
              << " unique) in "
              << std::chrono::duration<double>(t1 - t0).count() << "s\n";
    std::fflush(stderr);

    // ---- Pass 2: scan each read, accumulate per-(bubble,alt) hits,
    // assign read to its argmax-alt-per-bubble.
    //
    // For each read R:
    //   hits_per_bubble_alt: bubble_id → vector<uint32_t> sized (n_alts)
    //   For each k-mer in R that's in kmer_index:
    //     For each (bubble_id, alt_idx) entry:
    //       hits[bubble_id][alt_idx]++
    //   For each bubble with hits:
    //     If total_hits >= min_alt_hits:
    //       sort alts by hits descending
    //       if best * (1 - margin) > second_best: assign R to best alt
    //       (else: ambiguous, skip this read for this bubble)
    //
    // R can be assigned to multiple bubbles (one per bubble it overlaps).
    const std::size_t n_reads = reads.n_reads();
    out_read_outcomes.assign(n_reads, VoterReadOutcome{});
    for (ReadId rid = 0; rid < n_reads; ++rid) {
        out_read_outcomes[rid].read_id = rid;
        out_read_outcomes[rid].kind = VoterReadOutcome::ABSENT;
    }

    // Parallel pass 2 with thread-local accumulators (merged at end).
    // The naive approach (mutex per write) would serialize the inner
    // loop; we instead let each thread accumulate into its own
    // tl_votes vector and merge once when the parallel section
    // closes.
    //
    // FastqIndex is not thread-safe for random read (shared FILE*
    // position) — each thread opens its own file handle via
    // open_thread_reader(). Memory cost: 8 threads × O(n_bubbles ×
    // n_alts × 4B) = ~250 KB per thread + per-thread read_per_alt
    // vectors.
    std::atomic<std::size_t> n_assigned{0}, n_bridge{0},
                             n_novel{0}, n_absent{0};
    std::atomic<std::size_t> reads_done{0};
    auto t_voting_start = std::chrono::steady_clock::now();

    // Per-thread vote accumulators — sized to match out_votes layout.
    std::vector<std::vector<BubbleVote>> tl_votes(n_threads);
    for (auto& v : tl_votes) {
        v.assign(bubbles.size(), {});
        for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
            v[bi].read_counts.assign(bubbles[bi].alt_paths.size(), 0);
            v[bi].reads_per_alt.assign(bubbles[bi].alt_paths.size(), {});
        }
    }

    #pragma omp parallel
    {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        auto thread_reader = reads.open_thread_reader();
        auto& my_votes = tl_votes[tid];
        std::string read_seq;
        std::vector<std::uint64_t> read_kmers;
        std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> hits;
        hits.reserve(1024);

        #pragma omp for schedule(dynamic, 4096)
        for (ReadId rid = 0; rid < n_reads; ++rid) {
        reads.read_seq_via(thread_reader, rid, read_seq);
        if (read_seq.empty()) {
            ++n_absent;
            reads_done.fetch_add(1, std::memory_order_relaxed);
            continue;
        }
        scan_canonical_kmers(read_seq, opts.k, read_kmers);
        const std::size_t n_total_kmers = read_kmers.size();
        hits.clear();
        std::uint32_t total_alt_hits = 0;
        for (auto kmer : read_kmers) {
            auto it = kmer_index.find(kmer);
            if (it == kmer_index.end()) continue;
            for (const auto& [bid, alt_idx] : it->second) {
                auto& v = hits[bid];
                if (v.size() <= alt_idx) v.resize(alt_idx + 1, 0);
                ++v[alt_idx];
                ++total_alt_hits;
            }
        }

        // Decision tree (preserving every read with any bubble overlap):
        //
        //   total_alt_hits < min_alt_hits → ABSENT (truly out-of-bubble)
        //   else:
        //     For each bubble in hits:
        //       argmax alt + margin test:
        //         clear → assign to alt
        //         ambiguous → mark this bubble as bridge for this read
        //     If any clear assignment → ASSIGNED (primary = first clear)
        //     Else if any bridge → BRIDGE
        //     Else: read overlapped bubbles but no decision → check
        //       match-fraction: if alt-hits / total_kmers < threshold →
        //       NOVEL (overlaps bubble region with mostly-new sequence);
        //       else: ABSENT (sparse but unhelpful).
        if (total_alt_hits < opts.min_alt_hits) {
            out_read_outcomes[rid].kind = VoterReadOutcome::ABSENT;
            ++n_absent;
        } else {
            bool any_assigned = false;
            bool any_bridge = false;
            std::uint32_t first_assigned_bid = ~0u;
            std::uint8_t  first_assigned_alt = 0xFF;
            for (auto& [bid, v] : hits) {
                if (v.empty()) continue;
                // argmax + runner-up
                std::uint32_t best_idx = 0;
                std::uint32_t best = 0, second = 0;
                for (std::uint32_t i = 0; i < v.size(); ++i) {
                    if (v[i] > best) {
                        second = best;
                        best = v[i]; best_idx = i;
                    } else if (v[i] > second) {
                        second = v[i];
                    }
                }
                if (best == 0) continue;
                if (second > 0 &&
                    static_cast<float>(best) <
                    static_cast<float>(second) * (1.0f + opts.assign_margin)) {
                    // Ambiguous → bridge (CNV-junction evidence)
                    my_votes[bid].bridge_reads.push_back(rid);
                    any_bridge = true;
                } else {
                    ++my_votes[bid].read_counts[best_idx];
                    my_votes[bid].reads_per_alt[best_idx].push_back(rid);
                    if (!any_assigned) {
                        first_assigned_bid = bid;
                        first_assigned_alt = static_cast<std::uint8_t>(best_idx);
                    }
                    any_assigned = true;
                }
            }
            if (any_assigned) {
                out_read_outcomes[rid].kind = VoterReadOutcome::ASSIGNED;
                out_read_outcomes[rid].primary_bubble = first_assigned_bid;
                out_read_outcomes[rid].primary_alt    = first_assigned_alt;
                ++n_assigned;
            } else if (any_bridge) {
                out_read_outcomes[rid].kind = VoterReadOutcome::BRIDGE;
                ++n_bridge;
            } else {
                // Read had alt-hits but no per-bubble decision was made.
                // Check match-fraction for NOVEL classification.
                const float match_frac =
                    n_total_kmers > 0
                        ? static_cast<float>(total_alt_hits) /
                              static_cast<float>(n_total_kmers)
                        : 0.0f;
                if (match_frac < opts.novel_match_frac_threshold) {
                    out_read_outcomes[rid].kind = VoterReadOutcome::NOVEL;
                    ++n_novel;
                } else {
                    out_read_outcomes[rid].kind = VoterReadOutcome::ABSENT;
                    ++n_absent;
                }
            }
        }
        const std::size_t r_done =
            reads_done.fetch_add(1, std::memory_order_relaxed) + 1;
        if (opts.log_every > 0 &&
            r_done % opts.log_every == 0 &&
            tid == 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_voting_start).count();
            std::fprintf(stderr,
                "[phaser/voter] voting %zu/%zu (%.1f%%) "
                "assigned=%zu bridge=%zu novel=%zu absent=%zu elapsed=%.0fs\n",
                r_done, n_reads,
                100.0 * static_cast<double>(r_done) /
                static_cast<double>(n_reads),
                n_assigned.load(), n_bridge.load(),
                n_novel.load(), n_absent.load(), el);
            std::fflush(stderr);
        }
        }  // omp for
    }  // omp parallel

    // Merge thread-local accumulators into out_votes.
    auto t_merge_start = std::chrono::steady_clock::now();
    for (const auto& thread_v : tl_votes) {
        for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
            if (bi >= thread_v.size()) continue;
            const auto& src = thread_v[bi];
            auto& dst = out_votes[bi];
            for (std::size_t a = 0; a < src.read_counts.size(); ++a) {
                if (a >= dst.read_counts.size()) continue;
                dst.read_counts[a] += src.read_counts[a];
                dst.reads_per_alt[a].insert(
                    dst.reads_per_alt[a].end(),
                    src.reads_per_alt[a].begin(),
                    src.reads_per_alt[a].end());
            }
            dst.bridge_reads.insert(dst.bridge_reads.end(),
                                    src.bridge_reads.begin(),
                                    src.bridge_reads.end());
        }
    }
    auto t_merge_end = std::chrono::steady_clock::now();
    std::cerr << "[phaser/voter] tl_votes merge done in "
              << std::chrono::duration<double>(t_merge_end - t_merge_start).count()
              << "s\n";
    std::fflush(stderr);
    {
        // Empty the "extra" reads-processed loop placeholder removed
        // below — closing scope was repositioned for omp_parallel.
    }

    // ---- Pass 3: compute VAFs from read_counts.
    for (auto& v : out_votes) {
        std::uint32_t total = 0;
        for (auto c : v.read_counts) total += c;
        if (total == 0) {
            std::fill(v.vafs.begin(), v.vafs.end(), 0.0f);
            continue;
        }
        for (std::size_t i = 0; i < v.read_counts.size(); ++i) {
            v.vafs[i] = static_cast<float>(v.read_counts[i]) /
                        static_cast<float>(total);
        }
    }

    auto t_end = std::chrono::steady_clock::now();
    const std::size_t reads_processed = reads_done.load();
    const std::size_t na = n_assigned.load(), nb = n_bridge.load(),
                      nn = n_novel.load(), nab = n_absent.load();
    const double pct_assigned = 100.0 *
        static_cast<double>(na) /
        static_cast<double>(std::max<std::size_t>(1, reads_processed));
    const double pct_kept = 100.0 *
        static_cast<double>(na + nb + nn) /
        static_cast<double>(std::max<std::size_t>(1, reads_processed));
    std::cerr << "[phaser/voter] done in "
              << std::chrono::duration<double>(t_end - t0).count() << "s — "
              << reads_processed << " reads scanned | "
              << "assigned=" << na << " (" << pct_assigned << "%) "
              << "bridge=" << nb << " "
              << "novel=" << nn << " "
              << "absent=" << nab << " | "
              << "captured=" << pct_kept << "%\n";
    std::fflush(stderr);
}

void BubbleVoter::populate_branch_trees(
    std::vector<Bubble>& bubbles,
    const std::vector<BubbleVote>& votes)
{
    for (std::size_t bi = 0; bi < bubbles.size(); ++bi) {
        auto& b = bubbles[bi];
        if (bi >= votes.size()) continue;
        const auto& v = votes[bi];
        if (v.read_counts.empty()) continue;

        // Mirror VAFs into the Bubble (downstream consumers expect this).
        b.read_counts = v.read_counts;
        b.vafs = v.vafs;

        // Build a flat tree: root founder + one child per alt.
        b.tree.nodes.clear();
        b.tree.root_node_idx = 0;

        BranchNode root;
        root.depth = 0;
        root.vaf = 1.0f;
        root.alt_idx = ~0u;
        // Children indices filled below.
        b.tree.nodes.push_back(root);

        for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
            BranchNode child;
            child.depth = 1;
            child.vaf = (a < v.vafs.size()) ? v.vafs[a] : 0.0f;
            child.alt_idx = static_cast<std::uint32_t>(a);
            if (a < v.reads_per_alt.size()) {
                child.supporting_reads = v.reads_per_alt[a];
            }
            b.tree.nodes.push_back(child);
            b.tree.nodes[0].children.push_back(
                static_cast<std::uint32_t>(b.tree.nodes.size() - 1));
        }
    }
}

}  // namespace branch::wg::phaser
