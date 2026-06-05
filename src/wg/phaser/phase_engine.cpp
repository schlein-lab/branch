#include "phase_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <unordered_map>

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

// Emit each k-mer as its canonical 2-bit-packed form (NOT a hash). Must
// stay byte-identical with the router's `canonical_kmer_bits()` since
// the het_pair table on the consumer side keys on those raw bits — any
// extra hashing here breaks lookup (v13 measurement: votes 0/0).
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

}  // namespace

void PhaseEngine::prepare(
    FastqIndex& reads,
    const std::vector<HetPair>& het_pairs)
{
    auto t_start = std::chrono::steady_clock::now();
    const std::size_t N = reads.n_reads();
    std::cerr << "[phaser/pe] prepare: indexing " << het_pairs.size()
              << " het-pairs over " << N << " reads\n";
    std::fflush(stderr);

    // Map allele-A/B kmer hashes → (hetpair_idx, allele)
    std::unordered_map<std::uint64_t, std::pair<std::uint32_t, std::uint8_t>>
        allele_lut;
    allele_lut.reserve(het_pairs.size() * 2 + 16);
    for (std::size_t i = 0; i < het_pairs.size(); ++i) {
        allele_lut[het_pairs[i].kmer_a] = {static_cast<std::uint32_t>(i), 0};
        allele_lut[het_pairs[i].kmer_b] = {static_cast<std::uint32_t>(i), 1};
    }
    // Snapshot anchor nodes — needed in run_iteration for anchor-based
    // voting after het_pairs goes out of scope.
    hetpair_anchors_.clear();
    hetpair_anchors_.reserve(het_pairs.size());
    for (const auto& hp : het_pairs) {
        hetpair_anchors_.push_back(hp.anchor_node);
    }

    // ReadId == sequential 0..N-1 from FastqIndex; signatures aligned.
    read_signatures_.assign(N, {});

    std::vector<std::uint64_t> hashes;
    const int k = 21;  // het-pair k-mer size, must match Phase 0 settings
    std::size_t processed = 0;
    reads.for_each([&](ReadId rid, const std::string& seq) {
        if (processed % 200'000 == 0 && processed > 0) {
            const double el = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t_start).count();
            std::fprintf(stderr,
                "[phaser/pe] prepare %zu/%zu (%.1f%%) elapsed=%.0fs\n",
                processed, N,
                100.0 * static_cast<double>(processed) /
                static_cast<double>(N), el);
            std::fflush(stderr);
        }
        scan_canonical_kmers(seq, k, hashes);
        if (rid >= read_signatures_.size()) return;
        auto& sig = read_signatures_[rid];
        sig.clear();
        for (auto h : hashes) {
            auto it = allele_lut.find(h);
            if (it != allele_lut.end()) {
                sig.emplace_back(it->second.first, it->second.second);
            }
        }
        ++processed;
    });
}

PhaseEngineProgress PhaseEngine::run_iteration(
    std::vector<UnitigNode>& unitigs,
    std::vector<Bubble>& bubbles,
    std::vector<ReadAssignment>& assignments,
    const PhaseEngineOpts& opts)
{
    PhaseEngineProgress p;

    // Build per-bubble vote tally. For each unitig in a bubble, iterate
    // its member reads; each read brings allele 0/1 evidence per
    // het-pair index; those votes determine which alt path is h1 vs h2.

    // First pass: tag SHARED for every unitig outside a bubble.
    for (auto& u : unitigs) {
        if (u.bubble == ~0u && u.kind == NodeKind::UNCLASSIFIED) {
            u.kind = NodeKind::SHARED;
        }
    }

    // Build node -> bubble alt-path index lookup.
    // For each bubble, alt_paths holds the head NodeId of each path.
    // We classify the head's full chain as one side.
    std::unordered_map<NodeId, std::pair<BubbleId, std::uint8_t>> node_to_alt;
    for (auto& b : bubbles) {
        for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
            // Walk forward from alt_paths[a] until reaching the sink.
            NodeId cur = b.alt_paths[a];
            int guard = 0;
            while (cur != b.sink_anchor && cur < unitigs.size() && guard++ < 64) {
                node_to_alt[cur] = {b.id, static_cast<std::uint8_t>(a)};
                if (unitigs[cur].next_nodes.empty()) break;
                cur = unitigs[cur].next_nodes[0];
            }
        }
    }

    // v0.10 Variante B: bubble_voter has already done the read→alt
    // assignment by sequence similarity. We consume Bubble.read_counts
    // directly — no per-iteration re-counting from home_node.
    //
    // For each bubble:
    //   - Total reads supporting this bubble = sum(read_counts)
    //   - Per-alt support = read_counts[alt_idx]
    //   - We decide using these directly; het-pair signatures stay
    //     available only for B-3 (molecular-clock cross-check), not
    //     used here.
    std::uint64_t total_supporting = 0;
    for (const auto& b : bubbles) {
        for (auto c : b.read_counts) total_supporting += c;
    }
    std::cerr << "[phaser/pe] voter input: " << total_supporting
              << " reads attributed across " << bubbles.size() << " bubbles\n";
    std::fflush(stderr);

    // Helper: walk an alt-path and tag its interior nodes (skipping
    // anchors, which stay UNCLASSIFIED → SHARED).
    auto tag_alt = [&](const Bubble& b, NodeId start, NodeKind kind) {
        NodeId cur = start;
        int guard = 0;
        while (cur != b.sink_anchor && cur < unitigs.size() && guard++ < 64) {
            if (unitigs[cur].kind == NodeKind::UNCLASSIFIED) {
                unitigs[cur].kind = kind;
            }
            if (unitigs[cur].next_nodes.empty()) break;
            cur = unitigs[cur].next_nodes[0];
        }
    };

    // Decide each bubble: dispatch on BubbleType (set by classifier from
    // voter VAFs). Voter+classifier already answered "what is this
    // bubble?" via sequence-similarity read attribution; phase_engine
    // only does the final node-kind tagging.
    //
    //   HET (D0)        [50, 50]            → alt 0 → H1, alt 1 → H2
    //   BRANCH (D1)     [50, 25, 25]        → parent → H1, children → BRANCH
    //   BIB (D2)        [50, 25, 12.5, 12.5]→ parent → H1, 3 children → BRANCH
    //   TRIALLELIC / NOISY / UNKNOWN        → skip
    //
    // H1/H2 labels are per-bubble arbitrary (consistent only by parent_alt
    // = highest VAF). The anchored phaser flips globally when needed.
    for (auto& b : bubbles) {
        if (b.phased) continue;

        std::uint64_t total = 0;
        for (auto c : b.read_counts) total += c;
        if (total < static_cast<std::uint64_t>(opts.min_supporting_overlaps)) {
            continue;
        }

        if (b.type == BubbleType::HET && b.alt_paths.size() == 2) {
            tag_alt(b, b.alt_paths[0], NodeKind::BUBBLE_H1);
            tag_alt(b, b.alt_paths[1], NodeKind::BUBBLE_H2);
            b.phased = true;
            ++p.bubbles_newly_phased;
        } else if ((b.type == BubbleType::BRANCH ||
                    b.type == BubbleType::BRANCH_IN_BRANCH)
                   && b.parent_alt < b.alt_paths.size()) {
            for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
                NodeKind k = (a == b.parent_alt) ? NodeKind::BUBBLE_H1
                                                 : NodeKind::BRANCH;
                tag_alt(b, b.alt_paths[a], k);
            }
            b.phased = true;
            ++p.bubbles_newly_phased;
        }
    }

    // Build node → branch_idx map for BRANCH-classified bubbles. Each
    // non-parent alt of every BRANCH / BRANCH_IN_BRANCH bubble gets a
    // unique branch_idx (globally across the run) so that downstream
    // output can emit branch_0.fa, branch_1.fa, …
    std::unordered_map<NodeId, std::uint32_t> node_to_branch_idx;
    std::uint32_t next_branch_idx = 0;
    for (auto& b : bubbles) {
        if (b.type != BubbleType::BRANCH &&
            b.type != BubbleType::BRANCH_IN_BRANCH) continue;
        for (std::size_t a = 0; a < b.alt_paths.size(); ++a) {
            if (a == b.parent_alt) continue;
            NodeId cur = b.alt_paths[a];
            int guard = 0;
            while (cur != b.sink_anchor && cur < unitigs.size() && guard++ < 64) {
                if (unitigs[cur].kind == NodeKind::BRANCH) {
                    node_to_branch_idx[cur] = next_branch_idx;
                }
                if (unitigs[cur].next_nodes.empty()) break;
                cur = unitigs[cur].next_nodes[0];
            }
            ++next_branch_idx;
        }
    }

    // Per-read tag: based on home_node kind.
    for (auto& a : assignments) {
        if (a.tag != ReadTag::UNTAGGED) continue;
        if (a.home_node >= unitigs.size()) continue;
        const auto k = unitigs[a.home_node].kind;
        ReadTag t = ReadTag::UNTAGGED;
        switch (k) {
            case NodeKind::SHARED:    t = ReadTag::SHARED; break;
            case NodeKind::BUBBLE_H1: t = ReadTag::H1_ONLY; break;
            case NodeKind::BUBBLE_H2: t = ReadTag::H2_ONLY; break;
            case NodeKind::BRANCH: {
                t = ReadTag::BRANCH;
                auto it = node_to_branch_idx.find(a.home_node);
                if (it != node_to_branch_idx.end()) {
                    a.branch_idx = it->second;
                }
                break;
            }
            default: break;
        }
        if (t != ReadTag::UNTAGGED) {
            a.tag = t;
            a.confidence = 1.0f;  // tag derived from phased graph
            ++p.reads_newly_tagged;
        }
    }

    // Force-assign remainder if requested.
    if (opts.force_assign_remainder) {
        for (auto& a : assignments) {
            if (a.tag == ReadTag::UNTAGGED) {
                a.tag = ReadTag::UNCERTAIN;
                a.confidence = 0.0f;
                ++p.reads_newly_tagged;
            }
        }
    }
    for (auto& a : assignments) {
        if (a.tag == ReadTag::UNTAGGED) ++p.reads_still_untagged;
    }
    return p;
}

void PhaseEngine::apply_voter_outcomes(
    const std::vector<UnitigNode>& unitigs,
    const std::vector<Bubble>& bubbles,
    const std::vector<VoterReadOutcome>& outcomes,
    std::vector<ReadAssignment>& assignments)
{
    // Map alt-head NodeId → NodeKind for fast lookup. (Bubble alt-paths
    // hold the head NodeId of each path; its kind was stamped by
    // run_iteration based on classifier output.)
    std::unordered_map<NodeId, NodeKind> alt_head_kind;
    for (const auto& b : bubbles) {
        for (NodeId head : b.alt_paths) {
            if (head < unitigs.size()) {
                alt_head_kind[head] = unitigs[head].kind;
            }
        }
    }

    std::size_t n_overrides_h1 = 0, n_overrides_h2 = 0,
                n_overrides_branch = 0, n_bridge = 0, n_novel = 0;
    for (auto& a : assignments) {
        if (a.read_id >= outcomes.size()) continue;
        const auto& o = outcomes[a.read_id];
        switch (o.kind) {
            case VoterReadOutcome::ASSIGNED: {
                if (o.primary_bubble >= bubbles.size()) break;
                const auto& b = bubbles[o.primary_bubble];
                if (o.primary_alt >= b.alt_paths.size()) break;
                NodeId head = b.alt_paths[o.primary_alt];
                auto it = alt_head_kind.find(head);
                if (it == alt_head_kind.end()) break;
                switch (it->second) {
                    case NodeKind::BUBBLE_H1:
                        a.tag = ReadTag::H1_ONLY; ++n_overrides_h1; break;
                    case NodeKind::BUBBLE_H2:
                        a.tag = ReadTag::H2_ONLY; ++n_overrides_h2; break;
                    case NodeKind::BRANCH:
                        a.tag = ReadTag::BRANCH; ++n_overrides_branch; break;
                    default: break;  // bubble didn't phase → keep prior tag
                }
                a.branch_path = {o.primary_bubble, o.primary_alt};
                a.confidence = 1.0f;
                break;
            }
            case VoterReadOutcome::BRIDGE:
                a.tag = ReadTag::BRANCH_BRIDGE;
                a.confidence = 0.5f;  // by definition not committed to one alt
                ++n_bridge;
                break;
            case VoterReadOutcome::NOVEL:
                a.tag = ReadTag::BRANCH_NOVEL;
                a.confidence = 0.5f;
                ++n_novel;
                break;
            case VoterReadOutcome::ABSENT:
                break;  // leave whatever home_node-based tagging produced
        }
    }
    std::cerr << "[phaser/pe] voter→tag bridge: "
              << "h1=" << n_overrides_h1
              << " h2=" << n_overrides_h2
              << " branch=" << n_overrides_branch
              << " bridge=" << n_bridge
              << " novel=" << n_novel << "\n";
    std::fflush(stderr);
}

}  // namespace branch::wg::phaser
