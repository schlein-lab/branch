#include "phase_engine.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <unordered_map>

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
        default:            return 4;
    }
}

void scan_kmer_hashes(const std::string& seq, int k,
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
            out.push_back(hash64(std::min(fwd, rev)));
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
}

}  // namespace

void PhaseEngine::prepare(
    const std::vector<std::pair<ReadId, std::string>>& reads,
    const std::vector<HetPair>& het_pairs)
{
    // Map allele-A/B kmer hashes → (hetpair_idx, allele)
    std::unordered_map<std::uint64_t, std::pair<std::uint32_t, std::uint8_t>>
        allele_lut;
    allele_lut.reserve(het_pairs.size() * 2 + 16);
    for (std::size_t i = 0; i < het_pairs.size(); ++i) {
        allele_lut[het_pairs[i].kmer_a] = {static_cast<std::uint32_t>(i), 0};
        allele_lut[het_pairs[i].kmer_b] = {static_cast<std::uint32_t>(i), 1};
    }

    // Find the largest read_id to size the signatures array.
    ReadId max_id = 0;
    for (auto& [rid, _] : reads) max_id = std::max(max_id, rid);
    read_signatures_.assign(static_cast<std::size_t>(max_id) + 1, {});

    std::vector<std::uint64_t> hashes;
    int k = 21;  // het-pair k-mer size, must match Phase 0 settings
    for (const auto& [rid, seq] : reads) {
        scan_kmer_hashes(seq, k, hashes);
        auto& sig = read_signatures_[rid];
        sig.clear();
        for (auto h : hashes) {
            auto it = allele_lut.find(h);
            if (it != allele_lut.end()) {
                sig.emplace_back(it->second.first, it->second.second);
            }
        }
    }
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

    // Per-bubble per-alt vote totals: votes[b.id][alt] = total het-pair
    // evidence summed over all reads sitting on that alt's nodes.
    struct AltScore { std::uint64_t allele0 = 0; std::uint64_t allele1 = 0; };
    std::unordered_map<std::uint64_t, AltScore> bubble_alt_scores;
    auto key = [](BubbleId b, std::uint8_t a) {
        return (static_cast<std::uint64_t>(b) << 8) | a;
    };

    for (auto& a : assignments) {
        if (a.tag != ReadTag::UNTAGGED) continue;
        if (a.read_id >= read_signatures_.size()) continue;
        if (a.home_node >= unitigs.size()) continue;
        auto it = node_to_alt.find(a.home_node);
        if (it == node_to_alt.end()) continue;
        for (auto& [hp_idx, allele] : read_signatures_[a.read_id]) {
            auto& s = bubble_alt_scores[key(it->second.first,
                                            it->second.second)];
            if (allele == 0) ++s.allele0; else ++s.allele1;
        }
    }

    // Decide each bubble: alt with strongest skewed vote becomes
    // BUBBLE_H1; the other BUBBLE_H2. Untied/insufficient: leave unphased.
    for (auto& b : bubbles) {
        if (b.phased) continue;
        if (b.alt_paths.size() != 2) continue;  // skip multi-allelic for v1

        AltScore s0 = bubble_alt_scores[key(b.id, 0)];
        AltScore s1 = bubble_alt_scores[key(b.id, 1)];
        const std::uint64_t total0 = s0.allele0 + s0.allele1;
        const std::uint64_t total1 = s1.allele0 + s1.allele1;
        if (total0 < static_cast<std::uint64_t>(opts.min_supporting_overlaps)
         || total1 < static_cast<std::uint64_t>(opts.min_supporting_overlaps)) {
            continue;
        }
        // alt 0 prefers allele X; alt 1 prefers allele Y. If both prefer
        // the same allele, the bubble isn't really diploid → skip.
        const auto winner0 = (s0.allele0 >= s0.allele1) ? 0 : 1;
        const auto winner1 = (s1.allele0 >= s1.allele1) ? 0 : 1;
        if (winner0 == winner1) continue;
        const double r0 = static_cast<double>(std::max(s0.allele0, s0.allele1)) /
                          static_cast<double>(total0);
        const double r1 = static_cast<double>(std::max(s1.allele0, s1.allele1)) /
                          static_cast<double>(total1);
        if (r0 < opts.min_ratio || r1 < opts.min_ratio) continue;

        // Phase: alt 0 → BUBBLE_H1, alt 1 → BUBBLE_H2.
        // Walk each alt path, tag intermediate unitigs.
        auto tag_alt = [&](NodeId start, NodeKind kind) {
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
        tag_alt(b.alt_paths[0], NodeKind::BUBBLE_H1);
        tag_alt(b.alt_paths[1], NodeKind::BUBBLE_H2);
        b.phased = true;
        ++p.bubbles_newly_phased;
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
            case NodeKind::BRANCH:    t = ReadTag::BRANCH; break;
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

}  // namespace branch::wg::phaser
