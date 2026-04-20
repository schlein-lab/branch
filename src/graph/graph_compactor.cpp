// BRANCH v0.2 — Unitig compaction implementation with consensus building.
//
// Algorithm (now four passes):
//   1) Build in-degree / out-degree arrays and per-node adjacency indices
//      into the input edge vector.
//   2) Identify chain starts: a node N starts a unitig if either
//        - in_degree(N) != 1, or
//        - the unique predecessor P has out_degree(P) != 1.
//      From each start, walk forward while the current node has out_degree
//      exactly 1 AND the successor has in_degree exactly 1. Collect all
//      visited nodes into the same unitig.
//   3) Any node that is not yet assigned to a unitig must be on a pure
//      cycle (every node in the cycle has in=out=1). Break such a cycle
//      arbitrarily at an unvisited node and walk it as its own unitig.
//   4) Closed-twin merge. Two distinct singleton unitigs {u}, {v} are
//      "closed twins" iff (N(u) ∪ {u}) == (N(v) ∪ {v}) as undirected
//      sets AND length_bp(u) == length_bp(v). In a read-overlap graph
//      this means u and v represent the SAME genomic region covered by
//      reads of the SAME length — i.e. reads drawn from the same allele
//      that all mutually overlap. Such reads form a tournament/clique
//      that the linear-chain pass CANNOT collapse (every clique member
//      has in>1 and out>1, so every member is its own unitig). Merging
//      closed twins is the missing structural rule. Length equality is
//      the "different-allele firewall" — reads of unequal length never
//      merge even if they share the closed neighborhood.
//   5) Isolated-node bundle-merge. Nodes with in_degree == 0 AND
//      out_degree == 0 contribute no structure (no edges touch them).
//      These are the residue of graph_filter's containment drop — the
//      filter removes every edge incident on a contained node but
//      cannot erase the node entry itself (that would renumber
//      NodeIds). Keeping them as singleton unitigs bloats node_count
//      1:1 with dropped reads; dropping them outright destroys the
//      coverage signal length-bucket VAF aggregation reads out of the
//      compacted graph. Compromise: collapse isolated singletons into
//      one orphan unitig PER DISTINCT length_bp. One bookkeeping node
//      per read-length class preserves the signal without the blow-up.
//   6) Build the output graph: one new node per SURVIVING unitig with
//        length_bp = sum of member lengths (chains with mixed lengths)
//                    or the shared representative length (closed-twin
//                    clusters, orphan-per-length bundles, equal-length
//                    chains)
//        copy_count = chain start's copy_count
//        consensus = MSA-derived consensus from member sequences
//      Copy only inter-unitig edges (from_unitig != to_unitig), keeping
//      the original read_support of the boundary-crossing edge.

#include "graph/graph_compactor.hpp"
#include "consensus/majority_voter.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace branch::graph {

namespace {

constexpr NodeId kUnassigned = std::numeric_limits<NodeId>::max();

// Position-wise majority voting for unaligned sequences.
// Sequences of different lengths are padded with conceptual gaps
// (positions past the end contribute nothing to the base count).
// Adequate for well-overlapping unitig members where indels are rare;
// a POA-based aligner (abPOA or equivalent) would improve quality on
// pathological regions but is not vendored yet — revisit when the
// assembler starts hitting high-indel targets.
std::string simple_majority_consensus(const std::vector<std::string>& sequences) {
    if (sequences.empty()) {
        return {};
    }
    if (sequences.size() == 1) {
        return sequences[0];
    }

    // Find max length
    std::size_t max_len = 0;
    for (const auto& seq : sequences) {
        max_len = std::max(max_len, seq.size());
    }

    std::string result;
    result.reserve(max_len);

    for (std::size_t pos = 0; pos < max_len; ++pos) {
        // Count bases at this position
        std::array<int, 5> counts{}; // A=0, C=1, G=2, T=3, other=4

        for (const auto& seq : sequences) {
            if (pos >= seq.size()) continue;
            char c = seq[pos];
            switch (c) {
                case 'A': case 'a': ++counts[0]; break;
                case 'C': case 'c': ++counts[1]; break;
                case 'G': case 'g': ++counts[2]; break;
                case 'T': case 't': ++counts[3]; break;
                default: ++counts[4]; break;
            }
        }

        // Find majority base
        int max_count = 0;
        int max_idx = 4;
        for (int i = 0; i < 4; ++i) {
            if (counts[i] > max_count) {
                max_count = counts[i];
                max_idx = i;
            }
        }

        constexpr char bases[] = "ACGTN";
        result += bases[max_idx];
    }

    return result;
}

// Build consensus from member node sequences.
// Strategy:
//   - 0 sequences with content -> empty consensus
//   - 1 sequence -> use directly
//   - 2+ sequences -> simple_majority_consensus (position-wise majority
//     with length padding). A POA-based aligner would improve indel-
//     heavy regions but is not yet vendored.
std::string build_unitig_consensus(
    const std::vector<NodeId>& members,
    const LosslessGraph& input) {

    // Collect non-empty consensus sequences from member nodes
    std::vector<std::string> sequences;
    sequences.reserve(members.size());

    for (NodeId m : members) {
        const auto& node = input.node(m);
        if (!node.consensus.empty()) {
            sequences.push_back(node.consensus);
        }
    }

    if (sequences.empty()) {
        return {};
    }

    if (sequences.size() == 1) {
        return sequences[0];
    }

    return simple_majority_consensus(sequences);
}

// Degree arrays + unique-edge lookups for the linear-chain pass.
struct DegreeTables {
    std::vector<std::uint32_t> in_degree;
    std::vector<std::uint32_t> out_degree;
    std::vector<std::size_t>   unique_out_edge;
    std::vector<std::size_t>   unique_in_edge;
};

DegreeTables compute_degree_tables(std::size_t n, const std::vector<Edge>& edges) {
    DegreeTables t;
    t.in_degree.assign(n, 0);
    t.out_degree.assign(n, 0);
    t.unique_out_edge.assign(n, std::numeric_limits<std::size_t>::max());
    t.unique_in_edge.assign(n, std::numeric_limits<std::size_t>::max());

    for (const auto& e : edges) {
        if (e.from < n) ++t.out_degree[e.from];
        if (e.to < n)   ++t.in_degree[e.to];
    }
    for (std::size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        if (e.from < n && t.out_degree[e.from] == 1) {
            t.unique_out_edge[e.from] = i;
        }
        if (e.to < n && t.in_degree[e.to] == 1) {
            t.unique_in_edge[e.to] = i;
        }
    }
    return t;
}

// Walk linear chains + pure cycles from `input` into `unitigs`, writing
// each member's unitig index into `old_to_new`. Returns the unitigs
// vector unchanged (members list order matches chain-walk order).
void walk_linear_chains(std::size_t n,
                        const std::vector<Edge>& edges,
                        const DegreeTables& t,
                        std::vector<std::vector<NodeId>>& unitigs,
                        std::vector<NodeId>& old_to_new) {
    auto is_chain_start = [&](NodeId v) -> bool {
        if (t.in_degree[v] != 1) return true;  // source or branching
        // Safe: in_degree[v] == 1 guarantees unique_in_edge[v] was set.
        const Edge& pe = edges[t.unique_in_edge[v]];
        return t.out_degree[pe.from] != 1;
    };

    for (NodeId v = 0; v < n; ++v) {
        if (old_to_new[v] != kUnassigned) continue;
        if (!is_chain_start(v)) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            members.push_back(cur);
            old_to_new[cur] = unitig_id;

            if (t.out_degree[cur] != 1) break;
            NodeId next = edges[t.unique_out_edge[cur]].to;
            if (next >= n) break;                         // defensive
            if (t.in_degree[next] != 1) break;
            if (old_to_new[next] != kUnassigned) break;   // cycle-safety
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }

    // Pure cycles: every node in-out==1 and no chain start exists.
    for (NodeId v = 0; v < n; ++v) {
        if (old_to_new[v] != kUnassigned) continue;

        NodeId unitig_id = static_cast<NodeId>(unitigs.size());
        std::vector<NodeId> members;
        NodeId cur = v;
        while (true) {
            if (old_to_new[cur] != kUnassigned) break;
            members.push_back(cur);
            old_to_new[cur] = unitig_id;
            if (t.out_degree[cur] != 1) break;
            NodeId next = edges[t.unique_out_edge[cur]].to;
            if (next >= n) break;
            if (t.in_degree[next] != 1) break;
            cur = next;
        }
        unitigs.push_back(std::move(members));
    }
}

// Build, per node, the sorted *closed* undirected neighborhood
// (N(v) ∪ {v}). Self-loops are tolerated: the {v} term guarantees v
// itself is always present exactly once regardless of self-edges.
std::vector<std::vector<NodeId>>
build_closed_neighborhoods(std::size_t n, const std::vector<Edge>& edges) {
    std::vector<std::vector<NodeId>> nbr(n);
    for (const auto& e : edges) {
        if (e.from < n && e.to < n && e.from != e.to) {
            nbr[e.from].push_back(e.to);
            nbr[e.to].push_back(e.from);
        }
    }
    for (NodeId v = 0; v < static_cast<NodeId>(n); ++v) {
        nbr[v].push_back(v);                              // closed
        std::sort(nbr[v].begin(), nbr[v].end());
        nbr[v].erase(std::unique(nbr[v].begin(), nbr[v].end()), nbr[v].end());
    }
    return nbr;
}

// Merge unitigs that are "closed twins": both are singletons whose sole
// member has the same closed undirected neighborhood and the same
// length_bp. Returns a remap table: for each old unitig index, the
// (possibly merged) new unitig index. Also rewrites `unitigs` in place
// so the final unitig set reflects the merges, and rewrites
// `old_to_new` to track the new indices for every old NodeId.
//
// Why length_bp gates the merge: in a read-overlap graph, two reads
// covering the same genomic locus have (near-)identical closed
// neighborhoods (every third read overlapping one typically overlaps
// the other). But they may come from DIFFERENT alleles — e.g. a 9 kb
// k=2 read and a 19 kb k=6 read from the synthetic 4-allele fixture.
// Collapsing those loses the allele distinction. Gating on exact
// length_bp is a conservative proxy for "same allele read length";
// real polymerase data needs ±ε tolerance which will be revisited
// when a CN-aware cluster step lands.
void merge_closed_twins(std::size_t n,
                        const LosslessGraph& input,
                        const std::vector<Edge>& edges,
                        std::vector<std::vector<NodeId>>& unitigs,
                        std::vector<NodeId>& old_to_new) {
    if (unitigs.size() < 2) return;

    // Index: unitig_id -> sole member if the unitig is a singleton, else
    // kUnassigned (signalling "not mergeable by this pass").
    std::vector<NodeId> singleton_member(unitigs.size(), kUnassigned);
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (unitigs[u].size() == 1) singleton_member[u] = unitigs[u][0];
    }

    auto closed_nbrs = build_closed_neighborhoods(n, edges);

    // Group singleton unitigs by (length_bp, closed_neighborhood).
    // Two unitigs are in the same group iff both keys agree exactly.
    //
    // Representation: we keep a vector of (length, neighborhood_vec,
    // unitig_id) and sort lexicographically. Adjacent entries that
    // match on (length, neighborhood_vec) form a class.
    struct GroupKey {
        std::uint32_t length_bp;
        const std::vector<NodeId>* nbr;
        NodeId unitig_id;
    };
    std::vector<GroupKey> keys;
    keys.reserve(unitigs.size());
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (singleton_member[u] == kUnassigned) continue;
        NodeId m = singleton_member[u];
        keys.push_back(GroupKey{
            .length_bp = input.node(m).length_bp,
            .nbr = &closed_nbrs[m],
            .unitig_id = static_cast<NodeId>(u),
        });
    }

    std::sort(keys.begin(), keys.end(),
              [](const GroupKey& a, const GroupKey& b) {
                  if (a.length_bp != b.length_bp) return a.length_bp < b.length_bp;
                  return *a.nbr < *b.nbr;  // lexicographic on sorted vectors
              });

    // Canonical unitig index for every old unitig. Initially identity.
    std::vector<NodeId> canonical(unitigs.size());
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        canonical[u] = static_cast<NodeId>(u);
    }

    for (std::size_t i = 0; i < keys.size();) {
        std::size_t j = i + 1;
        while (j < keys.size() &&
               keys[j].length_bp == keys[i].length_bp &&
               *keys[j].nbr == *keys[i].nbr) {
            ++j;
        }
        // [i, j) is a twin-class. If size >= 2, merge into the lowest
        // unitig_id in the class (stable, reproducible).
        if (j - i >= 2) {
            NodeId rep = keys[i].unitig_id;
            for (std::size_t k = i; k < j; ++k) {
                if (keys[k].unitig_id < rep) rep = keys[k].unitig_id;
            }
            for (std::size_t k = i; k < j; ++k) {
                canonical[keys[k].unitig_id] = rep;
            }
        }
        i = j;
    }

    // Determine which of the original unitigs survive (those whose
    // canonical index maps to themselves). Assign new contiguous ids.
    std::vector<NodeId> old_unitig_to_new(unitigs.size(), kUnassigned);
    std::vector<std::vector<NodeId>> new_unitigs;
    new_unitigs.reserve(unitigs.size());
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (canonical[u] == static_cast<NodeId>(u)) {
            old_unitig_to_new[u] = static_cast<NodeId>(new_unitigs.size());
            new_unitigs.push_back(std::move(unitigs[u]));
        }
    }
    // Merge non-representative unitigs' members into the rep's member list.
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (canonical[u] == static_cast<NodeId>(u)) continue;
        NodeId rep = canonical[u];
        NodeId new_rep = old_unitig_to_new[rep];
        for (NodeId m : unitigs[u]) {
            new_unitigs[new_rep].push_back(m);
        }
    }

    // Rewrite old_to_new: every old NodeId now points to the new unitig
    // index (after merge, after compaction).
    for (NodeId v = 0; v < static_cast<NodeId>(n); ++v) {
        if (old_to_new[v] == kUnassigned) continue;
        NodeId old_u = old_to_new[v];
        NodeId canon_u = canonical[old_u];
        old_to_new[v] = old_unitig_to_new[canon_u];
    }

    unitigs = std::move(new_unitigs);
}

// Collapse every isolated unitig (members all have in=out=0) into a
// one orphan unitig PER distinct length_bp. Rationale: these nodes
// are the residue of graph_filter's containment-drop pass — the filter
// erases their edges but cannot renumber NodeIds, so the node entries
// linger. They carry no topology, but dropping them entirely would
// destroy coverage-like signals that downstream consumers
// (e.g. length-bucket VAF aggregation) read out of the compacted
// graph. Grouping by length_bp preserves those signals: isolated
// reads from the same allele cluster into one node; reads from
// different alleles stay distinguishable by length.
//
// Multi-member unitigs (linear chains and closed-twin clusters) can
// never be fully isolated: chains have internal edges; twin clusters
// have the very neighborhood-agreement edge set that gave them their
// closed-twin property. So only singleton unitigs are eligible.
void merge_isolated_unitigs(std::size_t n,
                            const LosslessGraph& input,
                            const DegreeTables& t,
                            std::vector<std::vector<NodeId>>& unitigs,
                            std::vector<NodeId>& old_to_new) {
    // Identify isolated singletons and their length.
    std::vector<bool> is_isolated(unitigs.size(), false);
    std::vector<std::uint32_t> iso_length(unitigs.size(), 0);
    std::size_t isolated_count = 0;
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (unitigs[u].size() != 1) continue;
        NodeId m = unitigs[u].front();
        if (t.in_degree[m] == 0 && t.out_degree[m] == 0) {
            is_isolated[u] = true;
            iso_length[u] = input.node(m).length_bp;
            ++isolated_count;
        }
    }
    if (isolated_count < 2) return;

    std::vector<NodeId> old_unitig_to_new(unitigs.size(), kUnassigned);
    std::vector<std::vector<NodeId>> kept;
    kept.reserve(unitigs.size() - isolated_count + 16);

    // First pass: keep non-isolated unitigs.
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (is_isolated[u]) continue;
        old_unitig_to_new[u] = static_cast<NodeId>(kept.size());
        kept.push_back(std::move(unitigs[u]));
    }

    // Second pass: fold isolated singletons into one orphan per length.
    // Linear search over len_to_bundle is fine — the number of
    // distinct isolated-node lengths is bounded by the number of
    // distinct read lengths in the input, which is typically ≤ a few
    // dozen (polymerase run-classes). A map is overkill.
    struct LenBundle {
        std::uint32_t length_bp;
        NodeId        new_unitig_id;
    };
    std::vector<LenBundle> len_to_bundle;
    for (std::size_t u = 0; u < unitigs.size(); ++u) {
        if (!is_isolated[u]) continue;
        NodeId bundle_id = kUnassigned;
        for (const auto& lb : len_to_bundle) {
            if (lb.length_bp == iso_length[u]) {
                bundle_id = lb.new_unitig_id;
                break;
            }
        }
        if (bundle_id == kUnassigned) {
            bundle_id = static_cast<NodeId>(kept.size());
            len_to_bundle.push_back({iso_length[u], bundle_id});
            kept.push_back(std::vector<NodeId>{});
        }
        old_unitig_to_new[u] = bundle_id;
        kept[bundle_id].push_back(unitigs[u].front());
    }

    // Rewrite old_to_new.
    for (NodeId v = 0; v < static_cast<NodeId>(n); ++v) {
        if (old_to_new[v] == kUnassigned) continue;
        old_to_new[v] = old_unitig_to_new[old_to_new[v]];
    }
    unitigs = std::move(kept);
}

// Pick a single length_bp value for a unitig with `members`. Rules:
//   - Singleton: the member's length.
//   - Multi-member, all-same length: that shared length (closed-twin
//     cluster, equal-length linear chain, or orphan-per-length bundle
//     — all three want the representative length, NOT the sum).
//   - Multi-member, mixed lengths: linear chain with varying read
//     lengths. Sum lengths as before. (Twin clusters and orphan
//     bundles are same-length by construction, so this branch only
//     fires for true chains.)
std::uint32_t pick_unitig_length(const std::vector<NodeId>& members,
                                 const LosslessGraph& input) {
    if (members.size() == 1) {
        return input.node(members.front()).length_bp;
    }

    const std::uint32_t first_len = input.node(members.front()).length_bp;
    bool all_same_length = true;
    std::uint64_t total_length = 0;
    for (NodeId m : members) {
        const std::uint32_t len = input.node(m).length_bp;
        if (len != first_len) all_same_length = false;
        total_length += len;
    }
    if (all_same_length) return first_len;
    return (total_length > std::numeric_limits<std::uint32_t>::max())
        ? std::numeric_limits<std::uint32_t>::max()
        : static_cast<std::uint32_t>(total_length);
}

}  // namespace

CompactionResult compact_unitigs(const LosslessGraph& input) {
    CompactionResult result;
    const std::size_t n = input.node_count();
    result.old_to_new.assign(n, kUnassigned);

    if (n == 0) {
        return result;
    }

    const auto& edges = input.edges();
    auto t = compute_degree_tables(n, edges);

    std::vector<std::vector<NodeId>> unitigs;
    unitigs.reserve(n);
    walk_linear_chains(n, edges, t, unitigs, result.old_to_new);

    merge_closed_twins(n, input, edges, unitigs, result.old_to_new);
    merge_isolated_unitigs(n, input, t, unitigs, result.old_to_new);

    // Build the compacted graph with consensus.
    LosslessGraph& out = result.compacted;
    for (const auto& members : unitigs) {
        const std::uint32_t len32 = pick_unitig_length(members, input);
        const auto& start_node = input.node(members.front());
        NodeId new_id = out.add_node(len32, start_node.copy_count);

        // RC:i on the unitig S-line is the sum of member-node read_support.
        // Always summed, regardless of how pick_unitig_length chose the length
        // (chain-sum vs representative-for-twins) — coverage is additive.
        std::uint64_t total_support = 0;
        for (NodeId m : members) total_support += input.node(m).read_support;
        out.node(new_id).read_support =
            (total_support > std::numeric_limits<std::uint32_t>::max())
                ? std::numeric_limits<std::uint32_t>::max()
                : static_cast<std::uint32_t>(total_support);

        std::string consensus = build_unitig_consensus(members, input);
        if (!consensus.empty()) {
            out.node(new_id).consensus = std::move(consensus);
        }
    }

    // Emit inter-unitig edges (skip edges touching dropped/isolated nodes
    // — they have old_to_new == kUnassigned, which can only happen when
    // both endpoints are isolated, but defensive-guard anyway).
    for (const auto& e : edges) {
        if (e.from >= n || e.to >= n) continue;
        NodeId uf = result.old_to_new[e.from];
        NodeId ut = result.old_to_new[e.to];
        if (uf == kUnassigned || ut == kUnassigned) continue;
        if (uf == ut) continue;
        out.add_edge(uf, ut, e.read_support);
    }

    return result;
}

CompactionResult compact_unitigs_with_sequences(
    const LosslessGraph& input,
    const std::vector<std::string>& node_sequences) {

    CompactionResult result;
    const std::size_t n = input.node_count();
    result.old_to_new.assign(n, kUnassigned);

    if (n == 0) {
        return result;
    }

    const auto& edges = input.edges();
    auto t = compute_degree_tables(n, edges);

    std::vector<std::vector<NodeId>> unitigs;
    unitigs.reserve(n);
    walk_linear_chains(n, edges, t, unitigs, result.old_to_new);

    merge_closed_twins(n, input, edges, unitigs, result.old_to_new);
    merge_isolated_unitigs(n, input, t, unitigs, result.old_to_new);

    // ---- Build compacted graph with consensus from sequences ----
    // Split into three sub-passes so the expensive per-unitig consensus
    // computation can run in parallel:
    //   (a) Sequential: allocate a new_id per unitig and gather its member
    //       sequences. add_node mutates `out` and must stay single-threaded.
    //   (b) Parallel: compute simple_majority_consensus(seqs[i]) for every
    //       unitig. Pure per-unitig, no shared state.
    //   (c) Sequential: assign each consensus into out.node(new_ids[i]).
    LosslessGraph& out = result.compacted;
    const std::size_t u_count = unitigs.size();

    std::vector<NodeId> new_ids(u_count);
    std::vector<std::vector<std::string>> seqs_per_unitig(u_count);

    // (a)
    for (std::size_t i = 0; i < u_count; ++i) {
        const auto& members = unitigs[i];
        const std::uint32_t len32 = pick_unitig_length(members, input);
        std::uint64_t total_support = 0;
        for (NodeId m : members) total_support += input.node(m).read_support;
        const auto& start_node = input.node(members.front());
        new_ids[i] = out.add_node(len32, start_node.copy_count);
        out.node(new_ids[i]).read_support =
            (total_support > std::numeric_limits<std::uint32_t>::max())
                ? std::numeric_limits<std::uint32_t>::max()
                : static_cast<std::uint32_t>(total_support);

        auto& seqs = seqs_per_unitig[i];
        seqs.reserve(members.size());
        for (NodeId m : members) {
            if (m < node_sequences.size() && !node_sequences[m].empty()) {
                seqs.push_back(node_sequences[m]);
            }
        }
    }

    // (b) — the expensive bit, now parallel
    std::vector<std::string> consensuses(u_count);
    #pragma omp parallel for schedule(dynamic, 16)
    for (std::size_t i = 0; i < u_count; ++i) {
        if (!seqs_per_unitig[i].empty()) {
            consensuses[i] = simple_majority_consensus(seqs_per_unitig[i]);
        }
    }

    // (c)
    std::size_t empty_consensus_unitigs = 0;
    for (std::size_t i = 0; i < u_count; ++i) {
        if (!consensuses[i].empty()) {
            out.node(new_ids[i]).consensus = std::move(consensuses[i]);
        } else {
            // Finding 4 cascade: empty consensus here means the downstream
            // write_bed_with_refs path cannot FASTA-dump this node, so
            // minimap2 never sees it and the BED row falls back to chrom=NA.
            // Upstream fix is tracked under Finding 1 (unitig collapse);
            // see feat/compactor-debug. Count-and-warn until then.
            ++empty_consensus_unitigs;
        }
    }
    if (empty_consensus_unitigs > 0) {
        std::cerr << "BRANCH: graph_compactor produced "
                  << empty_consensus_unitigs
                  << " unitigs with empty consensus (no member sequences "
                     "supplied). Downstream BED rows will have chrom=NA for "
                     "these. Root cause tracked in Finding 1 (unitig "
                     "collapse, feat/compactor-debug).\n";
    }

    // Emit inter-unitig edges
    for (const auto& e : edges) {
        if (e.from >= n || e.to >= n) continue;
        NodeId uf = result.old_to_new[e.from];
        NodeId ut = result.old_to_new[e.to];
        if (uf == kUnassigned || ut == kUnassigned) continue;
        if (uf == ut) continue;
        out.add_edge(uf, ut, e.read_support);
    }

    return result;
}

}  // namespace branch::graph
