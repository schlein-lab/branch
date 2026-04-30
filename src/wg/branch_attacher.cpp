// BRANCH v0.5 — Phase 2 implementation.
//
// Sketch-based transitive attachment. For an ambig read, "anchor"
// means we found ≥ jaccard_lo * sketch_size shared minimizers with
// either:
//   - one of the curated haplotypes (direct anchor, depth = 1)
//   - a previously-attached branch (transitive anchor via that branch)
//   - another already-anchored ambig read (transitive depth +1)
//
// We bound transitive depth at profile_.max_transitive_depth (HiFi=3,
// ONT=5). Reads that fail to anchor within the depth budget go to the
// orphan pool.

#include "wg/branch_attacher.hpp"
#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
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

// Compact MinHash sketch: collect the smallest N canonical k-mer hashes.
constexpr std::size_t kSketchSize = 64;

std::vector<std::uint64_t> sketch_seq(std::string_view seq, std::size_t k) {
    std::vector<std::uint64_t> kmers;
    if (seq.size() < k) return {};
    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t kmer = 0;
    std::size_t filled = 0;
    kmers.reserve(seq.size());
    for (char c : seq) {
        std::uint8_t b = encode_base(c);
        if (b > 3) { filled = 0; kmer = 0; continue; }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;
        kmers.push_back(canonical_kmer_bits(kmer, k));
    }
    if (kmers.empty()) return {};
    std::sort(kmers.begin(), kmers.end());
    kmers.erase(std::unique(kmers.begin(), kmers.end()), kmers.end());
    if (kmers.size() > kSketchSize) kmers.resize(kSketchSize);
    return kmers;
}

std::size_t shared_count(const std::vector<std::uint64_t>& a,
                         const std::vector<std::uint64_t>& b) {
    std::size_t i = 0, j = 0, n = 0;
    while (i < a.size() && j < b.size()) {
        if (a[i] == b[j]) { ++n; ++i; ++j; }
        else if (a[i] < b[j]) ++i;
        else ++j;
    }
    return n;
}

}  // namespace

BranchAttacher::BranchAttacher(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void BranchAttacher::attach(
    const std::vector<std::string>& hap_seqs,
    const std::vector<BranchCandidate>& candidates,
    const std::vector<std::pair<std::string, std::string>>& ambig_reads,
    std::vector<AttachedBranch>& out_branches,
    std::vector<OrphanRead>& out_orphans) const {
    out_branches.clear();
    out_orphans.clear();

    // 1. Direct branches (curator output).
    std::uint32_t bidx = 0;
    for (const auto& c : candidates) {
        AttachedBranch b{};
        std::ostringstream name;
        name << "branch_" << ++bidx;
        b.branch_name = name.str();
        b.seq = c.alt_seq;
        b.anchor_hap = c.hap_idx;
        b.anchor_pos = c.pos_in_hap;
        b.transitive_depth = 0;
        b.confidence = c.q_weighted_evidence;
        out_branches.push_back(std::move(b));
    }

    if (ambig_reads.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::size_t min_share = std::max<std::size_t>(
        4, static_cast<std::size_t>(kSketchSize * profile_.sketch_jaccard_lo));
    const int max_depth = std::max(1, profile_.max_transitive_depth);

    // 2. Sketch the haplotypes and attached branches.
    std::vector<std::vector<std::uint64_t>> hap_sketches;
    hap_sketches.reserve(hap_seqs.size());
    for (const auto& h : hap_seqs) hap_sketches.push_back(sketch_seq(h, k));

    std::vector<std::vector<std::uint64_t>> branch_sketches;
    branch_sketches.reserve(out_branches.size());
    for (const auto& b : out_branches) branch_sketches.push_back(sketch_seq(b.seq, k));

    // 3. Sketch each ambig read.
    std::vector<std::vector<std::uint64_t>> read_sketches(ambig_reads.size());
    for (std::size_t i = 0; i < ambig_reads.size(); ++i) {
        read_sketches[i] = sketch_seq(ambig_reads[i].second, k);
    }

    // 4. BFS over reads. Initial frontier: any read that shares ≥ min_share
    //    with a haplotype or a branch.
    struct AnchorInfo { std::uint32_t hap; std::uint32_t pos; std::uint32_t depth; };
    std::vector<AnchorInfo> anchor(ambig_reads.size(),
                                   AnchorInfo{0, 0, std::numeric_limits<std::uint32_t>::max()});
    std::queue<std::uint32_t> frontier;

    for (std::uint32_t i = 0; i < ambig_reads.size(); ++i) {
        if (read_sketches[i].empty()) continue;
        std::size_t best_share = 0;
        std::uint32_t best_hap = 0;
        for (std::size_t h = 0; h < hap_sketches.size(); ++h) {
            std::size_t s = shared_count(read_sketches[i], hap_sketches[h]);
            if (s > best_share) { best_share = s; best_hap = static_cast<std::uint32_t>(h); }
        }
        std::uint32_t best_depth = 1;  // direct attach to hap
        std::uint32_t best_pos = 0;
        for (std::size_t b = 0; b < branch_sketches.size(); ++b) {
            std::size_t s = shared_count(read_sketches[i], branch_sketches[b]);
            if (s > best_share) {
                best_share = s;
                best_hap = out_branches[b].anchor_hap;
                best_pos = out_branches[b].anchor_pos;
                best_depth = out_branches[b].transitive_depth + 1;
            }
        }
        if (best_share >= min_share) {
            anchor[i] = AnchorInfo{best_hap, best_pos, best_depth};
            frontier.push(i);
        }
    }

    // 5. Inverted minimizer index for fast read↔read shared-sketch lookup.
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> sketch_index;
    sketch_index.reserve(ambig_reads.size() * kSketchSize);
    for (std::uint32_t i = 0; i < read_sketches.size(); ++i) {
        for (std::uint64_t h : read_sketches[i]) sketch_index[h].push_back(i);
    }

    // 6. BFS to propagate anchors transitively.
    while (!frontier.empty()) {
        std::uint32_t i = frontier.front(); frontier.pop();
        if (anchor[i].depth >= static_cast<std::uint32_t>(max_depth)) continue;

        // Candidates: any read sharing a sketch hash with i.
        std::unordered_map<std::uint32_t, std::uint16_t> co_share;
        for (std::uint64_t h : read_sketches[i]) {
            auto it = sketch_index.find(h);
            if (it == sketch_index.end()) continue;
            for (std::uint32_t j : it->second) {
                if (j == i) continue;
                if (anchor[j].depth != std::numeric_limits<std::uint32_t>::max()) continue;
                ++co_share[j];
            }
        }
        for (auto& [j, s] : co_share) {
            if (s < min_share) continue;
            std::uint32_t new_depth = anchor[i].depth + 1;
            if (new_depth > static_cast<std::uint32_t>(max_depth)) continue;
            anchor[j] = AnchorInfo{anchor[i].hap, anchor[i].pos, new_depth};
            frontier.push(j);
        }
    }

    // 7. Emit AttachedBranch for each anchored ambig read; orphan the rest.
    for (std::uint32_t i = 0; i < ambig_reads.size(); ++i) {
        if (anchor[i].depth == std::numeric_limits<std::uint32_t>::max()) {
            out_orphans.push_back(OrphanRead{ambig_reads[i].first, ambig_reads[i].second});
            continue;
        }
        AttachedBranch b{};
        std::ostringstream name;
        name << "branch_" << ++bidx;
        b.branch_name = name.str();
        b.seq = ambig_reads[i].second;
        b.anchor_hap = anchor[i].hap;
        b.anchor_pos = anchor[i].pos;
        b.transitive_depth = anchor[i].depth;
        // Confidence ≈ 1 / depth (deeper = less certain), capped at 1.
        b.confidence = 1.0 / static_cast<double>(anchor[i].depth + 1);
        out_branches.push_back(std::move(b));
    }
}

}  // namespace branch::wg
