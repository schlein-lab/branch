// BRANCH v0.5 — Phase 5 implementation.
//
// Orphan re-clustering by MinHash sketches + union-find. Reads share
// a cluster if their sketches overlap by ≥ profile_.sketch_jaccard_lo
// fraction. Cluster consensus hint = the longest member sequence
// (placeholder for proper POA — sufficient for downstream visualisation
// and "give me a representative" use cases).

#include "wg/orphan_clusterer.hpp"
#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <cstdint>
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

// Union-find.
struct UF {
    std::vector<std::uint32_t> p;
    explicit UF(std::size_t n) : p(n) {
        for (std::size_t i = 0; i < n; ++i) p[i] = static_cast<std::uint32_t>(i);
    }
    std::uint32_t find(std::uint32_t x) {
        while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
        return x;
    }
    void unite(std::uint32_t a, std::uint32_t b) {
        a = find(a); b = find(b);
        if (a != b) p[a] = b;
    }
};

}  // namespace

OrphanClusterer::OrphanClusterer(const ::branch::graph::TechProfile& profile)
    : profile_(profile) {}

void OrphanClusterer::cluster(const std::vector<OrphanRead>& orphans,
                              std::vector<OrphanCluster>& out) const {
    out.clear();
    if (orphans.empty()) return;

    const std::size_t k = profile_.minimizer_k;
    const std::size_t min_share = std::max<std::size_t>(
        4, static_cast<std::size_t>(kSketchSize * profile_.sketch_jaccard_lo));

    std::vector<std::vector<std::uint64_t>> sketches(orphans.size());
    for (std::size_t i = 0; i < orphans.size(); ++i)
        sketches[i] = sketch_seq(orphans[i].seq, k);

    UF uf(orphans.size());
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> idx;
    idx.reserve(orphans.size() * kSketchSize);
    for (std::uint32_t i = 0; i < sketches.size(); ++i)
        for (std::uint64_t h : sketches[i]) idx[h].push_back(i);

    for (std::uint32_t i = 0; i < sketches.size(); ++i) {
        if (sketches[i].empty()) continue;
        std::unordered_map<std::uint32_t, std::uint16_t> co;
        for (std::uint64_t h : sketches[i]) {
            auto it = idx.find(h);
            if (it == idx.end()) continue;
            for (std::uint32_t j : it->second) {
                if (j <= i) continue;
                ++co[j];
            }
        }
        for (auto& [j, s] : co) {
            if (s >= min_share) uf.unite(i, j);
        }
    }

    std::unordered_map<std::uint32_t, std::vector<std::uint32_t>> groups;
    for (std::uint32_t i = 0; i < orphans.size(); ++i) {
        groups[uf.find(i)].push_back(i);
    }

    std::uint32_t cid = 0;
    for (auto& [_, members] : groups) {
        OrphanCluster c{};
        c.cluster_id = cid++;
        // Pick longest member as consensus hint.
        std::size_t longest_idx = members[0];
        for (std::uint32_t m : members) {
            if (orphans[m].seq.size() > orphans[longest_idx].seq.size()) longest_idx = m;
        }
        c.consensus_hint = orphans[longest_idx].seq;
        c.ref_attempt_chrom = "unmapped";
        c.ref_attempt_score = 0.0;
        c.member_read_ids.reserve(members.size());
        for (std::uint32_t m : members) c.member_read_ids.push_back(orphans[m].read_id);
        out.push_back(std::move(c));
    }
}

}  // namespace branch::wg
