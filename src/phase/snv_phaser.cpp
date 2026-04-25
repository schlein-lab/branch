// BRANCH — SNV phaser implementation.

#include "phase/snv_phaser.hpp"

#include <algorithm>
#include <unordered_set>
#include <utility>

#include "graph/read_db.hpp"

namespace branch::phase {

namespace {

// Map (unitig, pos) -> aggregator {ref_count, alt_count, base bag}.
struct PosKey {
    branch::graph::NodeId unitig;
    std::uint32_t pos;
    bool operator==(const PosKey& o) const noexcept {
        return unitig == o.unitig && pos == o.pos;
    }
};
struct PosKeyHash {
    std::size_t operator()(const PosKey& k) const noexcept {
        return (static_cast<std::size_t>(k.unitig) << 32) ^ k.pos;
    }
};

struct AlleleBag {
    // base -> reads carrying that base at the position
    std::unordered_map<char, std::vector<branch::graph::ReadId>> reads_by_base;
};

}  // namespace

std::vector<InformativeSNV>
find_informative_snvs(const branch::graph::ReadDB& db,
                      const SnvPhaserConfig& cfg) {
    // First pass: aggregate alt-reads per position from per-read deltas.
    std::unordered_map<PosKey, AlleleBag, PosKeyHash> agg;
    for (const auto& e : db.all()) {
        if (e.bucket != branch::graph::ReadBucket::Mapped) continue;
        if (e.path.empty()) continue;
        const auto& head = e.path.front();
        for (const auto& d : e.deltas) {
            if (d.op != 'S') continue;
            PosKey k{head.unitig, head.start + d.pos};
            agg[k].reads_by_base[d.base].push_back(e.id);
        }
    }

    // Second pass: for each candidate position, also collect the reads
    // that *span* it via their path step but carry NO delta there —
    // those are the reference-allele observations we'd otherwise miss.
    std::unordered_map<PosKey, std::vector<branch::graph::ReadId>, PosKeyHash>
        ref_reads_by_pos;
    for (const auto& [key, bag] : agg) {
        std::unordered_set<branch::graph::ReadId> alt_carriers;
        for (const auto& [_, reads] : bag.reads_by_base) {
            alt_carriers.insert(reads.begin(), reads.end());
        }
        for (const auto& e : db.all()) {
            if (e.bucket != branch::graph::ReadBucket::Mapped) continue;
            if (alt_carriers.count(e.id)) continue;
            // Does any path step of e cover this position?
            for (const auto& s : e.path) {
                if (s.unitig == key.unitig &&
                    s.start <= key.pos && s.end > key.pos) {
                    ref_reads_by_pos[key].push_back(e.id);
                    break;
                }
            }
        }
    }

    std::vector<InformativeSNV> out;
    for (auto& [key, bag] : agg) {
        const auto ref_it = ref_reads_by_pos.find(key);
        const std::size_t ref_count = (ref_it == ref_reads_by_pos.end()) ? 0 : ref_it->second.size();
        if (ref_count < cfg.min_ref_reads) continue;

        for (auto& [base, reads] : bag.reads_by_base) {
            if (reads.size() < cfg.min_alt_reads) continue;
            InformativeSNV snv;
            snv.unitig = key.unitig;
            snv.pos = key.pos;
            snv.ref_base = '.';   // base of the consensus; would come
                                   // from graph.node(unitig).consensus[pos]
                                   // once that wiring is in.
            snv.alt_base = base;
            snv.alt_reads = std::move(reads);
            if (ref_it != ref_reads_by_pos.end()) snv.ref_reads = ref_it->second;
            out.push_back(std::move(snv));
        }
    }

    std::sort(out.begin(), out.end(),
              [](const InformativeSNV& a, const InformativeSNV& b) {
                  if (a.unitig != b.unitig) return a.unitig < b.unitig;
                  return a.pos < b.pos;
              });
    return out;
}

namespace {

// Count reads that carry consistent alleles at both SNVs (either both
// ref or both alt). When this count is >= min_link_reads, the two SNVs
// belong to the same phase block.
std::size_t count_linked_reads(const InformativeSNV& a,
                                const InformativeSNV& b) {
    std::unordered_set<branch::graph::ReadId> a_ref(a.ref_reads.begin(), a.ref_reads.end());
    std::unordered_set<branch::graph::ReadId> a_alt(a.alt_reads.begin(), a.alt_reads.end());
    std::size_t linked = 0;
    for (auto rid : b.ref_reads) {
        if (a_ref.count(rid)) ++linked;
    }
    for (auto rid : b.alt_reads) {
        if (a_alt.count(rid)) ++linked;
    }
    return linked;
}

}  // namespace

std::vector<PhaseBlock>
build_phase_blocks(const std::vector<InformativeSNV>& snvs,
                   const SnvPhaserConfig& cfg) {
    std::vector<PhaseBlock> blocks;
    if (snvs.empty()) return blocks;

    PhaseBlock current;
    current.snvs.push_back(snvs.front());

    for (std::size_t i = 1; i < snvs.size(); ++i) {
        const auto& prev = current.snvs.back();
        const auto& cur = snvs[i];
        // Different unitig => start new block.
        // Otherwise, link only if consistent reads cross both SNVs.
        const bool same_unitig = (prev.unitig == cur.unitig);
        const std::size_t links = same_unitig ? count_linked_reads(prev, cur) : 0;
        if (same_unitig && links >= cfg.min_link_reads) {
            current.snvs.push_back(cur);
        } else {
            blocks.push_back(std::move(current));
            current = PhaseBlock{};
            current.snvs.push_back(cur);
        }
    }
    if (!current.snvs.empty()) blocks.push_back(std::move(current));

    // Within each block, assign reads to haplotypes by majority vote
    // over their per-SNV calls. Hap 0 = ref of first SNV, Hap 1 = alt.
    for (auto& block : blocks) {
        const auto& anchor = block.snvs.front();
        std::unordered_set<branch::graph::ReadId> hap0(anchor.ref_reads.begin(),
                                                         anchor.ref_reads.end());
        std::unordered_set<branch::graph::ReadId> hap1(anchor.alt_reads.begin(),
                                                         anchor.alt_reads.end());
        std::unordered_map<branch::graph::ReadId, std::pair<int,int>> votes;
        for (const auto& snv : block.snvs) {
            for (auto rid : snv.ref_reads) {
                if (hap0.count(rid)) ++votes[rid].first;
                else if (hap1.count(rid)) ++votes[rid].second;
            }
            for (auto rid : snv.alt_reads) {
                if (hap1.count(rid)) ++votes[rid].first;  // consistent with hap1
                else if (hap0.count(rid)) ++votes[rid].second;  // anti-vote
            }
        }
        for (const auto& [rid, vc] : votes) {
            const auto [pro, anti] = vc;
            if (pro > anti) {
                block.read_haplotype[rid] = (hap0.count(rid) ? 0 : 1);
            } else if (anti > pro) {
                block.read_haplotype[rid] = (hap0.count(rid) ? 1 : 0);
            } else {
                block.read_haplotype[rid] = 255;
            }
        }
    }
    return blocks;
}

}  // namespace branch::phase
