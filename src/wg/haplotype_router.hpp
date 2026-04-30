// BRANCH v0.5 — Phase 1.1: read → haplotype routing.
//
// Loads the Phase 0 (canonical_kmer_bits, count) table, finds het-SNP
// candidate pairs (1-bp-edit-distance k-mers, both with count near
// expected_coverage/2 and combined count near expected_coverage), assigns
// a deterministic local allele label per pair, then runs an anchor pass
// (highest-het-density read seeds hap1 allele pattern) and a voting pass
// (each read → match/mismatch count → hap1, hap2, or ambig).

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "graph/kmer_sketch.hpp"

namespace branch::wg {

enum class Haplotype : std::uint8_t { Hap1 = 0, Hap2 = 1, Ambig = 2 };

struct ReadRouting {
    Haplotype hap = Haplotype::Ambig;
    std::uint32_t votes_hap1 = 0;  // matches against hap1 allele pattern
    std::uint32_t votes_hap2 = 0;  // mismatches (= matches against hap2)
    [[nodiscard]] double confidence() const noexcept;
};

class HaplotypeRouter {
public:
    /// `expected_coverage` is the single-copy peak (typically the
    /// detected_coverage from Phase 0, or user-overridden via CLI).
    HaplotypeRouter(const ::branch::graph::TechProfile& profile,
                    const std::string& kmer_hist_path,
                    int expected_coverage);

    /// Build the het-pair table from the loaded counter. `low_frac` and
    /// `high_frac` define the per-allele count window relative to the
    /// expected single-copy coverage; `sum_frac_lo`/`sum_frac_hi` define
    /// how tight the combined-count constraint is.
    void find_het_pairs(double low_frac = 0.30, double high_frac = 0.70,
                        double sum_frac_lo = 0.70, double sum_frac_hi = 1.40);

    /// Anchor pass: scan a sample of reads, pick the one with the highest
    /// het-pair hit count, and lock its observed allele pattern as hap1.
    /// Empty `reads` is a no-op (route_read will then return Ambig for
    /// every input).
    void set_anchor_from_reads(
        const std::vector<std::pair<std::string, std::string>>& reads,
        std::size_t max_scan = 200'000);

    /// Route one read.
    [[nodiscard]] ReadRouting route_read(std::string_view seq) const noexcept;

    /// Bulk write a TSV with one row per read: read_id, hap, votes_hap1,
    /// votes_hap2, confidence. Used downstream for the #disp block.
    void write_assignments_tsv(
        const std::vector<std::pair<std::string, std::string>>& reads,
        const std::string& out_path,
        double min_confidence = 0.6) const;

    /// Diagnostics.
    [[nodiscard]] std::size_t n_recurrent_kmers() const noexcept { return counter_.size(); }
    [[nodiscard]] std::size_t n_het_pairs() const noexcept { return n_het_pairs_; }
    [[nodiscard]] bool has_anchor() const noexcept { return !hap1_pattern_.empty(); }
    [[nodiscard]] std::size_t bytes() const noexcept;

private:
    ::branch::graph::TechProfile profile_;
    int expected_coverage_;

    /// Recurrent k-mer counter loaded from Phase 0 output. Key =
    /// canonical kmer-bits (decodable to ATCG); value = total count.
    std::unordered_map<std::uint64_t, std::uint32_t> counter_;

    /// Het-pair table: canonical_bits → (pair_id, allele_label).
    /// Pair_id is a small dense integer; allele_label is 0 or 1, picked
    /// deterministically (the lex-smaller bits get label 0).
    struct HetEntry {
        std::uint32_t pair_id;
        std::uint8_t  allele;     // 0 or 1
    };
    std::unordered_map<std::uint64_t, HetEntry> het_;
    std::size_t n_het_pairs_ = 0;

    /// hap1_pattern_[pair_id] = expected allele label (0 or 1) on hap1.
    /// Filled by set_anchor_from_reads.
    std::unordered_map<std::uint32_t, std::uint8_t> hap1_pattern_;

    void load_counter_(const std::string& path);
};

}  // namespace branch::wg
