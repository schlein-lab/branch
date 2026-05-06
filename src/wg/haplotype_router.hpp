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
    /// If `lazy_load_counter` is true (the default has changed to true),
    /// the counter_ map is NOT loaded eagerly — it is only loaded if a
    /// caller invokes `find_het_pairs()` (CPU fallback). The GPU het-pair
    /// path skips counter_ entirely (fed via install_het_pairs_from_pairs)
    /// so the unordered_map's ~48 B/entry overhead — ~95 GB for 2 Gbil
    /// recurrent k-mers — is avoided. v08 OOMed at 144 GB MaxRSS exactly
    /// because of this eager load.
    HaplotypeRouter(const ::branch::graph::TechProfile& profile,
                    const std::string& kmer_hist_path,
                    int expected_coverage,
                    bool lazy_load_counter = true);

    /// Build the het-pair table from the loaded counter. `low_frac` and
    /// `high_frac` define the per-allele count window relative to the
    /// expected single-copy coverage; `sum_frac_lo`/`sum_frac_hi` define
    /// how tight the combined-count constraint is.
    void find_het_pairs(double low_frac = 0.30, double high_frac = 0.70,
                        double sum_frac_lo = 0.70, double sum_frac_hi = 1.40);

    /// Populate the het-pair table directly from a list of pre-computed
    /// (kmer_a, kmer_b) tuples — used by the Phase 1.1 GPU launcher to
    /// install the GPU-built het table without re-running the CPU loop.
    /// Each input tuple becomes two `het_` entries (one per side, with
    /// allele label 0 = lex-smaller bits).
    void install_het_pairs_from_pairs(
        const std::vector<std::uint64_t>& a,
        const std::vector<std::uint64_t>& b);

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

    /// Export het_pair allele tuples in canonical (kmer_a, kmer_b) form,
    /// where allele 0 maps to kmer_a and allele 1 to kmer_b. Used by the
    /// v0.6 phaser to seed its phase_engine signature index.
    [[nodiscard]] std::vector<std::pair<std::uint64_t, std::uint64_t>>
        export_het_pair_kmers() const;

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

    // Path retained for lazy loading on first find_het_pairs() call.
    std::string counter_path_;
    bool counter_loaded_ = false;
};

}  // namespace branch::wg
