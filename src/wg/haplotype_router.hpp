// BRANCH v0.5 — Phase 1.1: read → haplotype routing.
//
// Loads the Phase 0 k-mer counter, derives het-SNP candidate pairs
// (1-bp-edit-distance k-mers with counts each ≈ expected_coverage/2
// summing to ≈ expected_coverage), then assigns each read to hap1,
// hap2, or ambig by majority vote across the het-SNP k-mers it
// contains.

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
    std::uint32_t votes_hap1 = 0;
    std::uint32_t votes_hap2 = 0;
    /// Confidence = max(v1, v2) / (v1 + v2 + 1).
    /// Below this threshold → Ambig regardless of majority.
    [[nodiscard]] double confidence() const noexcept;
};

class HaplotypeRouter {
public:
    /// `kmer_hist_path` is the binary Phase 0 output. `expected_coverage`
    /// is the single-copy peak (typically detected_coverage from Phase 0
    /// result, or user-overridden via CLI).
    HaplotypeRouter(const ::branch::graph::TechProfile& profile,
                    const std::string& kmer_hist_path,
                    int expected_coverage);

    /// Builds the het-pair table. Must be called before route_read().
    /// `low_frac` and `high_frac` define the window
    /// [expected_coverage * low_frac, expected_coverage * high_frac]
    /// that het-allele k-mers fall in. Defaults 0.3 .. 0.7.
    void find_het_pairs(double low_frac = 0.3, double high_frac = 0.7);

    /// Route one read.
    [[nodiscard]] ReadRouting route_read(std::string_view seq) const noexcept;

    /// Write a TSV with one row per read: read_id, hap, votes_hap1, votes_hap2.
    /// Used by Phase 1.2 to bulk-load and split the input into per-hap
    /// FASTQ files.
    void write_assignments_tsv(
        const std::vector<std::pair<std::string, std::string>>& reads,
        const std::string& out_path,
        double min_confidence = 0.6) const;

    /// Diagnostics.
    [[nodiscard]] std::size_t n_het_pairs() const noexcept { return het_pair_count_ / 2; }
    [[nodiscard]] std::size_t bytes() const noexcept;

private:
    ::branch::graph::TechProfile profile_;
    int expected_coverage_;

    /// Recurrent k-mer counter loaded from Phase 0 output.
    /// Memory-efficient layout: parallel arrays of (key, count).
    std::vector<std::uint64_t> ckeys_;
    std::vector<std::uint32_t> ccounts_;

    /// Het-pair table: kmer_hash → (paired_kmer_hash, allele_id 0 or 1).
    /// Size = 2 × n_het_pairs (each pair has two entries pointing to each other).
    struct HetEntry {
        std::uint64_t partner_hash;
        std::uint8_t  allele;     // 0 or 1
    };
    std::unordered_map<std::uint64_t, HetEntry> het_;
    std::size_t het_pair_count_ = 0;

    void load_counter_(const std::string& path);
    [[nodiscard]] std::uint32_t lookup_count_(std::uint64_t hash) const noexcept;
};

}  // namespace branch::wg
