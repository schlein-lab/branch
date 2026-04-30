// BRANCH v0.5 — Phase 0: k-mer histogram + bimodality detection.
//
// Reads the input FASTQ stream twice. First pass: bloom filter
// records "seen at least once". Second pass: only k-mers that the
// bloom says recurred get added to the exact counter. Final step:
// histogram of counts → bimodality detection identifies expected_coverage
// (single-copy peak) and expected_coverage/2 (het-allele peak), which
// drives Phase 1 SNP discovery.

#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "graph/kmer_sketch.hpp"
#include "wg/bloom_filter.hpp"
#include "wg/kmer_counter.hpp"

namespace branch::wg {

/// Result of Phase 0.
struct KmerHistResult {
    /// Count distribution: histogram[i] = number of k-mers that occurred i times.
    /// histogram.back() is the sum for everything ≥ histogram.size()-1.
    std::vector<std::uint64_t> histogram;

    /// Detected single-copy peak (= expected coverage).
    /// 0 if no peak detected (input too small).
    int detected_coverage = 0;

    /// Detected het-allele peak (≈ expected_coverage/2).
    /// 0 if no second peak (haploid input or mono-allelic).
    int detected_het_coverage = 0;

    /// Z-score separation of the two peaks. Higher = cleaner bimodality.
    double bimodality_z = 0.0;

    /// Total number of distinct k-mers that occurred ≥ 2 times.
    std::uint64_t n_recurrent_kmers = 0;

    /// Total number of input reads + bp processed.
    std::uint64_t n_reads = 0;
    std::uint64_t n_bases = 0;
};

class KmerHistBuilder {
public:
    /// `profile` provides minimizer_k. `expected_n_distinct_kmers` is
    /// used to size the bloom filter; over-estimate is cheap (more
    /// bits) under-estimate is expensive (high FP rate). Default
    /// 30 G is sized for a human whole-genome ONT input at 30× coverage.
    explicit KmerHistBuilder(const ::branch::graph::TechProfile& profile,
                             std::size_t expected_n_distinct_kmers = 30'000'000'000ULL);

    /// Pass 1: feed every k-mer into the bloom filter.
    void pass1_add_read(std::string_view seq) noexcept;

    /// Pass 2: feed every k-mer; if bloom says "maybe seen", increment
    /// exact counter. Must come after all pass1_add_read calls and a
    /// call to switch_to_pass2().
    void pass2_add_read(std::string_view seq) noexcept;

    /// Transition between phases. Internally just freezes the bloom
    /// filter (no further adds) and prepares the counter.
    void switch_to_pass2();

    /// Run final histogram + bimodality detection. Caller hands in
    /// `expected_coverage` from the TechProfile or CLI override; the
    /// detector confirms or corrects it.
    [[nodiscard]] KmerHistResult finalize(int expected_coverage,
                                          double bimodality_z_threshold);

    /// Write counter table + histogram + result block to a binary
    /// intermediate file. Phase 1 reads this back.
    void write_to(const std::string& path,
                  const KmerHistResult& result) const;

    /// Diagnostics.
    [[nodiscard]] std::size_t bloom_bytes() const noexcept { return bloom_.bytes(); }
    [[nodiscard]] std::size_t counter_bytes() const noexcept { return counter_.bytes(); }
    [[nodiscard]] double bloom_saturation() const noexcept { return bloom_.saturation(); }

private:
    const ::branch::graph::TechProfile profile_;
    BloomFilter bloom_;
    KmerCounter counter_;
    bool in_pass2_ = false;
    std::uint64_t n_reads_ = 0;
    std::uint64_t n_bases_ = 0;

    void scan_kmers_(std::string_view seq, bool pass2) noexcept;
};

/// Compute canonical k-mer hash (forward / reverse complement, smaller).
[[nodiscard]] std::uint64_t canonical_kmer_hash(std::uint64_t kmer, std::size_t k) noexcept;

}  // namespace branch::wg
