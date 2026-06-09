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
#include <unordered_map>
#include <vector>

#include "graph/kmer_sketch.hpp"
#include "wg/bloom_filter.hpp"
#include "wg/external_kmer_counter.hpp"

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
    /// `max_ram_bytes` caps the Phase-0 second-pass occurrence buffer (the
    /// former exact counter was the whole-genome RAM wall); beyond it, runs
    /// spill to `scratch_dir` (BeeGFS) and are merged exactly at finalize.
    explicit KmerHistBuilder(const ::branch::graph::TechProfile& profile,
                             std::size_t expected_n_distinct_kmers = 30'000'000'000ULL,
                             std::size_t max_ram_bytes = 16ULL << 30,
                             std::string scratch_dir = "/tmp");

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
    [[nodiscard]] std::size_t counter_bytes() const noexcept { return ext_counter_.peak_buffer_bytes(); }
    [[nodiscard]] double bloom_saturation() const noexcept { return bloom_.saturation(); }

private:
    const ::branch::graph::TechProfile profile_;
    BloomFilter bloom_;
    ExternalKmerCounter ext_counter_;
    // Sorted (key, count) arrays materialized by finalize(); reused by write_to().
    std::vector<std::uint64_t> keys_;
    std::vector<std::uint32_t> counts_;
    bool in_pass2_ = false;
    std::uint64_t n_reads_ = 0;
    std::uint64_t n_bases_ = 0;

    void scan_kmers_(std::string_view seq, bool pass2) noexcept;
};

/// Compute canonical k-mer hash (forward / reverse complement, smaller hash).
/// Used by tests + bloom filter; not used as counter key (see canonical_kmer_bits).
[[nodiscard]] std::uint64_t canonical_kmer_hash(std::uint64_t kmer, std::size_t k) noexcept;

/// Canonical k-mer bits = min(forward_bits, rc_bits). Unique per canonical
/// k-mer (no hash collision) so it's safe to use directly as the counter
/// key. Phase 1.1 reads the bits straight back from the dumped counter
/// table to enumerate 1-bp variants for het-pair detection.
[[nodiscard]] std::uint64_t canonical_kmer_bits(std::uint64_t kmer, std::size_t k) noexcept;

/// Inverse of the bit-packing used in scan_kmers_: rebuild the ATCG sequence
/// from k-mer bits. ('A','C','G','T' for codes 0..3.)
void kmer_bits_to_seq(std::uint64_t bits, std::size_t k, char* out) noexcept;

/// Load the (canonical_kmer_bits, count) table from a Phase 0 .bin dump.
/// Used by Phase 3 Stage 1 (cheap whole-genome screen) without going
/// through the HaplotypeRouter — a lightweight reader.
void load_kmer_counter_dump(
    const std::string& path,
    std::unordered_map<std::uint64_t, std::uint32_t>& out_counter);

/// Same as load_kmer_counter_dump but loads directly into two parallel
/// sorted vectors (key, count) — bypasses the std::unordered_map detour
/// that costs ~48 B/entry of bucket overhead. For ~2 Gbil recurrent
/// k-mers (whole-genome HiFi) this saves ~70 GB of host RAM and is the
/// only viable load path under the 140 GB SLURM cap.
///
/// Output arrays are sorted ascending by key. Caller can then binary-
/// search them or pass directly to the GPU launchers (which already
/// expect sorted parallel arrays).
void load_kmer_counter_sorted_arrays(
    const std::string& path,
    std::vector<std::uint64_t>& out_keys,
    std::vector<std::uint32_t>& out_counts);

/// Build a KmerHistResult (histogram + bimodality) from (keys, counts)
/// pairs as returned by the Phase 0 GPU launcher. Mirrors the bimodality
/// logic of KmerHistBuilder::finalize.
[[nodiscard]] KmerHistResult finalize_from_keys_counts(
    int expected_coverage,
    double bimodality_z_threshold,
    const std::vector<std::uint64_t>& keys,
    const std::vector<std::uint32_t>& counts,
    std::uint64_t n_reads,
    std::uint64_t n_bases);

/// Write the same .bin dump format that KmerHistBuilder::write_to emits,
/// but from (keys, counts) pairs instead of a builder. Used by the
/// Phase 0 GPU launcher to persist the counter for Phase 1.1 + Stage 1.
void write_phase0_dump(
    const std::string& path,
    std::size_t k,
    const KmerHistResult& result,
    const std::vector<std::uint64_t>& keys,
    const std::vector<std::uint32_t>& counts);

}  // namespace branch::wg
