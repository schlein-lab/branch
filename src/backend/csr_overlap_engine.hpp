#pragma once

// BRANCH — CSR overlap engine (WGS-scale replacement for the
// unordered_map-based path in cpu_backend.cpp).
//
// Motivation (measured on a 768 GB SLURM node, 2026-08-08): the legacy
// compute_overlaps_impl holds three full-batch structures at once —
// all_hits, buckets (unordered_map<hash, vector<BucketEntry>>) and
// pair_matches (unordered_map<PairKey, vector<MatchInfo>>). On 5k dense
// HiFi reads that is a 37.7 GB peak; projected ~2.9 TB at 30x WGS.
//
// This engine follows the hifiasm layout instead:
//   * two-pass CSR index: count occurrences per hash into an
//     open-addressing table, prefix-sum into one flat entry array,
//     then scatter. One allocation for all positions, no per-key
//     vectors, no rehash spikes.
//   * per-query streaming match: each query read gathers its candidate
//     hits into a thread-local scratch vector, sorts it, and runs the
//     exact legacy clustering logic per candidate pair. Nothing about
//     the pair survives the query iteration except emitted OverlapPairs.
//   * seeding cap: hashes occurring more often than `seed_cap` are not
//     indexed. This caps the O(bucket²) candidate blowup on repeats.
//     LOSSLESS NOTE: the cap only affects candidate SEEDING — reads
//     remain anchored through their other minimizers, and capped hashes
//     are reported (CsrIndexStats) instead of dropped silently. Regions
//     whose every minimizer is capped surface in the unresolved report
//     rather than vanishing.
//
// Determinism: table fill order is the (fixed) read order, the entry
// array is laid out by table slot order, per-query scratch is sorted,
// and the final output is canonically sorted by (read_a, read_b,
// offset_a). Results are bit-identical for any thread count.

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "graph/kmer_sketch.hpp"

namespace branch::backend {

// Aggregate statistics from CSR index construction. `capped_hashes`
// holds up to `kMaxReportedCapped` of the largest capped buckets for
// the overcap report (TSV written by the CLI layer).
struct CsrIndexStats {
    static constexpr std::size_t kMaxReportedCapped = 20;

    std::uint64_t n_reads{0};
    std::uint64_t n_hits_total{0};       // minimizers sketched over all reads
    std::uint64_t n_distinct_hashes{0};
    std::uint64_t n_entries_indexed{0};  // hits that made it into the index
    std::uint64_t n_hashes_capped{0};    // distinct hashes over seed_cap
    std::uint64_t n_entries_dropped{0};  // hits belonging to capped hashes

    struct CappedHash {
        std::uint64_t hash;
        std::uint64_t count;
    };
    std::vector<CappedHash> capped_hashes;  // largest first, <= kMaxReportedCapped
};

// One index entry, 8 bytes: read (32 bit) | pos (24 bit) | strand (8 bit).
struct CsrEntry {
    std::uint64_t packed;

    [[nodiscard]] std::uint32_t read_id() const noexcept {
        return static_cast<std::uint32_t>(packed >> 32);
    }
    [[nodiscard]] std::uint32_t pos() const noexcept {
        return static_cast<std::uint32_t>((packed >> 8) & 0xFFFFFFULL);
    }
    [[nodiscard]] std::uint8_t strand() const noexcept {
        return static_cast<std::uint8_t>(packed & 0xFFULL);
    }

    static CsrEntry make(std::uint32_t read_id, std::uint32_t pos,
                         std::uint8_t strand) noexcept {
        return CsrEntry{(static_cast<std::uint64_t>(read_id) << 32) |
                        (static_cast<std::uint64_t>(pos & 0xFFFFFFu) << 8) |
                        strand};
    }
};

static_assert(sizeof(CsrEntry) == 8, "CsrEntry must stay 8 bytes");

// Open-addressing hash table (power-of-two capacity, linear probing)
// mapping minimizer hash -> [count during build, entry offset after
// finalise]. Deliberately hand-rolled: std::unordered_map's node-based
// layout is what we are escaping from.
class CsrMinimizerIndex {
public:
    // Build from a read batch. Sketching runs on `threads` workers
    // (chunked, deterministic order); table insert and entry scatter
    // are sequential consumers, so the layout is thread-independent.
    // seed_cap == 0 means "no cap".
    void build(const ReadBatch& batch,
               const graph::MinimizerProfile& profile,
               std::size_t seed_cap,
               unsigned int threads,
               CsrIndexStats* stats);

    // Lookup all indexed entries for `hash`. Empty span when the hash
    // is unknown or was capped.
    [[nodiscard]] std::span<const CsrEntry> lookup(std::uint64_t hash) const noexcept;

    [[nodiscard]] std::size_t entry_count() const noexcept { return entries_.size(); }
    [[nodiscard]] std::size_t distinct_hashes() const noexcept { return n_keys_; }

private:
    struct Slot {
        std::uint64_t key;
        std::uint32_t count;   // occurrences (build phase)
        std::uint32_t offset;  // start into entries_ (finalised), or kCapped
    };
    static constexpr std::uint64_t kEmptyKey = 0xFFFFFFFFFFFFFFFFULL;
    static constexpr std::uint32_t kCapped = 0xFFFFFFFFu;

    [[nodiscard]] std::size_t find_slot(std::uint64_t key) const noexcept;
    std::size_t find_or_insert_slot(std::uint64_t key) noexcept;

    std::vector<Slot> slots_;
    std::vector<CsrEntry> entries_;
    std::size_t n_keys_{0};
    std::size_t mask_{0};
};

// Tuning knobs for the CSR overlap pass. Field semantics and defaults
// match the legacy OverlapConfig in cpu_backend.cpp exactly, so both
// engines produce identical pairs when seed_cap == legacy bucket cap.
struct CsrOverlapConfig {
    std::size_t min_matches{5};
    std::int32_t offset_tolerance{100};
    float strand_majority{0.80f};
    unsigned int threads{1};
    std::size_t seed_cap{500};  // HiFi WGS default; 0 = uncapped
};

// Compute overlaps for the batch using a freshly built CSR index.
// Appends canonical-sorted OverlapPairs to `out`. Returns the number
// of pairs appended. `stats` (optional) receives index statistics.
std::size_t compute_overlaps_csr(const ReadBatch& batch,
                                 const CsrOverlapConfig& cfg,
                                 std::vector<OverlapPair>& out,
                                 CsrIndexStats* stats = nullptr);

}  // namespace branch::backend
