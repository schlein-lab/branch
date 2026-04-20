#pragma once

// BRANCH P1.3 — VAF-band benchmark statistics.
//
// Reusable helpers to bin called bubbles/nodes by the *true* (target)
// variant-allele-frequency they originated from and to compute
// per-band retrieval and accuracy metrics.
//
// Band edges follow the mutation-frequency buckets used by GIAB and
// ClinVar for low-VAF-variant benchmarking:
//
//     (0.00, 0.05]    ultra-low      — clonal mosaicism, FFPE artefacts
//     (0.05, 0.10]    low            — sub-clonal CNV
//     (0.10, 0.25]    moderate       — minor-allele CNV
//     (0.25, 0.50]    balanced-low   — het-ish
//     (0.50, 1.00]    dominant       — major allele
//
// These are hard-coded as named constants so the harness output is
// reproducible across BRANCH versions. If you retune the edges, bump
// kVafBandSchemaVersion in the JSON report so downstream consumers can
// detect the change.
//
// This file is NOT test-only — it is intended to be linked by future
// P1.1 (`analyze`) and P1.2 (`classify`) consumers that want to
// stratify precision/recall against truth-mixture fixtures. Keep it
// dependency-free beyond the STL.

#include <array>
#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace branch::analysis {

// ---------------------------------------------------------------------
// Band edges (right-closed intervals: (lo, hi]).
// ---------------------------------------------------------------------

struct VafBandEdges {
    double lo;         // exclusive lower bound
    double hi;         // inclusive upper bound
    const char* name;  // stable, JSON-safe identifier
};

// Fixed, documented band table. Keep synced with kVafBandSchemaVersion.
inline constexpr std::array<VafBandEdges, 5> kVafBands = {{
    {0.00, 0.05, "0.00-0.05"},
    {0.05, 0.10, "0.05-0.10"},
    {0.10, 0.25, "0.10-0.25"},
    {0.25, 0.50, "0.25-0.50"},
    {0.50, 1.00, "0.50-1.00"},
}};

inline constexpr const char* kVafBandSchemaVersion = "p1.3-vaf-band-1";

// Default recall gate for the low-VAF regime (vaf <= 0.10) used by the
// benchmark test. Set permissively — the goal of P1.3 is to ship the
// harness, not to hit a bar. Tune upwards as the pipeline improves.
inline constexpr double kLowBandRecallFloor = 0.5;

// Return the index of the band that a true VAF falls into. Returns
// kVafBands.size() (i.e. out-of-range) for VAFs <= 0.0 or > 1.0.
[[nodiscard]] std::size_t assign_band(double true_vaf) noexcept;

// ---------------------------------------------------------------------
// Truth + call records.
// ---------------------------------------------------------------------

// One row from the fixture meta TSV (columns:
// allele, copy_count, vaf_target, n_reads, cassette_len_bp).
struct TruthAllele {
    std::string allele_id;
    int copy_count = 0;
    double vaf_target = 0.0;
    std::uint64_t n_reads = 0;
    std::uint64_t cassette_len_bp = 0;
};

// One BRANCH node (post-compaction) promoted to a "called bubble".
// Identified by node_id; copy_count comes from CN:i, read_support
// from RC:i, est_vaf from read_support / sum(read_support).
struct CalledBubble {
    std::uint32_t node_id = 0;
    int est_copy_count = 0;
    double est_vaf = 0.0;
    std::uint32_t read_support = 0;
    std::uint32_t length_bp = 0;
};

// ---------------------------------------------------------------------
// Metrics.
// ---------------------------------------------------------------------

// Per-band counters. Populated by compute_band_stats().
struct BandStats {
    std::string band_name;
    double lo = 0.0;
    double hi = 0.0;

    std::size_t n_truth = 0;              // truth alleles that fall in this band
    std::size_t n_called = 0;             // calls whose matched truth is in this band
    std::size_t n_true_positive = 0;      // calls matched to a truth allele in-band
    std::size_t n_false_positive = 0;     // calls that matched nothing in-band
    std::size_t n_copy_count_correct = 0; // TP subset with est_copy_count == truth

    // |est_vaf - truth.vaf_target| residuals for all TPs in this band.
    std::vector<double> vaf_abs_residuals;

    // Derived metrics; NaN if undefined (division by zero).
    [[nodiscard]] double recall() const noexcept;
    [[nodiscard]] double precision() const noexcept;
    [[nodiscard]] double median_abs_vaf_residual() const noexcept;
    [[nodiscard]] double copy_count_accuracy() const noexcept;
};

struct BandBenchmarkReport {
    std::string schema_version = kVafBandSchemaVersion;
    std::string sample_name;
    std::size_t n_truth_total = 0;
    std::size_t n_called_total = 0;
    std::array<BandStats, kVafBands.size()> bands{};
};

// Maximum absolute VAF distance at which a call is considered to have
// "detected" a truth allele. Picked loosely because BRANCH v0.2 does
// not infer CN and assigns est_vaf from read-support fractions of
// unitig nodes — so any call whose est_vaf lands within
// kVafMatchTolerance of the truth target is treated as a detection.
//
// Tighten this (e.g. to 0.03) once the classifier + VAF estimator
// start producing confidence-weighted per-allele calls.
inline constexpr double kVafMatchTolerance = 0.10;

// Compute per-band statistics.
//
// Matching strategy (documented for reproducibility):
//   1. Each truth allele is assigned to its band by its `vaf_target`.
//   2. For each truth allele we look for the call whose est_vaf is
//      nearest the truth's vaf_target. If that distance is
//      <= kVafMatchTolerance we count the truth as DETECTED
//      (-> TP in the truth's band, with residual = |est - truth|).
//      Otherwise the truth is MISSED (contributes only to the band's
//      n_truth denominator, lowering recall).
//   3. Any call left unmatched after step 2 is a FALSE POSITIVE
//      attributed to the band its own est_vaf falls into.
//   4. Copy-count accuracy per band = fraction of TP matches where
//      call.est_copy_count == truth.copy_count. BRANCH v0.2 emits
//      CN:i:1 for every node pre-classifier, so this metric is
//      expected to be low (often 0) until the CN-inference pass
//      lands; it is still reported so the harness doubles as a
//      regression canary for CN rollout.
//   5. Recall per band = n_TP_in_band / n_truth_in_band.
//   6. Precision per band = n_TP_in_band / (n_TP + n_FP)_in_band.
//
// One-to-one matching: each call can only be claimed by a single
// truth (the first truth to pick it, scanning in truth-ID order).
// This prevents one over-represented call from inflating recall
// across multiple bands.
[[nodiscard]] BandBenchmarkReport compute_band_stats(
    const std::vector<TruthAllele>& truth,
    const std::vector<CalledBubble>& calls,
    const std::string& sample_name);

// ---------------------------------------------------------------------
// Report I/O.
// ---------------------------------------------------------------------

// Parse meta TSV emitted by tests/fixtures/generate_branch_fixture.
// Returns empty vector and sets `err` on failure.
[[nodiscard]] std::vector<TruthAllele> parse_fixture_meta_tsv(
    const std::string& path, std::string* err = nullptr);

// Parse BRANCH's GFA output to build CalledBubble records. Uses the
// S-line LN:i (length), RC:i (read support), and CN:i (copy count).
// est_vaf is set to read_support / sum(read_support across S-lines).
// If sum(RC) == 0 (v0.2 pre-RC path) est_vaf is left at 0.0; callers
// should prefer aggregate_nodes_by_length() in that regime.
// Returns empty vector and sets `err` on failure.
[[nodiscard]] std::vector<CalledBubble> parse_gfa_nodes(
    const std::string& path, std::string* err = nullptr);

// Aggregate a per-node call list into one CalledBubble per
// length-bucket. Nodes whose length_bp rounds to the same bucket
// (bucket = round(length_bp / bucket_size)) are merged into a single
// CalledBubble with:
//   node_id         = lowest node_id in the bucket
//   est_copy_count  = mode of per-node est_copy_count
//   length_bp       = bucket_size * bucket_index (center of bucket)
//   read_support    = sum(read_support)
//   est_vaf         = bucket.size() / total_nodes
//
// Rationale: BRANCH v0.2 emits one node per read in the 4-way bubble
// fixture because unitig compaction can't collapse reads from
// distinct haplotypes (different cassette copy counts -> different
// middle sequence). Grouping nodes by length_bp reconstructs the
// haplotype counts from the graph output — a reasonable "called
// bubble" proxy until BRANCH acquires explicit allele merging.
//
// bucket_size in bp. 500 is a sane default: small enough to keep
// distinct cassette-copy-count alleles apart, large enough to
// tolerate minor length drift from overlap trimming.
//
// Kept as a documented fallback for GFAs where RC:i is missing / all
// zero (pre-P2.2 outputs). The preferred path is build_calls_from_gfa
// below, which returns RC:i-driven est_vaf when RC:i is populated.
[[nodiscard]] std::vector<CalledBubble> aggregate_nodes_by_length(
    const std::vector<CalledBubble>& nodes,
    std::uint32_t bucket_size_bp = 500);

// Which path build_calls_from_gfa took. Reported for test diagnostics
// so we can assert the RC:i-driven path was chosen in regression tests.
enum class CalledBubbleSource {
    kReadSupport,   // RC:i sum > 0 -> per-node VAF derived from RC:i
    kLengthBucket,  // RC:i all zero / missing -> length-bucket fallback
};

// Build CalledBubble records from a parsed GFA node list, preferring
// RC:i (per-node read support) when populated. Logic:
//
//   - If any node has RC:i > 0 (and the total RC:i sum > 0), emit one
//     CalledBubble per node with est_vaf = RC:i / sum(RC:i). This is
//     the P2.2+ path: every node now carries real read-support so the
//     per-node VAF is meaningful.
//   - Otherwise fall back to aggregate_nodes_by_length(bucket_size_bp)
//     (pre-P2.2 behaviour). Length-bucket est_vaf is "fraction of
//     graph nodes in this bucket" — a structural proxy that was the
//     only sane option when every S-line had RC:i:0.
//
// `source_out`, if non-null, receives which branch ran.
[[nodiscard]] std::vector<CalledBubble> build_calls_from_gfa(
    const std::vector<CalledBubble>& nodes,
    std::uint32_t bucket_size_bp = 500,
    CalledBubbleSource* source_out = nullptr);

// Write BandBenchmarkReport as JSON. Hand-rolled, no external deps.
//
// JSON schema (stable within kVafBandSchemaVersion):
// {
//   "schema_version": "p1.3-vaf-band-1",
//   "sample": "<name>",
//   "n_truth": N, "n_called": M,
//   "bands": [
//     {
//       "name": "0.00-0.05", "lo": 0.0, "hi": 0.05,
//       "n_truth": K, "n_called": C,
//       "tp": T, "fp": F,
//       "recall": R, "precision": P,
//       "median_abs_vaf_residual": X,
//       "copy_count_accuracy": Y
//     }, ...
//   ]
// }
//
// NaN metrics are serialised as JSON null.
[[nodiscard]] bool write_report_json(const BandBenchmarkReport& report,
                                     const std::string& path,
                                     std::string* err = nullptr);

// Write a human-readable fixed-width table to `os`.
void write_report_table(const BandBenchmarkReport& report, std::ostream& os);

}  // namespace branch::analysis
