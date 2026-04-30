// BRANCH v0.5 — Phase 3 Stage 1: cheap whole-genome k-mer-coverage screen.
//
// Purpose
//   Pre-screen every base of every haplotype with cheap, parallelisable
//   signals that reveal CN-anomalous regions WITHOUT designing primers
//   anywhere. Then Phase 3 Stage 2 designs vPCR amplicons only in the
//   regions Stage 1 flagged. This makes it tractable to run a precise
//   counter on a 3 Gbp human haplotype where naive design would emit
//   billions of amplicon candidates.
//
// No coding-vs-non-coding bias
//   Stage 1 is intentionally annotation-blind. It does not consult any
//   GTF/BED/GENCODE input, does not exclude or downweight regions based
//   on coding status. Annotation is a post-hoc interpretation step
//   over the Stage 1 + 2 output.
//
// Sensitivity-first
//   We OR four independent signals so a window can flag for any of:
//     - k-mer-coverage anomaly (z ≥ kZThreshold) — CN gain or loss
//     - dimer-entropy ≤ kRepEntropy — repetitive (satellite, STR, retro)
//     - distinct-kmer count anomaly — sequence-complexity outlier
//     - minimizer-density anomaly — coverage drop-out / breakpoint
//   We also force-flag known paralog hotspots (IGH, KIR, MHC) even when
//   no signal fires, so users never miss the canonical CN-rich loci.
//   We sweep three window sizes {300, 1000, 3000} bp at 50 % overlap so
//   focal CNVs don't fall through the raster.
//
// Output
//   - Stage1Window[] per haplotype window (full track, BED-emittable)
//   - CandidateRegion[] with adjacent-flag merging; this feeds Stage 2.

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "graph/kmer_sketch.hpp"

namespace branch::wg {

/// Window class flags (OR-able into a single byte).
enum Stage1Flag : std::uint8_t {
    kStage1Elevated  = 1u << 0,  // |z(mean_kmer_count)| ≥ threshold, mean above
    kStage1Depleted  = 1u << 1,  // |z(mean_kmer_count)| ≥ threshold, mean below
    kStage1Repeat    = 1u << 2,  // dimer-entropy ≤ kRepEntropy
    kStage1ComplexOut= 1u << 3,  // distinct-kmer count outlier
    kStage1MinDrop   = 1u << 4,  // minimizer density drop
    kStage1Paralog   = 1u << 5,  // pinned paralog locus (IGH/KIR/MHC etc)
};

struct Stage1Window {
    std::uint32_t hap_idx;
    std::uint32_t start_bp;
    std::uint32_t length_bp;
    // Raw signals
    float         mean_kmer_count = 0;
    std::uint32_t distinct_kmer_count = 0;
    float         dimer_entropy = 0;
    float         minimizer_density = 0;  // minimizers / bp
    // Standardised against per-haplotype median + MAD
    float         z_mean = 0;
    float         z_distinct = 0;
    // OR'd flags
    std::uint8_t  flags = 0;
};

struct CandidateRegion {
    std::uint32_t hap_idx;
    std::uint32_t start_bp;
    std::uint32_t end_bp;
    std::uint8_t  flags;          // OR of contributing windows
    float         max_abs_z;      // strongest signal in this region
    std::uint32_t n_windows;
};

/// Default thresholds — sensitive (over-flag) by design. The paper-rate
/// can be tuned by CLI overrides.
inline constexpr double  kStage1ZThreshold       = 1.5;
inline constexpr double  kStage1RepEntropy       = 0.40;   // ≤ → repetitive
inline constexpr double  kStage1MinDensityZ      = 1.5;
inline constexpr std::uint32_t kStage1WindowSizes[] = { 300, 1000, 3000 };
inline constexpr double  kStage1OverlapFraction  = 0.50;
inline constexpr std::uint32_t kStage1MaxRegionGap = 2000;  // bp; merge windows within this

/// Pinned paralog loci on GRCh38 (chr index → (start_bp, end_bp)).
/// These are mapped onto haplotype-coordinate space heuristically:
/// when a haplotype contains paralog-name k-mers (built once as part of
/// design), the loci are flagged automatically. For now we expose this
/// list as data so callers can pin specific (hap_idx, start, end) ranges.
struct PinnedLocus { std::string name; std::uint32_t start_bp; std::uint32_t end_bp; };

class Stage1Screener {
public:
    explicit Stage1Screener(const ::branch::graph::TechProfile& profile);

    /// Compute Stage1Window[] for one haplotype.
    ///   - hap_seq: the haplotype sequence
    ///   - counter: canonical_kmer_bits → count from Phase 0
    ///   - hap_idx: haplotype id (0 or 1) stamped into output
    ///   - pinned: optional list of (hap_idx-relative) ranges to pin
    void screen_haplotype(
        std::uint32_t hap_idx,
        const std::string& hap_seq,
        const std::unordered_map<std::uint64_t, std::uint32_t>& counter,
        std::vector<Stage1Window>& out_windows,
        const std::vector<PinnedLocus>& pinned = {}) const;

    /// Merge adjacent flagged windows into CandidateRegion[].
    void merge_into_regions(const std::vector<Stage1Window>& windows,
                            std::vector<CandidateRegion>& out_regions) const;

    /// CLI-overridable thresholds.
    double z_threshold = kStage1ZThreshold;
    double rep_entropy = kStage1RepEntropy;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
