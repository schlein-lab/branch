// BRANCH — LLmap alignment backend (public interface)
//
// This is the single seam between BRANCH and the LLmap alignment engine.
// It replaces the former minimap2/ksw2 dependency entirely:
//
//   * align_pair()      — pairwise base-level alignment of two given
//                         sequences (query vs target) → CIGAR + edit
//                         distance. Former ksw_extz2_sse() use.
//   * ReferenceMapper   — seed-chain-extend mapping of a query against an
//                         indexed reference set → best hit(s) with CIGAR
//                         and identity. Former mm_idx + mm_map use.
//
// DESIGN: every LLmap (`llmap::classical::*`) type is confined to the
// .cpp. This header is plain C++20 and exposes only std/POD types, so the
// rest of BRANCH neither sees LLmap headers nor needs C++23. Only
// llmap_align_backend.cpp is compiled as C++23 and links the LLmap
// static libraries. This keeps the LLmap coupling auditable in one place
// and lets the LLmap version be bumped without touching callers.

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace branch::align::llmap_backend {

// ---------------------------------------------------------------------------
// Pairwise base-level alignment (ksw2 replacement)
// ---------------------------------------------------------------------------

// Penalty/scope knobs for align_pair(). Defaults mirror LLmap's WFA2
// gap-affine model (penalties, not scores) and are sane for the
// high-identity branch-vs-pangenome-path comparison somatic_delta needs.
struct PairwiseOptions {
    int  mismatch_penalty = 4;
    int  gap_open         = 6;   // affine: cost = gap_open + (len-1)*gap_extend
    int  gap_extend       = 2;
    bool global           = true;  // true = end-to-end (edit distance semantics)
    // Abort guards — keep a stray megabase contig from blowing up the DP.
    std::size_t max_query_len = 200000;
    std::size_t max_ref_len   = 200000;
};

struct PairwiseAlignment {
    std::string cigar;            // BAM-style CIGAR (e.g. "13=1X19=")
    int  edit_distance = 0;       // mismatches + inserted + deleted bases
    int  score         = 0;       // LLmap alignment score (lower = better)
    float identity     = 0.0f;    // matches / aligned length
    int  query_start = 0, query_end = 0;
    int  ref_start   = 0, ref_end   = 0;
    std::size_t num_matches = 0, num_mismatches = 0;
    std::size_t num_insertions = 0, num_deletions = 0;
};

// Align `query` against `target`. Returns nullopt only when alignment is
// impossible (empty input is handled by the caller) or LLmap aborts on a
// limit. Callers that need a guaranteed number (e.g. somatic_delta) should
// fall back to a direct difference count on nullopt.
[[nodiscard]] std::optional<PairwiseAlignment>
align_pair(std::string_view query,
           std::string_view target,
           const PairwiseOptions& opts = {});

// Convenience: edit distance only, with `cigar_out` optionally filled.
// Returns max(len) on failure so it never throws — a drop-in for the old
// ksw2_edit_distance() contract.
[[nodiscard]] int
pairwise_edit_distance(std::string_view query,
                       std::string_view target,
                       const PairwiseOptions& opts = {},
                       std::string* cigar_out = nullptr);

// ---------------------------------------------------------------------------
// Sequence-vs-reference mapping (minimap2 mm_idx + mm_map replacement)
// ---------------------------------------------------------------------------

// A single mapping hit, flattened from llmap::classical::ClassicalAlignment.
struct RefHit {
    std::string ref_name;
    std::uint32_t ref_id = 0;
    int  ref_start = 0, ref_end = 0;       // 0-based, end exclusive
    int  query_start = 0, query_end = 0;
    bool is_forward = true;
    bool is_primary = true;
    std::string cigar;
    int   score       = 0;
    float identity    = 0.0f;
    int   aligned_len = 0;   // query bases spanned by the alignment
    int   mismatches  = 0;   // substitutions + indel bases within the span
    std::uint32_t mapq = 0;
};

// Maps queries against a fixed reference set. Build once, query many —
// the minimizer index and reference sequences are held for the mapper's
// lifetime. Thread-compatible for concurrent const queries.
class ReferenceMapper {
public:
    // Preset selects k/w and identity floor:
    //   "asm5"/"map-hifi" → k=19,w=19, identity≥0.85 (high-identity assembly/HiFi)
    //   "map-ont"         → k=15,w=10, identity≥0.70 (noisier ONT)
    // Unknown presets fall back to map-hifi.
    ReferenceMapper(const std::string& ref_fasta_path,
                    const std::string& preset = "map-hifi");
    ReferenceMapper(const std::vector<std::string>& ref_names,
                    const std::vector<std::string>& ref_seqs,
                    const std::string& preset = "map-hifi");
    ~ReferenceMapper();

    ReferenceMapper(const ReferenceMapper&) = delete;
    ReferenceMapper& operator=(const ReferenceMapper&) = delete;
    ReferenceMapper(ReferenceMapper&&) noexcept;
    ReferenceMapper& operator=(ReferenceMapper&&) noexcept;

    [[nodiscard]] bool valid() const noexcept;

    // Primary (best) hit, or nullopt if the query maps nowhere above the
    // identity/length floor.
    [[nodiscard]] std::optional<RefHit> best(std::string_view query) const;

    // All reported hits (primary first), up to the preset's max_alignments.
    [[nodiscard]] std::vector<RefHit> all(std::string_view query) const;

    // Batch primary mapping: map every sequence in `seqs` and return its
    // best hit (nullopt where nothing maps). Parallelized across `n_threads`
    // via LLmap's own thread pool — the right entry point for WG-scale read
    // streams (one minimizer index, many reads). Order matches `seqs`.
    [[nodiscard]] std::vector<std::optional<RefHit>>
    best_batch(const std::vector<std::string>& seqs, int n_threads) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace branch::align::llmap_backend
