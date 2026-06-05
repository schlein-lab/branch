// BRANCH v0.4.3 — Somatic delta computation.
//
// Computes edit distance between branch consensus sequences and their
// closest HPRC haplotype path using ksw2 alignment.

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace branch::project {

/// Result of somatic delta computation for a single branch vs reference.
struct SomaticDelta {
    std::string branch_id;          ///< Branch/query name from FASTA
    std::string ref_name;           ///< Reference name
    int edit_distance{0};           ///< Total edit distance
    std::string cigar;              ///< CIGAR string from alignment
    std::int64_t aligned_query_len{0};  ///< Aligned query length
    std::int64_t aligned_ref_len{0};    ///< Aligned reference length
    double identity{0.0};           ///< Sequence identity (matches / alignment_len)
};

/// Configuration for somatic delta computation.
struct SomaticDeltaOptions {
    int match{1};           ///< Match score
    int mismatch{-2};       ///< Mismatch penalty (negative)
    int gap_open{4};        ///< Gap open penalty
    int gap_extend{2};      ///< Gap extend penalty
    int threads{4};         ///< Number of threads (for future parallelization)
    /// Skip branches longer than this against ref windows that are also
    /// long. ksw2 in plain DP mode allocates O(query × target) memory;
    /// running it on multi-megabase pairs blows past 50 GiB RSS without
    /// returning anything biologically meaningful (the cost is dominated
    /// by gaps, identity collapses to ~length-ratio). Default 5 000 bp
    /// keeps the typical bubble-resolution branches in scope and skips
    /// the long span-contigs — those are already covered by the
    /// linear+pangenome layers above.
    /// 0 disables the cap.
    int max_branch_len_bp{5000};
};

/// Compute somatic deltas for all branches against reference paths.
///
/// @param branches_fasta  Path to branch consensus FASTA
/// @param ref_paths       Vector of (ref_name, ref_fasta_path) pairs
/// @param opts            Alignment options
/// @param err_out         Optional error message output
/// @return Vector of somatic delta results
std::vector<SomaticDelta> compute_somatic_deltas(
    const std::string& branches_fasta,
    const std::vector<std::pair<std::string, std::string>>& ref_paths,
    const SomaticDeltaOptions& opts,
    std::string* err_out);

/// Compute edit distance between two sequences via the LLmap classical
/// engine (gap-affine WFA2/Gotoh). This is the low-level alignment function.
///
/// @param query    Query sequence
/// @param target   Target/reference sequence
/// @param opts     Alignment options
/// @param cigar_out  Output CIGAR string (optional)
/// @return Edit distance
int aligned_edit_distance(
    const std::string& query,
    const std::string& target,
    const SomaticDeltaOptions& opts,
    std::string* cigar_out);

}  // namespace branch::project
