// BRANCH v0.4 — Somatic delta computation.
//
// Computes edit distance between branch consensus sequences and their
// closest HPRC haplotype path. Edits represent candidate somatic variants.

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace branch::project {

/// Type of somatic edit.
enum class EditType {
    SNV,    ///< Single nucleotide variant
    INS,    ///< Insertion
    DEL,    ///< Deletion
    SUB,    ///< Complex substitution (MNV)
};

/// A single somatic edit (variant call).
struct SomaticEdit {
    EditType type;
    std::uint32_t pos;         ///< Position in branch coordinates
    std::string ref_seq;       ///< Reference sequence
    std::string alt_seq;       ///< Alternate sequence
    std::uint32_t ref_pos{0};  ///< Position in reference coordinates
};

/// Result of somatic delta computation for a single branch.
struct SomaticDeltaResult {
    std::string branch_id;
    std::string closest_path;       ///< HPRC path used for comparison
    std::vector<SomaticEdit> edits;
    double vaf{0.0};                ///< Variant allele frequency
    double coverage{0.0};           ///< Branch coverage
    std::uint32_t edit_distance{0}; ///< Total edit distance
};

/// Configuration for somatic delta computation.
struct SomaticDeltaConfig {
    std::string branch_fasta;    ///< Branch consensus FASTA
    std::string hprc_sequences;  ///< Path to HPRC path sequences (extracted from GAF)
    int threads{1};
};

/// Compute edit distance between two sequences using ksw2.
///
/// @param seq1  First sequence
/// @param seq2  Second sequence
/// @return Edit distance
///
/// TODO: Implement using existing ksw2 wrapper
std::uint32_t edit_distance(const std::string& seq1, const std::string& seq2);

/// Compute somatic deltas for all branches.
///
/// @param config  Delta computation configuration
/// @return Vector of delta results
///
/// TODO: Implement in v0.4
std::vector<SomaticDeltaResult> compute_somatic_deltas(const SomaticDeltaConfig& config);

/// Extract edits from alignment CIGAR.
///
/// @param cigar     CIGAR string from alignment
/// @param ref_seq   Reference sequence
/// @param query_seq Query sequence
/// @return Vector of somatic edits
///
/// TODO: Implement in v0.4
std::vector<SomaticEdit> extract_edits_from_cigar(
    const std::string& cigar,
    const std::string& ref_seq,
    const std::string& query_seq);

}  // namespace branch::project
