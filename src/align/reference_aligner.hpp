// BRANCH — Reference Aligner
// Wrapper around minimap2/ksw2 for query-vs-reference alignment.

#pragma once

#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace branch::align {

/// Result of aligning a query sequence against a reference.
struct AlignmentResult {
    int32_t score{0};           ///< Alignment score (ksw2 score)
    float identity{0.0f};       ///< Sequence identity (matches / alignment_length)
    std::string cigar;          ///< CIGAR string
    int32_t ref_start{0};       ///< Reference start position (0-based)
    int32_t ref_end{0};         ///< Reference end position (exclusive)
    int32_t query_start{0};     ///< Query start position (0-based)
    int32_t query_end{0};       ///< Query end position (exclusive)
};

/// Reference aligner using minimap2/ksw2.
///
/// Usage:
///   ReferenceAligner aligner("reference.fa");
///   auto result = aligner.align(query_seq);
///   if (result) {
///       float identity = result->identity;
///   }
class ReferenceAligner {
public:
    /// Construct aligner with reference FASTA path.
    /// @param ref_path Path to reference FASTA file.
    /// @throws std::runtime_error if reference cannot be loaded.
    explicit ReferenceAligner(const std::string& ref_path);
    
    ~ReferenceAligner();
    
    // Non-copyable, movable
    ReferenceAligner(const ReferenceAligner&) = delete;
    ReferenceAligner& operator=(const ReferenceAligner&) = delete;
    ReferenceAligner(ReferenceAligner&&) noexcept;
    ReferenceAligner& operator=(ReferenceAligner&&) noexcept;
    
    /// Align a query sequence against the reference.
    /// @param query Query sequence (DNA, uppercase ACGTN).
    /// @return AlignmentResult if alignment found, nullopt otherwise.
    [[nodiscard]] std::optional<AlignmentResult> align(std::string_view query) const;
    
    /// Align query and return identity score (0.0-1.0).
    /// Convenience method for feature extraction.
    /// @return Identity if aligned, 0.0f otherwise.
    [[nodiscard]] float align_identity(std::string_view query) const;
    
    /// Check if aligner is valid (reference loaded successfully).
    [[nodiscard]] bool is_valid() const noexcept;
    
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace branch::align
