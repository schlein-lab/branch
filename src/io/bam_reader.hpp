#pragma once

// BRANCH v0.2 — Streaming BAM reader.
//
// Yields one (name, sequence) tuple per read using htslib. Secondary and
// supplementary alignments are skipped so that each read is emitted exactly
// once regardless of how many alignments it produced. Reverse-strand reads
// are reverse-complemented so that the caller (minimizer sketcher, overlap
// backend) sees the original sequencing orientation rather than the
// reference-forward orientation baked into the BAM record.
//
// BAMs from HiFi / linked-read pipelines frequently carry unmapped records
// (PacBio CCS with `pbmm2 --unmapped`) — those are yielded normally with
// their stored sequence.
//
// NOTE: This is read-only; we do not inspect CIGAR, MD, or auxiliary tags.
// BRANCH treats the alignment purely as a transport format for sequences.

#include <cstddef>
#include <string>

// We deliberately do NOT forward-declare htslib's types here: htslib uses
// typedefs (e.g. `typedef htsFile samFile`), not struct tags, so a
// forward `struct samFile` is a distinct incomplete type and produces
// confusing compile errors. BamReader keeps the handles as `void*` in
// this header and only reinterprets them in bam_reader.cpp.

namespace branch::io {

// Minimal sequential reader. Not copyable, not thread-safe.
class BamReader {
public:
    BamReader() = default;
    ~BamReader();

    BamReader(const BamReader&) = delete;
    BamReader& operator=(const BamReader&) = delete;
    BamReader(BamReader&&) noexcept;
    BamReader& operator=(BamReader&&) noexcept;

    // Open a BAM/CRAM/SAM file. Returns false if the file cannot be opened
    // or its header cannot be parsed.
    bool open(const std::string& path);

    // Fetch the next read. Secondary (0x100) and supplementary (0x800)
    // alignments are silently skipped. Reverse-strand (0x10) reads are
    // reverse-complemented in-place so `seq` is always in the read's
    // original sequencing orientation.
    //
    // Returns false at EOF or on a non-recoverable decoding error.
    bool next(std::string& name, std::string& seq);

    // Close handles and free the record buffer. Idempotent; implicit on dtor.
    void close();

    // Number of reads successfully yielded so far.
    [[nodiscard]] std::size_t reads_read() const noexcept { return n_read_; }

private:
    // Opaque to callers; real htslib types live in the .cpp.
    void* file_ = nullptr;     // htsFile* (samFile*).
    void* header_ = nullptr;   // sam_hdr_t*.
    void* record_ = nullptr;   // bam1_t*.
    std::size_t n_read_{0};
};

}  // namespace branch::io
