#pragma once

// BRANCH v0.2 — Streaming FASTQ reader.
//
// Reads a plain-text (uncompressed) FASTQ file record-by-record.
// gz-support is intentionally deferred to keep v0.2 self-contained
// (no zlib dep yet); consumers that need .gz input can pipe through
// `zcat` into stdin.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <istream>
#include <memory>
#include <string>

namespace branch::io {

struct FastqRecord {
    std::string name;
    std::string sequence;
    std::string qualities;
};

class FastqReader {
public:
    // Open a file on disk. File may be plain FASTQ only in v0.2.
    explicit FastqReader(const std::string& path);

    // Wrap an existing stream (e.g. std::cin). Reader does not own it.
    explicit FastqReader(std::istream& stream);

    // Disallow copy, allow move.
    FastqReader(const FastqReader&) = delete;
    FastqReader& operator=(const FastqReader&) = delete;
    FastqReader(FastqReader&&) noexcept = default;
    FastqReader& operator=(FastqReader&&) noexcept = default;

    // True if the underlying stream opened successfully.
    [[nodiscard]] bool ok() const noexcept;

    // Read the next record. Returns false on EOF or malformed input.
    bool next_record(FastqRecord& out);

    // Count of records successfully returned so far.
    [[nodiscard]] std::size_t records_read() const noexcept { return n_read_; }

private:
    std::unique_ptr<std::ifstream> owned_file_;
    std::istream* stream_;
    std::size_t n_read_{0};
};

}  // namespace branch::io
