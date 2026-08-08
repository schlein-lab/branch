#pragma once

// BRANCH v0.1 — Overlap spill store.
//
// WGS-scale overlap detection produces far more OverlapPair records than
// fit in RAM at once. This module gives the overlap phase a disk-backed
// staging area: OverlapSpillWriter appends unsorted batches to a spill
// file as they're produced, merge_spills folds any number of spill files
// into one sorted, deduplicated store (falling back to an external
// chunked sort + k-way merge once the input exceeds a caller-supplied RAM
// budget), and OverlapStoreReader streams the merged result back with
// O(buffer) memory regardless of total file size.
//
// File format ("BROV1"): a 5-byte magic, an 8-byte little-endian-host
// record count (patched in place when the writer closes), followed by
// that many OverlapPair records written verbatim.
// branch::backend::OverlapPair is a 24-byte POD with a static_assert-
// enforced layout, so a raw byte copy is all serialization needs to do —
// spill files are read back on the same build that wrote them. Spill
// files and the final merged store share this exact format, so
// OverlapStoreReader reads both.
//
// Error handling: OverlapSpillWriter/OverlapStoreReader::open() return
// bool as elsewhere in src/io/; OverlapStoreReader additionally exposes
// last_error() with a specific message (bad magic vs. truncated vs.
// cannot-open) since a plain bool can't carry that detail and the caller
// needs to know which. merge_spills orchestrates several of these
// open/read steps across possibly many files, so — matching the
// std::runtime_error convention used for file problems elsewhere in this
// codebase (e.g. align/reference_aligner.cpp) — it throws
// std::runtime_error with a descriptive message instead of threading an
// error string through every intermediate step.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <span>
#include <string>
#include <vector>

#include "backend/backend_vtable.hpp"

namespace branch::io {

using OverlapPair = branch::backend::OverlapPair;

// Appends OverlapPair records to a BROV1 spill file. Records may be
// appended in any order — sorting and dedup happen later in
// merge_spills. Not copyable, not thread-safe.
class OverlapSpillWriter {
public:
    OverlapSpillWriter() = default;
    ~OverlapSpillWriter();

    OverlapSpillWriter(const OverlapSpillWriter&) = delete;
    OverlapSpillWriter& operator=(const OverlapSpillWriter&) = delete;
    OverlapSpillWriter(OverlapSpillWriter&&) noexcept;
    OverlapSpillWriter& operator=(OverlapSpillWriter&&) noexcept;

    // Creates (truncating) `path` and writes a placeholder header.
    // Returns false if the file cannot be opened for writing.
    bool open(const std::string& path);

    // Buffers `pair`/`pairs` for write; buffered internally (>=1 MiB)
    // rather than issuing a write syscall per record.
    void append(const OverlapPair& pair);
    void append(std::span<const OverlapPair> pairs);

    // Flushes buffered records, patches the header's record count via a
    // seek, and closes the file. Idempotent (a second call is a no-op
    // returning true). Also invoked from the destructor; call close()
    // explicitly if you need to observe failure.
    bool close();

    [[nodiscard]] bool is_open() const noexcept { return file_.is_open(); }
    [[nodiscard]] std::uint64_t records_written() const noexcept { return record_count_; }

private:
    void flush_buffer();

    std::ofstream file_;
    std::vector<OverlapPair> buffer_;
    std::uint64_t record_count_ = 0;
};

// Merges `spill_paths` (each a BROV1 file, records in any internal
// order) into one BROV1 file at `out_path`, ordered by
// (read_a, read_b, offset_a). Records that share (read_a, read_b) are
// collapsed into a single record, keeping the larger overlap_len; any
// remaining tie is broken by a fixed total order over the rest of the
// fields (see overlap_less in the .cpp), so merging the same inputs
// twice — in any order, via either the in-RAM or external path —
// always produces byte-identical output.
//
// If the total input fits within max_ram_bytes, the merge runs as one
// in-memory sort. Otherwise the input is split into sorted runs of up
// to max_ram_bytes each (spilled as temp files under `tmp_dir`, removed
// before returning), which are then combined with a k-way merge; at no
// point does materially more than max_ram_bytes of record data sit in
// memory.
//
// Throws std::runtime_error if a spill file cannot be opened or is
// corrupt/truncated, or if a temp/output file cannot be created.
void merge_spills(const std::vector<std::string>& spill_paths,
                   const std::string& out_path,
                   const std::string& tmp_dir,
                   std::uint64_t max_ram_bytes = std::uint64_t{1} << 30);

// Sequential reader over a BROV1 file (a raw spill, or a merge_spills
// output — same format). open() reads only the fixed-size header;
// next() streams one record at a time, so peak memory is O(1) (a single
// record) regardless of record_count(). Not copyable, not thread-safe.
class OverlapStoreReader {
public:
    OverlapStoreReader() = default;
    ~OverlapStoreReader() = default;

    OverlapStoreReader(const OverlapStoreReader&) = delete;
    OverlapStoreReader& operator=(const OverlapStoreReader&) = delete;
    OverlapStoreReader(OverlapStoreReader&&) noexcept;
    OverlapStoreReader& operator=(OverlapStoreReader&&) noexcept;

    // Opens `path` and validates the header. Returns false if the file
    // cannot be opened, the magic doesn't match, or the file is shorter
    // than its header claims — call last_error() for which one.
    bool open(const std::string& path);
    void close();

    [[nodiscard]] bool is_open() const noexcept { return file_.is_open(); }
    [[nodiscard]] std::uint64_t record_count() const noexcept { return record_count_; }
    [[nodiscard]] std::uint64_t records_read() const noexcept { return records_read_; }
    [[nodiscard]] const std::string& last_error() const noexcept { return last_error_; }

    // Reads the next record into `out`. Returns false at end of stream
    // (not an error — check last_error() to distinguish EOF from a
    // mid-stream read failure, which also sets it).
    bool next(OverlapPair& out);

private:
    std::ifstream file_;
    std::uint64_t record_count_ = 0;
    std::uint64_t records_read_ = 0;
    std::string last_error_;
};

}  // namespace branch::io
