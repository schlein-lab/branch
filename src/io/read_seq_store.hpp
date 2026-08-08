#pragma once

// BRANCH v0.1 — 2-bit-packed read sequence store.
//
// WGS-scale read sets don't fit in RAM as std::string per read. This
// module packs each read's sequence to 2 bits/base on disk (.rsq) and
// gives random-access lookup (fetch/fetch_window) that touches only the
// index at open() and one seek+read per call — never the whole file.
//
// File layout ("BRSQ1"), written in this order:
//   [sequence data]  packed 2-bit bases, one contiguous span per read,
//                     in read_id order (read_id 0, 1, 2, ... — implicit
//                     from write order, never stored explicitly)
//   [exceptions]      (u32 read_id, u32 pos) per non-ACGT base, ascending
//                     by (read_id, pos) — the order ReadSeqStoreWriter
//                     naturally produces them in, so no sort is needed
//   [index]           (u64 offset, u32 length) per read, read_id order —
//                     offset/length address the sequence-data region;
//                     length is in BASES, not packed bytes
//   [footer]          fixed-size, at absolute EOF: magic "BRSQ1" (5B) +
//                     index_offset (u64) + exceptions_offset (u64) +
//                     exception_count (u64) + read_count (u64)
//
// Only exactly 'A'/'C'/'G'/'T' (uppercase) pack as a base; every other
// character — N, IUPAC ambiguity codes, lowercase — is recorded as an
// exception and packed as a placeholder 'A'. Exceptions carry no code,
// only a position: re-applying one always yields 'N'. This matches how
// raw reads are actually encoded (basecallers emit uppercase ACGTN;
// the wider IUPAC alphabet is a reference/consensus concept, not a raw
// read one) and the task's own (read_id, pos) exception schema, which
// has no field for the original character.
//
// IndexEntry is 12 bytes on disk (u64 + u32, no padding) but 16 bytes
// once loaded into a std::vector<IndexEntry> in RAM: alignof(u64) == 8
// forces the in-memory struct's size up to a multiple of 8. That 16B/
// read figure is exactly what the loaded index is expected to cost.
//
// fetch_window's out-of-range behavior mirrors std::string::substr:
// `start > sequence_length` throws std::out_of_range (an invalid
// starting point), but `start + len` beyond the end just clips the
// result to whatever remains — never throws for that case alone.
// fetch(read_id) is fetch_window(read_id, 0, full_length).
//
// Error handling: mirrors src/io/overlap_store.hpp — open() returns
// bool with last_error() carrying detail (bad magic / truncated /
// cannot-open); invalid read_id or fetch_window start is programmer
// error, reported via std::out_of_range rather than a bool, since
// there's no caller-recoverable "try the next file" scenario for it.

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace branch::io {

// Streams read sequences to a BRSQ1 file. Read IDs are implicit and
// sequential, assigned in add() call order starting at 0. Not copyable,
// not thread-safe.
class ReadSeqStoreWriter {
public:
    ReadSeqStoreWriter() = default;
    ~ReadSeqStoreWriter();

    ReadSeqStoreWriter(const ReadSeqStoreWriter&) = delete;
    ReadSeqStoreWriter& operator=(const ReadSeqStoreWriter&) = delete;
    ReadSeqStoreWriter(ReadSeqStoreWriter&&) noexcept;
    ReadSeqStoreWriter& operator=(ReadSeqStoreWriter&&) noexcept;

    // Creates (truncating) `path`. Returns false if it cannot be opened
    // for writing.
    bool open(const std::string& path);

    // Packs and streams `seq` to disk immediately; only the (tiny) index
    // and exception lists accumulate in memory, not sequence data.
    // Returns the assigned read_id.
    std::uint32_t add(const std::string& seq);

    // Writes the exceptions/index/footer and closes. Idempotent (a
    // second call is a no-op returning true). Also invoked from the
    // destructor; call close() explicitly if you need to observe
    // failure.
    bool close();

    [[nodiscard]] bool is_open() const noexcept { return file_.is_open(); }
    [[nodiscard]] std::uint32_t reads_written() const noexcept {
        return static_cast<std::uint32_t>(index_.size());
    }

private:
    struct PendingIndexEntry {
        std::uint64_t offset;
        std::uint32_t length;
    };
    struct PendingException {
        std::uint32_t read_id;
        std::uint32_t pos;
    };

    std::ofstream file_;
    std::vector<PendingIndexEntry> index_;
    std::vector<PendingException> exceptions_;
    std::uint64_t write_offset_ = 0;
};

// Random-access reader over a BRSQ1 file. open() loads only the index
// (record_count() * 16B in RAM) and the exception list (typically
// sparse); it never loads sequence data. fetch()/fetch_window() each
// perform exactly one seek + read against the sequence-data region.
// Not copyable, not thread-safe.
class ReadSeqStoreReader {
public:
    ReadSeqStoreReader() = default;
    ~ReadSeqStoreReader() = default;

    ReadSeqStoreReader(const ReadSeqStoreReader&) = delete;
    ReadSeqStoreReader& operator=(const ReadSeqStoreReader&) = delete;
    ReadSeqStoreReader(ReadSeqStoreReader&&) noexcept;
    ReadSeqStoreReader& operator=(ReadSeqStoreReader&&) noexcept;

    // Opens `path`, validates the footer, and loads the index and
    // exception list. Returns false if the file cannot be opened, the
    // magic doesn't match, or the layout is internally inconsistent
    // (truncated) — call last_error() for which one.
    bool open(const std::string& path);
    void close();

    [[nodiscard]] bool is_open() const noexcept { return file_.is_open(); }
    [[nodiscard]] std::uint32_t read_count() const noexcept {
        return static_cast<std::uint32_t>(index_.size());
    }
    [[nodiscard]] const std::string& last_error() const noexcept { return last_error_; }

    // Length in bases of `read_id`. Throws std::out_of_range if
    // read_id >= read_count().
    [[nodiscard]] std::uint32_t sequence_length(std::uint32_t read_id) const;

    // Fetches the full sequence for `read_id`. Throws std::out_of_range
    // if read_id >= read_count().
    std::string fetch(std::uint32_t read_id);

    // Fetches up to `len` bases starting at `start` (both in bases).
    // Throws std::out_of_range if read_id >= read_count() or
    // start > sequence_length(read_id). If start + len runs past the
    // end of the sequence, the result is silently clipped to whatever
    // remains (may be shorter than `len`, possibly empty) — see the
    // file-level comment for the rationale.
    std::string fetch_window(std::uint32_t read_id, std::uint32_t start, std::uint32_t len);

private:
    struct IndexEntry {
        std::uint64_t offset;
        std::uint32_t length;
    };
    struct ExceptionEntry {
        std::uint32_t read_id;
        std::uint32_t pos;
    };

    std::ifstream file_;
    std::vector<IndexEntry> index_;
    std::vector<ExceptionEntry> exceptions_;
    std::string last_error_;
};

}  // namespace branch::io
