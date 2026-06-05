#pragma once
//
// FastqIndex — random-access reader for (gz)FASTQ files.
//
// Built once at pipeline start; sequences are loaded on demand instead of
// being held in RAM for the whole run. This is the central piece of the
// memory-budget refactor: the previous design held every read's full
// sequence in `vector<pair<ReadId, string>>` for the lifetime of the
// pipeline (60 GB on HiFi-WG). With this index the only persistent state
// is the offset table (~50 MB for 4M reads).
//
// Two access modes:
//   1. for_each(F): sequential streaming — cheapest, used by sketcher.
//   2. read_seq(rid, out): random access by ReadId — used by gfa_writer.
//
// Both modes return raw nucleotide strings (no headers, no quals).
// gzipped FASTQ is handled transparently. For gzip we still seek by
// uncompressed offset because we maintain a separate read-by-read scan
// position; this is fine because random access in v1 is sequential
// across the unitig list (gfa_writer walks unitigs in id order).

#include "types.hpp"

#include <cstdio>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace branch::wg::phaser {

class FastqIndex {
public:
    FastqIndex() = default;
    ~FastqIndex();

    FastqIndex(const FastqIndex&) = delete;
    FastqIndex& operator=(const FastqIndex&) = delete;

    // Scan the FASTQ once and record (rid, file_offset, seq_length) per
    // record. Throws std::runtime_error on parse failure.
    //
    // `assign_ids` controls whether ReadIds are assigned sequentially
    // 0..N-1 (true) or parsed from the FASTQ header (false). The phaser
    // pipeline currently uses sequential ids (callers map back to header
    // strings via read_id_strings).
    void build(const std::string& fastq_path, bool assign_ids = true);

    // Alternative builder for tests: stores the sequences in memory.
    // No disk I/O; read_seq() / for_each() use the stored copies.
    // ReadId is taken from each pair's first element. id_strings can be
    // empty (then numeric strings are used).
    void build_from_memory(
        std::vector<std::pair<ReadId, std::string>> reads_vec,
        std::vector<std::string> id_strings = {});

    // Load one read's sequence into `out`. Resizes `out` to seq length.
    void read_seq(ReadId rid, std::string& out);

    // Iterate every read sequentially: F(rid, seq_view).
    // Cheaper than calling read_seq(rid) N times because there's no seek.
    template <class F>
    void for_each(F&& fn) {
        if (in_memory_) {
            for (const auto& [rid, seq] : in_memory_reads_) {
                fn(rid, seq);
            }
            return;
        }
        rewind();
        std::string seq;
        ReadId rid = 0;
        while (read_next(seq)) {
            fn(rid, seq);
            ++rid;
        }
    }

    std::size_t n_reads() const {
        return in_memory_ ? in_memory_reads_.size() : entries_.size();
    }
    Position    seq_length(ReadId rid) const {
        if (in_memory_) {
            return rid < in_memory_reads_.size()
                   ? static_cast<Position>(in_memory_reads_[rid].second.size())
                   : 0;
        }
        return rid < entries_.size() ? entries_[rid].seq_length : 0;
    }

    // Approximate size of the in-memory index in bytes (not counting
    // header strings).
    std::size_t index_bytes() const {
        return entries_.size() * sizeof(Entry);
    }

    // Original FASTQ headers (for round-tripping read names to outputs).
    const std::vector<std::string>& read_id_strings() const {
        return id_strings_;
    }

    // Path to the (possibly decompressed) FASTQ file currently backing
    // this index. Empty for in-memory mode. Used by AnchoredPhaser to
    // hand the file to minimap2.
    const std::string& fastq_path() const { return fastq_path_; }
    bool is_in_memory() const { return in_memory_; }

    // Thread-safe random read access.
    //
    // The standard read_seq() uses a single shared FILE* held in
    // `reader_` and is therefore NOT thread-safe — two threads calling
    // it concurrently would race on the file position. For OpenMP
    // parallel voting we need independent file cursors per thread.
    //
    // Usage:
    //   #pragma omp parallel
    //   {
    //       auto handle = idx.open_thread_reader();
    //       #pragma omp for
    //       for (...) idx.read_seq_via(handle, rid, out);
    //   }
    // Each ThreadReader owns its own FILE* / in-memory ref.
    struct ThreadReader {
        std::FILE* fp = nullptr;
        bool in_memory = false;
        const FastqIndex* parent = nullptr;
        ~ThreadReader();
        ThreadReader() = default;
        ThreadReader(const ThreadReader&) = delete;
        ThreadReader& operator=(const ThreadReader&) = delete;
        ThreadReader(ThreadReader&& o) noexcept;
        ThreadReader& operator=(ThreadReader&& o) noexcept;
    };

    ThreadReader open_thread_reader() const;
    void read_seq_via(ThreadReader& h, ReadId rid, std::string& out) const;

private:
    struct Entry {
        std::uint64_t seq_offset;   // byte offset of seq line start
        std::uint32_t seq_length;
    };

    std::vector<Entry>       entries_;
    std::vector<std::string> id_strings_;
    std::string              fastq_path_;
    bool                     is_gzipped_ = false;

    // In-memory mode (test path).
    bool in_memory_ = false;
    std::vector<std::pair<ReadId, std::string>> in_memory_reads_;

    // Streaming reader state. We keep a single reader open to amortize
    // gzip-stream startup cost; for random access we may have to reopen.
    struct Reader;
    Reader* reader_ = nullptr;

    void open_reader();
    void close_reader();
    void rewind();
    bool read_next(std::string& seq);  // returns false at EOF
    void seek_to(std::uint64_t offset);
};

}  // namespace branch::wg::phaser
