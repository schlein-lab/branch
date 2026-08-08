#include "io/overlap_store.hpp"

#include <algorithm>
#include <array>
#include <filesystem>
#include <queue>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace branch::io {

namespace {

constexpr std::array<char, 5> kMagic = {'B', 'R', 'O', 'V', '1'};
constexpr std::streamoff kRecordCountOffset = static_cast<std::streamoff>(kMagic.size());
constexpr std::streamoff kHeaderSize =
    kRecordCountOffset + static_cast<std::streamoff>(sizeof(std::uint64_t));
constexpr std::size_t kMinBufferBytes = 1u << 20;  // >=1 MiB, per task spec.
constexpr std::size_t kMinBufferRecords = kMinBufferBytes / sizeof(OverlapPair);

void write_u64(std::ostream& os, std::uint64_t v) {
    os.write(reinterpret_cast<const char*>(&v), sizeof(v));
}

bool read_u64(std::istream& is, std::uint64_t& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(v));
    return static_cast<bool>(is);
}

// Total order over OverlapPair used by both the in-RAM sort and the
// k-way merge heap. The task only requires ordering by
// (read_a, read_b, offset_a); the remaining fields are appended so two
// records can never tie unless they're fully identical. That's what
// makes "keep the larger overlap_len on a (read_a, read_b) duplicate,
// first wins on an exact tie" reduction give the same answer regardless
// of which sort/heap implementation walks the data — required for
// merge-twice-is-byte-identical determinism.
bool overlap_less(const OverlapPair& lhs, const OverlapPair& rhs) noexcept {
    if (lhs.read_a != rhs.read_a) return lhs.read_a < rhs.read_a;
    if (lhs.read_b != rhs.read_b) return lhs.read_b < rhs.read_b;
    if (lhs.offset_a != rhs.offset_a) return lhs.offset_a < rhs.offset_a;
    if (lhs.overlap_len != rhs.overlap_len) return lhs.overlap_len < rhs.overlap_len;
    if (lhs.offset_b != rhs.offset_b) return lhs.offset_b < rhs.offset_b;
    if (lhs.diff_count != rhs.diff_count) return lhs.diff_count < rhs.diff_count;
    return lhs.strand < rhs.strand;
}

bool same_pair(const OverlapPair& lhs, const OverlapPair& rhs) noexcept {
    return lhs.read_a == rhs.read_a && lhs.read_b == rhs.read_b;
}

// Collapses a (read_a, read_b, offset_a)-sorted run of records down to
// one record per distinct (read_a, read_b), keeping the largest
// overlap_len (first-encountered wins on an exact tie — well-defined
// because `records` is already totally ordered by overlap_less).
void reduce_sorted_run(const std::vector<OverlapPair>& records, std::vector<OverlapPair>& out) {
    std::size_t i = 0;
    while (i < records.size()) {
        OverlapPair best = records[i];
        std::size_t j = i + 1;
        while (j < records.size() && same_pair(records[j], best)) {
            if (records[j].overlap_len > best.overlap_len) best = records[j];
            ++j;
        }
        out.push_back(best);
        i = j;
    }
}

// Writes `records` (assumed already sorted+deduped) to a fresh BROV1
// file. Returns false if the file cannot be created.
bool write_store(const std::string& path, const std::vector<OverlapPair>& records) {
    OverlapSpillWriter w;
    if (!w.open(path)) return false;
    w.append(std::span<const OverlapPair>(records));
    return w.close();
}

// Deletes a set of paths when it goes out of scope, success or
// exception — used to clean up merge_spills' temporary sorted runs.
class ScopedTempFiles {
public:
    void add(std::string path) { paths_.push_back(std::move(path)); }
    [[nodiscard]] const std::vector<std::string>& paths() const noexcept { return paths_; }

    ~ScopedTempFiles() {
        for (const auto& p : paths_) {
            std::error_code ec;
            std::filesystem::remove(p, ec);
        }
    }

private:
    std::vector<std::string> paths_;
};

// One "current front" record of a k-way merge input, tagged with which
// reader it came from so the merge loop knows what to refill.
struct HeapEntry {
    OverlapPair record;
    std::size_t source_index;
};

struct HeapEntryGreater {
    bool operator()(const HeapEntry& lhs, const HeapEntry& rhs) const noexcept {
        return overlap_less(rhs.record, lhs.record);  // inverted: min-heap
    }
};

// Streams every record out of `readers` in ascending overlap_less order
// (a k-way merge over inputs that are each already internally sorted),
// applying the (read_a, read_b)-dedup reduction as records are popped,
// and appends the result to `writer`. Because the merge visits records
// in the same total order a single global sort would, this produces the
// exact same output as reduce_sorted_run(sort(everything)) would — the
// in-RAM and external paths are cross-consistent, not just
// self-consistent.
void k_way_merge_into(std::vector<OverlapStoreReader>& readers, OverlapSpillWriter& writer) {
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapEntryGreater> heap;
    for (std::size_t i = 0; i < readers.size(); ++i) {
        OverlapPair rec;
        if (readers[i].next(rec)) heap.push(HeapEntry{rec, i});
    }

    bool have_pending = false;
    OverlapPair pending{};
    while (!heap.empty()) {
        HeapEntry top = heap.top();
        heap.pop();

        OverlapPair refill;
        if (readers[top.source_index].next(refill)) {
            heap.push(HeapEntry{refill, top.source_index});
        }

        if (have_pending && same_pair(pending, top.record)) {
            if (top.record.overlap_len > pending.overlap_len) pending = top.record;
        } else {
            if (have_pending) writer.append(pending);
            pending = top.record;
            have_pending = true;
        }
    }
    if (have_pending) writer.append(pending);
}

}  // namespace

// ---------------------------------------------------------------------
// OverlapSpillWriter
// ---------------------------------------------------------------------

OverlapSpillWriter::~OverlapSpillWriter() { close(); }

OverlapSpillWriter::OverlapSpillWriter(OverlapSpillWriter&& other) noexcept
    : file_(std::move(other.file_)),
      buffer_(std::move(other.buffer_)),
      record_count_(other.record_count_) {
    other.record_count_ = 0;
}

OverlapSpillWriter& OverlapSpillWriter::operator=(OverlapSpillWriter&& other) noexcept {
    if (this != &other) {
        close();
        file_ = std::move(other.file_);
        buffer_ = std::move(other.buffer_);
        record_count_ = other.record_count_;
        other.record_count_ = 0;
    }
    return *this;
}

bool OverlapSpillWriter::open(const std::string& path) {
    close();
    file_.open(path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!file_.is_open()) return false;
    file_.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
    write_u64(file_, 0);  // placeholder record count, patched at close()
    if (!file_) {
        file_.close();
        return false;
    }
    buffer_.clear();
    buffer_.reserve(kMinBufferRecords);
    record_count_ = 0;
    return true;
}

void OverlapSpillWriter::append(const OverlapPair& pair) {
    buffer_.push_back(pair);
    ++record_count_;
    if (buffer_.size() >= kMinBufferRecords) flush_buffer();
}

void OverlapSpillWriter::append(std::span<const OverlapPair> pairs) {
    for (const auto& p : pairs) append(p);
}

void OverlapSpillWriter::flush_buffer() {
    if (buffer_.empty()) return;
    file_.write(reinterpret_cast<const char*>(buffer_.data()),
                static_cast<std::streamsize>(buffer_.size() * sizeof(OverlapPair)));
    buffer_.clear();
}

bool OverlapSpillWriter::close() {
    if (!file_.is_open()) return true;
    flush_buffer();
    file_.seekp(kRecordCountOffset);
    write_u64(file_, record_count_);
    file_.flush();
    const bool ok = static_cast<bool>(file_);
    file_.close();
    return ok;
}

// ---------------------------------------------------------------------
// OverlapStoreReader
// ---------------------------------------------------------------------

OverlapStoreReader::OverlapStoreReader(OverlapStoreReader&& other) noexcept
    : file_(std::move(other.file_)),
      record_count_(other.record_count_),
      records_read_(other.records_read_),
      last_error_(std::move(other.last_error_)) {
    other.record_count_ = 0;
    other.records_read_ = 0;
}

OverlapStoreReader& OverlapStoreReader::operator=(OverlapStoreReader&& other) noexcept {
    if (this != &other) {
        close();
        file_ = std::move(other.file_);
        record_count_ = other.record_count_;
        records_read_ = other.records_read_;
        last_error_ = std::move(other.last_error_);
        other.record_count_ = 0;
        other.records_read_ = 0;
    }
    return *this;
}

bool OverlapStoreReader::open(const std::string& path) {
    close();
    file_.open(path, std::ios::binary | std::ios::in);
    if (!file_.is_open()) {
        last_error_ = "OverlapStoreReader: cannot open file: " + path;
        return false;
    }

    std::array<char, kMagic.size()> magic{};
    file_.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    std::uint64_t declared_count = 0;
    if (!file_ || !read_u64(file_, declared_count)) {
        last_error_ = "OverlapStoreReader: truncated header in " + path;
        file_.close();
        return false;
    }
    if (!std::equal(magic.begin(), magic.end(), kMagic.begin())) {
        last_error_ = "OverlapStoreReader: bad magic in " + path + " (not a BROV1 file)";
        file_.close();
        return false;
    }

    std::error_code ec;
    const auto file_size = std::filesystem::file_size(path, ec);
    if (!ec) {
        const auto expected = static_cast<std::uintmax_t>(kHeaderSize)
            + static_cast<std::uintmax_t>(declared_count) * sizeof(OverlapPair);
        if (file_size < expected) {
            last_error_ = "OverlapStoreReader: truncated file " + path + ": header claims "
                + std::to_string(declared_count) + " records but file is too short";
            file_.close();
            return false;
        }
    }

    record_count_ = declared_count;
    records_read_ = 0;
    last_error_.clear();
    return true;
}

void OverlapStoreReader::close() {
    file_.close();
    file_.clear();
    record_count_ = 0;
    records_read_ = 0;
}

bool OverlapStoreReader::next(OverlapPair& out) {
    if (!file_.is_open() || records_read_ >= record_count_) return false;
    file_.read(reinterpret_cast<char*>(&out), sizeof(OverlapPair));
    if (!file_) {
        last_error_ = "OverlapStoreReader: truncated record data (expected "
            + std::to_string(record_count_) + " records, got "
            + std::to_string(records_read_) + ")";
        return false;
    }
    ++records_read_;
    return true;
}

// ---------------------------------------------------------------------
// merge_spills
// ---------------------------------------------------------------------

void merge_spills(const std::vector<std::string>& spill_paths,
                   const std::string& out_path,
                   const std::string& tmp_dir,
                   std::uint64_t max_ram_bytes) {
    // Opening each spill is O(1) (header-only) and gives an exact total
    // record count without a preliminary full scan.
    std::vector<OverlapStoreReader> spills(spill_paths.size());
    std::uint64_t total_records = 0;
    for (std::size_t i = 0; i < spill_paths.size(); ++i) {
        if (!spills[i].open(spill_paths[i])) {
            throw std::runtime_error(spills[i].last_error());
        }
        total_records += spills[i].record_count();
    }
    const std::uint64_t total_bytes = total_records * sizeof(OverlapPair);

    if (total_bytes <= max_ram_bytes) {
        std::vector<OverlapPair> all;
        all.reserve(static_cast<std::size_t>(total_records));
        OverlapPair rec;
        for (auto& reader : spills) {
            while (reader.next(rec)) all.push_back(rec);
            if (!reader.last_error().empty()) throw std::runtime_error(reader.last_error());
        }
        std::sort(all.begin(), all.end(), overlap_less);
        std::vector<OverlapPair> deduped;
        deduped.reserve(all.size());
        reduce_sorted_run(all, deduped);
        if (!write_store(out_path, deduped)) {
            throw std::runtime_error("merge_spills: cannot create output file: " + out_path);
        }
        return;
    }

    // External path: split the concatenated input into sorted runs of
    // up to max_ram_bytes each (spilled under tmp_dir), then k-way merge
    // those runs. Dedup happens once, during the k-way merge itself —
    // a (read_a, read_b) pair split across two runs would survive an
    // in-run-only dedup, so it can't happen per-chunk.
    std::error_code dir_ec;
    std::filesystem::create_directories(tmp_dir, dir_ec);

    const std::size_t chunk_capacity =
        std::max<std::size_t>(1, static_cast<std::size_t>(max_ram_bytes / sizeof(OverlapPair)));

    ScopedTempFiles temp_runs;
    std::vector<OverlapPair> chunk;
    chunk.reserve(chunk_capacity);
    std::size_t run_index = 0;
    auto flush_chunk = [&]() {
        if (chunk.empty()) return;
        std::sort(chunk.begin(), chunk.end(), overlap_less);
        std::string run_path =
            tmp_dir + "/overlap_merge_run_" + std::to_string(run_index++) + ".BROV1.tmp";
        if (!write_store(run_path, chunk)) {
            throw std::runtime_error("merge_spills: cannot create temp run file: " + run_path);
        }
        temp_runs.add(run_path);
        chunk.clear();
    };

    OverlapPair rec;
    for (auto& reader : spills) {
        while (reader.next(rec)) {
            chunk.push_back(rec);
            if (chunk.size() >= chunk_capacity) flush_chunk();
        }
        if (!reader.last_error().empty()) throw std::runtime_error(reader.last_error());
    }
    flush_chunk();

    std::vector<OverlapStoreReader> runs(temp_runs.paths().size());
    for (std::size_t i = 0; i < runs.size(); ++i) {
        if (!runs[i].open(temp_runs.paths()[i])) {
            throw std::runtime_error(runs[i].last_error());
        }
    }

    OverlapSpillWriter writer;
    if (!writer.open(out_path)) {
        throw std::runtime_error("merge_spills: cannot create output file: " + out_path);
    }
    k_way_merge_into(runs, writer);
    if (!writer.close()) {
        throw std::runtime_error("merge_spills: failed writing output file: " + out_path);
    }
}

}  // namespace branch::io
