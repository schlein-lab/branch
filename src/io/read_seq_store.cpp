#include "io/read_seq_store.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace branch::io {

namespace {

constexpr std::array<char, 5> kMagic = {'B', 'R', 'S', 'Q', '1'};
constexpr std::array<char, 4> kBaseChars = {'A', 'C', 'G', 'T'};

// On-disk record sizes (deliberately unpadded — these describe bytes
// written/read one field at a time, not a struct layout).
constexpr std::uint64_t kIndexRecordSize = sizeof(std::uint64_t) + sizeof(std::uint32_t);  // 12
constexpr std::uint64_t kExceptionRecordSize = 2 * sizeof(std::uint32_t);                  // 8
constexpr std::uint64_t kFooterSize =
    static_cast<std::uint64_t>(kMagic.size()) + 4 * sizeof(std::uint64_t);  // 37

void write_u32(std::ostream& os, std::uint32_t v) {
    os.write(reinterpret_cast<const char*>(&v), sizeof(v));
}

void write_u64(std::ostream& os, std::uint64_t v) {
    os.write(reinterpret_cast<const char*>(&v), sizeof(v));
}

bool read_u32(std::istream& is, std::uint32_t& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(v));
    return static_cast<bool>(is);
}

bool read_u64(std::istream& is, std::uint64_t& v) {
    is.read(reinterpret_cast<char*>(&v), sizeof(v));
    return static_cast<bool>(is);
}

}  // namespace

// ---------------------------------------------------------------------
// ReadSeqStoreWriter
// ---------------------------------------------------------------------

ReadSeqStoreWriter::~ReadSeqStoreWriter() { close(); }

ReadSeqStoreWriter::ReadSeqStoreWriter(ReadSeqStoreWriter&& other) noexcept
    : file_(std::move(other.file_)),
      index_(std::move(other.index_)),
      exceptions_(std::move(other.exceptions_)),
      write_offset_(other.write_offset_) {
    other.write_offset_ = 0;
}

ReadSeqStoreWriter& ReadSeqStoreWriter::operator=(ReadSeqStoreWriter&& other) noexcept {
    if (this != &other) {
        close();
        file_ = std::move(other.file_);
        index_ = std::move(other.index_);
        exceptions_ = std::move(other.exceptions_);
        write_offset_ = other.write_offset_;
        other.write_offset_ = 0;
    }
    return *this;
}

bool ReadSeqStoreWriter::open(const std::string& path) {
    close();
    file_.open(path, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!file_.is_open()) return false;
    index_.clear();
    exceptions_.clear();
    write_offset_ = 0;
    return true;
}

std::uint32_t ReadSeqStoreWriter::add(const std::string& seq) {
    const auto read_id = static_cast<std::uint32_t>(index_.size());
    const auto length = static_cast<std::uint32_t>(seq.size());
    const auto packed_bytes = static_cast<std::uint32_t>((length + 3) / 4);

    if (packed_bytes > 0) {
        std::vector<std::uint8_t> packed(packed_bytes, 0);
        for (std::uint32_t i = 0; i < length; ++i) {
            std::uint8_t code = 0;
            switch (seq[i]) {
                case 'A': code = 0; break;
                case 'C': code = 1; break;
                case 'G': code = 2; break;
                case 'T': code = 3; break;
                default: exceptions_.push_back({read_id, i}); break;
            }
            const auto shift = static_cast<unsigned>((i % 4) * 2);
            packed[i / 4] =
                static_cast<std::uint8_t>(packed[i / 4] | (static_cast<unsigned>(code) << shift));
        }
        file_.write(reinterpret_cast<const char*>(packed.data()),
                    static_cast<std::streamsize>(packed_bytes));
    }

    index_.push_back({write_offset_, length});
    write_offset_ += packed_bytes;
    return read_id;
}

bool ReadSeqStoreWriter::close() {
    if (!file_.is_open()) return true;

    const std::uint64_t exceptions_offset = write_offset_;
    for (const auto& e : exceptions_) {
        write_u32(file_, e.read_id);
        write_u32(file_, e.pos);
    }

    const std::uint64_t index_offset =
        exceptions_offset + static_cast<std::uint64_t>(exceptions_.size()) * kExceptionRecordSize;
    for (const auto& entry : index_) {
        write_u64(file_, entry.offset);
        write_u32(file_, entry.length);
    }

    file_.write(kMagic.data(), static_cast<std::streamsize>(kMagic.size()));
    write_u64(file_, index_offset);
    write_u64(file_, exceptions_offset);
    write_u64(file_, static_cast<std::uint64_t>(exceptions_.size()));
    write_u64(file_, static_cast<std::uint64_t>(index_.size()));

    file_.flush();
    const bool ok = static_cast<bool>(file_);
    file_.close();
    return ok;
}

// ---------------------------------------------------------------------
// ReadSeqStoreReader
// ---------------------------------------------------------------------

ReadSeqStoreReader::ReadSeqStoreReader(ReadSeqStoreReader&& other) noexcept
    : file_(std::move(other.file_)),
      index_(std::move(other.index_)),
      exceptions_(std::move(other.exceptions_)),
      last_error_(std::move(other.last_error_)) {}

ReadSeqStoreReader& ReadSeqStoreReader::operator=(ReadSeqStoreReader&& other) noexcept {
    if (this != &other) {
        close();
        file_ = std::move(other.file_);
        index_ = std::move(other.index_);
        exceptions_ = std::move(other.exceptions_);
        last_error_ = std::move(other.last_error_);
    }
    return *this;
}

bool ReadSeqStoreReader::open(const std::string& path) {
    close();
    file_.open(path, std::ios::binary | std::ios::in);
    if (!file_.is_open()) {
        last_error_ = "ReadSeqStoreReader: cannot open file: " + path;
        return false;
    }

    file_.seekg(0, std::ios::end);
    const std::streamoff file_size = file_.tellg();
    if (file_size < 0 || static_cast<std::uint64_t>(file_size) < kFooterSize) {
        last_error_ = "ReadSeqStoreReader: truncated file (smaller than footer): " + path;
        file_.close();
        return false;
    }
    const auto ufile_size = static_cast<std::uint64_t>(file_size);

    file_.seekg(file_size - static_cast<std::streamoff>(kFooterSize));
    std::array<char, kMagic.size()> magic{};
    file_.read(magic.data(), static_cast<std::streamsize>(magic.size()));
    std::uint64_t index_offset = 0;
    std::uint64_t exceptions_offset = 0;
    std::uint64_t exception_count = 0;
    std::uint64_t declared_read_count = 0;
    if (!file_ || !read_u64(file_, index_offset) || !read_u64(file_, exceptions_offset) ||
        !read_u64(file_, exception_count) || !read_u64(file_, declared_read_count)) {
        last_error_ = "ReadSeqStoreReader: truncated footer in " + path;
        file_.close();
        return false;
    }
    if (!std::equal(magic.begin(), magic.end(), kMagic.begin())) {
        last_error_ = "ReadSeqStoreReader: bad magic in " + path + " (not a BRSQ1 file)";
        file_.close();
        return false;
    }

    const std::uint64_t expected_index_offset =
        exceptions_offset + exception_count * kExceptionRecordSize;
    const std::uint64_t expected_file_size =
        index_offset + declared_read_count * kIndexRecordSize + kFooterSize;
    if (exceptions_offset > ufile_size || index_offset != expected_index_offset ||
        expected_file_size != ufile_size) {
        last_error_ = "ReadSeqStoreReader: corrupt or truncated layout in " + path;
        file_.close();
        return false;
    }

    index_.clear();
    index_.reserve(static_cast<std::size_t>(declared_read_count));
    file_.seekg(static_cast<std::streamoff>(index_offset));
    for (std::uint64_t i = 0; i < declared_read_count; ++i) {
        std::uint64_t offset = 0;
        std::uint32_t length = 0;
        if (!read_u64(file_, offset) || !read_u32(file_, length)) {
            last_error_ = "ReadSeqStoreReader: truncated index in " + path;
            file_.close();
            index_.clear();
            return false;
        }
        index_.push_back({offset, length});
    }

    exceptions_.clear();
    exceptions_.reserve(static_cast<std::size_t>(exception_count));
    file_.seekg(static_cast<std::streamoff>(exceptions_offset));
    for (std::uint64_t i = 0; i < exception_count; ++i) {
        std::uint32_t rid = 0;
        std::uint32_t pos = 0;
        if (!read_u32(file_, rid) || !read_u32(file_, pos)) {
            last_error_ = "ReadSeqStoreReader: truncated exceptions in " + path;
            file_.close();
            index_.clear();
            exceptions_.clear();
            return false;
        }
        exceptions_.push_back({rid, pos});
    }

    last_error_.clear();
    return true;
}

void ReadSeqStoreReader::close() {
    file_.close();
    file_.clear();
    index_.clear();
    exceptions_.clear();
}

std::uint32_t ReadSeqStoreReader::sequence_length(std::uint32_t read_id) const {
    if (read_id >= index_.size()) {
        throw std::out_of_range("ReadSeqStoreReader::sequence_length: read_id " +
                                 std::to_string(read_id) + " out of range (have " +
                                 std::to_string(index_.size()) + " reads)");
    }
    return index_[read_id].length;
}

std::string ReadSeqStoreReader::fetch(std::uint32_t read_id) {
    return fetch_window(read_id, 0, sequence_length(read_id));
}

std::string ReadSeqStoreReader::fetch_window(std::uint32_t read_id, std::uint32_t start,
                                              std::uint32_t len) {
    if (read_id >= index_.size()) {
        throw std::out_of_range("ReadSeqStoreReader::fetch_window: read_id " +
                                 std::to_string(read_id) + " out of range (have " +
                                 std::to_string(index_.size()) + " reads)");
    }
    const IndexEntry entry = index_[read_id];
    if (start > entry.length) {
        throw std::out_of_range("ReadSeqStoreReader::fetch_window: start " +
                                 std::to_string(start) + " past sequence length " +
                                 std::to_string(entry.length) + " for read_id " +
                                 std::to_string(read_id));
    }
    const std::uint32_t avail = entry.length - start;
    const std::uint32_t clipped_len = std::min(len, avail);

    std::string result(clipped_len, 'A');
    if (clipped_len == 0) return result;

    const std::uint32_t first_base_byte = start / 4;
    const std::uint32_t last_base_index = start + clipped_len - 1;
    const std::uint32_t last_base_byte = last_base_index / 4;
    const std::uint32_t bytes_to_read = last_base_byte - first_base_byte + 1;

    std::vector<std::uint8_t> packed(bytes_to_read);
    file_.seekg(static_cast<std::streamoff>(entry.offset + first_base_byte));
    file_.read(reinterpret_cast<char*>(packed.data()), static_cast<std::streamsize>(bytes_to_read));
    if (!file_) {
        throw std::runtime_error("ReadSeqStoreReader::fetch_window: read failure for read_id " +
                                  std::to_string(read_id));
    }

    for (std::uint32_t i = 0; i < clipped_len; ++i) {
        const std::uint32_t global_base = start + i;
        const std::uint32_t byte_index = global_base / 4 - first_base_byte;
        const auto bit_shift = static_cast<unsigned>((global_base % 4) * 2);
        const unsigned code = (static_cast<unsigned>(packed[byte_index]) >> bit_shift) & 0x3u;
        result[i] = kBaseChars[code];
    }

    const auto lo = std::lower_bound(
        exceptions_.begin(), exceptions_.end(), read_id,
        [](const ExceptionEntry& e, std::uint32_t rid) { return e.read_id < rid; });
    const auto hi = std::upper_bound(
        exceptions_.begin(), exceptions_.end(), read_id,
        [](std::uint32_t rid, const ExceptionEntry& e) { return rid < e.read_id; });
    for (auto it = lo; it != hi; ++it) {
        if (it->pos >= start && it->pos < start + clipped_len) {
            result[it->pos - start] = 'N';
        }
    }

    return result;
}

}  // namespace branch::io
