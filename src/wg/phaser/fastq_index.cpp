#include "fastq_index.hpp"

#include <zlib.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <sys/types.h>
#include <unistd.h>

namespace branch::wg::phaser {

struct FastqIndex::Reader {
    std::FILE* fp = nullptr;
};

FastqIndex::~FastqIndex() {
    close_reader();
}

void FastqIndex::open_reader() {
    if (!reader_) reader_ = new Reader();
    if (!reader_->fp) {
        reader_->fp = std::fopen(fastq_path_.c_str(), "rb");
        if (!reader_->fp) {
            throw std::runtime_error("FastqIndex: cannot open " + fastq_path_);
        }
    }
}

void FastqIndex::close_reader() {
    if (reader_) {
        if (reader_->fp) { std::fclose(reader_->fp); reader_->fp = nullptr; }
        delete reader_;
        reader_ = nullptr;
    }
}

void FastqIndex::rewind() {
    open_reader();
    std::fseek(reader_->fp, 0, SEEK_SET);
}

void FastqIndex::seek_to(std::uint64_t offset) {
    open_reader();
    std::fseek(reader_->fp, static_cast<long>(offset), SEEK_SET);
}

namespace {

bool ends_with(const std::string& s, const char* suffix) {
    std::size_t n = std::strlen(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

// Decompress gz to scratch tempfile. Returns path of the decompressed
// file. Caller is responsible for cleanup (or rely on /tmp evictio).
// Probe a directory by attempting to create a 0-byte file there; returns
// true if writable. Used to skip $TMPDIR when SLURM points it at an
// unmounted /beegfs/tmp/... path that the compute node can't access.
bool dir_is_writable(const std::string& dir) {
    if (dir.empty()) return false;
    const std::string probe = dir + "/.fqx_probe_" +
                              std::to_string(::getpid());
    std::FILE* f = std::fopen(probe.c_str(), "wb");
    if (!f) return false;
    std::fclose(f);
    std::remove(probe.c_str());
    return true;
}

std::string gunzip_to_scratch(const std::string& gz_path) {
    // Resolve scratch dir: explicit env > SLURM TMPDIR > /tmp.
    // Each candidate is probed for writability before use; SLURM at our
    // HPC site sets TMPDIR to a /beegfs/tmp/... path that is sometimes
    // unwritable from the compute node, which cost v11 a 13-min run.
    std::vector<std::string> candidates;
    if (const char* p = std::getenv("FASTQ_INDEX_SCRATCH")) candidates.push_back(p);
    if (const char* p = std::getenv("TMPDIR")) candidates.push_back(p);
    candidates.push_back("/tmp");
    std::string base;
    for (const auto& c : candidates) {
        if (dir_is_writable(c)) { base = c; break; }
    }
    if (base.empty()) {
        throw std::runtime_error(
            "FastqIndex: no writable scratch dir found (tried "
            "FASTQ_INDEX_SCRATCH, TMPDIR, /tmp)");
    }
    std::string out_path = base + "/fastq_index_" +
                           std::to_string(::getpid()) + ".fastq";
    std::cerr << "[fastq_index] using scratch " << base << "\n";
    std::fflush(stderr);
    std::cerr << "[fastq_index] decompressing " << gz_path
              << " -> " << out_path << "\n";
    std::fflush(stderr);

    auto t0 = std::chrono::steady_clock::now();
    gzFile gz = gzopen(gz_path.c_str(), "rb");
    if (!gz) throw std::runtime_error("gzopen failed: " + gz_path);
    std::FILE* out = std::fopen(out_path.c_str(), "wb");
    if (!out) {
        gzclose(gz);
        throw std::runtime_error("fopen failed: " + out_path);
    }
    constexpr std::size_t BUFSIZE = 1 << 20;
    std::vector<char> buf(BUFSIZE);
    int n;
    std::uint64_t total = 0;
    while ((n = gzread(gz, buf.data(), BUFSIZE)) > 0) {
        if (std::fwrite(buf.data(), 1, static_cast<std::size_t>(n), out) !=
            static_cast<std::size_t>(n)) {
            std::fclose(out); gzclose(gz);
            throw std::runtime_error("fwrite failed during decompress");
        }
        total += static_cast<std::uint64_t>(n);
    }
    if (n < 0) {
        std::fclose(out); gzclose(gz);
        throw std::runtime_error("gzread error during decompress");
    }
    std::fclose(out);
    gzclose(gz);
    auto t1 = std::chrono::steady_clock::now();
    const double sec = std::chrono::duration<double>(t1 - t0).count();
    std::cerr << "[fastq_index] decompressed "
              << static_cast<double>(total) / (1024.0 * 1024.0 * 1024.0)
              << " GiB in " << sec << "s\n";
    std::fflush(stderr);
    return out_path;
}

}  // namespace

void FastqIndex::build(const std::string& fastq_path, bool assign_ids) {
    (void)assign_ids;  // sequential ids only in v1; reserved for future
    auto t0 = std::chrono::steady_clock::now();

    if (ends_with(fastq_path, ".gz")) {
        fastq_path_ = gunzip_to_scratch(fastq_path);
        is_gzipped_ = true;
    } else {
        fastq_path_ = fastq_path;
    }

    open_reader();
    std::FILE* fp = reader_->fp;
    std::fseek(fp, 0, SEEK_SET);

    constexpr std::size_t LBUF = 1 << 20;
    std::vector<char> line(LBUF);
    int line_idx = 0;
    std::string current_id;

    auto t_hb = t0;
    while (true) {
        std::uint64_t before = static_cast<std::uint64_t>(std::ftell(fp));
        if (!std::fgets(line.data(), LBUF, fp)) break;
        std::size_t len = std::strlen(line.data());
        if (len > 0 && line[len - 1] == '\n') --len;

        switch (line_idx) {
        case 0:  // header
            if (len < 1 || line[0] != '@') {
                throw std::runtime_error(
                    "FastqIndex: expected header at offset " +
                    std::to_string(before));
            }
            current_id.assign(line.data() + 1, len - 1);
            for (std::size_t i = 0; i < current_id.size(); ++i) {
                if (current_id[i] == ' ' || current_id[i] == '\t') {
                    current_id.resize(i);
                    break;
                }
            }
            break;
        case 1: {  // seq
            Entry e;
            e.seq_offset = before;
            e.seq_length = static_cast<std::uint32_t>(len);
            entries_.push_back(e);
            id_strings_.push_back(current_id);

            if (entries_.size() % 500'000 == 0) {
                auto now = std::chrono::steady_clock::now();
                const double el = std::chrono::duration<double>(now - t0).count();
                std::cerr << "[fastq_index] indexed " << entries_.size()
                          << " reads (elapsed=" << el << "s)\n";
                std::fflush(stderr);
                t_hb = now;
            }
            break;
        }
        case 2:  // '+'
        case 3:  // quals
        default:
            break;
        }
        line_idx = (line_idx + 1) % 4;
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[fastq_index] done: " << entries_.size() << " reads in "
              << std::chrono::duration<double>(t1 - t0).count() << "s ("
              << index_bytes() / (1024 * 1024) << " MB index)\n";
    std::fflush(stderr);
}

void FastqIndex::read_seq(ReadId rid, std::string& out) {
    if (in_memory_) {
        if (rid >= in_memory_reads_.size()) { out.clear(); return; }
        out = in_memory_reads_[rid].second;
        return;
    }
    if (rid >= entries_.size()) {
        out.clear();
        return;
    }
    open_reader();
    seek_to(entries_[rid].seq_offset);
    out.resize(entries_[rid].seq_length);
    std::size_t got = std::fread(out.data(), 1, out.size(), reader_->fp);
    if (got != out.size()) out.resize(got);
}

void FastqIndex::build_from_memory(
    std::vector<std::pair<ReadId, std::string>> reads_vec,
    std::vector<std::string> id_strings)
{
    in_memory_ = true;
    in_memory_reads_ = std::move(reads_vec);
    if (id_strings.empty()) {
        id_strings_.resize(in_memory_reads_.size());
        for (std::size_t i = 0; i < in_memory_reads_.size(); ++i) {
            id_strings_[i] = std::to_string(in_memory_reads_[i].first);
        }
    } else {
        id_strings_ = std::move(id_strings);
    }
}

FastqIndex::ThreadReader::~ThreadReader() {
    if (fp) std::fclose(fp);
}

FastqIndex::ThreadReader::ThreadReader(ThreadReader&& o) noexcept
    : fp(o.fp), in_memory(o.in_memory), parent(o.parent) {
    o.fp = nullptr; o.parent = nullptr;
}

FastqIndex::ThreadReader& FastqIndex::ThreadReader::operator=(
    ThreadReader&& o) noexcept
{
    if (this == &o) return *this;
    if (fp) std::fclose(fp);
    fp = o.fp;        o.fp = nullptr;
    in_memory = o.in_memory;
    parent = o.parent; o.parent = nullptr;
    return *this;
}

FastqIndex::ThreadReader FastqIndex::open_thread_reader() const {
    ThreadReader r;
    r.parent = this;
    r.in_memory = in_memory_;
    if (!in_memory_) {
        r.fp = std::fopen(fastq_path_.c_str(), "rb");
        if (!r.fp) {
            throw std::runtime_error(
                "FastqIndex::open_thread_reader: cannot open " + fastq_path_);
        }
    }
    return r;
}

void FastqIndex::read_seq_via(ThreadReader& h, ReadId rid,
                              std::string& out) const
{
    if (h.in_memory) {
        if (rid >= in_memory_reads_.size()) { out.clear(); return; }
        out = in_memory_reads_[rid].second;
        return;
    }
    if (rid >= entries_.size() || !h.fp) { out.clear(); return; }
    std::fseek(h.fp, static_cast<long>(entries_[rid].seq_offset), SEEK_SET);
    out.resize(entries_[rid].seq_length);
    std::size_t got = std::fread(out.data(), 1, out.size(), h.fp);
    if (got != out.size()) out.resize(got);
}

bool FastqIndex::read_next(std::string& seq) {
    open_reader();
    constexpr std::size_t LBUF = 1 << 20;
    static thread_local std::vector<char> line(LBUF);

    if (!std::fgets(line.data(), LBUF, reader_->fp)) return false;  // hdr
    if (!std::fgets(line.data(), LBUF, reader_->fp)) return false;  // seq
    std::size_t len = std::strlen(line.data());
    if (len > 0 && line[len - 1] == '\n') --len;
    seq.assign(line.data(), len);
    if (!std::fgets(line.data(), LBUF, reader_->fp)) return false;  // +
    if (!std::fgets(line.data(), LBUF, reader_->fp)) return false;  // qual
    return true;
}

}  // namespace branch::wg::phaser
