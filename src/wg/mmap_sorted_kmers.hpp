// BRANCH — memory-mapped view over the sorted (key, count) records of a
// Phase-0 .bin dump.
//
// Phase 1.1 het-pair detection needs only (a) sequential iteration over the
// recurrent k-mers and (b) count lookup of arbitrary 1-bp-neighbour k-mers.
// Both are served by the on-disk sorted-by-key records: iterate the span, or
// binary-search it. Replacing the former in-RAM std::unordered_map (~48 B per
// entry → ~95 GB for whole-genome, the v08 OOM) with this mmap view drops the
// process-heap footprint to ~0 — the 12 B/entry file stays in the OS page
// cache, shared and reclaimable.
//
// .bin layout (see write_phase0_dump): MAGIC[8] | n_recurrent u64 | n_reads u64
// | n_bases u64 | cov i32 | het i32 | z f64 | k u32 | hist_n u32 |
// histogram[hist_n] u64 | records[n_recurrent] {key u64, count u32} (packed 12 B).

#pragma once

#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>

namespace branch::wg {

class MmapSortedKmers {
public:
    MmapSortedKmers() = default;
    ~MmapSortedKmers() { close(); }

    MmapSortedKmers(const MmapSortedKmers&) = delete;
    MmapSortedKmers& operator=(const MmapSortedKmers&) = delete;

    static constexpr std::size_t kRecBytes = 12;  // u64 key + u32 count, packed

    // mmap the records region of `path`. Returns false on any error.
    bool open(const std::string& path) {
        close();
        fd_ = ::open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return false;
        struct stat st {};
        if (::fstat(fd_, &st) != 0) { close(); return false; }
        map_len_ = static_cast<std::size_t>(st.st_size);
        if (map_len_ < 56) { close(); return false; }
        void* m = ::mmap(nullptr, map_len_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (m == MAP_FAILED) { map_ = nullptr; close(); return false; }
        map_ = m;
        const char* base = static_cast<const char*>(map_);

        if (std::memcmp(base, "KHIST05", 7) != 0) { close(); return false; }
        std::uint64_t n_recurrent = 0;
        std::memcpy(&n_recurrent, base + 8, sizeof(n_recurrent));
        std::uint32_t hist_n = 0;
        std::memcpy(&hist_n, base + 52, sizeof(hist_n));  // 8 + 44 header

        const std::size_t recs_off =
            56 + static_cast<std::size_t>(hist_n) * sizeof(std::uint64_t);
        if (recs_off + n_recurrent * kRecBytes > map_len_) { close(); return false; }
        recs_ = base + recs_off;
        n_ = static_cast<std::size_t>(n_recurrent);
        // The records are scanned (iterate) and randomly probed (binary
        // search); hint the kernel accordingly.
        ::madvise(const_cast<char*>(recs_), n_ * kRecBytes, MADV_WILLNEED);
        return true;
    }

    void close() {
        if (map_) { ::munmap(map_, map_len_); map_ = nullptr; }
        if (fd_ >= 0) { ::close(fd_); fd_ = -1; }
        recs_ = nullptr; n_ = 0; map_len_ = 0;
    }

    [[nodiscard]] bool empty() const noexcept { return n_ == 0; }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }

    [[nodiscard]] std::uint64_t key_at(std::size_t i) const noexcept {
        std::uint64_t k;
        std::memcpy(&k, recs_ + i * kRecBytes, sizeof(k));
        return k;
    }
    [[nodiscard]] std::uint32_t count_at(std::size_t i) const noexcept {
        std::uint32_t c;
        std::memcpy(&c, recs_ + i * kRecBytes + sizeof(std::uint64_t), sizeof(c));
        return c;
    }

    // Count for `key`, or 0 if absent (recurrent counts are always ≥ 2, so 0
    // is an unambiguous "not present"). Binary search on the sorted keys.
    [[nodiscard]] std::uint32_t find(std::uint64_t key) const noexcept {
        std::size_t lo = 0, hi = n_;
        while (lo < hi) {
            std::size_t mid = lo + (hi - lo) / 2;
            std::uint64_t mk = key_at(mid);
            if (mk < key) lo = mid + 1;
            else if (mk > key) hi = mid;
            else return count_at(mid);
        }
        return 0;
    }

private:
    int         fd_ = -1;
    void*       map_ = nullptr;
    std::size_t map_len_ = 0;
    const char* recs_ = nullptr;
    std::size_t n_ = 0;
};

}  // namespace branch::wg
