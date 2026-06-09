#include "wg/external_kmer_counter.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <queue>
#include <stdexcept>

namespace branch::wg {

namespace {
constexpr std::uint32_t kCountMax = std::numeric_limits<std::uint32_t>::max();

struct Record {
    std::uint64_t key;
    std::uint32_t count;
};
}  // namespace

ExternalKmerCounter::ExternalKmerCounter(std::size_t max_ram_bytes,
                                         std::string scratch_dir)
    : cap_(std::max<std::size_t>(max_ram_bytes / sizeof(std::uint64_t),
                                 1u << 20)),  // floor 1M keys (8 MB)
      scratch_(std::move(scratch_dir)) {
    buf_.reserve(cap_);
}

ExternalKmerCounter::~ExternalKmerCounter() {
    for (const auto& p : runs_) std::remove(p.c_str());
}

void ExternalKmerCounter::add(std::uint64_t key) {
    buf_.push_back(key);
    if (buf_.size() >= cap_) spill();
}

void ExternalKmerCounter::spill() {
    if (buf_.empty()) return;
    std::sort(buf_.begin(), buf_.end());

    std::string path = scratch_ + "/branch_kmc_run_" +
                       std::to_string(run_seq_++) + ".bin";
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("ExternalKmerCounter: cannot open run " + path);

    // Run-length encode the sorted buffer into (key, count) records.
    std::size_t i = 0;
    const std::size_t n = buf_.size();
    while (i < n) {
        const std::uint64_t k = buf_[i];
        std::size_t j = i + 1;
        while (j < n && buf_[j] == k) ++j;
        Record r;
        r.key = k;
        r.count = static_cast<std::uint32_t>(
            std::min<std::size_t>(j - i, kCountMax));
        f.write(reinterpret_cast<const char*>(&r), sizeof(r));
        i = j;
    }
    f.close();
    runs_.push_back(std::move(path));
    buf_.clear();
}

void ExternalKmerCounter::finalize(std::vector<std::uint64_t>& out_keys,
                                   std::vector<std::uint32_t>& out_counts) {
    out_keys.clear();
    out_counts.clear();

    // Fast path: everything fits in one buffer, never spilled → sort + RLE.
    if (runs_.empty()) {
        std::sort(buf_.begin(), buf_.end());
        std::size_t i = 0;
        const std::size_t n = buf_.size();
        while (i < n) {
            const std::uint64_t k = buf_[i];
            std::size_t j = i + 1;
            while (j < n && buf_[j] == k) ++j;
            out_keys.push_back(k);
            out_counts.push_back(static_cast<std::uint32_t>(
                std::min<std::size_t>(j - i, kCountMax)));
            i = j;
        }
        std::vector<std::uint64_t>().swap(buf_);
        return;
    }

    // Spill any residual buffer so the merge sees a uniform set of runs.
    spill();

    // k-way merge: each run is sorted ascending by key; sum counts for equal
    // keys across runs. Min-heap keyed by the runs' current record key.
    const std::size_t R = runs_.size();
    std::vector<std::ifstream> in(R);
    std::vector<Record> cur(R);
    std::vector<bool> live(R, false);

    auto load = [&](std::size_t r) -> bool {
        Record rec;
        if (in[r].read(reinterpret_cast<char*>(&rec), sizeof(rec))) {
            cur[r] = rec;
            return true;
        }
        return false;
    };

    using HeapItem = std::pair<std::uint64_t, std::size_t>;  // (key, run)
    std::priority_queue<HeapItem, std::vector<HeapItem>,
                        std::greater<HeapItem>> heap;

    for (std::size_t r = 0; r < R; ++r) {
        in[r].open(runs_[r], std::ios::binary);
        if (in[r] && load(r)) {
            live[r] = true;
            heap.push({cur[r].key, r});
        }
    }

    while (!heap.empty()) {
        const std::uint64_t key = heap.top().first;
        std::uint64_t total = 0;
        // Drain every run currently at this key.
        while (!heap.empty() && heap.top().first == key) {
            const std::size_t r = heap.top().second;
            heap.pop();
            total += cur[r].count;
            if (load(r)) heap.push({cur[r].key, r});
            else live[r] = false;
        }
        out_keys.push_back(key);
        out_counts.push_back(static_cast<std::uint32_t>(
            std::min<std::uint64_t>(total, kCountMax)));
    }

    for (std::size_t r = 0; r < R; ++r) {
        in[r].close();
        std::remove(runs_[r].c_str());
    }
    runs_.clear();
    std::vector<std::uint64_t>().swap(buf_);
}

}  // namespace branch::wg
