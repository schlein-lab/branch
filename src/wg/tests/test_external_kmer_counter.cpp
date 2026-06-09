// Unit test for ExternalKmerCounter — verifies the spill + k-way-merge path
// produces exactly the same sorted (key, count) result as an in-memory count,
// including when the buffer is tiny enough to force many spills. Standalone
// (no GoogleTest); assertions kept live via -UNDEBUG.

#include "wg/external_kmer_counter.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <string>
#include <vector>

using branch::wg::ExternalKmerCounter;

namespace {

std::string scratch() {
    if (const char* p = std::getenv("TMPDIR")) {
        if (*p) return p;
    }
    return "/tmp";
}

// Run the counter with a given RAM budget over `keys`, compare against an
// exact std::map reference.
void check(const std::vector<std::uint64_t>& keys, std::size_t max_ram_bytes) {
    std::map<std::uint64_t, std::uint64_t> ref;
    for (auto k : keys) ++ref[k];

    ExternalKmerCounter c(max_ram_bytes, scratch());
    for (auto k : keys) c.add(k);
    std::vector<std::uint64_t> out_keys;
    std::vector<std::uint32_t> out_counts;
    c.finalize(out_keys, out_counts);

    // Same number of distinct keys.
    assert(out_keys.size() == ref.size());
    // Sorted ascending, exact counts, matching the reference.
    std::uint64_t prev = 0;
    bool first = true;
    for (std::size_t i = 0; i < out_keys.size(); ++i) {
        if (!first) assert(out_keys[i] > prev);  // strictly increasing, deduped
        first = false;
        prev = out_keys[i];
        auto it = ref.find(out_keys[i]);
        assert(it != ref.end());
        assert(out_counts[i] == it->second);
    }
}

}  // namespace

int main() {
    // 1. Tiny input, no spill (buffer floor is 1M keys → single in-memory run).
    check({5, 1, 5, 3, 1, 1, 9}, 1ULL << 30);

    // 2. Many duplicates with a forced-small budget → many spills + merge.
    //    Budget floors at 1M keys internally, so use a large input to exceed it.
    {
        std::mt19937_64 rng(1234);
        std::vector<std::uint64_t> keys;
        keys.reserve(5'000'000);
        // 50k distinct keys, each repeated ~100x in random order.
        for (int rep = 0; rep < 100; ++rep)
            for (std::uint64_t k = 0; k < 50'000; ++k)
                keys.push_back(k * 2654435761ULL);  // spread the key space
        std::shuffle(keys.begin(), keys.end(), rng);
        check(keys, 8ULL << 20);  // 8 MB budget → ~5 spills, exercises merge
    }

    // 3. Empty input.
    check({}, 1ULL << 30);

    std::printf("test_external_kmer_counter: all tests passed\n");
    return 0;
}
