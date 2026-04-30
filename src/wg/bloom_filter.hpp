// BRANCH v0.5 — Bloom filter for two-pass k-mer counting.
//
// First pass over the read stream populates the bloom; second pass
// only adds to the exact counter k-mers that test positive (i.e.
// likely seen at least twice). This keeps the exact counter sized
// to "k-mers that recurred", not "all distinct k-mers", which for
// a 90 Gbp human ONT input is roughly a 10x reduction.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace branch::wg {

class BloomFilter {
public:
    /// `n_expected_keys` ≈ how many distinct k-mers we expect.
    /// `bits_per_key`   tunes false-positive rate (5 → ~10% FP, 7 → ~3%, 10 → ~1%).
    /// Memory used = n_expected_keys * bits_per_key / 8 bytes.
    BloomFilter(std::size_t n_expected_keys, std::size_t bits_per_key = 6);

    /// Add a key. Idempotent.
    void add(std::uint64_t key) noexcept;

    /// Returns true if `add(key)` was likely called. False = definitely not.
    [[nodiscard]] bool maybe_contains(std::uint64_t key) const noexcept;

    /// Total bits set / total bits — for diagnostics. >0.5 means we are
    /// over-saturated and should have used more bits_per_key.
    [[nodiscard]] double saturation() const noexcept;

    /// Memory footprint in bytes.
    [[nodiscard]] std::size_t bytes() const noexcept;

private:
    std::vector<std::uint64_t> bits_;
    std::size_t n_bits_;
    int n_hashes_;       // optimal: bits_per_key * ln(2) ≈ 0.69 * bits_per_key

    // Three independent splitmix-style hash variants of the same key,
    // covering up to n_hashes_ ≤ 7 hashes for 6-10 bits/key.
    std::uint64_t hash_i(std::uint64_t key, int i) const noexcept;
};

}  // namespace branch::wg
