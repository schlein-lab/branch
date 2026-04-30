// BRANCH v0.5 — Bloom filter implementation.

#include "wg/bloom_filter.hpp"
#include <cmath>

namespace branch::wg {

namespace {
constexpr std::uint64_t splitmix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
}  // namespace

BloomFilter::BloomFilter(std::size_t n_expected_keys, std::size_t bits_per_key) {
    if (n_expected_keys == 0) n_expected_keys = 1;
    n_bits_ = n_expected_keys * bits_per_key;
    // Round up to multiple of 64.
    n_bits_ = ((n_bits_ + 63) / 64) * 64;
    bits_.assign(n_bits_ / 64, 0);
    // Optimal k for given m/n is ln(2) * m/n ≈ 0.693 * bits_per_key.
    n_hashes_ = static_cast<int>(std::round(0.693 * static_cast<double>(bits_per_key)));
    if (n_hashes_ < 1) n_hashes_ = 1;
    if (n_hashes_ > 7) n_hashes_ = 7;
}

std::uint64_t BloomFilter::hash_i(std::uint64_t key, int i) const noexcept {
    // Different salts for independent hash families.
    static constexpr std::uint64_t salts[7] = {
        0x9e3779b97f4a7c15ULL, 0xbf58476d1ce4e5b9ULL, 0x94d049bb133111ebULL,
        0xff51afd7ed558ccdULL, 0xc4ceb9fe1a85ec53ULL, 0x85ebca77c2b2ae63ULL,
        0xcc9e2d51a4093822ULL,
    };
    return splitmix64(key ^ salts[i]);
}

void BloomFilter::add(std::uint64_t key) noexcept {
    for (int i = 0; i < n_hashes_; ++i) {
        std::uint64_t h = hash_i(key, i) % n_bits_;
        bits_[h / 64] |= (1ULL << (h % 64));
    }
}

bool BloomFilter::maybe_contains(std::uint64_t key) const noexcept {
    for (int i = 0; i < n_hashes_; ++i) {
        std::uint64_t h = hash_i(key, i) % n_bits_;
        if (!(bits_[h / 64] & (1ULL << (h % 64)))) return false;
    }
    return true;
}

double BloomFilter::saturation() const noexcept {
    std::size_t set = 0;
    for (auto w : bits_) set += static_cast<std::size_t>(__builtin_popcountll(w));
    return n_bits_ ? static_cast<double>(set) / static_cast<double>(n_bits_) : 0.0;
}

std::size_t BloomFilter::bytes() const noexcept { return bits_.size() * sizeof(std::uint64_t); }

}  // namespace branch::wg
