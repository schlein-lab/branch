#include <limits>
// BRANCH v0.5 — Exact counter implementation.

#include "wg/kmer_counter.hpp"

namespace branch::wg {

namespace {
constexpr std::uint64_t splitmix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
}

KmerCounter::KmerCounter(std::size_t n_expected_recurrent_keys) {
    if (n_expected_recurrent_keys == 0) n_expected_recurrent_keys = 1;
    // Round up to next power of two and aim for load factor 0.5.
    std::size_t need = n_expected_recurrent_keys * 2;
    std::size_t cap = 1;
    while (cap < need) cap <<= 1;
    if (cap < 16) cap = 16;
    slots_.assign(cap, Slot{});
    mask_ = cap - 1;
    n_ = 0;
}

std::size_t KmerCounter::probe_start(std::uint64_t key) const noexcept {
    return splitmix64(key) & mask_;
}

void KmerCounter::grow() {
    std::vector<Slot> old = std::move(slots_);
    std::size_t new_cap = (mask_ + 1) << 1;
    slots_.assign(new_cap, Slot{});
    mask_ = new_cap - 1;
    n_ = 0;
    for (auto& s : old) {
        if (s.count > 0) add(s.key, s.count);
    }
}

void KmerCounter::add(std::uint64_t key, std::uint32_t delta) noexcept {
    if (n_ * 2 >= slots_.size()) grow();
    std::size_t i = probe_start(key);
    while (true) {
        Slot& s = slots_[i];
        if (s.count == 0) {
            s.key = key;
            s.count = delta;
            ++n_;
            return;
        }
        if (s.key == key) {
            // Saturate at uint32 max
            std::uint64_t newc = static_cast<std::uint64_t>(s.count) + delta;
            s.count = newc > std::numeric_limits<std::uint32_t>::max()
                      ? std::numeric_limits<std::uint32_t>::max()
                      : static_cast<std::uint32_t>(newc);
            return;
        }
        i = (i + 1) & mask_;
    }
}

std::uint32_t KmerCounter::get(std::uint64_t key) const noexcept {
    std::size_t i = probe_start(key);
    while (true) {
        const Slot& s = slots_[i];
        if (s.count == 0) return 0;
        if (s.key == key) return s.count;
        i = (i + 1) & mask_;
    }
}

}  // namespace branch::wg
