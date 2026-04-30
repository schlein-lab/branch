// BRANCH v0.5 — Open-addressed exact k-mer counter.
//
// Used in the second pass of Phase 0: we iterate the read stream a
// second time, and for every k-mer that the bloom filter says is
// "maybe seen ≥ 2 times", we either insert it here with count=2 (we
// already know one prior occurrence existed because the bloom is
// populated) or increment its existing count. K-mers seen exactly once
// (singletons) never enter this table, which keeps memory bounded to
// the recurrence-relevant subset of the k-mer universe.
//
// Storage: 16 bytes per entry (8-byte key + 4-byte count + 4-byte pad).
// At 0.5 load factor the slot table doubles → 32 bytes per stored key.
// For a 90 Gbp human ONT input we expect ~3 G recurrent k-mers →
// ~96 GB exact-counter footprint, fits in std-partition (768 GB).

#pragma once

#include <cstdint>
#include <cstddef>
#include <functional>
#include <vector>

namespace branch::wg {

class KmerCounter {
public:
    /// `n_expected_recurrent_keys` is the bloom-survivor estimate,
    /// not the total k-mer universe. Pre-allocates 2 * that for
    /// load factor 0.5.
    explicit KmerCounter(std::size_t n_expected_recurrent_keys);

    /// Add a key with count delta (default +1). Idempotent on (key, 0).
    void add(std::uint64_t key, std::uint32_t delta = 1) noexcept;

    /// Returns count for key, or 0 if absent.
    [[nodiscard]] std::uint32_t get(std::uint64_t key) const noexcept;

    /// Iteration: visit every (key, count) pair via `f`. Used by the
    /// downstream histogram builder + bimodality scanner.
    template <typename F>
    void for_each(F&& f) const {
        for (std::size_t i = 0; i < slots_.size(); ++i) {
            if (slots_[i].count > 0) f(slots_[i].key, slots_[i].count);
        }
    }

    /// Number of stored entries (load).
    [[nodiscard]] std::size_t size() const noexcept { return n_; }

    /// Memory footprint in bytes.
    [[nodiscard]] std::size_t bytes() const noexcept {
        return slots_.size() * sizeof(Slot);
    }

private:
    struct Slot {
        std::uint64_t key   = 0;
        std::uint32_t count = 0;  // 0 = empty slot
        std::uint32_t pad   = 0;
    };
    std::vector<Slot> slots_;
    std::size_t mask_;          // n_slots - 1, n_slots is power of two
    std::size_t n_;             // load

    [[nodiscard]] std::size_t probe_start(std::uint64_t key) const noexcept;
    void grow();
};

}  // namespace branch::wg
