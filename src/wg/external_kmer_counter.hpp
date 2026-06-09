// BRANCH — RAM-budgeted exact k-mer counter (external merge-sort).
//
// Drop-in replacement for the in-RAM open-addressed KmerCounter in Phase 0's
// second pass. Instead of holding every recurrent k-mer in a hash table
// (~96 GB for a human whole genome — the Phase-0 RAM wall), it buffers k-mer
// occurrences up to a fixed budget, then spills a sorted, run-length-counted
// run to disk and clears the buffer. finalize() k-way-merges the runs, summing
// counts for equal keys, and produces the SAME sorted (keys, counts) arrays the
// downstream histogram + Phase-1.1 het-pair detection already consume.
//
// Peak RAM is bounded by the buffer budget (--max-ram-gb) regardless of input
// size; counts are EXACT (no approximation), so het-pair detection quality is
// unchanged. Scratch I/O lives on BeeGFS (TMPDIR).

#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace branch::wg {

class ExternalKmerCounter {
public:
    // `max_ram_bytes` caps the in-memory occurrence buffer; `scratch_dir` holds
    // spilled runs (should be on BeeGFS). A budget large enough to hold all
    // occurrences means no spill (single in-memory sort) — identical result.
    ExternalKmerCounter(std::size_t max_ram_bytes, std::string scratch_dir);
    ~ExternalKmerCounter();

    ExternalKmerCounter(const ExternalKmerCounter&) = delete;
    ExternalKmerCounter& operator=(const ExternalKmerCounter&) = delete;

    // Record one occurrence of `key` (called once per recurrent k-mer hit).
    void add(std::uint64_t key);

    // Merge all buffered + spilled occurrences into sorted distinct keys with
    // summed counts. Consumes the counter (deletes run files). Counts are
    // saturated at UINT32_MAX.
    void finalize(std::vector<std::uint64_t>& out_keys,
                  std::vector<std::uint32_t>& out_counts);

    // Peak in-memory buffer footprint in bytes (for diagnostics).
    [[nodiscard]] std::size_t peak_buffer_bytes() const noexcept {
        return cap_ * sizeof(std::uint64_t);
    }
    [[nodiscard]] std::size_t n_spilled_runs() const noexcept {
        return runs_.size();
    }

private:
    void spill();  // sort buffer, run-length count, write a run file

    std::vector<std::uint64_t> buf_;
    std::size_t                cap_;       // max keys held before a spill
    std::string                scratch_;
    std::vector<std::string>   runs_;      // spilled run file paths
    std::size_t                run_seq_ = 0;
};

}  // namespace branch::wg
