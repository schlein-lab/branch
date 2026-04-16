#pragma once

// BRANCH v0.1 - Reverse Index: Node -> Reads traversing it.
//
// Per long-read-assembly expert consultation: without this index,
// traversal lookups degenerate to O(Reads) per node. With it, O(1)
// amortized at the cost of approximately 4.8 GB RAM for a human WGS
// graph at 30x HiFi coverage.
//
// Layout is CSR-like: a single offsets array keyed by NodeId, a
// single flat read_ids array. Append-only build; a build()-to-freeze
// step flips from mutable per-node vectors to the frozen CSR form.

#include <cstddef>
#include <cstdint>
#include <vector>
#include <cassert>
#include <span>
#include <stdexcept>

#include "graph/delta_read.hpp"

namespace branch::graph {

class ReverseIndex {
public:
    ReverseIndex() = default;

    // Mutable build phase: append (node, read) pairs.
    void reserve_nodes(std::size_t n_nodes) {
        buckets_.resize(n_nodes);
    }

    void add(NodeId node, ReadId read) {
        if (frozen_) {
            throw std::logic_error("add() called after freeze()");
        }
        if (node >= buckets_.size()) {
            buckets_.resize(static_cast<std::size_t>(node) + 1);
        }
        buckets_[node].push_back(read);
    }

    // Freeze mutable buckets into CSR form: offsets_ + read_ids_.
    // After freeze(), add() must not be called.
    void freeze() {
        assert(!frozen_);
        std::size_t total = 0;
        for (const auto& b : buckets_) total += b.size();

        offsets_.clear();
        offsets_.reserve(buckets_.size() + 1);
        read_ids_.clear();
        read_ids_.reserve(total);

        offsets_.push_back(0);
        for (auto& b : buckets_) {
            for (auto r : b) read_ids_.push_back(r);
            offsets_.push_back(static_cast<std::uint32_t>(read_ids_.size()));
        }
        buckets_.clear();
        buckets_.shrink_to_fit();
        frozen_ = true;
    }

    [[nodiscard]] bool is_frozen() const noexcept { return frozen_; }
    [[nodiscard]] std::size_t node_count() const noexcept {
        return frozen_ ? (offsets_.empty() ? 0 : offsets_.size() - 1)
                       : buckets_.size();
    }
    [[nodiscard]] std::size_t total_entries() const noexcept {
        return frozen_ ? read_ids_.size() : count_mutable_entries();
    }

    // Read lookup (frozen form only).
    [[nodiscard]] std::span<const ReadId> reads_for(NodeId node) const {
        assert(frozen_);
        assert(static_cast<std::size_t>(node) < node_count());
        auto start = offsets_[node];
        auto end = offsets_[static_cast<std::size_t>(node) + 1];
        return std::span<const ReadId>(
            read_ids_.data() + start, end - start);
    }

private:
    [[nodiscard]] std::size_t count_mutable_entries() const {
        std::size_t n = 0;
        for (const auto& b : buckets_) n += b.size();
        return n;
    }

    // Build-time mutable storage. Discarded after freeze().
    std::vector<std::vector<ReadId>> buckets_;

    // Frozen CSR form.
    std::vector<std::uint32_t> offsets_;
    std::vector<ReadId> read_ids_;
    bool frozen_{false};
};

}  // namespace branch::graph
