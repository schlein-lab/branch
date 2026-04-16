// BRANCH v0.1 — Minimizer sketcher implementation.
//
// Straightforward rolling-window implementation. Performance is
// adequate for v0.1 development; a SIMD version lands in v0.2.

#include "graph/minimizer_sketcher.hpp"

#include <algorithm>
#include <limits>

namespace branch::graph {

void sketch_read(std::string_view seq,
                 ReadId read_id,
                 std::vector<MinimizerHit>& out) noexcept {
    const std::size_t k = kMinimizerK;
    const std::size_t w = kMinimizerW;

    if (seq.size() < k) {
        return;
    }

    // Build initial k-mer, resetting on any ambiguous base.
    std::size_t filled = 0;
    std::uint64_t kmer = 0;
    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);

    // Rolling window of (hash, position) for the last w k-mers.
    struct WinEntry {
        std::uint64_t hash;
        std::uint32_t pos;
    };
    std::vector<WinEntry> window;
    window.reserve(w);

    std::uint64_t last_emitted_hash = std::numeric_limits<std::uint64_t>::max();
    std::uint32_t last_emitted_pos = std::numeric_limits<std::uint32_t>::max();

    auto flush_if_new_minimum = [&](std::uint32_t end_pos) {
        // Pick the smallest-hash entry in the window; emit it if it
        // differs from the previous emission (avoid consecutive
        // duplicate emissions).
        if (window.empty()) return;
        auto best = window.front();
        for (const auto& e : window) {
            if (e.hash < best.hash) best = e;
        }
        if (best.hash != last_emitted_hash || best.pos != last_emitted_pos) {
            out.push_back(MinimizerHit{
                .hash = best.hash,
                .read_id = read_id,
                .pos = best.pos,
            });
            last_emitted_hash = best.hash;
            last_emitted_pos = best.pos;
        }
        (void)end_pos;
    };

    for (std::size_t i = 0; i < seq.size(); ++i) {
        std::uint8_t b = encode_base(seq[i]);
        if (b > 3) {
            // Ambiguous base — reset
            filled = 0;
            kmer = 0;
            window.clear();
            continue;
        }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;

        if (filled < k) continue;

        std::uint32_t kmer_start = static_cast<std::uint32_t>(i + 1 - k);
        std::uint64_t h = canonical_hash(kmer, k);

        window.push_back(WinEntry{.hash = h, .pos = kmer_start});
        // Keep window at most w k-mers: drop the oldest.
        if (window.size() > w) {
            window.erase(window.begin());
        }

        // Once we have a full window, emit the minimum.
        if (window.size() == w) {
            flush_if_new_minimum(kmer_start);
        }
    }
}

}  // namespace branch::graph
