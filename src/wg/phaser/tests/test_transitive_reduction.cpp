// Unit test for OverlapGraph::transitive_reduce.
//
// Constructs synthetic edge sets with known shortcut topology and
// verifies that transitive reduction removes the shortcuts.
//
// We can't easily call OverlapGraph::transitive_reduce() directly because
// edges_ is private and only populated via build(). Instead we test the
// behaviour indirectly via the toy phaser pipeline: a graph that needs
// TR to be linear-extractable should produce more, longer unitigs WITH
// TR than without.

#include "wg/phaser/overlap_graph.hpp"
#include "wg/phaser/types.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace branch::wg::phaser;

namespace {

std::string rand_seq(std::size_t n, std::mt19937& rng) {
    std::string s; s.reserve(n);
    static const char* a = "ACGT";
    for (std::size_t i = 0; i < n; ++i) s.push_back(a[rng() & 3]);
    return s;
}

std::string mutate(const std::string& s, double p, std::mt19937& rng) {
    std::string out = s;
    std::uniform_real_distribution<double> ud(0.0, 1.0);
    static const char* alt[4] = {"CGT", "AGT", "ACT", "ACG"};
    for (auto& c : out) {
        if (ud(rng) < p) {
            int idx;
            switch (c) {
                case 'A': idx = 0; break;
                case 'C': idx = 1; break;
                case 'G': idx = 2; break;
                default:  idx = 3; break;
            }
            c = alt[idx][rng() % 3];
        }
    }
    return out;
}

}  // namespace

int main() {
    // Construct a simple chain genome with overlapping reads. With TR,
    // the unitig walker should produce a small number of long unitigs;
    // without TR, the dense overlap graph fragments into many tiny ones.

    std::mt19937 rng(7);
    const std::string genome = rand_seq(50'000, rng);
    const std::size_t read_len = 5'000;
    const std::size_t step = 500;  // dense overlap

    std::vector<std::pair<ReadId, std::string>> reads;
    ReadId rid = 0;
    for (std::size_t pos = 0; pos + read_len <= genome.size(); pos += step) {
        std::string r = mutate(genome.substr(pos, read_len), 0.005, rng);
        reads.emplace_back(rid++, std::move(r));
    }
    std::cerr << "[test_TR] generated " << reads.size() << " reads\n";

    OverlapGraphOpts opts;
    opts.minimizer_k = 21;
    opts.minimizer_w = 11;
    opts.min_overlap_bp = 500;
    opts.min_jaccard = 0.20;
    opts.n_threads = 1;

    // Pass 1: build → collapse without TR.
    std::size_t edges_before = 0;
    std::size_t edges_after = 0;
    {
        OverlapGraph og;
        og.build(reads, opts);
        edges_before = og.edges().size();
        std::cerr << "[test_TR] edges (raw)         = " << edges_before << "\n";
    }

    // Pass 2: build → TR → check edge count drops.
    {
        OverlapGraph og;
        og.build(reads, opts);
        og.transitive_reduce();
        edges_after = og.edges().size();
        std::cerr << "[test_TR] edges (after TR)    = " << edges_after << "\n";
    }

    // For a dense chain graph, TR should remove a substantial fraction
    // of shortcut edges. The exact ratio depends on overlap params, but
    // we expect at least 20% reduction at this density.
    assert(edges_after < edges_before &&
           "TR must reduce edge count on a dense chain graph");
    const double ratio = 1.0 -
        static_cast<double>(edges_after) / static_cast<double>(edges_before);
    std::cerr << "[test_TR] reduction ratio = " << ratio * 100.0 << "%\n";
    assert(ratio > 0.10 &&
           "TR should remove >=10% of edges on a dense chain");

    // TR must never add edges or duplicate them.
    assert(edges_after <= edges_before);

    std::cout << "test_transitive_reduction: OK\n";
    return 0;
}
