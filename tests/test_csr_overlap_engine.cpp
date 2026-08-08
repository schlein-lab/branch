// Tests for the CSR overlap engine (backend/csr_overlap_engine.hpp).
//
// The engine must be a drop-in replacement for the legacy
// compute_overlaps path: with seed_cap == 0 (matching the HiFi
// profile's uncapped buckets) both engines emit the identical set of
// OverlapPairs on the same batch. On top of equivalence, the CSR
// engine adds a seeding cap with overcap statistics and bit-identical
// output across thread counts.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "backend/backend_factory.hpp"
#include "backend/cpu_backend.hpp"
#include "backend/csr_overlap_engine.hpp"
#include "graph/minimizer_sketcher.hpp"

namespace bb = branch::backend;
namespace bg = branch::graph;

namespace {

std::string random_genome(std::size_t len, std::uint32_t seed) {
    static const char* bases = "ACGT";
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> d(0, 3);
    std::string g;
    g.reserve(len);
    for (std::size_t i = 0; i < len; ++i) g.push_back(bases[d(rng)]);
    return g;
}

std::string revcomp(std::string_view s) {
    std::string r;
    r.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        switch (*it) {
            case 'A': r.push_back('T'); break;
            case 'C': r.push_back('G'); break;
            case 'G': r.push_back('C'); break;
            case 'T': r.push_back('A'); break;
            default: r.push_back('N'); break;
        }
    }
    return r;
}

// Overlapping read tiling of a genome, every 4th read reverse-
// complemented so the strand path is exercised.
std::vector<std::string> tile_reads(const std::string& genome,
                                    std::size_t read_len, std::size_t stride) {
    std::vector<std::string> reads;
    for (std::size_t start = 0; start + read_len <= genome.size();
         start += stride) {
        std::string r = genome.substr(start, read_len);
        if (reads.size() % 4 == 3) r = revcomp(r);
        reads.push_back(std::move(r));
    }
    return reads;
}

bb::ReadBatch make_batch(const std::vector<std::string>& seqs) {
    bb::ReadBatch batch;
    batch.reads.reserve(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        batch.reads.push_back(bb::ReadBatch::Read{
            .seq = seqs[i],
            .id = static_cast<bg::ReadId>(i),
        });
    }
    return batch;
}

std::vector<bb::OverlapPair> run_legacy(const bb::ReadBatch& batch,
                                        unsigned int threads) {
    bb::set_cpu_overlap_threads(threads);
    auto backend = bb::make_cpu_backend();
    std::vector<bb::OverlapPair> pairs(2'000'000);
    std::size_t n = 0;
    backend.compute_overlaps(&batch, pairs, &n);
    pairs.resize(n);
    std::sort(pairs.begin(), pairs.end(),
              [](const bb::OverlapPair& x, const bb::OverlapPair& y) {
                  if (x.read_a != y.read_a) return x.read_a < y.read_a;
                  if (x.read_b != y.read_b) return x.read_b < y.read_b;
                  return x.offset_a < y.offset_a;
              });
    return pairs;
}

std::vector<bb::OverlapPair> run_csr(const bb::ReadBatch& batch,
                                     unsigned int threads,
                                     std::size_t seed_cap,
                                     bb::CsrIndexStats* stats = nullptr) {
    bb::CsrOverlapConfig cfg;
    cfg.threads = threads;
    cfg.seed_cap = seed_cap;
    std::vector<bb::OverlapPair> out;
    bb::compute_overlaps_csr(batch, cfg, out, stats);
    return out;
}

void expect_identical(const std::vector<bb::OverlapPair>& a,
                      const std::vector<bb::OverlapPair>& b) {
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        EXPECT_EQ(0, std::memcmp(&a[i], &b[i], sizeof(bb::OverlapPair)))
            << "pair " << i << " differs: legacy (" << a[i].read_a << ","
            << a[i].read_b << ",off_a=" << a[i].offset_a
            << ") vs csr (" << b[i].read_a << "," << b[i].read_b
            << ",off_a=" << b[i].offset_a << ")";
    }
}

class CsrOverlapEngine : public ::testing::Test {
protected:
    void SetUp() override { bg::set_current_profile(bg::kProfileHiFi); }
};

}  // namespace

TEST_F(CsrOverlapEngine, EquivalentToLegacyOnTiledGenome) {
    const auto genome = random_genome(50'000, /*seed=*/42);
    const auto reads = tile_reads(genome, 3'000, 500);
    ASSERT_GT(reads.size(), 50u);
    const auto batch = make_batch(reads);

    const auto legacy = run_legacy(batch, 1);
    ASSERT_GT(legacy.size(), 100u) << "test fixture must produce overlaps";
    const auto csr = run_csr(batch, 1, /*seed_cap=*/0);
    expect_identical(legacy, csr);
}

TEST_F(CsrOverlapEngine, DeterministicAcrossThreadCounts) {
    const auto genome = random_genome(80'000, /*seed=*/7);
    const auto reads = tile_reads(genome, 4'000, 700);
    const auto batch = make_batch(reads);

    const auto t1 = run_csr(batch, 1, 0);
    const auto t4 = run_csr(batch, 4, 0);
    expect_identical(t1, t4);
}

TEST_F(CsrOverlapEngine, SeedCapReducesRepeatCandidatesAndReports) {
    // 120 identical copies of one repeat read: every minimizer of the
    // repeat occurs 120x. With seed_cap=16 those hashes must be capped
    // and reported; the distinct tail (unique tiled reads) keeps its
    // overlaps.
    const auto repeat_unit = random_genome(2'000, /*seed=*/99);
    const auto genome = random_genome(30'000, /*seed=*/123);
    std::vector<std::string> reads(120, repeat_unit);
    const auto unique_reads = tile_reads(genome, 3'000, 600);
    reads.insert(reads.end(), unique_reads.begin(), unique_reads.end());
    const auto batch = make_batch(reads);

    bb::CsrIndexStats uncapped_stats;
    const auto uncapped = run_csr(batch, 2, 0, &uncapped_stats);
    EXPECT_EQ(uncapped_stats.n_hashes_capped, 0u);

    bb::CsrIndexStats capped_stats;
    const auto capped = run_csr(batch, 2, 16, &capped_stats);

    EXPECT_GT(capped_stats.n_hashes_capped, 0u);
    EXPECT_GT(capped_stats.n_entries_dropped, 0u);
    EXPECT_FALSE(capped_stats.capped_hashes.empty());
    EXPECT_GE(capped_stats.capped_hashes.front().count, 120u);
    EXPECT_LT(capped.size(), uncapped.size());

    // Unique-tail overlaps survive the cap: count pairs whose both ids
    // are in the unique range [120, ...).
    const auto count_unique_pairs = [](const std::vector<bb::OverlapPair>& v) {
        std::size_t n = 0;
        for (const auto& p : v) {
            if (p.read_a >= 120 && p.read_b >= 120) ++n;
        }
        return n;
    };
    EXPECT_EQ(count_unique_pairs(capped), count_unique_pairs(uncapped));
}

TEST_F(CsrOverlapEngine, EmptyAndSingleReadBatches) {
    bb::ReadBatch empty;
    std::vector<bb::OverlapPair> out;
    bb::CsrOverlapConfig cfg;
    EXPECT_EQ(0u, bb::compute_overlaps_csr(empty, cfg, out));
    EXPECT_TRUE(out.empty());

    const std::string one = random_genome(3'000, 5);
    std::vector<std::string> seqs{one};
    const auto batch = make_batch(seqs);
    EXPECT_EQ(0u, bb::compute_overlaps_csr(batch, cfg, out));
    EXPECT_TRUE(out.empty());
}

TEST_F(CsrOverlapEngine, NoSharedMinimizersMeansNoPairs) {
    // Two unrelated random sequences share (practically) no 21-mers.
    std::vector<std::string> seqs{random_genome(3'000, 1000),
                                  random_genome(3'000, 2000)};
    const auto batch = make_batch(seqs);
    std::vector<bb::OverlapPair> out;
    bb::CsrOverlapConfig cfg;
    bb::compute_overlaps_csr(batch, cfg, out);
    EXPECT_TRUE(out.empty());
}
