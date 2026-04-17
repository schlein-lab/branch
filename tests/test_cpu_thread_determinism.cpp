// Determinism test for cpu_backend: same input must produce identical
// OverlapPair outputs across runs and across thread counts, now that
// bucket + pair_matches iteration is driven by sorted key vectors.

#include <gtest/gtest.h>

#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"

namespace {

// Generate a deterministic set of "reads" that share k-mers: each read
// is a random 500 bp sequence, and reads overlap by appending shared
// prefixes/suffixes. 20 reads give the backend enough buckets to
// exercise iteration order.
std::vector<std::string> make_read_set(std::mt19937& rng, std::size_t n_reads,
                                       std::size_t read_len) {
    static const char bases[] = "ACGT";
    std::uniform_int_distribution<int> pick(0, 3);

    // A shared stem every read contains somewhere — forces overlaps.
    std::string stem(200, 'A');
    for (char& c : stem) c = bases[pick(rng)];

    std::vector<std::string> reads;
    reads.reserve(n_reads);
    for (std::size_t i = 0; i < n_reads; ++i) {
        std::string r;
        r.reserve(read_len);
        while (r.size() + stem.size() < read_len) r.push_back(bases[pick(rng)]);
        r.insert(i % (r.size() + 1), stem);
        while (r.size() < read_len) r.push_back(bases[pick(rng)]);
        reads.push_back(std::move(r));
    }
    return reads;
}

std::vector<branch::backend::OverlapPair> run_overlap(
    const std::vector<std::string>& reads, unsigned int threads) {
    branch::backend::set_cpu_overlap_threads(threads);
    auto bk = branch::backend::make_cpu_backend();

    branch::backend::ReadBatch batch;
    batch.reads.reserve(reads.size());
    for (std::size_t i = 0; i < reads.size(); ++i) {
        batch.reads.push_back({
            .seq = std::string_view(reads[i]),
            .id = static_cast<branch::graph::ReadId>(i),
        });
    }

    std::vector<branch::backend::OverlapPair> out(200'000);
    std::size_t n_out = 0;
    bk.compute_overlaps(&batch, out, &n_out);
    out.resize(n_out);
    return out;
}

bool pairs_equal(const std::vector<branch::backend::OverlapPair>& a,
                 const std::vector<branch::backend::OverlapPair>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::memcmp(&a[i], &b[i], sizeof(branch::backend::OverlapPair)) != 0) {
            return false;
        }
    }
    return true;
}

}  // namespace

TEST(CpuThreadDeterminism, SameInputSameOutputAcrossRuns) {
    std::mt19937 rng(42);
    auto reads = make_read_set(rng, 20, 500);

    auto run_a = run_overlap(reads, 1);
    auto run_b = run_overlap(reads, 1);
    EXPECT_TRUE(pairs_equal(run_a, run_b))
        << "Single-threaded run produced different outputs on repeated calls — "
           "iteration order is not deterministic";
}

TEST(CpuThreadDeterminism, SingleVsMultiThreadAgree) {
    std::mt19937 rng(12345);
    auto reads = make_read_set(rng, 20, 500);

    auto single = run_overlap(reads, 1);
    auto multi4 = run_overlap(reads, 4);
    EXPECT_TRUE(pairs_equal(single, multi4))
        << "4-thread run produced different OverlapPair sequence than "
           "single-threaded run — output ordering must be thread-count-"
           "invariant for downstream reproducibility.";
}

TEST(CpuThreadDeterminism, EmptyBatchReturnsNothing) {
    branch::backend::set_cpu_overlap_threads(2);
    auto bk = branch::backend::make_cpu_backend();
    branch::backend::ReadBatch batch;
    std::vector<branch::backend::OverlapPair> out(10);
    std::size_t n_out = 99;
    bk.compute_overlaps(&batch, out, &n_out);
    EXPECT_EQ(n_out, 0u);
}
