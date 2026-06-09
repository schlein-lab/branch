// BRANCH v0.5 — Phase 0 implementation.

#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace branch::wg {

namespace {

constexpr std::uint64_t splitmix64(std::uint64_t z) noexcept {
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

inline std::uint8_t encode_base(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // ambiguous = reset window
    }
}

inline std::uint64_t rc_kmer(std::uint64_t kmer, std::size_t k) noexcept {
    std::uint64_t rc = 0;
    for (std::size_t i = 0; i < k; ++i) {
        rc = (rc << 2) | ((~kmer) & 3ULL);
        kmer >>= 2;
    }
    return rc;
}

}  // namespace

std::uint64_t canonical_kmer_hash(std::uint64_t kmer, std::size_t k) noexcept {
    auto hf = splitmix64(kmer);
    auto hr = splitmix64(rc_kmer(kmer, k));
    return hf < hr ? hf : hr;
}

std::uint64_t canonical_kmer_bits(std::uint64_t kmer, std::size_t k) noexcept {
    auto rc = rc_kmer(kmer, k);
    return kmer < rc ? kmer : rc;
}

void kmer_bits_to_seq(std::uint64_t bits, std::size_t k, char* out) noexcept {
    static const char kBases[4] = {'A', 'C', 'G', 'T'};
    for (std::size_t i = 0; i < k; ++i) {
        out[k - 1 - i] = kBases[bits & 3ULL];
        bits >>= 2;
    }
}

KmerHistBuilder::KmerHistBuilder(const ::branch::graph::TechProfile& profile,
                                 std::size_t expected_n_distinct_kmers,
                                 std::size_t max_ram_bytes,
                                 std::string scratch_dir)
    : profile_(profile),
      bloom_(expected_n_distinct_kmers, /*bits_per_key=*/6),
      ext_counter_(max_ram_bytes, std::move(scratch_dir)) {
}

void KmerHistBuilder::scan_kmers_(std::string_view seq, bool pass2) noexcept {
    const std::size_t k = profile_.minimizer_k;
    if (seq.size() < k) return;

    const std::uint64_t mask = (k >= 32) ? ~0ULL : ((1ULL << (2 * k)) - 1ULL);
    std::uint64_t kmer = 0;
    std::size_t filled = 0;

    for (std::size_t i = 0; i < seq.size(); ++i) {
        std::uint8_t b = encode_base(seq[i]);
        if (b > 3) {
            filled = 0; kmer = 0;
            continue;
        }
        kmer = ((kmer << 2) | b) & mask;
        ++filled;
        if (filled < k) continue;

        // Counter key = canonical kmer-bits. Bloom uses the same key
        // and hashes it internally — no semantic difference, but downstream
        // consumers (haplotype_router) can decode the bits directly to
        // enumerate 1-bp variants for het-pair detection.
        std::uint64_t cbits = canonical_kmer_bits(kmer, k);
        if (pass2) {
            if (bloom_.maybe_contains(cbits)) ext_counter_.add(cbits);
        } else {
            bloom_.add(cbits);
        }
    }
}

void KmerHistBuilder::pass1_add_read(std::string_view seq) noexcept {
    if (in_pass2_) return;
    scan_kmers_(seq, /*pass2=*/false);
    ++n_reads_;
    n_bases_ += seq.size();
}

void KmerHistBuilder::switch_to_pass2() {
    in_pass2_ = true;
    n_reads_ = 0;
    n_bases_ = 0;
}

void KmerHistBuilder::pass2_add_read(std::string_view seq) noexcept {
    if (!in_pass2_) return;
    scan_kmers_(seq, /*pass2=*/true);
    ++n_reads_;
    n_bases_ += seq.size();
}

KmerHistResult KmerHistBuilder::finalize(int expected_coverage,
                                         double bimodality_z_threshold) {
    // External merge-sort the second-pass occurrences into the same sorted
    // (key, count) arrays the GPU path produces, then reuse the shared
    // histogram + bimodality logic. keys_/counts_ are retained for write_to().
    // Peak RAM is bounded by the counter's buffer budget, not the recurrent
    // k-mer count, so this no longer scales with genome size.
    ext_counter_.finalize(keys_, counts_);
    return finalize_from_keys_counts(expected_coverage, bimodality_z_threshold,
                                     keys_, counts_, n_reads_, n_bases_);
}

void KmerHistBuilder::write_to(const std::string& path,
                               const KmerHistResult& result) const {
    // Emit the identical .bin format via the shared keys/counts writer; the
    // sorted arrays were materialized by finalize(). (Same MAGIC, header,
    // histogram and (key,count) records the GPU path writes.)
    write_phase0_dump(path, profile_.minimizer_k, result, keys_, counts_);
}

KmerHistResult finalize_from_keys_counts(
    int expected_coverage,
    double bimodality_z_threshold,
    const std::vector<std::uint64_t>& keys,
    const std::vector<std::uint32_t>& counts,
    std::uint64_t n_reads,
    std::uint64_t n_bases) {
    KmerHistResult r;
    r.n_reads = n_reads;
    r.n_bases = n_bases;
    r.n_recurrent_kmers = keys.size();

    constexpr std::size_t HIST_MAX = 512;
    r.histogram.assign(HIST_MAX, 0);
    for (std::uint32_t c : counts) {
        std::size_t bin = std::min<std::size_t>(c, HIST_MAX - 1);
        ++r.histogram[bin];
    }

    int lo = std::max(4, expected_coverage / 3);
    int hi = std::min<int>(static_cast<int>(HIST_MAX) - 1, expected_coverage * 4);
    if (lo >= hi) return r;

    auto smooth = [&](int x) -> double {
        double s = 0; int n = 0;
        for (int d = -2; d <= 2; ++d) {
            int xx = x + d;
            if (xx >= 0 && xx < (int)r.histogram.size()) {
                s += static_cast<double>(r.histogram[xx]); ++n;
            }
        }
        return n ? s / n : 0.0;
    };

    int peak1 = lo;
    double v1 = smooth(lo);
    for (int x = lo + 1; x <= hi; ++x) {
        double v = smooth(x);
        if (v > v1) { v1 = v; peak1 = x; }
    }
    r.detected_coverage = peak1;

    int peak2 = 0;
    double v2 = 0.0;
    int hi2 = std::max(lo, static_cast<int>(peak1 * 0.7));
    for (int x = lo; x <= hi2; ++x) {
        double v = smooth(x);
        if (v > v2) { v2 = v; peak2 = x; }
    }

    if (peak2 > 0) {
        int valley = (peak1 + peak2) / 2;
        double vv = smooth(valley);
        if (vv > 0) {
            double z_low  = (v2 - vv) / std::sqrt(std::max(vv, 1.0));
            double z_high = (v1 - vv) / std::sqrt(std::max(vv, 1.0));
            r.bimodality_z = std::min(z_low, z_high);
            if (r.bimodality_z >= bimodality_z_threshold) {
                r.detected_het_coverage = peak2;
            }
        }
    }
    return r;
}

void write_phase0_dump(
    const std::string& path,
    std::size_t k,
    const KmerHistResult& result,
    const std::vector<std::uint64_t>& keys,
    const std::vector<std::uint32_t>& counts) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open: " + path);
    constexpr char MAGIC[8] = {'K','H','I','S','T','0','5','\0'};
    f.write(MAGIC, 8);

    std::uint64_t n_recurrent = result.n_recurrent_kmers;
    std::uint64_t n_reads = result.n_reads;
    std::uint64_t n_bases = result.n_bases;
    std::int32_t cov = result.detected_coverage;
    std::int32_t het = result.detected_het_coverage;
    double z = result.bimodality_z;
    std::uint32_t kk = static_cast<std::uint32_t>(k);
    f.write(reinterpret_cast<const char*>(&n_recurrent), sizeof(n_recurrent));
    f.write(reinterpret_cast<const char*>(&n_reads), sizeof(n_reads));
    f.write(reinterpret_cast<const char*>(&n_bases), sizeof(n_bases));
    f.write(reinterpret_cast<const char*>(&cov), sizeof(cov));
    f.write(reinterpret_cast<const char*>(&het), sizeof(het));
    f.write(reinterpret_cast<const char*>(&z), sizeof(z));
    f.write(reinterpret_cast<const char*>(&kk), sizeof(kk));

    std::uint32_t hist_n = static_cast<std::uint32_t>(result.histogram.size());
    f.write(reinterpret_cast<const char*>(&hist_n), sizeof(hist_n));
    f.write(reinterpret_cast<const char*>(result.histogram.data()),
            result.histogram.size() * sizeof(std::uint64_t));

    for (std::size_t i = 0; i < keys.size(); ++i) {
        f.write(reinterpret_cast<const char*>(&keys[i]), sizeof(std::uint64_t));
        f.write(reinterpret_cast<const char*>(&counts[i]), sizeof(std::uint32_t));
    }
}

void load_kmer_counter_dump(
    const std::string& path,
    std::unordered_map<std::uint64_t, std::uint32_t>& out_counter) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open Phase 0 dump: " + path);
    char magic[8];
    f.read(magic, 8);
    if (std::memcmp(magic, "KHIST05", 7) != 0)
        throw std::runtime_error("bad magic in Phase 0 file: " + path);

    std::uint64_t n_recurrent = 0, n_reads = 0, n_bases = 0;
    std::int32_t cov = 0, het = 0;
    double z = 0;
    std::uint32_t k = 0;
    f.read(reinterpret_cast<char*>(&n_recurrent), sizeof(n_recurrent));
    f.read(reinterpret_cast<char*>(&n_reads), sizeof(n_reads));
    f.read(reinterpret_cast<char*>(&n_bases), sizeof(n_bases));
    f.read(reinterpret_cast<char*>(&cov), sizeof(cov));
    f.read(reinterpret_cast<char*>(&het), sizeof(het));
    f.read(reinterpret_cast<char*>(&z), sizeof(z));
    f.read(reinterpret_cast<char*>(&k), sizeof(k));

    std::uint32_t hist_n = 0;
    f.read(reinterpret_cast<char*>(&hist_n), sizeof(hist_n));
    f.seekg(sizeof(std::uint64_t) * hist_n, std::ios::cur);

    out_counter.reserve(static_cast<std::size_t>(n_recurrent));
    for (std::uint64_t i = 0; i < n_recurrent; ++i) {
        std::uint64_t key = 0;
        std::uint32_t cnt = 0;
        f.read(reinterpret_cast<char*>(&key), sizeof(key));
        f.read(reinterpret_cast<char*>(&cnt), sizeof(cnt));
        out_counter.emplace(key, cnt);
    }
}

void load_kmer_counter_sorted_arrays(
    const std::string& path,
    std::vector<std::uint64_t>& out_keys,
    std::vector<std::uint32_t>& out_counts) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open Phase 0 dump: " + path);
    char magic[8];
    f.read(magic, 8);
    if (std::memcmp(magic, "KHIST05", 7) != 0)
        throw std::runtime_error("bad magic in Phase 0 file: " + path);

    std::uint64_t n_recurrent = 0, n_reads = 0, n_bases = 0;
    std::int32_t cov = 0, het = 0;
    double z = 0;
    std::uint32_t k = 0;
    f.read(reinterpret_cast<char*>(&n_recurrent), sizeof(n_recurrent));
    f.read(reinterpret_cast<char*>(&n_reads), sizeof(n_reads));
    f.read(reinterpret_cast<char*>(&n_bases), sizeof(n_bases));
    f.read(reinterpret_cast<char*>(&cov), sizeof(cov));
    f.read(reinterpret_cast<char*>(&het), sizeof(het));
    f.read(reinterpret_cast<char*>(&z), sizeof(z));
    f.read(reinterpret_cast<char*>(&k), sizeof(k));

    std::uint32_t hist_n = 0;
    f.read(reinterpret_cast<char*>(&hist_n), sizeof(hist_n));
    f.seekg(sizeof(std::uint64_t) * hist_n, std::ios::cur);

    // Memory profile (the whole point of this function):
    //   keys:   8 B × n_recurrent
    //   counts: 4 B × n_recurrent
    //   Total:  12 B/entry — vs. ~48 B/entry in unordered_map. For
    //   2 Gbil recurrent k-mers (HG002 HiFi WG) this saves ~70 GB and
    //   is the only path under the 140 GB SLURM cap.
    out_keys.clear();
    out_counts.clear();
    out_keys.resize(static_cast<std::size_t>(n_recurrent));
    out_counts.resize(static_cast<std::size_t>(n_recurrent));
    // Stream-read in 1 M-entry chunks to keep buffer overhead small.
    constexpr std::uint64_t kChunk = 1'000'000ULL;
    std::vector<char> buf(static_cast<std::size_t>(kChunk) * 12);
    std::uint64_t done = 0;
    while (done < n_recurrent) {
        std::uint64_t take = std::min(kChunk, n_recurrent - done);
        f.read(buf.data(), static_cast<std::streamsize>(take * 12));
        for (std::uint64_t i = 0; i < take; ++i) {
            std::memcpy(&out_keys[done + i],
                        buf.data() + i * 12, 8);
            std::memcpy(&out_counts[done + i],
                        buf.data() + i * 12 + 8, 4);
        }
        done += take;
    }

    // Slot order from the GPU is hash-order, not key-order. Sort once
    // here so downstream consumers (binary-search-by-key, GPU launchers
    // that expect sorted input) need no further reorganization.
    //
    // Sort approach: pack (key, count) into a single u128-like pair
    // (we just use std::vector<std::uint64_t> as packed representation
    //  with key in upper 64 bits and count in lower 32 bits — but that
    //  loses precision; instead we sort a struct-of-arrays via a 12-byte
    //  packed entry in a single vector and split after).
    //
    // Practical: use std::sort with std::execution::par on a packed
    // 12-byte buffer. For 2 G entries on 8 cores: ~30-90 s.
    //
    // Memory profile:
    //   packed (12 B × N) replaces (keys 8 + counts 4 = 12 B × N).
    //   Total during sort: 24 GB packed (no extra buffer).
    std::fprintf(stderr,
        "[wg/p0-load] sorting %zu (key,count) entries…\n", out_keys.size());
    std::fflush(stderr);
    auto t_sort = std::chrono::steady_clock::now();

    struct Pkc { std::uint64_t key; std::uint32_t cnt; };
    std::vector<Pkc> packed(out_keys.size());
    for (std::size_t i = 0; i < out_keys.size(); ++i) {
        packed[i].key = out_keys[i];
        packed[i].cnt = out_counts[i];
    }
    // Free the unsorted parallel arrays now (we keep the packed copy);
    // peak memory is just packed (24 GB) instead of 36 GB.
    std::vector<std::uint64_t>().swap(out_keys);
    std::vector<std::uint32_t>().swap(out_counts);

#ifdef _OPENMP
    // Crude parallel quicksort via OpenMP tasks: split the array into
    // 8 chunks, sort each, then merge with std::inplace_merge.
    const int nT = std::min(8, omp_get_max_threads());
    const std::size_t chunk = packed.size() / nT;
    #pragma omp parallel num_threads(nT)
    {
        const int tid = omp_get_thread_num();
        const std::size_t lo = tid * chunk;
        const std::size_t hi = (tid == nT - 1) ? packed.size()
                                               : (tid + 1) * chunk;
        std::sort(packed.begin() + lo, packed.begin() + hi,
            [](const Pkc& a, const Pkc& b) { return a.key < b.key; });
    }
    // Sequential O(log_2 nT) merges: pairs of adjacent chunks.
    for (std::size_t step = 1; step < static_cast<std::size_t>(nT); step *= 2) {
        for (std::size_t i = 0; i + step < static_cast<std::size_t>(nT); i += 2 * step) {
            std::size_t lo = i * chunk;
            std::size_t mid = (i + step) * chunk;
            std::size_t hi = std::min((i + 2 * step) * chunk, packed.size());
            std::inplace_merge(packed.begin() + lo,
                               packed.begin() + mid,
                               packed.begin() + hi,
                [](const Pkc& a, const Pkc& b) { return a.key < b.key; });
        }
    }
#else
    std::sort(packed.begin(), packed.end(),
        [](const Pkc& a, const Pkc& b) { return a.key < b.key; });
#endif

    const double sort_el = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_sort).count();
    std::fprintf(stderr,
        "[wg/p0-load] sort done in %.1fs, splitting back to parallel arrays\n",
        sort_el);
    std::fflush(stderr);

    out_keys.resize(packed.size());
    out_counts.resize(packed.size());
    for (std::size_t i = 0; i < packed.size(); ++i) {
        out_keys[i] = packed[i].key;
        out_counts[i] = packed[i].cnt;
    }
    std::vector<Pkc>().swap(packed);
}

}  // namespace branch::wg
