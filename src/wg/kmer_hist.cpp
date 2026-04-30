// BRANCH v0.5 — Phase 0 implementation.

#include "wg/kmer_hist.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>

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

KmerHistBuilder::KmerHistBuilder(const ::branch::graph::TechProfile& profile,
                                 std::size_t expected_n_distinct_kmers)
    : profile_(profile),
      bloom_(expected_n_distinct_kmers, /*bits_per_key=*/6),
      counter_(expected_n_distinct_kmers / 10) {  // ~10x reduction post-bloom
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

        std::uint64_t h = canonical_kmer_hash(kmer, k);
        if (pass2) {
            if (bloom_.maybe_contains(h)) counter_.add(h, 1);
        } else {
            bloom_.add(h);
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
    KmerHistResult r;
    r.n_reads = n_reads_;
    r.n_bases = n_bases_;
    r.n_recurrent_kmers = counter_.size();

    // Build histogram[count] over recurrent k-mers.
    constexpr std::size_t HIST_MAX = 512;
    r.histogram.assign(HIST_MAX, 0);
    counter_.for_each([&](std::uint64_t /*key*/, std::uint32_t cnt) {
        std::size_t bin = std::min<std::size_t>(cnt, HIST_MAX - 1);
        ++r.histogram[bin];
    });

    // Find peaks. We only look in [4, 4*expected_coverage] because
    // counts of 2-3 are mostly Q-error survivors that the bloom let
    // through, not informative signal.
    int lo = 4;
    int hi = std::min<int>(static_cast<int>(HIST_MAX) - 1, expected_coverage * 4);
    if (lo >= hi) return r;

    auto smooth = [&](int x) -> double {
        // 5-bin moving average for noise suppression.
        double s = 0; int n = 0;
        for (int d = -2; d <= 2; ++d) {
            int xx = x + d;
            if (xx >= 0 && xx < (int)r.histogram.size()) {
                s += static_cast<double>(r.histogram[xx]); ++n;
            }
        }
        return n ? s / n : 0.0;
    };

    // First peak — highest in [lo, hi]
    int peak1 = lo;
    double v1 = smooth(lo);
    for (int x = lo + 1; x <= hi; ++x) {
        double v = smooth(x);
        if (v > v1) { v1 = v; peak1 = x; }
    }
    r.detected_coverage = peak1;

    // Second peak — highest in [lo, peak1 * 0.7] (lower than primary)
    int peak2 = 0;
    double v2 = 0.0;
    int hi2 = std::max(lo, static_cast<int>(peak1 * 0.7));
    for (int x = lo; x <= hi2; ++x) {
        double v = smooth(x);
        if (v > v2) { v2 = v; peak2 = x; }
    }

    // Validate the second peak via z-separation between the two peaks
    // and the valley between them.
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

void KmerHistBuilder::write_to(const std::string& path,
                               const KmerHistResult& result) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open: " + path);

    // Magic + version
    constexpr char MAGIC[8] = {'K','H','I','S','T','0','5','\0'};
    f.write(MAGIC, 8);

    // Header
    std::uint64_t n_recurrent = result.n_recurrent_kmers;
    std::uint64_t n_reads = result.n_reads;
    std::uint64_t n_bases = result.n_bases;
    std::int32_t  cov = result.detected_coverage;
    std::int32_t  het = result.detected_het_coverage;
    double        z   = result.bimodality_z;
    std::uint32_t k   = static_cast<std::uint32_t>(profile_.minimizer_k);
    f.write(reinterpret_cast<const char*>(&n_recurrent), sizeof(n_recurrent));
    f.write(reinterpret_cast<const char*>(&n_reads), sizeof(n_reads));
    f.write(reinterpret_cast<const char*>(&n_bases), sizeof(n_bases));
    f.write(reinterpret_cast<const char*>(&cov), sizeof(cov));
    f.write(reinterpret_cast<const char*>(&het), sizeof(het));
    f.write(reinterpret_cast<const char*>(&z), sizeof(z));
    f.write(reinterpret_cast<const char*>(&k), sizeof(k));

    // Histogram
    std::uint32_t hist_n = static_cast<std::uint32_t>(result.histogram.size());
    f.write(reinterpret_cast<const char*>(&hist_n), sizeof(hist_n));
    f.write(reinterpret_cast<const char*>(result.histogram.data()),
            result.histogram.size() * sizeof(std::uint64_t));

    // Counter dump: (key, count) records, one per recurrent k-mer.
    counter_.for_each([&](std::uint64_t key, std::uint32_t cnt) {
        f.write(reinterpret_cast<const char*>(&key), sizeof(key));
        f.write(reinterpret_cast<const char*>(&cnt), sizeof(cnt));
    });

    f.close();
}

}  // namespace branch::wg
