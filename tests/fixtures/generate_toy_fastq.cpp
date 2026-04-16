// BRANCH v0.2 — deterministic toy FASTQ generator for E2E tests.
//
// Standalone C++20 program. Emits `n_reads` HiFi-like reads, all drawn
// as overlapping windows of a single 5kb synthetic random reference
// (fixed by seed). Reads are staggered so adjacent reads share
// `mean_overlap` bases — enough for the minimizer sketcher + overlap
// detector to find ~n_reads-1 adjacency overlaps and build a near-
// linear graph.
//
// Usage:
//   generate_toy_fastq <output.fastq> <seed> <n_reads> <read_len> <mean_overlap>
//
// Output is plain-text FASTQ (no gzip). Quality = '?' for every base
// (Phred 30). No N's in the sequence — the reference is drawn from
// the ACGT alphabet only so the minimizer window never resets.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

namespace {

// splitmix64 — same mix function used inside BRANCH's hasher, so the
// fixture is reproducible with the existing test infra.
std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

std::string random_reference(std::uint64_t seed, std::size_t len) {
    std::mt19937_64 rng(splitmix64(seed));
    std::uniform_int_distribution<int> base(0, 3);
    const char alphabet[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.resize(len);
    for (std::size_t i = 0; i < len; ++i) {
        s[i] = alphabet[base(rng)];
    }
    return s;
}

void die(const char* msg) {
    std::cerr << "generate_toy_fastq: " << msg << "\n";
    std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: generate_toy_fastq <output.fastq> <seed> "
                     "<n_reads> <read_len> <mean_overlap>\n";
        return 2;
    }
    const std::string out_path = argv[1];
    const std::uint64_t seed = std::strtoull(argv[2], nullptr, 10);
    const std::size_t n_reads = std::strtoull(argv[3], nullptr, 10);
    const std::size_t read_len = std::strtoull(argv[4], nullptr, 10);
    const std::size_t mean_overlap = std::strtoull(argv[5], nullptr, 10);

    if (n_reads == 0) die("n_reads must be >= 1");
    if (read_len == 0) die("read_len must be >= 1");
    if (mean_overlap >= read_len) die("mean_overlap must be < read_len");

    // Reference length chosen so every read window fits:
    //   stride = read_len - mean_overlap
    //   ref_len >= read_len + (n_reads - 1) * stride
    // We pad to at least 5000 bp so ambient entropy is high.
    const std::size_t stride = read_len - mean_overlap;
    std::size_t ref_len = read_len + (n_reads > 0 ? (n_reads - 1) * stride : 0);
    if (ref_len < 5000) ref_len = 5000;

    const std::string ref = random_reference(seed, ref_len);

    std::ofstream ofs(out_path);
    if (!ofs) die("cannot open output");

    const std::string qual(read_len, '?');
    for (std::size_t i = 0; i < n_reads; ++i) {
        const std::size_t start = i * stride;
        // Paranoid bounds-check — ref is sized to fit by construction.
        if (start + read_len > ref.size()) die("internal: read extends past reference");
        ofs << '@' << "toy_read_" << i << '\n';
        ofs.write(ref.data() + start, static_cast<std::streamsize>(read_len));
        ofs << '\n' << "+\n" << qual << '\n';
    }
    if (!ofs) die("write error");
    return 0;
}
