// BRANCH — synthetic CNV-bubble fixture generator.
//
// BRANCH's core job is detecting copy-number variation in assembly
// graph branches. This fixture builds a MULTI-ALLELIC BUBBLE at a
// known locus, with one 2-3 kb cassette repeated a different number
// of times on each path:
//
//   Flank-L ─┬─► (no cassette, 0 copies) ─┐
//            ├─► cassette × 2            ─┼─► Flank-R
//            ├─► cassette × 4            ─┤
//            └─► cassette × 6            ─┘
//
// Each allele is sampled with a configurable VAF (the fraction of the
// read pool drawn from that path) so the test can assert BRANCH's
// per-branch VAF estimate and copy-count inference.
//
// If BRANCH cannot see a 4-way bubble and does not recover the copy
// counts {0, 2, 4, 6} with VAFs close to the input mix, then real
// low-frequency CNV data is out of reach.
//
// Usage:
//   generate_branch_fixture <out.fastq> <seed> <total_reads> <read_len>
//                           [--ref OUT.fa] [--meta OUT.tsv]
//                           [--cassette-bp N]            default 2500
//                           [--vaf-mix V0,V2,V4,V6]      default 0.10,0.50,0.25,0.15
//
// Meta TSV columns:
//   allele  copy_count  vaf_target  n_reads_generated  cassette_len_bp

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace {

std::uint64_t splitmix64(std::uint64_t x) noexcept {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

std::string random_dna(std::uint64_t seed, std::size_t len) {
    std::mt19937_64 rng(splitmix64(seed));
    std::uniform_int_distribution<int> base(0, 3);
    const char alphabet[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.resize(len);
    for (std::size_t i = 0; i < len; ++i) s[i] = alphabet[base(rng)];
    return s;
}

void die(const char* msg) {
    std::cerr << "generate_branch_fixture: " << msg << '\n';
    std::exit(2);
}

std::vector<double> parse_vafs(const std::string& s) {
    std::vector<double> out;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        out.push_back(std::atof(tok.c_str()));
    }
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr <<
            "usage: generate_branch_fixture <out.fastq> <seed> <total_reads> "
            "<read_len> [--ref OUT.fa] [--meta OUT.tsv]\n"
            "       [--cassette-bp N] [--vaf-mix V0,V2,V4,V6]\n";
        return 2;
    }
    const std::string out_path = argv[1];
    const std::uint64_t seed = std::strtoull(argv[2], nullptr, 10);
    const std::size_t total_reads = std::strtoull(argv[3], nullptr, 10);
    const std::size_t read_len = std::strtoull(argv[4], nullptr, 10);

    std::string ref_out, meta_out;
    std::size_t cassette_bp = 2500;
    std::vector<double> vaf_mix = {0.10, 0.50, 0.25, 0.15};  // 0, 2, 4, 6 copies
    const std::vector<int> copy_counts = {0, 2, 4, 6};

    for (int i = 5; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--ref" && i + 1 < argc) ref_out = argv[++i];
        else if (a == "--meta" && i + 1 < argc) meta_out = argv[++i];
        else if (a == "--cassette-bp" && i + 1 < argc) {
            cassette_bp = std::strtoull(argv[++i], nullptr, 10);
        } else if (a == "--vaf-mix" && i + 1 < argc) {
            vaf_mix = parse_vafs(argv[++i]);
            if (vaf_mix.size() != copy_counts.size())
                die("--vaf-mix needs exactly 4 comma-separated values");
        } else {
            std::cerr << "unknown arg: " << a << '\n';
            return 2;
        }
    }

    if (total_reads < 50) die("total_reads must be >= 50 for meaningful coverage");
    if (read_len < 5000) die("read_len should be >= 5000 to span the bubble");
    if (cassette_bp < 500) die("cassette-bp should be >= 500");

    double vaf_sum = 0.0;
    for (double v : vaf_mix) vaf_sum += v;
    if (vaf_sum <= 0.0) die("--vaf-mix must sum to > 0");

    std::mt19937_64 rng(splitmix64(seed));

    // Generate cassette + two flanks with enough length that every read
    // fully spans the bubble (flank_L + max_cassette_block + flank_R).
    const std::size_t max_copies = 6;
    const std::size_t max_allele_body = max_copies * cassette_bp;
    const std::size_t flank_bp = 2000;
    const std::size_t full_len = flank_bp + max_allele_body + flank_bp;
    if (read_len < full_len) {
        std::cerr << "generate_branch_fixture: read_len ("
                  << read_len << ") < full_len (" << full_len
                  << "), reads won't span the biggest bubble.\n";
    }

    const std::string cassette = random_dna(splitmix64(seed ^ 0xCACECACEULL), cassette_bp);
    const std::string flank_L  = random_dna(splitmix64(seed ^ 0xF1A4CE1EULL), flank_bp);
    const std::string flank_R  = random_dna(splitmix64(seed ^ 0xF1A4CE2EULL), flank_bp);

    // Build the four alleles.
    auto make_allele = [&](int k) {
        std::string s;
        s.reserve(flank_bp + k * cassette_bp + flank_bp);
        s += flank_L;
        for (int i = 0; i < k; ++i) s += cassette;
        s += flank_R;
        return s;
    };
    std::vector<std::string> alleles;
    alleles.reserve(copy_counts.size());
    for (int k : copy_counts) alleles.push_back(make_allele(k));

    // Sample reads from allele pool weighted by vaf_mix.
    std::discrete_distribution<int> allele_pick(vaf_mix.begin(), vaf_mix.end());

    std::ofstream ofs(out_path);
    if (!ofs) die("cannot open output fastq");

    std::vector<std::size_t> reads_per_allele(copy_counts.size(), 0);
    for (std::size_t r = 0; r < total_reads; ++r) {
        const int a_idx = allele_pick(rng);
        const auto& hap = alleles[a_idx];
        const std::size_t take = std::min(read_len, hap.size());
        const std::size_t max_start = hap.size() - take;
        std::uniform_int_distribution<std::size_t> start_dist(0, max_start);
        const std::size_t start = start_dist(rng);
        const std::string qual(take, '?');
        ofs << "@syn_" << r << "_k" << copy_counts[a_idx] << '\n';
        ofs.write(hap.data() + start, static_cast<std::streamsize>(take));
        ofs << "\n+\n" << qual << '\n';
        ++reads_per_allele[a_idx];
    }
    if (!ofs) die("write error on fastq");

    if (!ref_out.empty()) {
        // Emit the 2-copy allele as the "reference" — a sensible
        // diploid backbone. The fixture meta still names each path.
        std::ofstream rf(ref_out);
        if (!rf) die("cannot open --ref output");
        const auto& ref_allele = alleles[1];  // 2-copy
        rf << ">synChr\n";
        for (std::size_t p = 0; p < ref_allele.size(); p += 80) {
            rf.write(ref_allele.data() + p,
                     static_cast<std::streamsize>(
                         std::min<std::size_t>(80, ref_allele.size() - p)));
            rf << '\n';
        }
        if (!rf) die("write error on --ref");
    }

    if (!meta_out.empty()) {
        std::ofstream m(meta_out);
        if (!m) die("cannot open --meta output");
        m << "#allele\tcopy_count\tvaf_target\tn_reads\tcassette_len_bp\n";
        for (std::size_t i = 0; i < copy_counts.size(); ++i) {
            m << "allele_" << i << '\t'
              << copy_counts[i] << '\t'
              << (vaf_mix[i] / vaf_sum) << '\t'
              << reads_per_allele[i] << '\t'
              << cassette_bp << '\n';
        }
        if (!m) die("write error on --meta");
    }

    std::cerr << "generate_branch_fixture: " << total_reads
              << " reads, 4 alleles copies={0,2,4,6} cassette="
              << cassette_bp << " bp\n";
    for (std::size_t i = 0; i < copy_counts.size(); ++i) {
        std::cerr << "  allele k=" << copy_counts[i]
                  << " target_vaf=" << (vaf_mix[i] / vaf_sum)
                  << " reads=" << reads_per_allele[i] << '\n';
    }
    return 0;
}
