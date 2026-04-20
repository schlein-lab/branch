// BRANCH — multi-locus stress fixture.
//
// Four independent loci with distinct structural challenges, all in one
// FASTQ, to stress classifier + conservation + low-band VAF gate at once:
//
//   Locus A  — clean 4-way branch (copies {0,2,4,6}, vafs {0.10,0.50,0.25,0.15})
//   Locus B  — true duplication (identical flanks + twin cassette, 2x depth)
//   Locus C  — mixed: duplication backbone + one low-VAF (5%) true branch allele
//   Locus D  — nested bubble: outer 2-way, inner 3-way inside one arm
//
// HiFi-style substitution noise at 0.1% per base is applied to every read
// so overlap + classifier see realistic mismatch rates rather than
// synthetic perfection.
//
// Each locus uses a disjoint seed family so cassette and flank sequences
// do not collide between loci. Read names carry the locus tag so the
// fixture meta TSV can be diff'd against BED output.
//
// Usage:
//   generate_branch_complex_fixture <out.fastq> <seed> <reads_per_locus> <read_len>
//                                   [--ref OUT.fa] [--meta OUT.tsv]
//                                   [--err-rate R]     default 0.001
//                                   [--low-vaf V]      default 0.05 (for Locus C)
//
// Meta TSV columns:
//   locus  structure  allele  copy_count  vaf_target  n_reads

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

constexpr std::size_t kFlankBp    = 2000;
constexpr std::size_t kCassetteBp = 2500;

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
    std::string s(len, 'A');
    for (auto& c : s) c = alphabet[base(rng)];
    return s;
}

void mutate_in_place(std::string& seq, std::mt19937_64& rng, double err_rate) {
    if (err_rate <= 0.0) return;
    std::uniform_real_distribution<double> u(0.0, 1.0);
    std::uniform_int_distribution<int> base(0, 3);
    const char alphabet[4] = {'A', 'C', 'G', 'T'};
    for (auto& c : seq) {
        if (u(rng) < err_rate) {
            char pick = alphabet[base(rng)];
            while (pick == c) pick = alphabet[base(rng)];
            c = pick;
        }
    }
}

void die(const char* msg) {
    std::cerr << "generate_branch_complex_fixture: " << msg << '\n';
    std::exit(2);
}

struct Allele {
    std::string locus;
    std::string structure;  // branch | duplication | mixed | nested
    std::string name;
    int copy_count;
    double vaf_target;
    std::string sequence;
    std::size_t n_reads;  // populated after sampling
};

// Build locus A: clean 4-way branch, copies {0,2,4,6}, vafs {0.10,0.50,0.25,0.15}.
std::vector<Allele> build_locus_A(std::uint64_t seed) {
    const std::string cassette = random_dna(splitmix64(seed ^ 0xA1ULL), kCassetteBp);
    const std::string flank_L  = random_dna(splitmix64(seed ^ 0xA2ULL), kFlankBp);
    const std::string flank_R  = random_dna(splitmix64(seed ^ 0xA3ULL), kFlankBp);
    const std::vector<int> copies = {0, 2, 4, 6};
    const std::vector<double> vafs = {0.10, 0.50, 0.25, 0.15};
    std::vector<Allele> out;
    for (std::size_t i = 0; i < copies.size(); ++i) {
        std::string body;
        for (int k = 0; k < copies[i]; ++k) body += cassette;
        out.push_back({"A", "branch", "A_k" + std::to_string(copies[i]),
                       copies[i], vafs[i], flank_L + body + flank_R, 0});
    }
    return out;
}

// Build locus B: duplication — two identical twin cassettes between shared flanks.
// Both "alleles" have the same sequence; reads are assigned to allele_0/allele_1
// by a coin flip to simulate 2x diploid depth with no structural divergence.
std::vector<Allele> build_locus_B(std::uint64_t seed) {
    const std::string cassette = random_dna(splitmix64(seed ^ 0xB1ULL), kCassetteBp);
    const std::string flank_L  = random_dna(splitmix64(seed ^ 0xB2ULL), kFlankBp);
    const std::string flank_R  = random_dna(splitmix64(seed ^ 0xB3ULL), kFlankBp);
    const std::string seq = flank_L + cassette + cassette + flank_R;
    return {
        {"B", "duplication", "B_hap0", 2, 0.50, seq, 0},
        {"B", "duplication", "B_hap1", 2, 0.50, seq, 0},
    };
}

// Build locus C: mixed — duplication backbone (both haps share 2-cassette body)
// plus one rare branch allele with 3 cassettes at configurable low VAF.
std::vector<Allele> build_locus_C(std::uint64_t seed, double low_vaf) {
    const std::string cassette = random_dna(splitmix64(seed ^ 0xC1ULL), kCassetteBp);
    const std::string flank_L  = random_dna(splitmix64(seed ^ 0xC2ULL), kFlankBp);
    const std::string flank_R  = random_dna(splitmix64(seed ^ 0xC3ULL), kFlankBp);
    const std::string backbone = flank_L + cassette + cassette + flank_R;
    const std::string rare     = flank_L + cassette + cassette + cassette + flank_R;
    const double backbone_vaf = (1.0 - low_vaf) / 2.0;
    return {
        {"C", "mixed", "C_hap0", 2, backbone_vaf, backbone, 0},
        {"C", "mixed", "C_hap1", 2, backbone_vaf, backbone, 0},
        {"C", "mixed", "C_rare", 3, low_vaf,      rare,     0},
    };
}

// Build locus D: nested bubble — outer 2-way (0 vs 1 cassette), and within the
// 1-cassette arm an inner 3-way that varies the cassette body SNP pattern.
std::vector<Allele> build_locus_D(std::uint64_t seed) {
    const std::string flank_L  = random_dna(splitmix64(seed ^ 0xD1ULL), kFlankBp);
    const std::string flank_R  = random_dna(splitmix64(seed ^ 0xD2ULL), kFlankBp);
    const std::string cassette0 = random_dna(splitmix64(seed ^ 0xD3ULL), kCassetteBp);
    std::string cassette_v1 = cassette0;
    std::string cassette_v2 = cassette0;
    std::string cassette_v3 = cassette0;
    cassette_v1[kCassetteBp / 3] = 'A';
    cassette_v2[kCassetteBp / 3] = 'C';
    cassette_v3[kCassetteBp / 3] = 'T';
    return {
        {"D", "nested", "D_outer_ref", 0, 0.40, flank_L + flank_R, 0},
        {"D", "nested", "D_inner_A",   1, 0.25, flank_L + cassette_v1 + flank_R, 0},
        {"D", "nested", "D_inner_C",   1, 0.20, flank_L + cassette_v2 + flank_R, 0},
        {"D", "nested", "D_inner_T",   1, 0.15, flank_L + cassette_v3 + flank_R, 0},
    };
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr <<
            "usage: generate_branch_complex_fixture <out.fastq> <seed> "
            "<reads_per_locus> <read_len> [--ref OUT.fa] [--meta OUT.tsv] "
            "[--err-rate R] [--low-vaf V]\n";
        return 2;
    }
    const std::string out_path   = argv[1];
    const std::uint64_t seed     = std::strtoull(argv[2], nullptr, 10);
    const std::size_t per_locus  = std::strtoull(argv[3], nullptr, 10);
    const std::size_t read_len   = std::strtoull(argv[4], nullptr, 10);

    std::string ref_out, meta_out;
    double err_rate = 0.001;
    double low_vaf  = 0.05;
    for (int i = 5; i < argc; ++i) {
        std::string_view a = argv[i];
        if (a == "--ref" && i + 1 < argc) ref_out  = argv[++i];
        else if (a == "--meta" && i + 1 < argc) meta_out = argv[++i];
        else if (a == "--err-rate" && i + 1 < argc) err_rate = std::atof(argv[++i]);
        else if (a == "--low-vaf"  && i + 1 < argc) low_vaf  = std::atof(argv[++i]);
        else { std::cerr << "unknown arg: " << a << '\n'; return 2; }
    }
    if (per_locus < 20) die("reads_per_locus must be >= 20");
    if (read_len  < 5000) die("read_len should be >= 5000 to span a bubble");
    if (err_rate < 0.0 || err_rate > 0.1) die("--err-rate out of [0, 0.1]");
    if (low_vaf  < 0.0 || low_vaf  > 0.5) die("--low-vaf out of [0, 0.5]");

    std::vector<Allele> all;
    auto extend = [&](std::vector<Allele> v) { for (auto& a : v) all.push_back(std::move(a)); };
    extend(build_locus_A(seed));
    extend(build_locus_B(seed + 1));
    extend(build_locus_C(seed + 2, low_vaf));
    extend(build_locus_D(seed + 3));

    std::mt19937_64 rng(splitmix64(seed ^ 0xFEEDFACEULL));
    std::ofstream ofs(out_path);
    if (!ofs) die("cannot open output fastq");

    std::size_t read_counter = 0;
    // Per locus, draw `per_locus` reads weighted by this locus's vaf_targets.
    // Alleles group by locus tag.
    std::vector<std::string> locus_tags = {"A", "B", "C", "D"};
    for (const auto& tag : locus_tags) {
        std::vector<std::size_t> idx;
        std::vector<double> w;
        for (std::size_t i = 0; i < all.size(); ++i) {
            if (all[i].locus == tag) { idx.push_back(i); w.push_back(all[i].vaf_target); }
        }
        std::discrete_distribution<int> pick(w.begin(), w.end());
        for (std::size_t r = 0; r < per_locus; ++r) {
            const std::size_t ai = idx[pick(rng)];
            auto& a = all[ai];
            if (a.sequence.size() < 100) continue;  // degenerate
            const std::size_t take = std::min(read_len, a.sequence.size());
            const std::size_t max_start = a.sequence.size() - take;
            std::uniform_int_distribution<std::size_t> start_dist(0, max_start);
            const std::size_t start = start_dist(rng);
            std::string slice = a.sequence.substr(start, take);
            mutate_in_place(slice, rng, err_rate);
            const std::string qual(slice.size(), '?');
            ofs << "@syn_" << read_counter++ << "_L" << a.locus
                << "_" << a.name << '\n'
                << slice << "\n+\n" << qual << '\n';
            ++a.n_reads;
        }
    }
    if (!ofs) die("write error on fastq");

    if (!ref_out.empty()) {
        std::ofstream rf(ref_out);
        if (!rf) die("cannot open --ref output");
        // Emit one representative contig per locus: A=k2 allele, B=dup, C=backbone, D=innerA.
        const auto emit = [&](const std::string& name, const std::string& seq) {
            rf << '>' << name << '\n';
            for (std::size_t p = 0; p < seq.size(); p += 80) {
                rf.write(seq.data() + p,
                         static_cast<std::streamsize>(
                             std::min<std::size_t>(80, seq.size() - p)));
                rf << '\n';
            }
        };
        for (const auto& a : all) {
            if ((a.locus == "A" && a.name == "A_k2") ||
                (a.locus == "B" && a.name == "B_hap0") ||
                (a.locus == "C" && a.name == "C_hap0") ||
                (a.locus == "D" && a.name == "D_inner_A")) {
                emit("ref_" + a.locus, a.sequence);
            }
        }
        if (!rf) die("write error on --ref");
    }

    if (!meta_out.empty()) {
        std::ofstream m(meta_out);
        if (!m) die("cannot open --meta output");
        m << "#locus\tstructure\tallele\tcopy_count\tvaf_target\tn_reads\n";
        for (const auto& a : all) {
            m << a.locus << '\t' << a.structure << '\t' << a.name << '\t'
              << a.copy_count << '\t' << a.vaf_target << '\t' << a.n_reads << '\n';
        }
        if (!m) die("write error on --meta");
    }

    std::cerr << "generate_branch_complex_fixture: wrote " << read_counter
              << " reads across " << locus_tags.size() << " loci "
              << "(err_rate=" << err_rate << ", low_vaf=" << low_vaf << ")\n";
    for (const auto& a : all) {
        std::cerr << "  " << a.locus << " " << a.structure << " " << a.name
                  << " cn=" << a.copy_count << " vaf=" << a.vaf_target
                  << " reads=" << a.n_reads << '\n';
    }
    return 0;
}
