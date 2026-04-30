// BRANCH v0.5 — `branch wg` orchestrator subcommand.
//
// Runs phases 0 → 5 in order and writes <out-prefix>.{fa,bed,gfa}.

#include "wg_cmd.hpp"

#include "wg/kmer_hist.hpp"
#include "wg/haplotype_router.hpp"
#include "wg/master_tiler.hpp"
#include "wg/master_curator.hpp"
#include "wg/branch_attacher.hpp"
#include "wg/vpcr_wg.hpp"
#include "wg/orphan_clusterer.hpp"
#include "wg/vpf_writer.hpp"
#include "graph/kmer_sketch.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#include <zlib.h>

namespace branch::cli {

namespace {

void print_wg_usage(std::ostream& os) {
    os << "branch wg — whole-genome lossless haplotype-master-tiling assembly\n"
       << "\n"
       << "Usage:\n"
       << "  branch wg --reads <fastq.gz> --read-tech {hifi,ont} \\\n"
       << "            --out-prefix <prefix> [--threads N]\n"
       << "            [--phase 0,1,2,3,4,5] [--no-vpcr] [--no-ref-annotation]\n"
       << "\n"
       << "Outputs:\n"
       << "  <prefix>.fa     FASTA: hap1 + hap2 + branches + orphans\n"
       << "  <prefix>.bed    BED: branch + vPCR positions + ref-class\n"
       << "  <prefix>.gfa    extended GFA1.2 with v0.5 metadata blocks\n";
}

struct WgArgs {
    std::string reads;
    std::string read_tech = "hifi";
    std::string out_prefix;
    int threads = 4;
    int expected_coverage = 0;       // 0 → auto-detect
    bool do_vpcr = true;
    bool do_ref = true;
    std::vector<int> phases = {0, 1, 2, 3, 4, 5};
    bool ok = false;
    std::string err;
};

WgArgs parse_wg_args(int argc, char** argv) {
    WgArgs a;
    for (int i = 1; i < argc; ++i) {
        std::string_view k = argv[i];
        auto need = [&](const char* n) -> const char* {
            if (i + 1 >= argc) { a.err = std::string("missing value for ") + n; return nullptr; }
            return argv[++i];
        };
        if (k == "--reads")           { auto v = need("--reads");      if (!v) return a; a.reads = v; }
        else if (k == "--read-tech")  { auto v = need("--read-tech");  if (!v) return a; a.read_tech = v; }
        else if (k == "--out-prefix") { auto v = need("--out-prefix"); if (!v) return a; a.out_prefix = v; }
        else if (k == "--threads")    { auto v = need("--threads");    if (!v) return a; a.threads = std::atoi(v); }
        else if (k == "--expected-coverage") { auto v = need("--expected-coverage"); if (!v) return a; a.expected_coverage = std::atoi(v); }
        else if (k == "--no-vpcr")    { a.do_vpcr = false; }
        else if (k == "--no-ref-annotation") { a.do_ref = false; }
        else if (k == "--help" || k == "-h") { a.err = "HELP"; return a; }
        else { a.err = std::string("unknown arg: ") + std::string(k); return a; }
    }
    if (a.reads.empty()) { a.err = "--reads is required"; return a; }
    if (a.out_prefix.empty()) { a.err = "--out-prefix is required"; return a; }
    a.ok = true;
    return a;
}

bool stream_fastq_gz(const std::string& path,
                     std::vector<std::pair<std::string, std::string>>& out) {
    gzFile f = gzopen(path.c_str(), "rb");
    if (!f) return false;
    constexpr int BUF = 1 << 20;
    std::string buf(BUF, '\0');
    std::string carry;
    int line_kind = 0;  // 0=name, 1=seq, 2=plus, 3=qual
    std::string cur_name;
    std::string cur_seq;
    while (true) {
        int n = gzread(f, buf.data(), BUF);
        if (n <= 0) break;
        carry.append(buf.data(), buf.data() + n);
        std::size_t pos = 0;
        while (true) {
            auto nl = carry.find('\n', pos);
            if (nl == std::string::npos) { carry.erase(0, pos); break; }
            std::string_view line(carry.data() + pos, nl - pos);
            switch (line_kind) {
                case 0: if (!line.empty() && line[0] == '@') {
                    cur_name.assign(line.data() + 1, line.size() - 1);
                    auto sp = cur_name.find_first_of(" \t");
                    if (sp != std::string::npos) cur_name.resize(sp);
                } break;
                case 1: cur_seq.assign(line.data(), line.size()); break;
                case 2: break;
                case 3:
                    out.emplace_back(std::move(cur_name), std::move(cur_seq));
                    cur_name.clear(); cur_seq.clear();
                    break;
            }
            line_kind = (line_kind + 1) & 3;
            pos = nl + 1;
        }
    }
    gzclose(f);
    return true;
}

}  // namespace

int run_wg(int argc, char** argv) {
    auto args = parse_wg_args(argc, argv);
    if (!args.ok) {
        if (args.err == "HELP") { print_wg_usage(std::cout); return 0; }
        std::cerr << "branch wg: " << args.err << "\n\n";
        print_wg_usage(std::cerr);
        return 2;
    }

    using clock = std::chrono::steady_clock;
    branch::wg::WgOutputs out;
    out.meta.assembler_version = "branch-v0.5.0";
    out.meta.read_tech = args.read_tech;
    out.meta.sample = args.out_prefix.substr(args.out_prefix.find_last_of('/') + 1);

    const auto& profile = (args.read_tech == "ont")
        ? branch::graph::kTechProfileONT
        : branch::graph::kTechProfileHiFi;

    int expected_cov = args.expected_coverage > 0 ? args.expected_coverage
                                                  : profile.expected_coverage;

    // Stream reads into RAM (first cut). Future commit: streaming + per-phase shards.
    std::cerr << "[wg] loading reads from " << args.reads << "\n";
    std::vector<std::pair<std::string, std::string>> reads;
    if (!stream_fastq_gz(args.reads, reads)) {
        std::cerr << "[wg] cannot read FASTQ: " << args.reads << "\n";
        return 3;
    }
    out.meta.input_reads = reads.size();
    std::uint64_t total_bp = 0;
    for (const auto& r : reads) total_bp += r.second.size();
    out.meta.input_bases = total_bp;
    std::cerr << "[wg] loaded " << reads.size() << " reads, " << total_bp << " bp\n";

    // ---------- Phase 0 ----------
    auto t0 = clock::now();
    std::cerr << "[wg] Phase 0: k-mer histogram + bimodality detection\n";
    // Size the bloom filter for expected distinct k-mers. For chr14
    // (~100 Mbp) at k=15 ONT preset there are ~50M distinct k-mers; for
    // a full WG ~20G. Cap at 1G to keep bloom RSS under control on a
    // single node — over-saturation just costs FP rate, not correctness.
    std::size_t expected_distinct = std::min<std::size_t>(
        std::max<std::size_t>(total_bp / 4, 100'000'000ULL),
        1'000'000'000ULL);
    branch::wg::KmerHistBuilder kh(profile, expected_distinct);
    for (const auto& [_, seq] : reads) kh.pass1_add_read(seq);
    kh.switch_to_pass2();
    for (const auto& [_, seq] : reads) kh.pass2_add_read(seq);
    auto khr = kh.finalize(expected_cov, profile.bimodality_z);
    out.meta.detected_coverage = khr.detected_coverage;
    out.meta.detected_het_coverage = khr.detected_het_coverage;
    out.meta.bimodality_z = khr.bimodality_z;
    std::string p0_path = args.out_prefix + ".phase0.bin";
    kh.write_to(p0_path, khr);
    out.phase_log.emplace_back("0",
        std::chrono::duration<double>(clock::now() - t0).count());
    std::cerr << "[wg] Phase 0 done. recurrent_kmers=" << khr.n_recurrent_kmers
              << " detected_cov=" << khr.detected_coverage
              << " detected_het=" << khr.detected_het_coverage
              << " bimodality_z=" << khr.bimodality_z << "\n";

    // ---------- Phase 1.1: routing ----------
    t0 = clock::now();
    std::cerr << "[wg] Phase 1.1: haplotype routing\n";
    int route_cov = khr.detected_coverage > 0 ? khr.detected_coverage : expected_cov;
    branch::wg::HaplotypeRouter router(profile, p0_path, route_cov);
    router.find_het_pairs();
    router.set_anchor_from_reads(reads);
    std::cerr << "[wg]   het_pairs=" << router.n_het_pairs()
              << " anchor=" << (router.has_anchor() ? "ok" : "none") << "\n";

    // Bin reads to hap1/hap2/ambig and stamp #disp records.
    std::vector<std::pair<std::string, std::string>> hap1_reads, hap2_reads, ambig_reads;
    hap1_reads.reserve(reads.size() / 2);
    hap2_reads.reserve(reads.size() / 2);
    for (const auto& [id, seq] : reads) {
        auto rr = router.route_read(seq);
        const char* label = "ambig";
        bool routed = false;
        if (rr.confidence() >= 0.6) {
            switch (rr.hap) {
                case branch::wg::Haplotype::Hap1:
                    label = "hap1"; hap1_reads.emplace_back(id, seq); routed = true; break;
                case branch::wg::Haplotype::Hap2:
                    label = "hap2"; hap2_reads.emplace_back(id, seq); routed = true; break;
                default: break;
            }
        }
        if (!routed) ambig_reads.emplace_back(id, seq);
        std::ostringstream info;
        info << "phase=1 outcome=" << label
             << " votes_h1=" << rr.votes_hap1
             << " votes_h2=" << rr.votes_hap2;
        out.read_disposition.emplace_back(id, info.str());
    }
    std::cerr << "[wg]   bin: hap1=" << hap1_reads.size()
              << " hap2=" << hap2_reads.size()
              << " ambig=" << ambig_reads.size() << "\n";

    // ---------- Phase 1.3: master tiling ----------
    std::cerr << "[wg] Phase 1.3: master tiling per haplotype\n";
    branch::wg::MasterTiler tiler(profile);
    out.haplotype_tiles.assign(2, {});
    out.haplotype_seqs.assign(2, std::string());
    tiler.build_tiling(hap1_reads, out.haplotype_tiles[0]);
    tiler.build_tiling(hap2_reads, out.haplotype_tiles[1]);
    out.haplotype_seqs[0] = tiler.render_haplotype_seq(out.haplotype_tiles[0], hap1_reads);
    out.haplotype_seqs[1] = tiler.render_haplotype_seq(out.haplotype_tiles[1], hap2_reads);
    std::cerr << "[wg]   hap1: tiles=" << out.haplotype_tiles[0].size()
              << " bp=" << out.haplotype_seqs[0].size() << "\n";
    std::cerr << "[wg]   hap2: tiles=" << out.haplotype_tiles[1].size()
              << " bp=" << out.haplotype_seqs[1].size() << "\n";

    // ---------- Phase 1.5: curation ----------
    std::cerr << "[wg] Phase 1.5: master curation\n";
    branch::wg::MasterCurator curator(profile);
    std::vector<branch::wg::CurationEvent> ev1, ev2;
    std::vector<branch::wg::BranchCandidate> bc1, bc2;
    curator.curate(0, out.haplotype_tiles[0], out.haplotype_seqs[0],
                   hap1_reads, ev1, bc1);
    curator.curate(1, out.haplotype_tiles[1], out.haplotype_seqs[1],
                   hap2_reads, ev2, bc2);
    out.curation_events = std::move(ev1);
    out.curation_events.insert(out.curation_events.end(),
                               ev2.begin(), ev2.end());
    std::vector<branch::wg::BranchCandidate> all_candidates;
    all_candidates.reserve(bc1.size() + bc2.size());
    all_candidates.insert(all_candidates.end(), bc1.begin(), bc1.end());
    all_candidates.insert(all_candidates.end(), bc2.begin(), bc2.end());
    std::cerr << "[wg]   curation_events=" << out.curation_events.size()
              << " branch_candidates=" << all_candidates.size() << "\n";
    out.phase_log.emplace_back("1",
        std::chrono::duration<double>(clock::now() - t0).count());

    // ---------- Phase 2: branch attachment ----------
    t0 = clock::now();
    std::cerr << "[wg] Phase 2: branch attachment\n";
    branch::wg::BranchAttacher att(profile);
    std::vector<branch::wg::OrphanRead> orphans;
    att.attach(out.haplotype_seqs, all_candidates, ambig_reads,
               out.branches, orphans);
    std::cerr << "[wg]   attached=" << out.branches.size()
              << " orphans=" << orphans.size() << "\n";
    out.phase_log.emplace_back("2",
        std::chrono::duration<double>(clock::now() - t0).count());

    // ---------- Phase 3: vPCR ----------
    if (args.do_vpcr) {
        t0 = clock::now();
        std::cerr << "[wg] Phase 3: vPCR auto-design + counting\n";
        branch::wg::VpcrAutoDesigner vp(profile);
        vp.design(out.haplotype_seqs, out.vpcr_amplicons);
        std::cerr << "[wg]   designed " << out.vpcr_amplicons.size()
                  << " amplicons\n";
        vp.count(out.vpcr_amplicons, reads, expected_cov);
        out.phase_log.emplace_back("3",
            std::chrono::duration<double>(clock::now() - t0).count());
    }

    // ---------- Phase 4: ref annotation ----------
    if (args.do_ref) {
        std::cerr << "[wg] Phase 4: ref annotation deferred to `branch project`\n";
        // The existing `branch project` subcommand can be invoked
        // separately on the FASTA output. Wiring it inline is a v0.5.1
        // follow-up — for now the emitted FASTA is project-ready.
    }

    // ---------- Phase 5: orphan clustering ----------
    t0 = clock::now();
    std::cerr << "[wg] Phase 5: orphan re-clustering\n";
    branch::wg::OrphanClusterer oc(profile);
    oc.cluster(orphans, out.orphan_clusters);
    std::cerr << "[wg]   orphan_clusters=" << out.orphan_clusters.size() << "\n";
    out.phase_log.emplace_back("5",
        std::chrono::duration<double>(clock::now() - t0).count());

    // ---------- Write outputs ----------
    branch::wg::VpfWriter writer(args.out_prefix);
    writer.write(out);
    std::cerr << "[wg] wrote " << args.out_prefix << ".{fa,bed,gfa}\n";
    return 0;
}

int run_wg_extract(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "branch extract — pull a section out of a v0.5 GFA file.\n"
                  << "Usage:\n"
                  << "  branch extract <prefix>.gfa <view>\n"
                  << "Views: hap1, hap2, branches, orphans, vpcr, disposition, curation, log\n";
        return 2;
    }
    std::string path(argv[1]);
    std::string view(argv[2]);
    std::ifstream f(path);
    if (!f) { std::cerr << "cannot open " << path << "\n"; return 3; }
    std::string line;
    auto starts_with = [](std::string_view s, std::string_view p) {
        return s.size() >= p.size() && std::memcmp(s.data(), p.data(), p.size()) == 0;
    };
    while (std::getline(f, line)) {
        std::string_view sv(line);
        if      (view == "hap1" && starts_with(sv, "S\thap1\t"))         { std::cout << line << "\n"; }
        else if (view == "hap2" && starts_with(sv, "S\thap2\t"))         { std::cout << line << "\n"; }
        else if (view == "branches" && starts_with(sv, "S\tbranch"))     { std::cout << line << "\n"; }
        else if (view == "orphans" && starts_with(sv, "S\torphan"))      { std::cout << line << "\n"; }
        else if (view == "vpcr" && starts_with(sv, "#vpcr"))             { std::cout << line << "\n"; }
        else if (view == "disposition" && starts_with(sv, "#disp"))      { std::cout << line << "\n"; }
        else if (view == "curation" && starts_with(sv, "#cur"))          { std::cout << line << "\n"; }
        else if (view == "log" && starts_with(sv, "#log"))               { std::cout << line << "\n"; }
    }
    return 0;
}

}  // namespace branch::cli
