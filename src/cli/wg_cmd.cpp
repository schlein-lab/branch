// BRANCH v0.5 — `branch wg` orchestrator subcommand.
//
// Runs phases 0 → 5 in order and writes <out-prefix>.{fa,bed,gfa}.

#include "wg_cmd.hpp"

#include "wg/kmer_hist.hpp"
#include "wg/haplotype_router.hpp"
#include "wg/master_tiler.hpp"
#include "wg/master_curator.hpp"
#include "wg/branch_attacher.hpp"
#include "wg/stage1_screen.hpp"
#include "wg/vpcr_wg.hpp"
#include "wg/orphan_clusterer.hpp"
#include "wg/vpf_writer.hpp"
#include "graph/kmer_sketch.hpp"
#include "gpu/wg_kernels.cuh"

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
    bool no_gpu = false;          // --no-gpu forces CPU-only Phase 3
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
        else if (k == "--no-gpu")     { a.no_gpu = true; }
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

    // ---------- Phase 0: k-mer histogram (GPU-first, CPU fallback) ----------
    auto t0 = clock::now();
    std::cerr << "[wg] Phase 0: k-mer histogram + bimodality detection\n";
    bool use_gpu_p0 = !args.no_gpu && branch::gpu::wg_kernels::is_gpu_available();
    branch::wg::KmerHistResult khr;
    std::string p0_path = args.out_prefix + ".phase0.bin";
    bool p0_done = false;
    if (use_gpu_p0) {
        std::vector<std::uint64_t> p0_keys;
        std::vector<std::uint32_t> p0_cnts;
        bool ok = branch::gpu::wg_kernels::launch_phase0_count(
            reads, profile.minimizer_k, p0_keys, p0_cnts);
        if (ok) {
            khr = branch::wg::finalize_from_keys_counts(
                expected_cov, profile.bimodality_z,
                p0_keys, p0_cnts,
                reads.size(), total_bp);
            branch::wg::write_phase0_dump(p0_path, profile.minimizer_k,
                                          khr, p0_keys, p0_cnts);
            std::cerr << "[wg]   phase0 on GPU: recurrent_kmers=" << p0_keys.size()
                      << "\n";
            p0_done = true;
        }
    }
    if (!p0_done) {
        std::size_t expected_distinct = std::min<std::size_t>(
            std::max<std::size_t>(total_bp / 4, 100'000'000ULL),
            1'000'000'000ULL);
        branch::wg::KmerHistBuilder kh(profile, expected_distinct);
        for (const auto& [_, seq] : reads) kh.pass1_add_read(seq);
        kh.switch_to_pass2();
        for (const auto& [_, seq] : reads) kh.pass2_add_read(seq);
        khr = kh.finalize(expected_cov, profile.bimodality_z);
        kh.write_to(p0_path, khr);
        std::cerr << "[wg]   phase0 on CPU: recurrent_kmers="
                  << khr.n_recurrent_kmers << "\n";
    }
    out.meta.detected_coverage = khr.detected_coverage;
    out.meta.detected_het_coverage = khr.detected_het_coverage;
    out.meta.bimodality_z = khr.bimodality_z;
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
    bool use_gpu_p11 = !args.no_gpu && branch::gpu::wg_kernels::is_gpu_available();
    bool p11_done = false;
    if (use_gpu_p11) {
        // Build sorted (key, count) arrays from the Phase 0 dump for the
        // GPU het-pair launcher; binary-searches partner kmers in-kernel.
        std::unordered_map<std::uint64_t, std::uint32_t> counter_map;
        branch::wg::load_kmer_counter_dump(p0_path, counter_map);
        std::vector<std::pair<std::uint64_t, std::uint32_t>> pairs(
            counter_map.begin(), counter_map.end());
        std::sort(pairs.begin(), pairs.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
        std::vector<std::uint64_t> p11_keys;
        std::vector<std::uint32_t> p11_cnts;
        p11_keys.reserve(pairs.size());
        p11_cnts.reserve(pairs.size());
        for (const auto& kv : pairs) {
            p11_keys.push_back(kv.first);
            p11_cnts.push_back(kv.second);
        }
        std::vector<std::uint64_t> ha, hb;
        bool ok = branch::gpu::wg_kernels::launch_phase11_het_pairs(
            p11_keys, p11_cnts, profile.minimizer_k, route_cov,
            0.30, 0.70, 0.70, 1.40, ha, hb);
        if (ok) {
            router.install_het_pairs_from_pairs(ha, hb);
            std::cerr << "[wg]   phase1.1 on GPU: het_pairs=" << ha.size() << "\n";
            p11_done = true;
        }
    }
    if (!p11_done) {
        router.find_het_pairs();
        std::cerr << "[wg]   phase1.1 on CPU: het_pairs=" << router.n_het_pairs() << "\n";
    }
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

    // ---------- Phase 1.5: curation (GPU-first, CPU fallback) ----------
    std::cerr << "[wg] Phase 1.5: master curation\n";
    bool use_gpu_p15 = !args.no_gpu && branch::gpu::wg_kernels::is_gpu_available();
    branch::wg::MasterCurator curator(profile);
    std::vector<branch::wg::CurationEvent> ev1, ev2;
    std::vector<branch::wg::BranchCandidate> bc1, bc2;
    int default_q = (profile.name && profile.name[0] == 'h') ? 30 : 18;
    auto curate_one = [&](std::uint32_t h_idx,
                          const std::string& seq,
                          const std::vector<std::pair<std::string, std::string>>& reads_h,
                          std::vector<branch::wg::CurationEvent>& evout,
                          std::vector<branch::wg::BranchCandidate>& bcout) {
        bool ok = false;
        if (use_gpu_p15) {
            ok = branch::gpu::wg_kernels::launch_phase15_curation(
                h_idx, seq, reads_h,
                profile.minimizer_k, profile.minimizer_w,
                profile.min_coverage_branch,
                profile.min_base_agreement,
                default_q, evout, bcout);
            if (ok) std::cerr << "[wg]   hap" << (h_idx+1)
                              << " curate on GPU: events=" << evout.size()
                              << " branches=" << bcout.size() << "\n";
        }
        if (!ok) {
            curator.curate(h_idx, out.haplotype_tiles[h_idx], seq, reads_h, evout, bcout);
            std::cerr << "[wg]   hap" << (h_idx+1)
                      << " curate on CPU: events=" << evout.size()
                      << " branches=" << bcout.size() << "\n";
        }
    };
    curate_one(0, out.haplotype_seqs[0], hap1_reads, ev1, bc1);
    curate_one(1, out.haplotype_seqs[1], hap2_reads, ev2, bc2);
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

    // ---------- Phase 2: branch attachment (GPU-first, CPU fallback) ----------
    t0 = clock::now();
    std::cerr << "[wg] Phase 2: branch attachment\n";
    bool use_gpu_p2 = !args.no_gpu && branch::gpu::wg_kernels::is_gpu_available();
    bool p2_done = false;
    branch::wg::BranchAttacher att(profile);
    std::vector<branch::wg::OrphanRead> orphans;

    if (use_gpu_p2) {
        // GPU direct-attach: per ambig read, sketch + best-shared-with-hap.
        // Ambig reads with shared >= jaccard_lo*kSketchSize get attached to
        // the best hap at depth=1; rest go to the orphan list. Curator's
        // direct branches are kept verbatim.
        std::vector<std::int32_t> anchored_hap;
        bool ok = branch::gpu::wg_kernels::launch_phase2_direct_attach(
            out.haplotype_seqs, ambig_reads, profile.minimizer_k,
            profile.sketch_jaccard_lo, anchored_hap);
        if (ok) {
            // Build out.branches: direct curator candidates + anchored ambigs.
            std::uint32_t bidx = 0;
            for (const auto& c : all_candidates) {
                branch::wg::AttachedBranch b{};
                std::ostringstream name; name << "branch_" << ++bidx;
                b.branch_name = name.str();
                b.seq = c.alt_seq;
                b.anchor_hap = c.hap_idx;
                b.anchor_pos = c.pos_in_hap;
                b.transitive_depth = 0;
                b.confidence = c.q_weighted_evidence;
                out.branches.push_back(std::move(b));
            }
            for (std::size_t r = 0; r < ambig_reads.size(); ++r) {
                if (anchored_hap[r] >= 0) {
                    branch::wg::AttachedBranch b{};
                    std::ostringstream name; name << "branch_" << ++bidx;
                    b.branch_name = name.str();
                    b.seq = ambig_reads[r].second;
                    b.anchor_hap = static_cast<std::uint32_t>(anchored_hap[r]);
                    b.anchor_pos = 0;
                    b.transitive_depth = 1;
                    b.confidence = 0.5;
                    out.branches.push_back(std::move(b));
                } else {
                    orphans.push_back(branch::wg::OrphanRead{
                        ambig_reads[r].first, ambig_reads[r].second});
                }
            }
            std::cerr << "[wg]   phase2 on GPU: attached=" << out.branches.size()
                      << " orphans=" << orphans.size() << "\n";
            p2_done = true;
        }
    }
    if (!p2_done) {
        att.attach(out.haplotype_seqs, all_candidates, ambig_reads,
                   out.branches, orphans);
        std::cerr << "[wg]   phase2 on CPU: attached=" << out.branches.size()
                  << " orphans=" << orphans.size() << "\n";
    }
    out.phase_log.emplace_back("2",
        std::chrono::duration<double>(clock::now() - t0).count());

    // ---------- Phase 3: two-stage vPCR (GPU-first, CPU fallback) ----------
    //   Stage 1 — cheap whole-genome k-mer-coverage screen (no coding
    //   bias, multi-window, OR'd signals). Output: candidate regions +
    //   per-window coverage track (emitted to <prefix>.cov.bed).
    //   Stage 2 — design 1-N primer pairs per candidate region, then
    //   count via k-mer-seed indexed scan.
    bool use_gpu = !args.no_gpu && branch::gpu::wg_kernels::is_gpu_available();
    std::cerr << "[wg] gpu=" << (use_gpu ? "available, will use" : "off") << "\n";
    if (args.do_vpcr) {
        t0 = clock::now();
        std::cerr << "[wg] Phase 3 Stage 1: cheap k-mer-coverage screen\n";
        // Reload Phase 0 counter (was freed by router after Phase 1.1).
        std::unordered_map<std::uint64_t, std::uint32_t> p0_counter;
        branch::wg::load_kmer_counter_dump(p0_path, p0_counter);
        std::cerr << "[wg]   loaded " << p0_counter.size()
                  << " recurrent k-mers from " << p0_path << "\n";

        // Pre-flatten counter for GPU launcher (sorted ascending).
        std::vector<std::uint64_t> p0_keys_sorted;
        std::vector<std::uint32_t> p0_cnts_sorted;
        if (use_gpu) {
            p0_keys_sorted.reserve(p0_counter.size());
            p0_cnts_sorted.reserve(p0_counter.size());
            std::vector<std::pair<std::uint64_t, std::uint32_t>> pairs(
                p0_counter.begin(), p0_counter.end());
            std::sort(pairs.begin(), pairs.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
            for (const auto& kv : pairs) {
                p0_keys_sorted.push_back(kv.first);
                p0_cnts_sorted.push_back(kv.second);
            }
            std::cerr << "[wg]   flattened counter for GPU ("
                      << p0_keys_sorted.size() << " entries)\n";
        }

        branch::wg::Stage1Screener screener(profile);
        std::vector<branch::wg::Stage1Window> all_windows;
        std::vector<branch::wg::CandidateRegion> all_regions;
        for (std::size_t h = 0; h < out.haplotype_seqs.size(); ++h) {
            std::vector<branch::wg::Stage1Window> wins;
            bool gpu_ok = false;
            if (use_gpu) {
                gpu_ok = branch::gpu::wg_kernels::launch_stage1(
                    static_cast<std::uint32_t>(h),
                    out.haplotype_seqs[h],
                    p0_keys_sorted, p0_cnts_sorted,
                    profile.minimizer_k, profile.minimizer_w,
                    branch::wg::kStage1ZThreshold,
                    branch::wg::kStage1RepEntropy,
                    branch::wg::kStage1MinDensityZ,
                    wins);
                if (gpu_ok) {
                    std::cerr << "[wg]   hap" << (h+1) << " stage1 on GPU: "
                              << wins.size() << " windows\n";
                }
            }
            if (!gpu_ok) {
                screener.screen_haplotype(static_cast<std::uint32_t>(h),
                                          out.haplotype_seqs[h], p0_counter, wins);
                std::cerr << "[wg]   hap" << (h+1) << " stage1 on CPU: "
                          << wins.size() << " windows\n";
            }
            all_windows.insert(all_windows.end(), wins.begin(), wins.end());
        }
        screener.merge_into_regions(all_windows, all_regions);
        std::cerr << "[wg]   stage1: " << all_windows.size()
                  << " windows, " << all_regions.size()
                  << " candidate regions flagged\n";

        // Emit Stage-1 coverage track as a separate BED next to the main one.
        std::string cov_bed_path = args.out_prefix + ".cov.bed";
        std::ofstream cov(cov_bed_path);
        if (cov) {
            for (const auto& w : all_windows) {
                cov << "hap" << (w.hap_idx + 1) << '\t'
                    << w.start_bp << '\t' << (w.start_bp + w.length_bp) << '\t'
                    << "win"
                    << '\t' << w.mean_kmer_count
                    << '\t' << w.distinct_kmer_count
                    << '\t' << w.dimer_entropy
                    << '\t' << w.minimizer_density
                    << '\t' << w.z_mean
                    << '\t' << static_cast<int>(w.flags) << '\n';
            }
            std::cerr << "[wg]   wrote " << cov_bed_path << "\n";
        }

        std::cerr << "[wg] Phase 3 Stage 2: design + count amplicons in "
                  << all_regions.size() << " regions\n";
        branch::wg::VpcrAutoDesigner vp(profile);
        vp.design(out.haplotype_seqs, all_regions, out.vpcr_amplicons,
                  /*n_per_region=*/3);
        std::cerr << "[wg]   designed " << out.vpcr_amplicons.size()
                  << " amplicons\n";

        bool stage2_gpu_ok = false;
        const int amp_min = 200;
        const int amp_max = (profile.name && profile.name[0] == 'h') ? 1500 : 3000;
        if (use_gpu) {
            stage2_gpu_ok = branch::gpu::wg_kernels::launch_stage2(
                out.vpcr_amplicons, reads,
                /*max_mm=*/2, amp_min, amp_max, expected_cov);
            if (stage2_gpu_ok) std::cerr << "[wg]   stage2 counted on GPU\n";
        }
        if (!stage2_gpu_ok) {
            vp.count(out.vpcr_amplicons, reads, expected_cov);
            std::cerr << "[wg]   stage2 counted on CPU\n";
        }
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
