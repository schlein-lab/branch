#include "phaser.hpp"

#include "balance_qc.hpp"
#include "bubble_finder.hpp"
#include "gfa_writer.hpp"
#include "overlap_graph.hpp"
#include "phase_engine.hpp"
#include "refiner.hpp"

#include <chrono>
#include <cstdio>
#include <iostream>
#include <unordered_map>

namespace branch::wg::phaser {

namespace {

void seed_assignments(
    const std::vector<std::pair<ReadId, std::string>>& reads,
    const std::vector<UnitigNode>& unitigs,
    std::vector<ReadAssignment>& out)
{
    out.clear();
    out.reserve(reads.size());
    // Map ReadId → home unitig
    std::unordered_map<ReadId, NodeId> home;
    home.reserve(reads.size() * 2);
    for (NodeId ni = 0; ni < unitigs.size(); ++ni) {
        for (ReadId rid : unitigs[ni].member_reads) home[rid] = ni;
    }
    for (const auto& [rid, _] : reads) {
        ReadAssignment a;
        a.read_id = rid;
        a.tag = ReadTag::UNTAGGED;
        a.confidence = 0.0f;
        auto it = home.find(rid);
        a.home_node = (it == home.end()) ? 0 : it->second;
        out.push_back(a);
    }
}

}  // namespace

PhaserResult run(const ReadsInput& reads_in, const PhaserOpts& opts) {
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();

    PhaserResult result;
    result.stats.n_reads_total = reads_in.reads.size();

    // ---- 1. Overlap graph ----
    std::cerr << "[phaser] step 1: overlap graph (" << reads_in.reads.size()
              << " reads)\n";
    std::fflush(stderr);
    OverlapGraph og;
    OverlapGraphOpts og_opts;
    og_opts.minimizer_k    = opts.minimizer_k;
    og_opts.minimizer_w    = opts.minimizer_w;
    og_opts.min_overlap_bp = opts.min_overlap_bp;
    og_opts.min_jaccard    = opts.min_jaccard;
    og.build(reads_in.reads, og_opts);
    og.collapse_to_unitigs(reads_in.reads, result.unitigs);

    // ---- 2. Bubbles ----
    std::cerr << "[phaser] step 2: bubble finder\n";
    std::fflush(stderr);
    BubbleFinder bf;
    bf.find(result.unitigs, result.bubbles);

    // ---- 3. Phase engine, iterative ----
    std::cerr << "[phaser] step 3: phase engine\n";
    std::fflush(stderr);
    seed_assignments(reads_in.reads, result.unitigs, result.assignments);

    PhaseEngine pe;
    pe.prepare(reads_in.reads, reads_in.het_pairs);
    if (reads_in.het_pairs.empty()) {
        std::cerr << "[phaser] WARNING: het_pairs empty — voting will not "
                     "find diploid signal; output will be SHARED-only.\n";
    } else {
        std::cerr << "[phaser] het_pairs=" << reads_in.het_pairs.size()
                  << " loaded into phase_engine signature index\n";
    }

    // Iterative build: start with strict thresholds, relax each pass.
    // Stop when newly-tagged drops below build_progress_floor of total
    // OR when max_build_iters is exhausted. Final pass is always
    // force-assign-remainder so no read exits UNTAGGED.
    int build_iter = 0;
    int min_supp = opts.iter1_min_supporting;
    double min_ratio = opts.iter1_min_ratio;
    bool finished = false;
    while (build_iter < opts.max_build_iters && !finished) {
        PhaseEngineOpts pe_opts;
        pe_opts.iter_n = build_iter + 1;
        pe_opts.min_supporting_overlaps = min_supp;
        pe_opts.min_ratio = min_ratio;
        // We do NOT force-assign in any non-final iteration; only the
        // explicit final sweep below should set UNCERTAIN tags.
        pe_opts.force_assign_remainder = false;

        auto progress = pe.run_iteration(
            result.unitigs, result.bubbles, result.assignments, pe_opts);

        IterationRecord rec;
        rec.iter_n = build_iter + 1;
        rec.min_supporting = min_supp;
        rec.min_ratio = min_ratio;
        rec.reads_newly_tagged   = progress.reads_newly_tagged;
        rec.bubbles_newly_phased = progress.bubbles_newly_phased;
        rec.reads_still_untagged = progress.reads_still_untagged;
        rec.force_assigned = false;
        result.build_trace.push_back(rec);

        std::cerr << "[phaser] build iter " << (build_iter + 1)
                  << ": min_supp=" << min_supp << " min_ratio=" << min_ratio
                  << " newly_tagged=" << progress.reads_newly_tagged
                  << " bubbles_phased=" << progress.bubbles_newly_phased
                  << " untagged_left=" << progress.reads_still_untagged << "\n";

        ++build_iter;

        // Convergence: <build_progress_floor of reads newly tagged AND
        // no untagged reads left → done. If untagged still left, force
        // assign sweep next.
        const std::size_t floor_count = static_cast<std::size_t>(
            opts.build_progress_floor *
            static_cast<double>(result.assignments.size()));
        if (progress.reads_newly_tagged < floor_count) {
            finished = true;
            break;
        }

        // Relax thresholds for next iteration.
        min_supp = std::max(1, min_supp / 2);
        min_ratio = std::max(0.51, min_ratio - 0.10);
    }

    // Final force-assign sweep — every still-untagged read gets UNCERTAIN.
    {
        PhaseEngineOpts final_opts;
        final_opts.iter_n = build_iter + 1;
        final_opts.min_supporting_overlaps = 1;
        final_opts.min_ratio = 0.51;
        final_opts.force_assign_remainder = true;
        auto progress = pe.run_iteration(
            result.unitigs, result.bubbles, result.assignments, final_opts);
        IterationRecord rec;
        rec.iter_n = build_iter + 1;
        rec.min_supporting = 1;
        rec.min_ratio = 0.51;
        rec.reads_newly_tagged = progress.reads_newly_tagged;
        rec.bubbles_newly_phased = progress.bubbles_newly_phased;
        rec.reads_still_untagged = progress.reads_still_untagged;
        rec.force_assigned = true;
        result.build_trace.push_back(rec);
        std::cerr << "[phaser] force-assign sweep: tagged_remainder="
                  << progress.reads_newly_tagged
                  << " untagged_left=" << progress.reads_still_untagged
                  << " (should be 0)\n";
    }
    result.stats.build_iterations = build_iter + 1;

    // ---- 4. Refinement (EM remap) ----
    std::cerr << "[phaser] step 4: refinement\n";
    std::fflush(stderr);
    Refiner rf;
    RefinerOpts rf_opts;
    rf_opts.change_floor = opts.refine_change_floor;
    rf_opts.max_iters    = opts.max_refine_iters;
    auto rp = rf.refine(result.unitigs, result.bubbles,
                        result.assignments, reads_in.reads, rf_opts);
    result.stats.refine_iterations = rp.iter_n;
    result.refine_trace = rp.per_iter_changes;

    // ---- 5. Balance QC ----
    std::cerr << "[phaser] step 5: balance QC\n";
    BalanceQc qc;
    BalanceQcOpts qc_opts;
    qc_opts.imbalance_warn = opts.imbalance_warn;
    qc_opts.imbalance_err  = opts.imbalance_err;
    qc.check(result.unitigs, result.bubbles, result.assignments,
             result.qc_flags, qc_opts);
    BalanceQc::summarize(result.unitigs,
                         result.stats.h1_bp, result.stats.h2_bp,
                         result.stats.shared_bp, result.stats.branch_bp);

    // Aggregate stats.
    for (const auto& u : result.unitigs) {
        switch (u.kind) {
            case NodeKind::SHARED:    ++result.stats.n_unitigs_shared; break;
            case NodeKind::BUBBLE_H1: ++result.stats.n_unitigs_h1;     break;
            case NodeKind::BUBBLE_H2: ++result.stats.n_unitigs_h2;     break;
            case NodeKind::BRANCH:    ++result.stats.n_unitigs_branch; break;
            default: break;
        }
    }
    for (const auto& b : result.bubbles) {
        ++result.stats.n_bubbles;
        if (b.phased) ++result.stats.n_bubbles_phased;
    }
    for (const auto& a : result.assignments) {
        switch (a.tag) {
            case ReadTag::H1_ONLY:   ++result.stats.n_reads_h1;        break;
            case ReadTag::H2_ONLY:   ++result.stats.n_reads_h2;        break;
            case ReadTag::SHARED:    ++result.stats.n_reads_shared;    break;
            case ReadTag::BRANCH:    ++result.stats.n_reads_branch;    break;
            case ReadTag::UNCERTAIN: ++result.stats.n_reads_uncertain; break;
            default: break;
        }
    }

    // ---- 6. Output ----
    std::cerr << "[phaser] step 6: output\n";
    GfaWriter gw;
    GfaWriterOpts gw_opts;
    gw_opts.out_gfa_path             = opts.out_gfa_path;
    gw_opts.out_h1_fa_path           = opts.out_h1_fa_path;
    gw_opts.out_h2_fa_path           = opts.out_h2_fa_path;
    gw_opts.out_assignments_tsv_path = opts.out_assignments_tsv_path;
    gw.write(result.unitigs, result.bubbles, result.assignments,
             reads_in.read_id_strings, gw_opts);

    auto t_end = clock::now();
    const double elapsed_s =
        std::chrono::duration<double>(t_end - t_start).count();

    std::cerr << "[phaser] done in " << elapsed_s << "s — "
              << "h1=" << result.stats.n_reads_h1
              << " h2=" << result.stats.n_reads_h2
              << " shared=" << result.stats.n_reads_shared
              << " branch=" << result.stats.n_reads_branch
              << " uncertain=" << result.stats.n_reads_uncertain
              << " | h1_bp=" << static_cast<long long>(result.stats.h1_bp)
              << " h2_bp=" << static_cast<long long>(result.stats.h2_bp)
              << " shared_bp=" << static_cast<long long>(result.stats.shared_bp)
              << " | qc_flags=" << result.qc_flags.size() << "\n";

    return result;
}

}  // namespace branch::wg::phaser
