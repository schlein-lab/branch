#include "phaser.hpp"

#include "anchored/anchored_phaser.hpp"
#include "balance_qc.hpp"
#include "bubble_classifier.hpp"
#include "bubble_finder.hpp"
#include "bubble_voter.hpp"
#include "fastq_index.hpp"
#include "gfa_writer.hpp"
#include "overlap_graph.hpp"
#include "phase_engine.hpp"
#include "reconcile/reconciler.hpp"
#include "refiner.hpp"

#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

namespace branch::wg::phaser {

namespace {

void seed_assignments(
    std::size_t n_reads,
    const std::vector<UnitigNode>& unitigs,
    std::vector<ReadAssignment>& out)
{
    out.clear();
    out.reserve(n_reads);

    // ReadId == 0..N-1 from FastqIndex; build a flat lookup table.
    std::vector<NodeId> home(n_reads, 0);
    for (NodeId ni = 0; ni < unitigs.size(); ++ni) {
        for (ReadId rid : unitigs[ni].member_reads) {
            if (rid < n_reads) home[rid] = ni;
        }
    }
    for (std::size_t rid = 0; rid < n_reads; ++rid) {
        ReadAssignment a;
        a.read_id = static_cast<ReadId>(rid);
        a.tag = ReadTag::UNTAGGED;
        a.confidence = 0.0f;
        a.home_node = home[rid];
        out.push_back(a);
    }
}

}  // namespace

PhaserResult run(const ReadsInput& reads_in, const PhaserOpts& opts) {
    using clock = std::chrono::steady_clock;
    auto t_start = clock::now();

    PhaserResult result;

    // ---- 0. Open the read source. Either a FASTQ on disk or an
    //         in-memory vector for tests; subsequent phases are agnostic.
    FastqIndex idx;
    if (!reads_in.fastq_path.empty()) {
        idx.build(reads_in.fastq_path);
    } else {
        idx.build_from_memory(reads_in.reads, reads_in.read_id_strings);
    }
    const std::size_t N = idx.n_reads();
    result.stats.n_reads_total = N;

    // ---- 1. Overlap graph ----
    std::cerr << "[phaser] step 1: overlap graph (" << N << " reads)\n";
    std::fflush(stderr);
    {
        OverlapGraph og;
        OverlapGraphOpts og_opts;
        og_opts.minimizer_k    = opts.minimizer_k;
        og_opts.minimizer_w    = opts.minimizer_w;
        og_opts.min_overlap_bp = opts.min_overlap_bp;
        og_opts.min_jaccard    = opts.min_jaccard;
        // Sketching needs the actual sequences. For WG-HiFi this is
        // ~60 GB. Materialize-once-then-drop keeps the floor low for the
        // rest of the pipeline (bubble_finder, balance_qc, etc.).
        {
            std::vector<std::pair<ReadId, std::string>> reads_vec;
            reads_vec.reserve(N);
            idx.for_each([&](ReadId rid, const std::string& seq) {
                reads_vec.emplace_back(rid, seq);
            });
            og.build(reads_vec, og_opts);
            // reads_vec destroyed here → 60 GB freed
        }
        og.transitive_reduce();
        og.collapse_to_unitigs(idx, result.unitigs);
    }
    // og destroyed → edges_ and adj freed before bubble_finder.

    // ---- 2. Bubbles ----
    std::cerr << "[phaser] step 2: bubble finder\n";
    std::fflush(stderr);
    BubbleFinder bf;
    bf.find(result.unitigs, result.bubbles);

    // ---- 3. Phase engine, iterative ----
    std::cerr << "[phaser] step 3: phase engine\n";
    std::fflush(stderr);
    seed_assignments(N, result.unitigs, result.assignments);

    // ---- 2b. Vote bubbles via sequence-similarity (v0.10 Variante B).
    //
    // bubble_voter materializes each alt's sequence, extracts canonical
    // k-mers, builds an inverted kmer→(bubble,alt) index, then scans
    // each read assigning it to its best-matching alt. The output is
    // per-bubble per-alt read counts that drive the VAF classification
    // — replaces the legacy home_node-based counting which starved
    // 99.7% of read evidence at WG scale (v16 measured 12.7K voting
    // reads / 4.3M total; voter should pull in 1-2M).
    //
    // The "clonally branched haplotype lineage" tree is populated here
    // at depth 1 (root + one child per alt). B-2 (per-position VAF
    // cascade for deeper sub-branch resolution within each alt) is a
    // follow-up pass.
    std::vector<BubbleVote> votes;
    std::vector<VoterReadOutcome> read_outcomes;
    {
        BubbleVoter bv;
        bv.vote(result.unitigs, result.bubbles, idx, votes, read_outcomes);
        BubbleVoter::populate_branch_trees(result.bubbles, votes);
    }

    // ---- 2c. Classify bubbles by VAF pattern (depth-0..3).
    //
    // Now consumes the voter's read_counts/vafs directly. The classifier
    // pattern-matches against canonical [50, 25, 25] / [50, 25, 12.5,
    // 12.5] / etc. families and stamps Bubble.type + branch_depth.
    {
        BubbleClassifier bc;
        BubbleClassifierOpts copts;
        copts.min_total_reads = 4;  // sane floor on noisy-pattern reject
        bc.classify_from_vafs(result.bubbles, copts);
    }

    PhaseEngine pe;
    pe.prepare(idx, reads_in.het_pairs);
    if (reads_in.het_pairs.empty()) {
        std::cerr << "[phaser] WARNING: het_pairs empty — voting will not "
                     "find diploid signal; output will be SHARED-only.\n";
    } else {
        std::cerr << "[phaser] het_pairs=" << reads_in.het_pairs.size()
                  << " loaded into phase_engine signature index\n";
    }

    int build_iter = 0;
    int min_supp = opts.iter1_min_supporting;
    double min_ratio = opts.iter1_min_ratio;
    bool finished = false;
    while (build_iter < opts.max_build_iters && !finished) {
        PhaseEngineOpts pe_opts;
        pe_opts.iter_n = build_iter + 1;
        pe_opts.min_supporting_overlaps = min_supp;
        pe_opts.min_ratio = min_ratio;
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

        const std::size_t floor_count = static_cast<std::size_t>(
            opts.build_progress_floor *
            static_cast<double>(result.assignments.size()));
        if (progress.reads_newly_tagged < floor_count) {
            finished = true;
            break;
        }

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

    // ---- 3b. Apply voter outcomes (v0.10 Variante B):
    //
    // Bridge the per-read voter classification into ReadAssignment.tag,
    // overriding home_node-based tagging. Reads the voter assigned to a
    // specific (bubble, alt) get H1/H2/BRANCH; reads spanning ≥2 alts
    // become BRANCH_BRIDGE; reads with novel sequence become
    // BRANCH_NOVEL. Reads the voter saw nothing of (out-of-bubble) keep
    // whatever home_node-based tag they got, which is SHARED for the
    // homozygous majority — correct.
    //
    // This is what closes the v16 gap: voter pulled in 1-2M reads but
    // their tags stayed SHARED because home_node was off-alt. Bridge
    // makes the voter-output authoritative.
    PhaseEngine::apply_voter_outcomes(result.unitigs, result.bubbles,
                                      read_outcomes, result.assignments);

    // ---- 4. Refinement (skipped at WG scale by Refiner itself) ----
    std::cerr << "[phaser] step 4: refinement\n";
    std::fflush(stderr);
    Refiner rf;
    RefinerOpts rf_opts;
    rf_opts.change_floor = opts.refine_change_floor;
    rf_opts.max_iters    = opts.max_refine_iters;
    auto rp = rf.refine(result.unitigs, result.bubbles,
                        result.assignments, idx, rf_opts);
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
        switch (b.type) {
            case BubbleType::HET:              ++result.stats.n_bubbles_het; break;
            case BubbleType::HOMOZ:            ++result.stats.n_bubbles_homoz; break;
            case BubbleType::BRANCH:           ++result.stats.n_bubbles_branch; break;
            case BubbleType::BRANCH_IN_BRANCH: ++result.stats.n_bubbles_branch_in_branch; break;
            case BubbleType::TRIALLELIC:       ++result.stats.n_bubbles_triallelic; break;
            case BubbleType::NOISY:            ++result.stats.n_bubbles_noisy; break;
            case BubbleType::EXTINCT_INTERMEDIATE:      ++result.stats.n_bubbles_noisy; break;
            case BubbleType::UNKNOWN:          break;
        }
    }
    for (const auto& a : result.assignments) {
        switch (a.tag) {
            case ReadTag::H1_ONLY:       ++result.stats.n_reads_h1;        break;
            case ReadTag::H2_ONLY:       ++result.stats.n_reads_h2;        break;
            case ReadTag::SHARED:        ++result.stats.n_reads_shared;    break;
            case ReadTag::BRANCH:        ++result.stats.n_reads_branch;    break;
            case ReadTag::UNCERTAIN:     ++result.stats.n_reads_uncertain; break;
            case ReadTag::BRANCH_BRIDGE: ++result.stats.n_reads_bridge;    break;
            case ReadTag::BRANCH_NOVEL:  ++result.stats.n_reads_novel;     break;
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
    gw_opts.out_branch_fa_path       = opts.out_branch_fa_path;
    gw_opts.out_bridge_fa_path       = opts.out_bridge_fa_path;
    gw_opts.out_novel_fa_path        = opts.out_novel_fa_path;
    gw_opts.out_assignments_tsv_path = opts.out_assignments_tsv_path;
    gw_opts.approx_overlap_bp        = opts.minimizer_w * 5;
    gw.write(result.unitigs, result.bubbles, result.assignments,
             idx.read_id_strings(), idx, gw_opts);

    // ---- 7. v0.8 dual-phaser: optional anchored track + reconciliation ----
    if (!opts.ref_hap1_fa_path.empty() && !opts.ref_hap2_fa_path.empty()) {
        std::cerr << "[phaser] step 7: anchored phaser (v0.8)\n";
        std::fflush(stderr);
        AnchoredPhaser ap;
        AnchoredPhaserOpts ap_opts;
        ap_opts.ref_hap1_fa_path  = opts.ref_hap1_fa_path;
        ap_opts.ref_hap2_fa_path  = opts.ref_hap2_fa_path;
        ap_opts.n_threads                = opts.anchored_n_threads;
        ap_opts.min_aligned_bp           = opts.anchored_min_aligned_bp;
        ap_opts.min_discriminating_sites = opts.anchored_min_discriminating_sites;
        ap_opts.preset                   = opts.anchored_preset;
        std::vector<ReadAssignment> assignments_b;
        ap.run(idx, assignments_b, ap_opts);

        std::cerr << "[phaser] step 8: reconciliation\n";
        std::fflush(stderr);
        Reconciler rc;
        ReconcilerOpts rc_opts;
        rc_opts.neural_model_path = opts.neural_voter_model_path;
        rc_opts.verbose = true;
        std::vector<ReconciledAssignment> reconciled;
        rc.reconcile(result.assignments, assignments_b, reconciled, rc_opts);

        if (!opts.out_reconciled_tsv_path.empty()) {
            std::ofstream r(opts.out_reconciled_tsv_path);
            if (!r) {
                std::cerr << "[phaser] cannot open "
                          << opts.out_reconciled_tsv_path << "\n";
            } else {
                r << "read_id\ttag_a\ttag_b\ttag_final\tsource\tconfidence"
                     "\tconf_a\tconf_b\thome_node_a\thome_node_b\n";
                const auto& names = idx.read_id_strings();
                for (const auto& x : reconciled) {
                    const std::string& nm =
                        (x.read_id < names.size())
                            ? names[x.read_id]
                            : std::to_string(x.read_id);
                    r << nm << '\t'
                      << read_tag_str(x.tag_a)   << '\t'
                      << read_tag_str(x.tag_b)   << '\t'
                      << read_tag_str(x.tag_final) << '\t'
                      << phase_source_str(x.source) << '\t'
                      << x.confidence << '\t'
                      << x.conf_a << '\t'
                      << x.conf_b << '\t'
                      << x.home_node_a << '\t'
                      << x.home_node_b << '\n';
                }
                std::cerr << "[phaser] wrote " << opts.out_reconciled_tsv_path
                          << " (" << reconciled.size() << " reads, "
                          << rc.stats().n_agreed << " agreed, "
                          << rc.stats().n_flagged << " flagged)\n";
            }
        }
    }

    auto t_end = clock::now();
    const double elapsed_s =
        std::chrono::duration<double>(t_end - t_start).count();

    std::cerr << "[phaser] done in " << elapsed_s << "s — "
              << "h1=" << result.stats.n_reads_h1
              << " h2=" << result.stats.n_reads_h2
              << " shared=" << result.stats.n_reads_shared
              << " branch=" << result.stats.n_reads_branch
              << " uncertain=" << result.stats.n_reads_uncertain
              << " bridge=" << result.stats.n_reads_bridge
              << " novel=" << result.stats.n_reads_novel
              << " | h1_bp=" << static_cast<long long>(result.stats.h1_bp)
              << " h2_bp=" << static_cast<long long>(result.stats.h2_bp)
              << " shared_bp=" << static_cast<long long>(result.stats.shared_bp)
              << " | qc_flags=" << result.qc_flags.size()
              << " | bubbles: het=" << result.stats.n_bubbles_het
              << " branch=" << result.stats.n_bubbles_branch
              << " bib=" << result.stats.n_bubbles_branch_in_branch
              << " tri=" << result.stats.n_bubbles_triallelic
              << " noisy=" << result.stats.n_bubbles_noisy << "\n";

    return result;
}

}  // namespace branch::wg::phaser
