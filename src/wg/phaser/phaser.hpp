#pragma once
//
// BRANCH v0.6 phaser — public API.
//
// One entry point: `run(reads, opts) -> PhaserResult`.
// Internally orchestrates overlap_graph → bubble_finder → phase_engine →
// refiner → balance_qc → gfa_writer. Each sub-module is independent and
// unit-testable in isolation.

#include "types.hpp"
#include "phase_engine.hpp"   // for HetPair
#include <string>
#include <vector>
#include <utility>

namespace branch::wg::phaser {

struct PhaserOpts {
    // Resource limits. The build loop terminates when either condition holds.
    std::size_t ram_budget_bytes = 100ULL << 30;   // 100 GB default
    int         max_build_iters  = 16;             // safety bound only

    // Termination criteria.
    double build_progress_floor = 0.05;   // <5% newly tagged → stop building
    double refine_change_floor  = 0.001;  // <0.1% reads moved → converged
    int    max_refine_iters     = 10;     // safety bound only

    // Overlap graph parameters.
    int         minimizer_k = 21;
    int         minimizer_w = 11;
    int         min_overlap_bp = 1000;
    double      min_jaccard    = 0.35;

    // Phasing thresholds — the FIRST iteration's strict thresholds.
    int         iter1_min_supporting = 10;
    double      iter1_min_ratio      = 0.95;

    // Branch (SV) detection.
    std::size_t branch_min_cluster_size = 5;
    double      branch_min_jaccard      = 0.60;

    // Imbalance QC.
    double      imbalance_warn = 0.05;   // 5%
    double      imbalance_err  = 0.20;   // 20%

    // Output paths. Empty disables that output.
    std::string out_gfa_path;
    std::string out_h1_fa_path;
    std::string out_h2_fa_path;
    std::string out_branch_fa_path;
    std::string out_bridge_fa_path;
    std::string out_novel_fa_path;
    std::string out_assignments_tsv_path;
    std::string out_qc_bed_path;

    // ---- v0.8 dual-phaser: optional pangenome-anchored second track ----
    //
    // When ref_hap1_fa_path AND ref_hap2_fa_path are both set, BRANCH
    // runs a second phaser that aligns reads against the pangenome
    // reference assemblies. The two assignment vectors then go into the
    // Reconciler, which agrees / neural-resolves / flags per read.
    //
    // HPRC is not ground-truth — it's a fast first model. The point of
    // running both is to expose disagreement, which is itself biological
    // signal.
    std::string ref_hap1_fa_path;
    std::string ref_hap2_fa_path;
    std::string out_reconciled_tsv_path;  // emits reconciled assignments
                                          // with provenance column
    int         anchored_n_threads = 8;
    int         anchored_min_align_score = 200;
    double      anchored_hap_call_margin = 0.10;
    std::string anchored_preset = "map-hifi";  // or "map-ont" for ONT
    std::string neural_voter_model_path;  // v0.9; empty for now
};

// Iteration log: one entry per build pass. Refinement gets its own
// entry list inside PhaserStats.
struct IterationRecord {
    int    iter_n              = 0;
    int    min_supporting      = 0;
    double min_ratio           = 0.0;
    std::size_t reads_newly_tagged   = 0;
    std::size_t bubbles_newly_phased = 0;
    std::size_t reads_still_untagged = 0;
    bool   force_assigned      = false;
};

struct PhaserResult {
    std::vector<UnitigNode>     unitigs;
    std::vector<Bubble>         bubbles;
    std::vector<ReadAssignment> assignments;
    std::vector<ImbalanceFlag>  qc_flags;
    std::vector<IterationRecord> build_trace;   // one entry per build iter
    std::vector<std::size_t>    refine_trace;    // reads remapped per pass
    PhaserStats                 stats;
};

// Reads passed in either as a path to a (gz)FASTQ on disk OR as an
// in-memory id/sequence vector for unit tests. The pipeline opens the
// FASTQ once, builds an offset index, and loads sequences on demand —
// the full dataset is never held in RAM simultaneously.
//
// At WG-HiFi scale this is critical: 4M reads × 15 KB ≈ 60 GB. The old
// vector-resident layout forced peak RSS up to the 140 GB SLURM cap.
struct ReadsInput {
    // Preferred for production runs. If non-empty, the in-memory vector
    // is ignored and the file is indexed.
    std::string fastq_path;

    // Test path: in-memory reads. Used when fastq_path is empty.
    std::vector<std::pair<ReadId, std::string>> reads;
    std::vector<std::string> read_id_strings;   // ReadId index → name

    // Het-pairs from the upstream k-mer histogram / haplotype router.
    // When empty, phase_engine cannot vote on bubbles → all unitigs end
    // up as SHARED. Caller is responsible for populating these with
    // genuine het signal (e.g. the het_pairs vector emitted by
    // find_het_pairs / set_anchor_from_reads in haplotype_router).
    std::vector<HetPair> het_pairs;
};

// Single entry point. Produces output files if paths are set in opts and
// always returns the in-memory PhaserResult.
PhaserResult run(const ReadsInput& reads, const PhaserOpts& opts);

}  // namespace branch::wg::phaser
