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
    std::string out_assignments_tsv_path;
    std::string out_qc_bed_path;
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

// Reads passed in as id-sequence pairs. Sequences may be released by the
// caller after run() returns; the PhaserResult holds only consensus
// sequences in unitigs, never raw read sequences.
struct ReadsInput {
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
