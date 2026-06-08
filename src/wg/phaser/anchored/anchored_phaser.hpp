#pragma once
//
// AnchoredPhaser — pangenome-anchored read phasing (v0.8 second track).
//
// Maps reads against a reference pangenome (two haplotype FASTAs) using the
// in-process LLmap classical engine and emits per-read haplotype assignments
// based on which haplotype each read aligns to best.
//
// This is NOT ground-truth — it's a fast first model with its own bias
// (HPRC sample composition, hifiasm assembly errors, annotation drift).
// Its output is one input to the Reconciliation Engine, where it gets
// confronted with the de-novo phaser's independent calls.
//
// v0.10: migrated off the external minimap2 binary entirely. Reads are
// streamed from FastqIndex in bounded batches and mapped in parallel via
// LLmap (branch::align::llmap_backend::ReferenceMapper) — no subprocess,
// no PAF temp files, no minimap2 dependency.

#include "../fastq_index.hpp"
#include "../types.hpp"
#include <string>
#include <vector>

namespace branch::wg::phaser {

struct AnchoredPhaserOpts {
    // Pangenome reference: two haplotype FASTA files (hap1, hap2) from a
    // representative HPRC sample. Future versions accept a true GFA graph.
    std::string ref_hap1_fa_path;
    std::string ref_hap2_fa_path;

    int n_threads = 8;
    // A hit counts only if it spans at least this many query bases (drops
    // tiny spurious alignments before the mismatch comparison).
    int min_aligned_bp = 500;
    // Require at least this net mismatch difference between the two haps to
    // call a haplotype; below it the read covers no discriminating site and
    // is SHARED (homozygous span).
    int min_discriminating_sites = 1;

    // LLmap mapping preset, matched to read tech:
    //   "map-hifi" (k=19,w=19, identity≥0.85) for PacBio HiFi/CCS,
    //   "map-ont"  (k=15,w=10, identity≥0.70) for ONT R10.
    std::string preset = "map-hifi";

    // Reads processed per parallel mapping batch. Bounds peak memory
    // (≈ batch × mean_read_len) while keeping the thread pool saturated.
    std::size_t batch_size = 4096;
};

class AnchoredPhaser {
public:
    // Run the full anchoring pass:
    //   1. Build an LLmap minimizer index for hap1 and hap2.
    //   2. Stream reads from FastqIndex in bounded batches; map each batch
    //      against both haplotypes in parallel.
    //   3. For each read: pick the winning haplotype by matched-base weight.
    //   4. Emit a ReadAssignment (H1_ONLY / H2_ONLY / SHARED / UNCERTAIN).
    void run(
        FastqIndex& reads,
        std::vector<ReadAssignment>& out_assignments,
        const AnchoredPhaserOpts& opts);
};

}  // namespace branch::wg::phaser
