# BRANCH — Lossless Assembly Graph for Low-Frequency CNV Discovery

## Overview
BRANCH is a HiFi-read genome assembler built to be state-of-the-art at low-frequency copy-number variants. It produces a lossless, CN-aware assembly graph where branches are graph bifurcations (not tumor clones). Every variant call carries VAF evidence from reads, in-silico PCR, and k-mer counts.

## Status
Phase 0 end-to-end pipeline is working on HiFi samples. Known gaps: unitig collapse not yet final, RAM consumption not fully deterministic, CPU utilisation low (single-threaded overlap), per-contig chromosome projection pending (see `branch project` below).

## Pipeline
```
reader → graph_build → graph_compactor → graph_filter → assemble
```

## Subcommands
- `branch assemble` — reads → minimizer overlap → lossless graph → GFA + FASTA + PAF + BED.
- `branch analyze` — mosdepth regions → copy-number inference with paralog awareness.
- `branch project` *(v0.4, scaffolded)* — three-layer reference projection: linear (CHM13 + GRCh38 via minimap2), pangenome (HPRC v1.1 via minigraph / GraphAligner), somatic delta vs. nearest pangenome path. See `docs/branch-project-design.md`.

## Output contract
- BED: branch intervals (chrom, start, end, branch_id, VAF, CN).
- Consensus FASTA per branch.
- VAF evidence channels: supporting reads, primer-bracketed in-silico PCR amplicons, k-mer counts on read sequence.
- Genome-wide repeat CN for main path and every branch, normalised against single-copy reference amplicons.

## Build
```
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```
Requires C++20, CMake ≥ 3.20, zlib, pthread. Vendored in `third_party/`: htslib, ksw2, abPOA. CUDA backend is opt-in via `-DBRANCH_BUILD_CUDA=ON`.

## Running on an HPC cluster
A `sbatch`-compatible driver lives in `workflow/`; adapt the partition, account, and paths to your own site. The pipeline is filesystem-agnostic — point `--fastq` / `--bam` at your reads and `--out` at a writable output directory.

## Tech stack
- C++20 core.
- htslib (BAM/CRAM I/O), ksw2 (affine-gap alignment), abPOA (partial-order consensus).
- CMake build; ASan + TSan required in CI.

## Repository layout
- `src/` — core C++ sources.
- `docs/architecture.md` — pipeline internals, graph data model, classification problem.
- `docs/branch-project-design.md` — reference-projection subcommand design.
- `docs/cnv_roadmap.md` — low-frequency CNV roadmap.
- `docs/graph-format-spec.md` — on-disk graph format.
- `workflow/` — SLURM / Snakemake driver scaffold.
- `tests/` — unit + integration tests (GoogleTest).

## Design principles
Lossless graph, CN-aware nodes, VAF-tagged edges, SV-first phasing, multi-allelic branches.
