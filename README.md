# BRANCH — Lossless Assembly Graph for Low-Frequency CNV Discovery

## Overview
BRANCH is a HiFi-read genome assembler built to be state-of-the-art at low-frequency copy-number variants. It produces a lossless, CN-aware assembly graph where branches are graph bifurcations (not tumor clones). Every variant call carries VAF evidence from reads, in-silico PCR, and k-mer counts.

## Status (Phase 0)
- End-to-end run on HG002 HiFi: Hummel-2 SLURM job 1823405 on node n093, 4:24 min wall time (2026-04-17).
- 8 HPRC genomes processed end-to-end.
- Known gaps: unitig collapse not yet final, RAM consumption non-deterministic, CPU utilisation ~6 %, per-contig chromosome = NA.

## Pipeline
```
reader → graph_build → graph_compactor → graph_filter → assemble
```

## Output contract
- BED: branch intervals (chrom, start, end, branch_id, VAF, CN).
- Consensus FASTA per branch.
- VAF evidence channels: supporting reads, primer-bracketed in-silico PCR amplicons, k-mer counts on read sequence.
- Genome-wide repeat CN for main path and every branch, normalised against single-copy reference amplicons.

## Build
```
cmake -S . -B build
cmake --build build -j
```

## Running Phase 0 on Hummel-2
```
sbatch workflow/slurm/phase0.slurm <hifi.fastq>
```
Outputs: `.gfa`, `.fasta`, `.bed`, `.paf`.

## Tech stack
- C++20 core (performance-driven decision).
- htslib (BAM/CRAM I/O), ksw2 (affine-gap alignment), abPOA (partial-order consensus).
- CMake build; ASan + TSan required in CI.

## Repository layout
- `src/` — core C++ sources.
- `docs/architecture.md` — pipeline internals, graph data model, classification problem.
- `docs/graph-format-spec.md` — on-disk graph format.
- `workflow/` — Snakemake / SLURM driver.
- `phase0/` — Phase 0 artefacts.
- `tests/` — unit + integration tests.

## Design principles
Lossless graph, CN-aware nodes, VAF-tagged edges, SV-first phasing, multi-allelic branches.
