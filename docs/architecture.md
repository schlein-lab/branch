# BRANCH Architecture

## Pipeline
```
HiFi FASTQ
   |
reader ---------> raw reads (streamed)
   |
graph_build ----> overlap graph (k-mer seeds + lossless edges)
   |
graph_compactor -> unitig compaction (collapse linear chains)
   |
graph_filter ---> transitive reduction + containment removal
   |
assemble -------> BED + consensus FASTA + VAF + repeat CN
```

## Graph data model
Nodes are CN-aware: each node carries a copy-number estimate from read depth plus optional k-mer support. Edges carry VAF tags computed from the reads that span the junction. The graph is lossless — reads are never discarded, only re-routed.

Branches are graph bifurcations. A branch is triggered when two parallel paths both have enough read support to be distinguished from sequencing error.

## Classification problem
Central question: is a pair of parallel paths a true branch, a duplication, or a mixed state? Discriminators, applied hierarchically when the call is mixed:
- Flanks: flanking sequence identity and coverage outside the branch region.
- Depth-sum: do the two paths together sum to the expected diploid depth?
- Read-span: how many reads fully span either path?
- Prior: known repeat content and CN from the reference (fall-back).

## Output contract
- BED columns: chrom, start, end, branch_id, VAF (0..1), CN (integer or fractional).
- Consensus FASTA per branch, one record per allele.
- VAF evidence: reads-supporting, in-silico PCR with primer-bracketed amplicons, k-mer count on read sequence.
- Repeat CN: genome-wide, main path and each branch, normalised by single-copy reference amplicons.

## Tech stack in practice
- C++20, CMake. Build: `cmake -S . -B build && cmake --build build -j`.
- htslib (BAM/CRAM), ksw2 (affine-gap alignment), abPOA (partial-order consensus).
- ASan + TSan required for every CI run.

## HPC
Hummel-2 cluster, SLURM, BeeGFS storage. Containers via Apptainer with `--nv` for CUDA. Jobs bind `/beegfs`; `TMPDIR=/tmp`; home is read-only inside jobs.
