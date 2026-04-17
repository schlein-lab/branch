# CNV Low-Frequency Roadmap

Long-term goal: be state-of-the-art at low-frequency copy-number variant
discovery from human HiFi bulk-blood sequencing. This document translates that
ambition into a concrete sequence of milestones anchored to measurable
benchmarks.

## What "low-frequency" means here

The target regime is clonal hematopoiesis and somatic mosaicism in bulk blood:
CNVs that are present in only a subset of circulating cells. Concretely:

- VAF range: 1 % – 20 %. Below 1 % is dominated by sequencing artefact in HiFi;
  above 20 % is solved by every diploid assembler.
- Size range: 300 bp – 50 kb. Below 300 bp overlaps with small-indel territory;
  above 50 kb can be recovered by read-depth tools.
- Region class: tandem repeats, segmental duplications, VNTRs, and any locus
  where parallel paths exist in the assembly graph. These are exactly the
  regions where read-mapping-only pipelines fail.

## Baselines we must beat

- **hifiasm-meta** on bulk blood — treats low-frequency branches as noise and
  collapses them.
- **Shasta** on HG002 somatic regions — solid on the haplotypes, but loses
  sub-diploid fractions.
- **Sniffles2 + pbsv on HiFi** — read-mapping tools; blind to CNVs inside
  complex repeats because reads cannot be uniquely placed.
- **Severus / SAVANA** — tumour-clone-oriented, fails the bulk-blood use case
  (trained for subclone structure, not sparse mosaicism).

Beating these means higher recall at VAF ≤ 5 % without trading precision.

## Truth sets

- **GIAB Challenging Medically Relevant Genes** (CMRG) — segmental duplications
  in the medically important region set; gold standard for SD accuracy.
- **HPRC year-1 pangenome variant set** — high-confidence CNV calls across 47
  assemblies; use the subset inside our 8 already-processed HPRC genomes.
- **Known CH loci** — DNMT3A, TET2, ASXL1, JAK2 CNV events from targeted
  long-read CH studies. Expected positives; drives the clinically relevant
  end of the ROC.
- **In-silico spike-in** — tile HPRC assemblies into synthetic mixtures at
  defined VAFs (1 %, 2 %, 5 %, 10 %) to measure recall as a function of VAF.
  This is the primary internal benchmark because real truth at VAF < 5 % does
  not exist at scale.

## Phase 0 → Phase 1

Phase 0 produced an end-to-end pipeline and surfaced four gaps that block any
meaningful CNV benchmark. They are prerequisites, not Phase 1 deliverables:

1. **Unitig collapse** — without real collapse the graph is dominated by
   sequencing-error bubbles, overwhelming the branch classifier.
2. **RAM non-determinism** — makes resource scaling unpredictable on SLURM
   and blocks reliable multi-genome runs.
3. **CPU ≈ 6 %** — 16× headroom is sitting on the table; until this is
   recovered, runtime forbids the spike-in grid we need.
4. **chrom = NA in BED** — output is unusable for any genomic benchmark;
   every truth set is chromosome-keyed.

These are tracked as separate diagnostic findings; see
`coding-zyrkel-1` report (WIP).

## Phase 1 — Measurable CNV discovery

Deliverables, in dependency order:

1. Phase 0 gaps closed (prerequisite).
2. **Spike-in benchmark harness**: a Snakemake workflow that constructs HPRC
   mixtures at VAF ∈ {1, 2, 5, 10, 20} %, runs BRANCH, and emits an ROC curve
   per VAF bin.
3. **HPRC self-benchmark**: recall/precision of BRANCH branches against the
   HPRC pangenome CNV set on the 8 genomes already processed.
4. **VAF precision benchmark**: for every discovered branch, compare the
   reads-supporting VAF, vPCR-amplicon VAF, and k-mer VAF to the true
   spike-in fraction. Report RMSE per evidence channel.
5. **First real CH run**: one bulk-blood HiFi sample (to be acquired) with
   orthogonal ddPCR validation of a handful of BRANCH-discovered CNVs.

Phase 1 is "done" when the spike-in grid shows recall ≥ 0.8 at VAF = 5 % and
precision ≥ 0.9 across the full grid, with per-evidence-channel VAF RMSE < 2 %.

## Phase 2 — Specialisation

Only after Phase 1 is locked.

- **Repeat-CN calibration**: genome-wide normalisation against single-copy
  reference amplicons; report integer + fractional CN jointly.
- **Branch classification refinement**: sharpen the branch vs. duplication vs.
  mixed discriminators using the Phase 1 truth-labelled training data.
- **GPU hotpath** (`src/gpu/`): accelerate the overlap and consensus passes if
  profiling after the CPU-6 % fix still shows wall-time dominated by those.
- **Submission to an external evaluator** (e.g. precisionFDA Truth Challenge
  sub-track) as the external SOTA claim.

## What stays out of scope

Per project scope: tumour subclones, plants, cancer heterogeneity. Everything
in this document assumes normal bulk blood. If a tumour use case arrives later,
it is a separate project, not an extension of this roadmap.
