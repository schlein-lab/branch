# BRANCH

**A HiFi long-read genome assembler aimed at low-frequency copy-number variation.**

BRANCH (Breakpoint-Resolved Assembly of Non-diploid Copy-number Heterogeneity) takes PacBio HiFi reads and produces a *lossless, copy-number-aware assembly graph* together with a per-branch evidence bundle: BED intervals, consensus FASTAs, variant-allele-frequency (VAF) channels, and genome-wide repeat copy numbers. It is being developed in the Schlein lab as a research-grade tool; Phase 0 has run end-to-end on HG002 and eight HPRC reference genomes.

This README is written for bioinformaticians who want to understand what BRANCH does and how to run it without first having to read the C++. The internals live in [`docs/architecture.md`](docs/architecture.md) and the on-disk graph format in [`docs/graph-format-spec.md`](docs/graph-format-spec.md).

## Why low-frequency CNVs are hard

Copy-number variants that segregate at high allele frequency (germline het/hom variants, large structural events) are reasonably well caught by current pipelines. The much harder regime is **CNVs present in a small fraction of the reads** — somatic mosaicism in bulk tissue, low-VAF repeat expansions, partial gains in segmental-duplication regions, mixed-population samples, etc. Standard short-read pipelines hit the wall at mappability and repeat collapse; standard long-read assemblers (hifiasm, Verkko, Flye) prioritise contiguous diploid haplotypes and either drop low-frequency alternative paths as noise or fold them into the consensus.

BRANCH inverts that tradeoff. The graph is treated as the primary output, not a scaffold to throw away. Any pair of parallel paths that has enough independent read support to be distinguished from sequencing error is preserved as a *branch*, tagged with its VAF, and propagated to all downstream evidence channels. The resulting bundle answers questions like "is this region heterozygous diploid, an expanded tandem repeat, a low-VAF mosaic event, or a paralog collapse?" with explicit numbers rather than a single best-guess contig.

## What "branch" means here

A **branch** in BRANCH is strictly a property of the assembly graph: a bifurcation where two (or more) parallel paths each carry independent read evidence. It is *not* a biological clone, lineage, or sub-population. The graph just tells you that the reads disagree at this locus in a non-stochastic way; the biological interpretation — heterozygosity, somatic mosaicism, repeat expansion, pooled haplotypes, segmental duplication — is left to downstream callers, which can use the VAF and copy-number annotations BRANCH provides.

Every branch is classified into one of four states (`B` = branch, `D` = duplication, `M` = mixed, `N` = non-separable) using flank identity, depth-sum, read-span, and reference prior, with the per-edge confidence written into the graph. Mixed states are decomposed hierarchically — see [`docs/architecture.md`](docs/architecture.md) for the discriminator details.

## Pipeline

```
HiFi FASTQ ──► reader ──► graph_build ──► graph_compactor ──► graph_filter ──► assemble
```

The five stages are deliberately small and individually testable:

- **reader** streams HiFi reads from FASTQ/BAM/CRAM via htslib, with quality and length filters configurable per-run.
- **graph_build** constructs the lossless overlap graph from minimizer seeds; reads are indexed but never discarded.
- **graph_compactor** collapses linear, unambiguous chains into unitigs (Phase 0: works but not yet final — see Status below).
- **graph_filter** applies transitive reduction and containment removal so only minimal, evidence-bearing edges survive.
- **assemble** walks the filtered graph to produce the output bundle: GFA, per-branch consensus FASTA via abPOA, BED intervals with VAF and CN, and a PAF of branch placements.

## Output bundle

A successful `branch assemble` run produces the following artefacts side-by-side, all keyed by `branch_id` so they can be cross-referenced:

- **`*.bed`** — branch intervals: `chrom, start, end, branch_id, VAF, CN`. Chromosome assignment requires a reference; otherwise `chrom = NA` (Phase 0 limitation, see roadmap).
- **`*.fa`** — per-branch consensus FASTA, one record per allele, built with abPOA partial-order consensus from the supporting reads.
- **VAF evidence**, three independent channels per branch:
  - reads spanning the branch (the direct count),
  - **primer-bracketed in-silico PCR amplicons** that simulate a wet-lab vPCR readout,
  - k-mer counts on the raw read sequence (mapping-free).
- **Repeat copy numbers**, genome-wide, for the main path and *every* branch, normalised against single-copy reference amplicons so absolute CN is comparable across runs.
- **`*.gfa`** — the full lossless graph in our [GFA-1.2 extension](docs/graph-format-spec.md), with CN-aware node tags (`CN`, `CV`, `CC`, `HP`) and VAF-tagged edges (`VF`, `VC`, `BT`, `CT`, `FS`).
- **`*.paf`** — branch-to-reference placement when a reference is supplied.

The three VAF channels are deliberately redundant: agreement across read-count, vPCR, and k-mer is a strong signal; disagreement flags a region for manual review.

## Status — Phase 0 versus the roadmap

**Working in Phase 0:** the full `reader → graph_build → graph_compactor → graph_filter → assemble` pipeline runs end-to-end. HG002 HiFi completed on Hummel-2 (SLURM job 1823405, node n093) in 4 minutes 24 seconds on 2026-04-17, and eight HPRC reference genomes have since been processed. GFA, FASTA, BED, and PAF are produced and pass the E2E test. Three subcommands are exposed via the `branch` CLI: `assemble` (full pipeline), `analyze` (mosdepth-based copy-number inference with paralog awareness), and `project` (scaffolded reference projection — linear via minimap2, pangenome via minigraph/GraphAligner, somatic delta vs. nearest pangenome path; see [`docs/branch-project-design.md`](docs/branch-project-design.md)).

**Known rough edges:** unitig compaction is wired in but not yet final; RAM consumption is non-deterministic across reruns; CPU utilisation is around 6 % (overlap stage is still effectively single-threaded); per-contig chromosome assignment is `NA` until `branch project` is fully wired. The CNV-discovery roadmap toward state-of-the-art low-frequency calling lives in [`docs/cnv_roadmap.md`](docs/cnv_roadmap.md).

## Build

BRANCH is C++20 with vendored dependencies, so the build is self-contained on a modern Linux toolchain:

```sh
# one-time: build the vendored htslib static lib
make -C third_party/htslib -j4 libhts.a

# configure + build + test
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Requirements: a C++20 compiler (GCC ≥ 13 or Clang ≥ 16), CMake ≥ 3.28, zlib, pthread. minimap2's ksw2 sources are pulled via `FetchContent` at configure time. The CUDA backend is opt-in (`-DBRANCH_BUILD_CUDA=ON`, fat-binary for sm_80 / sm_90). Sanitizers are wired but off by default; CI enables ASan and TSan on every run, and so should you before submitting changes (`-DBRANCH_USE_ASAN=ON -DBRANCH_USE_TSAN=ON`, in separate builds).

The compiled CLI lands at `build/src/cli/branch`. It is a single binary with subcommands: `branch assemble --fastq … --out …`, `branch analyze`, `branch project`.

## Running Phase 0 on Hummel-2

The Schlein-lab Phase-0 driver is a Snakemake workflow tailored to the Hummel-2 cluster (RRZ Hamburg). It expects BeeGFS storage (`/beegfs/u/$USER/humangenetik/ag/schlein/$USER`), Apptainer for tool containers, and the SLURM executor plugin. From a Hummel-2 login node:

```sh
# one-off: build the BRANCH binary on the cluster (BeeGFS — never on /home)
cmake -S . -B build && cmake --build build -j

# submit the workflow
./phase0/submit_hummel.sh                       # HG002 default
./phase0/submit_hummel.sh --config sample=NA24385
./phase0/submit_hummel.sh --dry-run             # show DAG without submitting
```

`submit_hummel.sh` enforces the cluster's quirks for you: no `--mem` (RRZ rejects it), `TMPDIR=/tmp` (home is read-only inside batch jobs), `--bind /beegfs` for Apptainer, the SLURM executor profile under `workflow/config/slurm/`, and a 50-job concurrent cap. Off-cluster, the same Snakefile runs via `snakemake -s workflow/Snakefile` once you point `branch_bin`, `hifi_dir`, and the reference paths at your own filesystem — the pipeline itself is filesystem-agnostic.

## Repository layout

- **`src/`** — core C++ sources, organised by stage: `io/`, `graph/`, `align/`, `consensus/`, `classify/`, `detect/`, `analysis/`, `vpcr/`, `project/`, `cli/`, plus `backend/` and `gpu/` for the optional CUDA path.
- **`tests/`** — GoogleTest unit + integration suite, including the end-to-end `test_e2e_assemble.cpp` that exercises the full pipeline.
- **`third_party/`** — vendored htslib (BAM/CRAM/VCF I/O). ksw2 (affine-gap alignment) is fetched from minimap2 at configure time; abPOA (partial-order consensus) is built in-tree.
- **`workflow/`** — Snakemake pipeline (`Snakefile`, `config/`) for cluster runs.
- **`phase0/`** — Hummel-2 submit helper plus Phase-0 task done-files.
- **`docs/`** — design and reference docs:
  - [`architecture.md`](docs/architecture.md) — pipeline internals, graph data model, classification.
  - [`graph-format-spec.md`](docs/graph-format-spec.md) — GFA-1.2 extension with CN/VAF tags.
  - [`branch-project-design.md`](docs/branch-project-design.md) — `branch project` subcommand design.
  - [`cnv_roadmap.md`](docs/cnv_roadmap.md) — low-frequency CNV roadmap.
- **`knowledge/`** — project knowledge catalog (SQLite).

## Design principles

The graph is **lossless**: reads are indexed and re-routed, never thrown away. Nodes are **CN-aware** rather than presumed-haploid, edges carry **VAF tags** computed from the reads that span the junction, **SVs are phased first** (small variants are resolved within the SV-defined blocks), and the entire model is **multi-allelic** — there is no implicit "the" alternative allele. Together these choices keep the graph honest about uncertainty, which is the prerequisite for catching low-frequency events that other assemblers smooth away.

## Scope

Primary target is **bulk-blood HiFi** sequencing data. **HPRC LCL reference genomes** serve as positive controls during development. A reference genome is optional for graph assembly but required for chromosome assignment and for the single-copy amplicons used in CN normalisation.

## License

See [`LICENSE`](LICENSE).
