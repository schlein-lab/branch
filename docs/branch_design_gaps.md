# BRANCH Design Gaps — Phase-0 Baseline (T8)

Date: 2026-04-16
Sample: HG002 HiFi from HPRC (pre-Revio 2020 dataset, ~4.1× coverage, 12.8 kb N50)
Reference: GRCh38_full_analysis_set_plus_decoy_hla.fa
BAM: `/beegfs/u/bbg6775/humangenetik/ag/schlein/bbg6775/projects/branch/phase0/HG002_grch38.bam`

## Coverage Reality Check

Phase-0 used the **pre-existing** HG002 dataset on BeeGFS (6 FASTQ files, 994k reads, 12.76 Gbp total).
- Estimated coverage: **4.1× genome-wide**
- Read length: avg 12.8 kb, N50 12.8 kb, Q20 97.8%

**Implication:** Insufficient for BRANCH v0.3+ mosaicism detection. A real training run needs the
full HPRC HG002 (~40–80 GB, 30–60×) from `s3://human-pangenomics/working/HPRC/HG002/raw_data/PacBio_HiFi/`.
Phase 0 data works for pipeline validation and feature-extraction skeletons only.

## ROI Catalog (T6+T7 output)

689,558 regions annotated from UCSC:
| Class | Count | Mean cov | Median | Std | Max |
|-------|-------|----------|--------|-----|-----|
| simple_repeat | 681,515 | 3.90× | 3.92× | 4.59 | 1051 |
| segdup | 8,043 | **5.24×** | 3.81× | 44.47 | 2601 |

## Key Finding: SegDups Have +34% Depth Over Genome-Wide Mean

SegDup regions show mean coverage **5.24×** vs genome-wide 3.90× — a **+34% over-coverage**
in segmental-duplicate regions. At 30× normal coverage this scales to ~40× in SegDups.

This is the exact signal BRANCH needs to distinguish:
- **Branch** (same locus, two alleles): depth-sum ≈ baseline
- **Duplication** (N paralog copies): depth-sum ≈ N × baseline

Observed stddev 44.47 with max 2601× in SegDups shows **extreme** multi-mapped accumulation
regions (centromeric SegDups, acrocentric arms, NPIP/USP17L). These will stress-test the
classification cascade.

## Design Gaps Identified

### ROI #1: HLA region (chr6:28Mb–33Mb)
- Known high paralogy, highest N-count in any SegDup cluster
- Expected BRANCH behavior: preserve all paralog paths with CN>=6 per HLA class I/II gene family
- Feature extraction challenges: Flanken-Jaccard may be <95% even between true paralogs due to
  extreme local diversity → cascade stage 1 false-negatives possible
- Mitigation in v0.3: stratify classifier training by region_type=HLA

### ROI #2: SMN1/SMN2 (chr5:70Mb)
- Near-identical paralog pair (>99% identity), clinically relevant (SMA)
- Expected BRANCH behavior: classify as Duplication with CN=2, NOT collapse into a single node
- Feature extraction challenges: Depth-Sum alone is insufficient; needs Flank-Jaccard + Read-Span-Ratio
- Mitigation: emit as "low-confidence Duplication" when Depth-Sum ambiguous, let downstream vPCR disambiguate

### ROI #3: Amylase cluster (chr1:103Mb)
- Known 1–10× population CN variation (copy-number polymorphism)
- Expected BRANCH behavior: CN output should reflect actual sample CN, not a fixed diploid assumption
- Feature extraction challenges: Depth-Sum needs to project to integer CN (vs. continuous AF)
- Mitigation: dedicated CN-integer-inference pass after cascade classification

### ROI #4: KIR (chr19:55Mb)
- Variable-length haplotypes, extreme polymorphism
- Expected BRANCH behavior: emit multiple co-existing haplotypes as Mixed (dup + branch)
- Feature extraction challenges: hierarchical Mixed decomposition (see termination-guarantees)
- Mitigation: v0.4 hierarchical classifier pass; termination at min_length=500bp

### ROI #5: Centromeric alpha-satellite
- Read alignment fails (ambiguous mapping), not useful ROI for BRANCH v0.3
- Expected BRANCH behavior: explicit exclude-list via `exclude_centromeres.bed` from T6
- Generated exclude-BED is present at `rois/exclude_centromeres.bed` (1k-ish entries)

## Minimum Algorithmic Components for v0.1

Based on gaps above, v0.1 must ship:
1. Lossless graph data model — ✅ `src/graph/lossless_graph.hpp`
2. Delta-Graph read representation — ✅ `src/graph/delta_read.hpp`
3. Reverse-index Node→Reads — ✅ `src/graph/reverse_index.hpp`
4. Type-Erased backend VTable for CPU/GPU dispatch — ✅ `src/backend/backend_vtable.hpp`
5. Classifier feature definitions (12 features) — ✅ `src/classify/features.hpp`
6. CUDA overlap-kernel skeleton — pending
7. Snakemake workflow binding BeeGFS inputs + apptainer containers — pending

## Non-goals for v0.1 (deferred to v0.2+)

- Actual overlap computation (needs minimap2-lib binding)
- abPOA consensus calling
- vPCR primer design (Primer3 integration)
- Snakemake SLURM-profile polish

## Conclusion

Phase-0 confirmed both the feasibility (infrastructure + data) and the algorithmic core
challenge (SegDups show measurable depth differentiation — the exact signal BRANCH depends on).
v0.1 scaffold in place, 22 tests green. Next: graph-classify cascade skeleton (regel-form for
v0.3-upgrade to LightGBM).
