# branch project — Three-Layer Reference Projection

## Overview

The `branch project` subcommand projects branches from the lossless graph (output of
`branch assemble`) against three reference layers to identify somatic and low-frequency
variants that diverge from known population haplotypes.

**Version target**: v0.4 (scaffold in v0.3)

## Motivation

Standard variant calling assumes diploid genomes with variants present at ~50% or ~100%
VAF. BRANCH targets somatic mosaicism and low-frequency structural variants (1-20% VAF)
which require:

1. Assembly of variant-supporting reads into "branches" (done by `branch assemble`)
2. Projection of branches against multiple reference layers to distinguish:
   - Known population variants (present in HPRC pangenome)
   - Reference-specific artifacts (CHM13 vs GRCh38 differences)
   - True somatic/mosaic variants (absent from HPRC, present in sample)

## Three Reference Layers

### 1. Linear Reference (CHM13 + GRCh38)

**Tool**: minimap2 with `asm20` preset (assembly-to-reference)

**Output**: `<sample>.linear.paf`

Each branch consensus sequence is aligned to both CHM13 and GRCh38. The PAF output
includes a custom `rf:Z:` tag indicating which reference (`CHM13` or `GRCh38`).

**Purpose**:
- Establish genomic coordinates for each branch
- Identify branches that align better to one reference (potential SV or assembly gap)
- Filter branches with MAPQ < 20 on both references as potentially novel

### 2. Pangenome Reference (HPRC v1.1)

**Tool**: minigraph `-c` (fast, CIGAR-free) or GraphAligner (HiFi precision)

**Output**: `<sample>.pangenome.gaf`

Branches are aligned to the HPRC v1.1 Minigraph-Cactus pangenome graph (GBZ format).
This identifies which known haplotype path each branch most closely matches.

**Purpose**:
- Identify the closest HPRC haplotype for each branch
- Branches matching a known haplotype path are "explained" variants
- Branches with no close match are candidate novel/somatic variants

### 3. Somatic Delta Layer

**Tool**: Internal C++ implementation using ksw2/WFA2

**Output**:
- `<sample>.somatic.vcf` — VCF with SNV/INDEL/SV calls
- `<sample>.branch-report.json` — Detailed per-branch analysis

For each branch, compute edit distance between:
- Branch consensus sequence
- Closest HPRC haplotype path sequence

Edits represent potential somatic variants. These are emitted as VCF records with
branch-derived VAF and coverage annotations.

## Input Files

| File | Source | Description |
|------|--------|-------------|
| `<sample>.gfa` | `branch assemble --out` | Lossless graph with branch nodes |
| `<sample>.fasta` | `branch assemble --fasta` | Branch consensus sequences |
| CHM13.mmi | External | minimap2 index for CHM13 |
| GRCh38.mmi | External | minimap2 index for GRCh38 |
| HPRC-v1.1.gbz | External | HPRC pangenome graph |

## Output Files

| File | Format | Description |
|------|--------|-------------|
| `<prefix>.linear.paf` | PAF + rf:Z: tag | Linear alignments to both references |
| `<prefix>.pangenome.gaf` | GAF | Pangenome graph alignments |
| `<prefix>.somatic.vcf` | VCF 4.2 | Somatic variant calls |
| `<prefix>.branch-report.json` | JSON | Per-branch detailed analysis |

### branch-report.json Schema

```json
{
  "sample": "HG002",
  "branches": [
    {
      "branch_id": "branch_001",
      "length_bp": 15234,
      "coverage": 45.2,
      "vaf": 0.12,
      "linear_alignments": {
        "CHM13": {"mapq": 60, "identity": 0.998},
        "GRCh38": {"mapq": 55, "identity": 0.995}
      },
      "closest_hprc_path": "HG00733#1#chr1:1000000-1015234",
      "somatic_edits": [
        {"type": "SNV", "pos": 5000, "ref": "A", "alt": "G"},
        {"type": "DEL", "pos": 8000, "ref": "ACGT", "alt": "A"}
      ],
      "unannotated": false
    }
  ],
  "summary": {
    "total_branches": 1234,
    "annotated": 1100,
    "unannotated": 134,
    "somatic_snvs": 5600,
    "somatic_indels": 890,
    "somatic_svs": 45
  }
}
```

## Module Structure

```
src/project/
├── project_cmd.hpp        # CLI interface (run_project declaration)
├── project_cmd.cpp        # CLI parsing + orchestration
├── linear_mapper.hpp      # minimap2 shell-out interface
├── linear_mapper.cpp      # minimap2 invocation + PAF merge
├── pangenome_mapper.hpp   # minigraph/GraphAligner interface
├── pangenome_mapper.cpp   # Pangenome alignment + GAF parse
├── somatic_delta.hpp      # Edit distance computation
├── somatic_delta.cpp      # ksw2-based delta calculation
├── project_report.hpp     # JSON/VCF output
├── project_report.cpp     # Report generation
└── CMakeLists.txt         # Module build config
```

## CLI Interface

```
branch project — three-layer reference projection

Usage:
  branch project --gfa <path.gfa> --fasta <path.fasta>
                 --ref-linear <path.mmi> [--ref-linear <path2.mmi>]
                 --ref-pangenome <path.gbz>
                 --out-prefix <path>
                 [--mapper linear=minimap2,pangenome=minigraph]
                 [--threads <N>]
                 [--emit-bam]

Required:
  --gfa <path>           GFA from branch assemble (nodes = branches)
  --fasta <path>         Branch consensus FASTA
  --ref-linear <path>    Linear reference (minimap2 .mmi or .fa), repeatable
                         Use name= prefix: --ref-linear name=CHM13:/path/chm13.mmi
  --ref-pangenome <path> Pangenome GBZ file, repeatable
  --out-prefix <path>    Output prefix for all result files

Optional:
  --mapper <spec>        Mapper selection (default: linear=minimap2,pangenome=minigraph)
                         Alternatives: pangenome=graphaligner
  --threads <N>          Thread count (default: available cores)
  --emit-bam             Also emit BAM (default: off, PAF/GAF only)
  --min-mapq <N>         Minimum MAPQ for annotation (default: 20)
  --help                 Show this help
```

## Unmapped / Unannotated Branches

Branches are marked as `UNANNOTATED` when:
- MAPQ < 20 on ALL linear references, OR
- No unique pangenome path match (multiple equally-close paths)

These represent the core value proposition: **novel sequences not present in HPRC**.

Unannotated branches should be:
1. Manually reviewed (IGV with BAM if `--emit-bam` was used)
2. Cross-referenced with known difficult regions (centromeres, segdups)
3. Validated with orthogonal methods (long-range PCR, Sanger)

## Parallelism Strategy

| Stage | Parallelism | Rationale |
|-------|-------------|-----------|
| Linear mapping | Per-branch | Independent alignments, N threads |
| Dual-ref | Sequential per branch | minimap2 instances are memory-heavy |
| Pangenome | Dual-graph parallel | CHM13-graph + GRCh38-graph simultaneously |
| Somatic delta | Per-branch | Edit distance is CPU-bound, parallelizes well |

## External Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >= 2.26 | Linear reference alignment |
| minigraph | >= 0.20 | Fast pangenome alignment |
| GraphAligner | >= 1.0.17 | High-precision pangenome alignment (optional) |

All external tools are invoked via shell-out (not linked). This provides:
- Version flexibility
- Easy debugging (intermediate files visible)
- No C++ ABI compatibility concerns

## Testing Strategy

### Unit Tests

- `somatic_delta::edit_distance()` against known ref/branch pairs
- PAF/GAF parser correctness
- JSON/VCF output format validation

### Integration Tests

- `branch assemble` -> `branch project` on HG002 smoke sample
- Verify: GRCh38 PAF has alignments
- Verify: HPRC GAF has alignments
- Verify: At least 1 somatic edit detected

### Smoke Tests

- `branch project --help` exits 0
- `branch project` (no args) exits 2 with usage

## Future Extensions (v0.5+)

- **Phasing integration**: Link branches to haplotype blocks
- **Population frequency**: Cross-reference with gnomAD SV
- **Methylation**: Integrate with ONT 5mC calls
- **Clinical annotation**: ClinVar/OMIM variant overlay

## References

1. Liao et al. (2023). A draft human pangenome reference. Nature.
2. Li (2023). Minigraph-Cactus pangenome. Nature Biotechnology.
3. Rautiainen & Marschall (2020). GraphAligner. Genome Biology.
