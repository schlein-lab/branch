# ONT R10 read support

`branch assemble` supports two read-technology presets:

| `--read-tech` | k  | w  | bucket cap | tuned for                                  |
|---------------|----|----|-----------:|--------------------------------------------|
| `hifi` (default) | 21 | 19 |        — | PacBio HiFi (~0.1 % residual error)        |
| `ont`            | 15 | 10 |       64 | ONT R10.4.1 SUP / DUPLEX (~1-2 % residual) |

The HiFi preset is what the assembler used historically; calling
`branch assemble` without `--read-tech` keeps that exact behaviour.

## Why ONT needs a separate preset

At HiFi parameters the minimizer sketch produces ≈ 1 minimizer per 10 bp
on noise-free input. ONT R10.4.1 SUP basecalls have ≈ 1-2 % residual
error after Dorado correction — at k=21 most error events split a
minimizer, so the same input region produces dramatically more
minimizer hits, *and* the hits are dominated by non-informative
multi-occurrence k-mers from basecaller-error patterns. The
all-pairs match table grows quadratically in bucket size and the
process OOMs.

Concrete failure mode observed on HG002 chr14 IGH region (1959 reads,
37 Mbp ONT R10.4.1 SUP basecalls):

| Cap     | Profile  | Outcome | Peak RSS | Time |
|---------|----------|---------|---------:|-----:|
| 8 GiB   | hifi     | OOM     | 6.8 GiB  | 4.6 s |
| 32 GiB  | hifi     | OOM     | 26.1 GiB | 12.1 s |
| 16 GiB  | **ont**  | success | 2.5 GiB  | 40.7 s |

The `ont` preset drops 15 799 high-multiplicity buckets (~2.9 M hits)
before the all-pairs phase, which is what keeps RAM bounded.

## Quickstart

```sh
# Pull a chr14 IGH slice from the ONT Open Data 2025.01 release.
# .bai index is downloaded once; samtools view streams the slice via
# HTTP range requests so we never store the full 172 GB BAM locally.
mkdir -p ont_test && cd ont_test
wget -q "https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam.bai" \
     -O HG002_PAW70337.bam.bai
samtools view -b "https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam" \
     chr14:105400000-106500000 -o HG002_igh.bam
samtools fastq HG002_igh.bam 2>/dev/null | gzip > HG002_igh.fastq.gz

# Assemble with the ONT preset.
zcat HG002_igh.fastq.gz | \
    branch assemble --fastq /dev/stdin --out HG002_ont.gfa \
                    --bed HG002_ont.bed --read-tech ont \
                    --max-memory 16G
```

## Reference

The ONT preset's k=15, w=10 mirrors minimap2's `map-ont` defaults; the
bucket cap of 64 was chosen so that 99th-percentile bucket sizes on the
HG002 IGH slice stay below cap (measured: P99 = 51 hits without cap).
A future tightening pass can move the cap into the per-locus repeat
masker once `branch project --pangenome` lands.
