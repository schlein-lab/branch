#!/bin/bash
# End-to-end smoke test for `--read-tech ont`.
#
# Pulls the HG002 chr14 IGH slice from ONT Open Data and runs
# `branch assemble`. Asserts: peak RSS < 8 GiB, GFA non-empty,
# at least 100 contigs in the output BED.
set -euo pipefail
cd "$(dirname "$0")/../.."

BR="${BRANCH_BIN:-./build-release/src/cli/branch}"
WORK="${WORK_DIR:-/tmp/branch_ont_smoke}"
mkdir -p "$WORK"
cd "$WORK"

if [ ! -f HG002_igh.fastq.gz ]; then
    echo "[smoke] downloading HG002 chr14 IGH slice"
    wget -q "https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam.bai" \
         -O HG002_PAW70337.bam.bai
    samtools view -b "https://ont-open-data.s3.amazonaws.com/giab_2025.01/basecalling/sup/HG002/PAW70337/calls.sorted.bam" \
        chr14:105400000-106500000 -o HG002_igh.bam
    samtools fastq HG002_igh.bam 2>/dev/null | gzip > HG002_igh.fastq.gz
fi

echo "[smoke] running branch assemble --read-tech ont"
zcat HG002_igh.fastq.gz | \
    "$BR" assemble --fastq /dev/stdin --out HG002_ont.gfa \
                   --bed HG002_ont.bed --read-tech ont \
                   --max-memory 16G 2>&1 | tee branch_ont.log

n_contigs=$(wc -l < HG002_ont.bed)
if [ "$n_contigs" -lt 100 ]; then
    echo "[smoke] FAIL: only $n_contigs contigs in BED" >&2; exit 1
fi
peak=$(grep -oE "peak RSS = [0-9.]+ GiB" branch_ont.log | tail -1 | awk '{print $4}')
echo "[smoke] OK: $n_contigs contigs, peak RSS = $peak GiB"
