#!/bin/bash
#============================================================================
# BRANCH Phase 0 Task T5: HG002 HiFi Reference Alignment
# Aligns 6 HiFi FASTQ files to GRCh38 using minimap2 map-hifi
#============================================================================
#SBATCH --job-name=branch_t5_align
#SBATCH --partition=std
#SBATCH --account=kubisch_std
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --export=NONE
#SBATCH --output=/beegfs/u/%u/humangenetik/ag/schlein/%u/logs/slurm_jobs/t5_align_%j.out
#SBATCH --error=/beegfs/u/%u/humangenetik/ag/schlein/%u/logs/slurm_jobs/t5_align_%j.err

set -eo pipefail  # no -u here: RRZ module init uses unbound vars

#--- Environment Setup (isolated from set -u) ---
set +u
source /sw/batch/init.sh
module load apptainer/1.2.5
set -u  # re-enable strict-unbound check after module init

#--- Paths ---
# Shared references (hardcoded - cluster-wide)
INPUT_DIR="/beegfs/u/bbg6775/humangenetik/shared/references/human_pangenome/hifi_reads/HG002/15kb"
REFERENCE="/beegfs/u/bbg6775/humangenetik/shared/references/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# User-specific paths (dynamic)
BEEGFS_BASE="/beegfs/u/${USER}/humangenetik/ag/schlein/${USER}"
CONTAINER_DIR="${BEEGFS_BASE}/containers"
OUTPUT_DIR="${BEEGFS_BASE}/projects/branch/phase0"
LOG_DIR="${BEEGFS_BASE}/logs/slurm_jobs"

# Containers
MINIMAP2_SIF="${CONTAINER_DIR}/minimap2_2.28.sif"
SAMTOOLS_SIF="${CONTAINER_DIR}/samtools_1.21.sif"

# Output files
OUTPUT_BAM="${OUTPUT_DIR}/HG002_grch38.bam"
OUTPUT_BAI="${OUTPUT_BAM}.bai"

#--- Pre-flight checks ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] T5 Alignment starting"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "CPUs: ${SLURM_CPUS_PER_TASK:-16}"
echo "Node: $(hostname)"

# Verify inputs exist
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR" >&2
    exit 1
fi

if [ ! -f "$REFERENCE" ]; then
    echo "ERROR: Reference not found: $REFERENCE" >&2
    exit 1
fi

if [ ! -f "$MINIMAP2_SIF" ]; then
    echo "ERROR: minimap2 container not found: $MINIMAP2_SIF" >&2
    exit 1
fi

if [ ! -f "$SAMTOOLS_SIF" ]; then
    echo "ERROR: samtools container not found: $SAMTOOLS_SIF" >&2
    exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

# Count input files
FASTQ_COUNT=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fastq" -type f 2>/dev/null | wc -l)
echo "Found $FASTQ_COUNT FASTQ files in $INPUT_DIR"

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "ERROR: No FASTQ files found" >&2
    exit 1
fi

#--- Alignment Pipeline ---
# Strategy: cat all FASTQ | minimap2 | samtools sort
# This avoids redundant reference loading (minimap2 loads ~8GB index once)

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting alignment pipeline"
echo "Input: $INPUT_DIR/*.fastq"
echo "Reference: $REFERENCE"
echo "Output: $OUTPUT_BAM"

# Single pipeline: cat → minimap2 → samtools sort
# minimap2 uses SLURM_CPUS_PER_TASK threads, samtools uses 4 for compression
cat "${INPUT_DIR}"/*.fastq | \
    apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$MINIMAP2_SIF" \
        minimap2 -ax map-hifi --secondary=no -t "${SLURM_CPUS_PER_TASK:-16}" \
        "$REFERENCE" - | \
    apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$SAMTOOLS_SIF" \
        samtools sort -@ 4 -O bam -o "$OUTPUT_BAM" -

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment complete, indexing BAM"

#--- Index BAM ---
apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$SAMTOOLS_SIF" \
    samtools index -@ 4 "$OUTPUT_BAM"

#--- Generate stats ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Generating alignment stats"

apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$SAMTOOLS_SIF" \
    samtools flagstat "$OUTPUT_BAM" > "${OUTPUT_BAM%.bam}.flagstat"

apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$SAMTOOLS_SIF" \
    samtools idxstats "$OUTPUT_BAM" > "${OUTPUT_BAM%.bam}.idxstats"

#--- Summary ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] T5 Alignment COMPLETE"
echo "Output files:"
ls -lh "${OUTPUT_DIR}/HG002_grch38"* 2>/dev/null || echo "(listing failed)"

echo ""
echo "=== Flagstat Summary ==="
head -5 "${OUTPUT_BAM%.bam}.flagstat" 2>/dev/null || echo "(flagstat not available)"
