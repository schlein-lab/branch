#!/bin/bash
#============================================================================
# BRANCH Phase 0 Task T6+T7: ROI annotation catalog + coverage profiles
# Uses PRE-EXISTING UCSC tracks on BeeGFS, no downloads
#============================================================================
#SBATCH --job-name=branch_t6_t7_roi
#SBATCH --partition=std
#SBATCH --account=kubisch_std
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --export=NONE
#SBATCH --output=/beegfs/u/%u/humangenetik/ag/schlein/%u/logs/slurm_jobs/t6_t7_roi_%j.out
#SBATCH --error=/beegfs/u/%u/humangenetik/ag/schlein/%u/logs/slurm_jobs/t6_t7_roi_%j.err

set -eo pipefail

# RRZ-init tolerant of unbound vars
set +u
source /sw/batch/init.sh
module load apptainer/1.2.5
set -u

# --- Paths (dynamic) ---
USER_BASE="/beegfs/u/${USER}/humangenetik/ag/schlein/${USER}"
PROJECT_DIR="${USER_BASE}/projects/branch/phase0"
ROI_DIR="${PROJECT_DIR}/rois"
CONTAINER_DIR="${USER_BASE}/containers"

# --- Inputs (shared references, cluster-wide) ---
BAM="${PROJECT_DIR}/HG002_grch38.bam"
SHARED_ANNOT="/beegfs/u/bbg6775/humangenetik/shared/annotations"
SEGDUP_UCSC="${SHARED_ANNOT}/segdups/hg38_genomicSuperDups.txt.gz"
STR_UCSC="${SHARED_ANNOT}/repeatmasker/hg38_simpleRepeat.txt.gz"
CENTRO_UCSC="${SHARED_ANNOT}/centromeres/hg38_centromeres.txt.gz"

MOSDEPTH_SIF="${CONTAINER_DIR}/mosdepth_0.3.13.sif"
BEDTOOLS_SIF="${CONTAINER_DIR}/bedtools_2.31.1.sif"

# --- Output files ---
SEGDUP_BED="${ROI_DIR}/roi_segdups.bed"
STR_BED="${ROI_DIR}/roi_simplerepeat.bed"
CENTRO_BED="${ROI_DIR}/exclude_centromeres.bed"
ROI_CATALOG="${ROI_DIR}/roi_catalog.bed"
COV_PREFIX="${ROI_DIR}/HG002_coverage"

# --- Preflight ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] T6+T7 starting — JobID ${SLURM_JOB_ID:-local}"
mkdir -p "$ROI_DIR"

for f in "$BAM" "$SEGDUP_UCSC" "$STR_UCSC" "$CENTRO_UCSC" "$MOSDEPTH_SIF" "$BEDTOOLS_SIF"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f" >&2; exit 1; }
done

# --- Step 1: UCSC → BED (chromStart is 0-based, chromEnd exclusive — matches UCSC table dump) ---
# UCSC genomicSuperDups columns: bin(1) chrom(2) chromStart(3) chromEnd(4) name(5) ... otherChrom(8)
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Converting UCSC tracks to BED"
zcat "$SEGDUP_UCSC" \
    | awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr/ {print $2, $3, $4, "segdup_"NR, "segdup"}' \
    | sort -k1,1 -k2,2n > "$SEGDUP_BED"

# UCSC simpleRepeat columns: bin(1) chrom(2) chromStart(3) chromEnd(4) name(5) period(6) copyNum(7) ...
zcat "$STR_UCSC" \
    | awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr/ {print $2, $3, $4, "trf_"$6"_"$7, "simple_repeat"}' \
    | sort -k1,1 -k2,2n > "$STR_BED"

# UCSC centromeres columns: bin(1) chrom(2) chromStart(3) chromEnd(4) name(5)
zcat "$CENTRO_UCSC" \
    | awk -F'\t' 'BEGIN{OFS="\t"} $2 ~ /^chr/ {print $2, $3, $4, "centromere", "centromere"}' \
    | sort -k1,1 -k2,2n > "$CENTRO_BED"

SEGDUP_N=$(wc -l < "$SEGDUP_BED")
STR_N=$(wc -l < "$STR_BED")
CENTRO_N=$(wc -l < "$CENTRO_BED")
echo "  SegDup entries: $SEGDUP_N"
echo "  SimpleRepeat entries: $STR_N"
echo "  Centromere entries: $CENTRO_N"

# --- Step 2: Merge ROI catalog (SegDup + STR, NOT centromeres) ---
# First column-3 bed for simple merge, drop overlapping intervals per-class
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Building merged ROI catalog"
cat "$SEGDUP_BED" "$STR_BED" \
    | sort -k1,1 -k2,2n \
    > "${ROI_CATALOG}.raw"

# Use bedtools merge per-class to collapse adjacents, then concat
apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$BEDTOOLS_SIF" \
    bedtools merge -i "$SEGDUP_BED" -c 5 -o distinct > "${ROI_DIR}/segdup_merged.bed"

apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$BEDTOOLS_SIF" \
    bedtools merge -i "$STR_BED" -c 5 -o distinct > "${ROI_DIR}/str_merged.bed"

cat "${ROI_DIR}/segdup_merged.bed" "${ROI_DIR}/str_merged.bed" \
    | sort -k1,1 -k2,2n > "$ROI_CATALOG"

ROI_N=$(wc -l < "$ROI_CATALOG")
echo "  Merged ROI catalog: $ROI_N regions"

# --- Step 3: mosdepth per-region coverage ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Running mosdepth"
cd "$ROI_DIR"
apptainer exec --bind /beegfs:/beegfs --env TMPDIR=/tmp "$MOSDEPTH_SIF" \
    mosdepth --by "$ROI_CATALOG" --no-per-base --threads 4 "$(basename "$COV_PREFIX")" "$BAM"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] mosdepth done, outputs:"
ls -lh "${COV_PREFIX}"* 2>/dev/null | head

# --- Step 4: Aggregate stats per ROI class ---
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Aggregating per-class stats"
export COV_PREFIX ROI_CATALOG
python3 << 'PYEOF'
import gzip, json, statistics
from collections import defaultdict

cov_file = f"{__import__('os').environ.get('COV_PREFIX','HG002_coverage')}.regions.bed.gz"
cat_file = f"{__import__('os').environ.get('ROI_CATALOG','roi_catalog.bed')}"

# Build catalog lookup: (chrom,start,end) -> class from the merged distinct class (col 4)
import os
ROI_CATALOG = os.environ['ROI_CATALOG']
COV_PREFIX  = os.environ['COV_PREFIX']

lookup = {}
with open(ROI_CATALOG) as f:
    for line in f:
        p = line.rstrip().split('\t')
        if len(p) >= 4:
            lookup[(p[0], int(p[1]), int(p[2]))] = p[3]

per_class = defaultdict(list)
total_regions = 0
with gzip.open(f"{COV_PREFIX}.regions.bed.gz", "rt") as f:
    for line in f:
        p = line.rstrip().split('\t')
        if len(p) < 4: continue
        chrom, start, end = p[0], int(p[1]), int(p[2])
        mean_cov = float(p[-1])
        cls = lookup.get((chrom,start,end), "unknown")
        per_class[cls].append(mean_cov)
        total_regions += 1

stats = {"total_regions": total_regions, "by_class": {}}
for cls, covs in per_class.items():
    n = len(covs)
    stats["by_class"][cls] = {
        "n_regions": n,
        "mean_coverage": round(statistics.mean(covs), 3) if n else 0,
        "median_coverage": round(statistics.median(covs), 3) if n else 0,
        "stdev_coverage": round(statistics.stdev(covs), 3) if n > 1 else 0,
        "min_coverage": round(min(covs), 3) if n else 0,
        "max_coverage": round(max(covs), 3) if n else 0,
    }

out = f"{os.path.dirname(COV_PREFIX)}/roi_stats.json"
with open(out, "w") as f:
    json.dump(stats, f, indent=2)
print(f"Stats written: {out}")
print(json.dumps(stats, indent=2))
PYEOF

echo "[$(date '+%Y-%m-%d %H:%M:%S')] T6+T7 COMPLETE"
ls -lh "$ROI_DIR"
