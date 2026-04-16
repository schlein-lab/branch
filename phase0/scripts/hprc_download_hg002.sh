#!/usr/bin/env bash
# BRANCH Phase-0 post-scaffold: download full HPRC HG002 HiFi raw data.
#
# Targets the public HPRC S3 bucket via its HTTPS endpoint so no awscli
# is required. Lists the raw_data/PacBio_HiFi/ folder via the S3 XML
# index and wgets every .bam (+ .bai where present).
#
# Runs idempotently: existing files with matching size are skipped.
# Safe on Hummel-2 frontnode (long-running background job pattern).

set -euo pipefail

USER_ID="${USER:-bbg6775}"
BEEGFS_BASE="/beegfs/u/${USER_ID}/humangenetik/shared/references/human_pangenome/hifi_reads/HG002"
TARGET_DIR="${BEEGFS_BASE}/raw_data_hprc"
LOG_DIR="${BEEGFS_BASE}"
LOG_FILE="${LOG_DIR}/hprc_download_$(date +%Y%m%d_%H%M%S).log"

BUCKET_BASE="https://human-pangenomics.s3.amazonaws.com"
PREFIX="working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/"

mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

echo "[$(date '+%F %T')] HPRC HG002 HiFi download starting"     | tee -a "$LOG_FILE"
echo "[$(date '+%F %T')] target=$TARGET_DIR"                    | tee -a "$LOG_FILE"

# Fetch S3 XML index recursively — S3 returns 1000 keys per page, so paginate
# with continuation-token. We iterate across the full prefix tree (no
# delimiter) and filter downsampled/ later.
INDEX_XML="${TARGET_DIR}/.s3_index.xml"
: > "$INDEX_XML"
CONT_TOKEN=""
PAGE=0
while : ; do
    PAGE=$((PAGE + 1))
    if [[ -z "$CONT_TOKEN" ]]; then
        URL="${BUCKET_BASE}/?prefix=${PREFIX}&list-type=2"
    else
        # URL-encode the token (it can contain /,+,=)
        TOKEN_ENC=$(python3 -c "import urllib.parse,sys; print(urllib.parse.quote(sys.argv[1], safe=''))" "$CONT_TOKEN")
        URL="${BUCKET_BASE}/?prefix=${PREFIX}&list-type=2&continuation-token=${TOKEN_ENC}"
    fi
    wget -q -O "${INDEX_XML}.page${PAGE}" "$URL"
    cat "${INDEX_XML}.page${PAGE}" >> "$INDEX_XML"
    NEXT=$(grep -oP '(?<=<NextContinuationToken>)[^<]+' "${INDEX_XML}.page${PAGE}" | head -1 || true)
    TRUNC=$(grep -oP '(?<=<IsTruncated>)[^<]+' "${INDEX_XML}.page${PAGE}" | head -1 || true)
    if [[ "$TRUNC" != "true" || -z "$NEXT" ]]; then
        break
    fi
    CONT_TOKEN="$NEXT"
done

# Extract keys + sizes; keep .bam and .pbi (HiFi index). Skip dirs and
# any downsampled copies (we want the full original reads).
python3 - "$INDEX_XML" > "${TARGET_DIR}/.to_fetch.tsv" <<'PYEOF'
import re, sys
path = sys.argv[1]
with open(path) as f:
    xml = f.read()
entries = re.findall(r"<Contents>(.*?)</Contents>", xml, flags=re.DOTALL)
for e in entries:
    k = re.search(r"<Key>(.*?)</Key>", e)
    s = re.search(r"<Size>(.*?)</Size>", e)
    if not (k and s):
        continue
    key, size = k.group(1), int(s.group(1))
    if size == 0:
        continue
    if "/downsampled/" in key:
        continue
    if key.endswith((".bam", ".bam.pbi", ".bam.bai",
                     ".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        print(f"{key}\t{size}")
PYEOF

N_FILES=$(wc -l < "${TARGET_DIR}/.to_fetch.tsv")
echo "[$(date '+%F %T')] index listed ${N_FILES} candidate files" | tee -a "$LOG_FILE"

# Pull each file. Idempotent: skip if size matches.
# Preserve the relative path under PREFIX to keep insert-size folders
# separate (15kb/, 19kb/, 20kb/, 25kb/, UW_supplement/).
while IFS=$'\t' read -r KEY EXPECTED_SIZE; do
    REL="${KEY#${PREFIX}}"
    LOCAL_PATH="${TARGET_DIR}/${REL}"
    mkdir -p "$(dirname "$LOCAL_PATH")"
    if [[ -f "$LOCAL_PATH" ]]; then
        ACTUAL=$(stat -c %s "$LOCAL_PATH")
        if [[ "$ACTUAL" == "$EXPECTED_SIZE" ]]; then
            echo "[$(date '+%F %T')] SKIP ${REL} (already present, size matches)" \
                | tee -a "$LOG_FILE"
            continue
        else
            echo "[$(date '+%F %T')] RESUME ${REL} (${ACTUAL}/${EXPECTED_SIZE})" \
                | tee -a "$LOG_FILE"
        fi
    fi
    URL="${BUCKET_BASE}/${KEY}"
    echo "[$(date '+%F %T')] GET ${URL}" | tee -a "$LOG_FILE"
    wget -q --continue --tries=5 --timeout=60 -O "$LOCAL_PATH" "$URL" \
        || echo "[$(date '+%F %T')] ERROR ${REL} wget failed" | tee -a "$LOG_FILE"
    if [[ -f "$LOCAL_PATH" ]]; then
        ACTUAL=$(stat -c %s "$LOCAL_PATH")
        echo "[$(date '+%F %T')] DONE ${REL} size=${ACTUAL}" | tee -a "$LOG_FILE"
    fi
done < "${TARGET_DIR}/.to_fetch.tsv"

echo "[$(date '+%F %T')] HPRC HG002 download complete" | tee -a "$LOG_FILE"
ls -lh "$TARGET_DIR" | tail -30 | tee -a "$LOG_FILE"
