#!/usr/bin/env bash
# ============================================================================
# T1: Container Pull Script for BRANCH Phase 0
# Pulls 5 utility containers (NOT assembly tools) to BeeGFS via /tmp staging
# ============================================================================
set -uo pipefail

# --- Paths (dynamic, based on $USER) ---
CONTAINER_BASE="/beegfs/u/${USER}/humangenetik/ag/schlein/${USER}/containers"
TMP_STAGE="/tmp/branch_container_pull_$$"
MANIFEST_FILE="${CONTAINER_BASE}/manifest.json"
REGISTRY="docker://quay.io/biocontainers"

# --- Tool definitions: name|version|tag ---
declare -a TOOLS=(
    "minimap2|2.28|2.28--he4a0461_0"
    "samtools|1.21|1.21--h50ea8bc_0"
    "mosdepth|0.3.9|0.3.9--hbce0af5_1"
    "bedtools|2.31.1|2.31.1--h13024bc_3"
    "seqkit|2.8.2|2.8.2--h9ee0642_0"
)

# --- Functions ---
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

cleanup() {
    [[ -d "$TMP_STAGE" ]] && rm -rf "$TMP_STAGE"
}
trap cleanup EXIT

# --- Setup ---
mkdir -p "$CONTAINER_BASE" "$TMP_STAGE"

# Initialize manifest array
declare -a MANIFEST_ENTRIES=()

log "Starting container pull to: $CONTAINER_BASE"
log "Staging directory: $TMP_STAGE"

# --- Pull loop ---
for tool_spec in "${TOOLS[@]}"; do
    IFS='|' read -r NAME VERSION TAG <<< "$tool_spec"
    SIF_NAME="${NAME}_${VERSION}.sif"
    SIF_PATH="${CONTAINER_BASE}/${SIF_NAME}"
    PULL_URI="${REGISTRY}/${NAME}:${TAG}"
    
    # Check idempotency
    if [[ -f "$SIF_PATH" ]]; then
        log "SKIP: $SIF_NAME already exists"
        EXISTING_SHA=$(sha256sum "$SIF_PATH" | awk '{print $1}')
        MANIFEST_ENTRIES+=("{\"tool\":\"$NAME\",\"version\":\"$VERSION\",\"tag\":\"$TAG\",\"sha256\":\"$EXISTING_SHA\",\"path\":\"$SIF_PATH\",\"status\":\"exists\"}")
        continue
    fi
    
    log "PULL: $NAME v$VERSION from $PULL_URI"
    
    TMP_SIF="${TMP_STAGE}/${SIF_NAME}"
    PULL_STDERR=""
    PULL_STATUS="success"
    
    # Pull to /tmp first (BeeGFS hates partial writes)
    if ! PULL_STDERR=$(apptainer pull --force "$TMP_SIF" "$PULL_URI" 2>&1); then
        err "Failed to pull $NAME: $PULL_STDERR"
        MANIFEST_ENTRIES+=("{\"tool\":\"$NAME\",\"version\":\"$VERSION\",\"tag\":\"$TAG\",\"sha256\":null,\"path\":null,\"status\":\"error\",\"error\":\"$PULL_STDERR\"}")
        continue
    fi
    
    # Compute SHA256
    SHA256=$(sha256sum "$TMP_SIF" | awk '{print $1}')
    log "SHA256: $SHA256"
    
    # Move to final location
    mv "$TMP_SIF" "$SIF_PATH"
    log "DONE: $SIF_NAME -> $SIF_PATH"
    
    MANIFEST_ENTRIES+=("{\"tool\":\"$NAME\",\"version\":\"$VERSION\",\"tag\":\"$TAG\",\"sha256\":\"$SHA256\",\"path\":\"$SIF_PATH\",\"status\":\"$PULL_STATUS\"}")
done

# --- Write JSON manifest ---
{
    echo "{"
    echo "  \"generated\": \"$(date -Iseconds)\","
    echo "  \"user\": \"$USER\","
    echo "  \"base_path\": \"$CONTAINER_BASE\","
    echo "  \"containers\": ["
    for i in "${!MANIFEST_ENTRIES[@]}"; do
        if [[ $i -lt $((${#MANIFEST_ENTRIES[@]} - 1)) ]]; then
            echo "    ${MANIFEST_ENTRIES[$i]},"
        else
            echo "    ${MANIFEST_ENTRIES[$i]}"
        fi
    done
    echo "  ]"
    echo "}"
} > "$MANIFEST_FILE"

log "Manifest written to: $MANIFEST_FILE"
log "Container pull complete."
