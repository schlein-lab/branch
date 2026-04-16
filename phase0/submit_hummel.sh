#!/usr/bin/env bash
# =============================================================================
# BRANCH v0.2 - Hummel-2 Snakemake submit helper
# =============================================================================
# Launches the BRANCH workflow on Hummel-2's SLURM fabric. This script is the
# canonical entry point for per-sample runs; it MUST be executed from a login
# node of the Hummel-2 cluster (not from your laptop, not from Hummel-1).
#
# Usage:
#   ./phase0/submit_hummel.sh                      # HG002 default
#   ./phase0/submit_hummel.sh --config sample=NA24385
#   ./phase0/submit_hummel.sh --dry-run            # show DAG without submitting
#
# All extra arguments are forwarded verbatim to snakemake.
#
# Per Hummel-2 rules:
#   - no --mem anywhere (RRZ rejects it)
#   - `set +u` around /sw/batch/init.sh (module init uses unbound vars)
#   - TMPDIR=/tmp (home is read-only inside batch jobs)
#   - apptainer --bind /beegfs so the SIFs can see project data
#   - SLURM executor plugin via workflow/config/slurm/config.yaml
# =============================================================================

set -euo pipefail

# --- BEFORE sourcing RRZ init: loosen unbound-var checking ------------------
# /sw/batch/init.sh and the module system touch $PS1, $BASH_ENV, $CC etc. with
# no default, which would trip `set -u`. We re-enable it after module load.
set +u

# Fast-fail if this script is running anywhere except Hummel-2.
# Login nodes match 'hummel*' or 'rrz*'; fall back to detecting BeeGFS.
HOSTNAME_SHORT="$(hostname -s 2>/dev/null || echo unknown)"
if [[ ! "${HOSTNAME_SHORT}" =~ ^(hummel|rrz|node|compute) ]] && [[ ! -d /beegfs/u ]]; then
    echo "ERROR: submit_hummel.sh must run on Hummel-2." >&2
    echo "  detected hostname: ${HOSTNAME_SHORT}" >&2
    echo "  /beegfs/u missing — this looks like a non-cluster host." >&2
    echo "  SSH to Hummel-2 first, then rerun." >&2
    exit 2
fi

# --- TMPDIR: /tmp (NOT /beegfs; NOT $HOME which is read-only in jobs) --------
export TMPDIR=/tmp

# --- RRZ module init ---------------------------------------------------------
# shellcheck disable=SC1091
source /sw/batch/init.sh
module load snakemake
module load apptainer

# Re-enable strict-unbound after the module machinery is set up.
set -u

# --- Paths (mirror Snakefile defaults) --------------------------------------
USER="${USER:-$(whoami)}"
BEEGFS_BASE="/beegfs/u/${USER}/humangenetik/ag/schlein/${USER}"
PROJECT_DIR="${PROJECT_DIR:-${BEEGFS_BASE}/projects/branch}"
CONTAINER_DIR="${CONTAINER_DIR:-${BEEGFS_BASE}/containers}"
LOG_DIR="${LOG_DIR:-${BEEGFS_BASE}/logs/slurm_jobs}"

mkdir -p "${LOG_DIR}" "${TMPDIR}"

# --- Sanity checks -----------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
CONFIG_FILE="${REPO_ROOT}/workflow/config/config.yaml"
SLURM_PROFILE="${REPO_ROOT}/workflow/config/slurm"

for path in "${SNAKEFILE}" "${CONFIG_FILE}" "${SLURM_PROFILE}/config.yaml"; do
    if [[ ! -e "${path}" ]]; then
        echo "ERROR: missing workflow asset: ${path}" >&2
        exit 3
    fi
done

# BRANCH binary: the workflow checks this too, but failing fast here gives a
# friendlier error than a Snakemake parse AssertionError.
BRANCH_BIN="${PROJECT_DIR}/build/src/cli/branch"
if [[ ! -x "${BRANCH_BIN}" ]]; then
    echo "WARNING: BRANCH binary not executable at ${BRANCH_BIN}" >&2
    echo "  Build it first with cmake -S ${REPO_ROOT} -B ${REPO_ROOT}/build && cmake --build ${REPO_ROOT}/build" >&2
    echo "  Or override via: ${0##*/} --config branch_bin=/path/to/branch" >&2
fi

# --- Snakemake invocation ----------------------------------------------------
# Notes:
#   - --profile points at the SLURM executor config (no --mem there)
#   - --software-deployment-method apptainer (replaces --use-singularity in Snakemake 9.x)
#   - --apptainer-args "--bind /beegfs" so SIFs see BeeGFS data
#   - --jobs 50 is the per-workflow concurrent-SLURM-job cap
#   - --rerun-incomplete + --keep-going: robust against single-job transient failures
cd "${REPO_ROOT}"

exec snakemake \
    -s "${SNAKEFILE}" \
    --profile "${SLURM_PROFILE}" \
    --configfile "${CONFIG_FILE}" \
    --software-deployment-method apptainer \
    --apptainer-args "--bind /beegfs" \
    --jobs 50 \
    --rerun-incomplete \
    --keep-going \
    "$@"
