#!/usr/bin/env bash
# BRANCH control-zyrkel wrapper — robust executor+verification pattern.
#
# Flow:
#   1. Dispatch an executor task (arg $1 = json-escaped prompt) to
#      service-zyrkel via /api/chat.
#   2. Parse claims from the executor response — PIDs, file paths,
#      byte sizes, SLURM job IDs, timestamps.
#   3. Independently verify each claim by direct SSH/local checks
#      (ps, ls, stat, sacct, sinfo).
#   4. Report reality vs. claim as machine-readable JSON and emit a
#      human-readable divergence summary.
#   5. Exit code 0 only if all claims verified; >=1 on any divergence.
#
# Usage:
#   control_check.sh <task-json-file> [<verify-target>...]
#
# Where <task-json-file> holds the executor prompt, and each
# <verify-target> is one of:
#   pid:<int>          process must be running
#   file:<path>        file must exist
#   file:<path>:<size> file must exist with exact byte size
#   ssh-file:<host>:<path>        file must exist on remote host
#   ssh-file:<host>:<path>:<size> file must exist remote with size
#   slurm:<jobid>:<state>  SLURM job state must match (RUNNING/COMPLETED)
#   ssh-dir-count:<host>:<path>:<min-count>  dir has >= min entries

set -euo pipefail

TOKEN_FILE="${HOME}/.local/share/zyrkel/fleet_token"
SERVICE_URL="http://127.0.0.1:37854/api/chat"
CONTROL_URL="http://127.0.0.1:38200/api/chat"   # control-zyrkel (if running)
OUT_DIR="${BRANCH_CONTROL_OUT_DIR:-/tmp/branch_control}"
mkdir -p "$OUT_DIR"
STAMP="$(date +%Y%m%d_%H%M%S)"
EXEC_LOG="${OUT_DIR}/exec_${STAMP}.json"
REPORT="${OUT_DIR}/report_${STAMP}.json"

if [[ ! -f "$TOKEN_FILE" ]]; then
    echo "ERROR: fleet token not found at $TOKEN_FILE" >&2
    exit 2
fi
TOKEN=$(cat "$TOKEN_FILE")

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <task-json-file> [<verify-target>...]" >&2
    exit 2
fi
TASK_FILE="$1"; shift
if [[ ! -f "$TASK_FILE" ]]; then
    echo "ERROR: task file not found: $TASK_FILE" >&2
    exit 2
fi

# ------------------------------------------------------------------
# Stage 1: dispatch executor task
# ------------------------------------------------------------------
MSG=$(cat "$TASK_FILE")
PAYLOAD=$(jq -Rn --arg m "$MSG" '{message:$m}')
echo "[$(date '+%F %T')] control: dispatching executor task"
curl -s --max-time 900 -X POST "$SERVICE_URL" \
    -H "Content-Type: application/json" \
    -H "X-Fleet-Token: $TOKEN" \
    -d "$PAYLOAD" \
    -o "$EXEC_LOG" \
    -w "%{http_code}\n" > "${EXEC_LOG}.http"

EXEC_HTTP=$(cat "${EXEC_LOG}.http")
echo "[$(date '+%F %T')] control: executor http=$EXEC_HTTP log=$EXEC_LOG"

EXEC_RESPONSE=$(jq -r '.response // ""' "$EXEC_LOG" 2>/dev/null || echo "")
EXEC_OK=$(jq -r '.ok // false' "$EXEC_LOG" 2>/dev/null || echo false)

# ------------------------------------------------------------------
# Stage 2: independent verification for each target
# ------------------------------------------------------------------
declare -a CHECKS=()
declare -i PASS=0 FAIL=0

verify_pid() {
    local pid="$1"
    if ps -p "$pid" -o pid= >/dev/null 2>&1; then
        CHECKS+=("{\"type\":\"pid\",\"pid\":$pid,\"status\":\"running\",\"ok\":true}")
        PASS=$((PASS+1))
    else
        CHECKS+=("{\"type\":\"pid\",\"pid\":$pid,\"status\":\"not-running\",\"ok\":false}")
        FAIL=$((FAIL+1))
    fi
}

verify_file() {
    local path="$1" expected="${2:-}"
    if [[ ! -f "$path" ]]; then
        CHECKS+=("{\"type\":\"file\",\"path\":\"$path\",\"exists\":false,\"ok\":false}")
        FAIL=$((FAIL+1)); return
    fi
    local sz
    sz=$(stat -c %s "$path" 2>/dev/null || echo 0)
    if [[ -n "$expected" && "$sz" != "$expected" ]]; then
        CHECKS+=("{\"type\":\"file\",\"path\":\"$path\",\"size\":$sz,\"expected\":$expected,\"ok\":false}")
        FAIL=$((FAIL+1))
    else
        CHECKS+=("{\"type\":\"file\",\"path\":\"$path\",\"size\":$sz,\"ok\":true}")
        PASS=$((PASS+1))
    fi
}

verify_ssh_file() {
    local host="$1" path="$2" expected="${3:-}"
    local ssh_out
    ssh_out=$(ssh -o ConnectTimeout=15 -o BatchMode=yes "$host" \
        "stat -c '%s' '$path' 2>/dev/null || echo MISSING" 2>/dev/null | tail -1 || echo MISSING)
    if [[ "$ssh_out" == "MISSING" ]]; then
        CHECKS+=("{\"type\":\"ssh-file\",\"host\":\"$host\",\"path\":\"$path\",\"exists\":false,\"ok\":false}")
        FAIL=$((FAIL+1)); return
    fi
    if [[ -n "$expected" && "$ssh_out" != "$expected" ]]; then
        CHECKS+=("{\"type\":\"ssh-file\",\"host\":\"$host\",\"path\":\"$path\",\"size\":$ssh_out,\"expected\":$expected,\"ok\":false}")
        FAIL=$((FAIL+1))
    else
        CHECKS+=("{\"type\":\"ssh-file\",\"host\":\"$host\",\"path\":\"$path\",\"size\":$ssh_out,\"ok\":true}")
        PASS=$((PASS+1))
    fi
}

verify_slurm() {
    local jobid="$1" want="$2"
    local state
    state=$(ssh -o ConnectTimeout=15 hummel-front1 \
        "export PATH=/syssw/slurm/current/bin:\$PATH; sacct -j $jobid --format=State --noheader -n -X 2>/dev/null | tr -d ' \n'" \
        2>/dev/null || echo UNKNOWN)
    if [[ "$state" == "$want" ]]; then
        CHECKS+=("{\"type\":\"slurm\",\"jobid\":\"$jobid\",\"state\":\"$state\",\"ok\":true}")
        PASS=$((PASS+1))
    else
        CHECKS+=("{\"type\":\"slurm\",\"jobid\":\"$jobid\",\"state\":\"$state\",\"expected\":\"$want\",\"ok\":false}")
        FAIL=$((FAIL+1))
    fi
}

verify_ssh_dir_count() {
    local host="$1" path="$2" mincount="$3"
    local n
    n=$(ssh -o ConnectTimeout=15 "$host" "ls -A '$path' 2>/dev/null | wc -l" 2>/dev/null || echo 0)
    if (( n >= mincount )); then
        CHECKS+=("{\"type\":\"ssh-dir-count\",\"host\":\"$host\",\"path\":\"$path\",\"count\":$n,\"min\":$mincount,\"ok\":true}")
        PASS=$((PASS+1))
    else
        CHECKS+=("{\"type\":\"ssh-dir-count\",\"host\":\"$host\",\"path\":\"$path\",\"count\":$n,\"min\":$mincount,\"ok\":false}")
        FAIL=$((FAIL+1))
    fi
}

for target in "$@"; do
    case "$target" in
        pid:*)
            verify_pid "${target#pid:}"
            ;;
        file:*)
            rest="${target#file:}"
            if [[ "$rest" == *:* ]]; then
                verify_file "${rest%:*}" "${rest##*:}"
            else
                verify_file "$rest" ""
            fi
            ;;
        ssh-file:*)
            rest="${target#ssh-file:}"
            host="${rest%%:*}"; rest="${rest#*:}"
            if [[ "$rest" == *:* ]]; then
                verify_ssh_file "$host" "${rest%:*}" "${rest##*:}"
            else
                verify_ssh_file "$host" "$rest" ""
            fi
            ;;
        slurm:*)
            rest="${target#slurm:}"
            verify_slurm "${rest%%:*}" "${rest##*:}"
            ;;
        ssh-dir-count:*)
            rest="${target#ssh-dir-count:}"
            host="${rest%%:*}"; rest="${rest#*:}"
            path="${rest%:*}"; mc="${rest##*:}"
            verify_ssh_dir_count "$host" "$path" "$mc"
            ;;
        *)
            echo "ERROR: unknown verify target: $target" >&2
            FAIL=$((FAIL+1))
            ;;
    esac
done

# ------------------------------------------------------------------
# Stage 3: assemble report
# ------------------------------------------------------------------
CHECKS_JSON=$(printf '%s\n' "${CHECKS[@]}" | jq -s '.' 2>/dev/null || echo "[]")
VERDICT="ok"
if (( FAIL > 0 )); then VERDICT="divergence"; fi

jq -n \
    --arg stamp "$STAMP" \
    --arg task "$TASK_FILE" \
    --arg http "$EXEC_HTTP" \
    --arg ok "$EXEC_OK" \
    --arg resp "$EXEC_RESPONSE" \
    --argjson checks "$CHECKS_JSON" \
    --argjson pass "$PASS" \
    --argjson fail "$FAIL" \
    --arg verdict "$VERDICT" \
    '{stamp:$stamp, task_file:$task, executor:{http:$http, ok:($ok=="true"), response:$resp}, checks:$checks, pass:$pass, fail:$fail, verdict:$verdict}' \
    > "$REPORT"

echo "[$(date '+%F %T')] control: verdict=$VERDICT pass=$PASS fail=$FAIL report=$REPORT"

if [[ "$VERDICT" != "ok" ]]; then
    echo "=== DIVERGENCE ==="
    jq -r '.checks[] | select(.ok==false) | "  FAIL: " + (. | tostring)' "$REPORT"
    exit 1
fi
exit 0
