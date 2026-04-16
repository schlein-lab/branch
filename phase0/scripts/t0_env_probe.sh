#!/usr/bin/env bash
set -u

# BRANCH Phase 0 Task T0 — Environment Probe for Hummel-2 Cluster
# Output: Structured JSON to stdout
# Safe: Read-only, no modifications, idempotent

export PATH=/syssw/slurm/current/bin:$PATH

# Helper: JSON-escape strings
json_escape() {
    printf '%s' "$1" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))'
}

# Helper: Run command with error capture
run_section() {
    local name="$1"
    shift
    local output
    local status
    output=$("$@" 2>&1) && status=0 || status=$?
    if [[ $status -eq 0 ]]; then
        printf '"%s": {"status": "ok", "data": %s}' "$name" "$(json_escape "$output")"
    else
        printf '"%s": {"status": "error", "exit_code": %d, "data": %s}' "$name" "$status" "$(json_escape "$output")"
    fi
}

# Start JSON output
printf '{\n'
printf '"probe_timestamp": "%s",\n' "$(date -Iseconds)"
printf '"probe_version": "t0_env_probe_v1",\n'

# 1. Hostname and uname
printf '%s,\n' "$(run_section "hostname" hostname)"
printf '%s,\n' "$(run_section "uname" uname -a)"

# 2. SLURM partitions
printf '%s,\n' "$(run_section "slurm_partitions" sinfo --format="%P %c %m %G %l")"

# 3. Account memberships
printf '%s,\n' "$(run_section "slurm_accounts" sacctmgr show assoc user=$USER format=account,partition --noheader --parsable2)"

# 4. Module availability (head -100)
module_avail_output=$(module avail 2>&1 | head -100) || module_avail_output="ERROR: module avail failed"
printf '"module_avail": {"status": "ok", "data": %s},\n' "$(json_escape "$module_avail_output")"

# 5. GCC versions (via module avail, not spider - spider needs initialized module system)
gcc_avail_output=$(module avail 2>&1 | grep -iE "gcc" | head -40) || gcc_avail_output="ERROR: module avail gcc failed"
printf '"gcc_versions": {"status": "ok", "data": %s},\n' "$(json_escape "$gcc_avail_output")"

# 6. CUDA versions (via module avail)
cuda_avail_output=$(module avail 2>&1 | grep -iE "cuda" | head -40) || cuda_avail_output="ERROR: module avail cuda failed"
printf '"cuda_versions": {"status": "ok", "data": %s},\n' "$(json_escape "$cuda_avail_output")"

# 7. BeeGFS quota and df
quota_output=$(quota -s 2>/dev/null) || quota_output="quota command not available or no quota set"
df_beegfs_output=$(df -h /beegfs 2>&1) || df_beegfs_output="ERROR: df /beegfs failed"
printf '"quota": {"status": "ok", "data": %s},\n' "$(json_escape "$quota_output")"
printf '"df_beegfs": {"status": "ok", "data": %s},\n' "$(json_escape "$df_beegfs_output")"

# 8. Hardware info
printf '%s,\n' "$(run_section "nproc" nproc)"
printf '%s,\n' "$(run_section "free" free -h)"
lscpu_output=$(lscpu | head -20) || lscpu_output="ERROR: lscpu failed"
printf '"lscpu": {"status": "ok", "data": %s},\n' "$(json_escape "$lscpu_output")"

# 9. GPU partition state (non-blocking, from frontnode — no srun job submission in T0)
gpu_partition_output=$(sinfo -p gpu --format="%P %c %m %G %t %D" 2>&1) || gpu_partition_output="ERROR: sinfo -p gpu failed"
printf '"gpu_partition": {"status": "ok", "data": %s}\n' "$(json_escape "$gpu_partition_output")"

printf '}\n'
