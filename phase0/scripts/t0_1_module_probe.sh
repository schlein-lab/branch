#!/usr/bin/env bash
set -u

# T0.1 Module Probe for Hummel-2
# Tries multiple module init paths and collects available modules as JSON

json_escape() {
    python3 -c "import json,sys; print(json.dumps(sys.stdin.read()))"
}

section_output() {
    local name="$1"
    local content="$2"
    local status="${3:-ok}"
    local escaped
    escaped=$(echo -n "$content" | json_escape)
    printf '"%s": {"status": "%s", "data": %s}' "$name" "$status" "$escaped"
}

try_module_init() {
    local init_path="$1"
    local init_name="$2"

    if [[ -f "$init_path" ]]; then
        # Subshell to isolate module environment
        (
            source "$init_path" 2>/dev/null
            if type module &>/dev/null; then
                echo "INIT_SUCCESS"
                echo "---AVAIL_START---"
                module avail 2>&1 | head -100
                echo "---AVAIL_END---"
                echo "---TERSE_START---"
                module --terse avail 2>&1 | head -200
                echo "---TERSE_END---"
                echo "---GCC_START---"
                module avail gcc 2>&1 | head -50
                echo "---GCC_END---"
                echo "---CUDA_START---"
                module avail cuda 2>&1 | head -50
                echo "---CUDA_END---"
            else
                echo "INIT_FAIL: module command not available after sourcing"
            fi
        )
    else
        echo "INIT_FAIL: $init_path does not exist"
    fi
}

echo "{"
printf '"timestamp": "%s",\n' "$(date -Iseconds)"
printf '"probe_version": "1.1",\n'
printf '"probe_type": "module_init",\n'

# Init path A: /etc/profile.d/modules.sh
init_a_out=$(try_module_init "/etc/profile.d/modules.sh" "profile.d" 2>&1) || init_a_out="ERROR: $init_a_out"
section_output "init_profile_d" "$init_a_out"
echo ","

# Init path B: /sw/batch/init.sh (Hummel-2 specific)
init_b_out=$(try_module_init "/sw/batch/init.sh" "sw_batch" 2>&1) || init_b_out="ERROR: $init_b_out"
section_output "init_sw_batch" "$init_b_out"
echo ","

# Init path C: /usr/share/Modules/init/bash
init_c_out=$(try_module_init "/usr/share/Modules/init/bash" "usr_share" 2>&1) || init_c_out="ERROR: $init_c_out"
section_output "init_usr_share" "$init_c_out"
echo ","

# Bonus: Check for Lmod specifically
lmod_out=""
if [[ -f "/usr/share/lmod/lmod/init/bash" ]]; then
    lmod_out=$(try_module_init "/usr/share/lmod/lmod/init/bash" "lmod" 2>&1) || lmod_out="ERROR: $lmod_out"
elif [[ -f "/opt/apps/lmod/lmod/init/bash" ]]; then
    lmod_out=$(try_module_init "/opt/apps/lmod/lmod/init/bash" "lmod" 2>&1) || lmod_out="ERROR: $lmod_out"
else
    lmod_out="INIT_FAIL: No Lmod init found at standard paths"
fi
section_output "init_lmod" "$lmod_out"
echo ","

# Environment info for debugging
env_info="MODULEPATH=${MODULEPATH:-unset}
LMOD_DIR=${LMOD_DIR:-unset}
MODULESHOME=${MODULESHOME:-unset}"
section_output "module_env_vars" "$env_info"
echo ""

echo "}"
