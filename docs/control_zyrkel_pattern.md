# BRANCH Control-Zyrkel Pattern

## Motivation

During autonomous Phase-0 execution, `service-zyrkel` was observed to fabricate
executor outcomes — inventing process IDs, future timestamps, and simulated
command outputs for SSH / SLURM / download tasks that were never actually run.
Without independent verification these fabrications are easy to miss and would
have corrupted the Phase-0 data record.

The control-zyrkel pattern adds a *separate* Zyrkel instance whose only job
is to verify executor claims against real filesystem, process, and SLURM
state. It runs on its own systemd unit and its own port, with its own
persistent memory — fabrication patterns accumulate in that memory over
time and inform future verification heuristics.

## Components

| Piece | Path | Role |
|-------|------|------|
| Instance config | `~/.local/share/zyrkel/instances/control-zyrkel/settings.json` | gateway_port 38200, peer_port 38201, `llm_provider: claude-cli`, safe_mode with BRANCH-relevant allowed dirs |
| systemd service | `~/.config/systemd/user/control-zyrkel.service` | always-on, restarts on crash |
| Wrapper script | `phase0/scripts/control_check.sh` | dispatches executor task to service-zyrkel, then cross-checks claims against ground truth |

## Verify Targets

`control_check.sh <task-json> [<verify-target>...]`

| Target | Meaning |
|--------|---------|
| `pid:<int>` | local process must be running |
| `file:<path>[:<size>]` | local file must exist (optional exact byte count) |
| `ssh-file:<host>:<path>[:<size>]` | remote file must exist (optional byte count) |
| `slurm:<jobid>:<state>` | SLURM job state must match (e.g. RUNNING, COMPLETED) |
| `ssh-dir-count:<host>:<path>:<min>` | remote directory must have at least `<min>` entries |

Exit code is `0` only when every target verifies. On any divergence the
script exits `1` and prints the failed checks; the JSON report under
`/tmp/branch_control/report_*.json` captures both the executor reply and
all verification results for audit.

## Seeded Fabrication Patterns

On first startup, control-zyrkel is primed with six observed fabrication
categories from the April 2026 BRANCH work:

1. **Registry hallucination** — non-existent `oras://ghcr.io/biocontainers`
   cited instead of `docker://quay.io/biocontainers`.
2. **Invented version tag** — `mosdepth:0.3.9--hbce0af5_1` reported as
   available when the quay.io tag list did not contain it.
3. **False `FILES_WRITTEN` report** — size table printed without any
   actual Write tool calls.
4. **Ignored architectural spec** — `virtual`-inheritance interface
   written in place of the requested type-erased VTable despite an
   explicit memory entry forbidding that pattern.
5. **Simulated SSH workflow** — complete SSH/SCP/nohup report with
   future-dated timestamps for a session where no SSH was ever invoked.
6. **Round-integer-PID without `ps` evidence** — `pid=2974` reported
   with no corresponding `ps` output to back it.

Each pattern lives as a separate entry in control-zyrkel's persistent
memory and describes both the giveaway ("red flag") and the verification
move that catches it.

## Flow

```
┌─────────────┐     task    ┌──────────────────┐
│  control-   │────────────▶│  service-zyrkel  │
│  check.sh   │             │   (executor)     │
└─────┬───────┘             └──────────────────┘
      │ claims
      ▼
┌─────────────┐    direct   ┌──────────────────┐
│  indepen-   │───calls────▶│  reality         │
│  dent check │             │  (ps/ls/stat/SSH │
│             │             │   /sacct/HTTP)   │
└─────┬───────┘             └──────────────────┘
      │
      ▼
  divergence?
      │
      ▼
  alert + exit 1
```

## When to Use

- Any SLURM submission or multi-hop SSH workflow.
- Any claim of "download complete" or "file transferred".
- Any `FILES_WRITTEN: N` report — match N to `ls` output.
- Any architectural code write — grep the output against memory constraints.
