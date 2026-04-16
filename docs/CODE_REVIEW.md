# BRANCH v0.1 — Code Review (Consolidated)

**Datum**: 2026-04-16
**Reviewer**: Zyrkel (service-zyrkel, 4 gesplittete Calls)
**Scope**: vollstaendiger C++20 Core + Tests + Snakemake-Workflow
**Commits reviewed**: `26b1727` (HEAD), `cab0ea8`, `745fbd9`

Quellen:
- [CODE_REVIEW_graph.md](CODE_REVIEW_graph.md)
- [CODE_REVIEW_classify.md](CODE_REVIEW_classify.md)
- [CODE_REVIEW_backend_gpu.md](CODE_REVIEW_backend_gpu.md)
- [CODE_REVIEW_tests_workflow.md](CODE_REVIEW_tests_workflow.md)

---

## Executive Summary

Das v0.1-Skelett steht architektonisch solide. Kern-Designprinzipien (lossless, CN-aware Nodes, VAF-tagged Edges, Backend-Vtable fuer CPU/GPU-Swap) sind korrekt umgesetzt. Zwei Ebenen-Gaps sind vor Phase-0-Run zu adressieren: **Memory-Safety** (keine Bounds-Checks im Graph-Layer) und **Classification-Logic** (Mixed-Klasse wird nie emittiert, Confidence nicht kalibriert). CUDA-Kernel ist bewusst als No-op-Skeleton angelegt — das ist OK fuer v0.1, aber Performance-Bugs im Launcher (per-call `cudaMalloc`, kein Stream) sollten vor v0.2 gefixt werden.

---

## Critical Issues (Priorisiert)

### HIGH — muss vor Phase-0 fixed werden
| # | Layer | Datei:Zeile | Problem |
|---|-------|-------------|---------|
| 1 | graph | [lossless_graph.hpp:70-71](../src/graph/lossless_graph.hpp#L70-L71) | `node(NodeId)` ohne Bounds-Check → UB |
| 2 | graph | [lossless_graph.hpp:83-85](../src/graph/lossless_graph.hpp#L83-L85) | `set_edge_vaf()` ohne Bounds-Check |
| 3 | graph | [reverse_index.hpp:33-37](../src/graph/reverse_index.hpp#L33-L37) | `add()` nach `freeze()` korrumpiert State silently |
| 4 | classify | [cascade.hpp:78-80](../src/classify/cascade.hpp#L78-L80) | Stage 2 setzt `confidence=1.0f` hart, Confidence unkalibriert |
| 5 | workflow | [Snakefile:~roi_catalog](../workflow/Snakefile) | `/tmp/segdup_merged.bed` → Daten nicht auf BeeGFS |

### MEDIUM — vor v0.2
| # | Layer | Problem |
|---|-------|---------|
| 6 | classify | `Mixed`-Klasse wird nie emittiert — rekursive Dekomposition fehlt komplett |
| 7 | classify | DepthRatio-Threshold 1.8 biologisch zu niedrig (Duplikation ~2.0x) |
| 8 | gpu | Per-call `cudaMalloc` in [overlap_kernel.cu:24-34](../src/gpu/overlap_kernel.cu#L24-L34) killt Batch-Perf |
| 9 | gpu | Keine CUDA-Error-Checks (`cudaGetLastError()` nach Launch) |
| 10 | gpu | Kein `cudaStream_t`-Parameter — blockiert Multi-GPU-Concurrency auf Hummel-2 H100s |
| 11 | tests | GPU-Kernel-Logik ungetestet (nur Stub) |
| 12 | ci | Keine `.github/workflows/` — ASan/TSan laufen nirgends automatisch |

### LOW
- 6 von 12 Features werden nicht in Cascade genutzt (Feature 3-7, 10-11) — fuer LightGBM in v0.3 dokumentieren
- `CascadeConfig::min_bubble_length_bp` / `min_coverage` / `max_recursion_depth` sind Dead-Config
- Hardcoded Username `bbg6775` in Snakefile (`SHARED`-Pfad)
- Node-Sequence-Storage fehlt (nur `length_bp`, keine tatsaechlichen Basen)
- Adjacency-Index (`edges_from_/edges_to_`) nur im Kommentar erwaehnt

---

## Design Feedback

**Gut:**
- Lossless-Property gewahrt: `PositionalDelta` (8B packed) + Read-as-Path, keine Daten-Collapse-Pfade
- CN-aware als Attribut (nicht Node-Duplikation) — korrekt
- VAF + `vaf_confidence` per Edge — QUASR-Integration vorbereitet
- Backend-VTable sauber typerased, move-only, CPU↔GPU-Swap ohne Recompile
- `.cuh` / `.cu`-Trennung korrekt: Header ohne `__global__` compilebar mit g++
- Stub-Fallback fuer `BRANCH_BUILD_CUDA=OFF` funktioniert
- Cascade-Pattern mit Early-Exit, Feature-Enum stabil fuer ML-Input-Index
- ASan/TSan/UBSan-CMake-Optionen vorhanden

**Fehlend:**
- Mixed-Klasse-Dekomposition (im Header-Kommentar beschrieben, nicht implementiert)
- Memory-Pool fuer GPU-Workspace (ersetzt per-call cudaMalloc)
- Fat-Binary-Runtime-Device-Selection (`cudaGetDeviceCount`/`cudaSetDevice`)
- Integration-Tests (End-to-End BAM → Graph → Classify)
- Fuzz-/Property-Tests (rapidcheck)
- Bounds-Error-Death-Tests (`EXPECT_DEATH` mit ASan)

---

## Phase-0-Readiness Checklist

Vor empirischem Baseline-Run (1 Genom, hifiasm auf Hummel-2) fixen:

- [ ] Bounds-Checks in `node()`, `set_edge_vaf()`, `add()` post-freeze
- [ ] `/tmp/segdup_merged.bed` → `{PROJECT_DIR}/tmp/` in Snakefile
- [ ] `.github/workflows/ci.yml` mit ASan/TSan-Matrix (CPU-only Target)
- [ ] Hardcoded `SHARED`-Pfad via `config.yaml` parametrisieren
- [ ] `min_bubble_length_bp` / `min_coverage` in `classify_one()` tatsaechlich pruefen

Non-blocker fuer Phase 0 (Paper-Review laut Projekt-Regel nicht Ziel):
- CUDA-Perf-Fixes (No-op-Kernel v0.1)
- Mixed-Klasse-Rekursion (nur Branch/Dup reicht fuer Baseline-Run)
- Integration-Tests (kann nach Phase-0 kommen)

---

## Konkrete Naechste Schritte (empfohlene Reihenfolge)

1. **PR 1**: HIGH-Issues 1-3 fixen + Death-Tests → Graph-Layer memory-safe
2. **PR 2**: Snakefile `/tmp`-Bug + `SHARED`-Parametrisierung → Hummel-2-ready
3. **PR 3**: `.github/workflows/ci.yml` mit ASan/TSan → CI aktiv
4. **Phase 0 Baseline-Run** auf Hummel-2
5. **PR 4** (nach Phase 0): Classification-Kalibrierung + Mixed-Klasse
6. **PR 5** (v0.2): CUDA-Memory-Pool + Stream-API + Multi-GPU-Dispatch
