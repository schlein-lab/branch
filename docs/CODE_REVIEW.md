# BRANCH v0.1 — Code Review (Consolidated)

**Datum**: 2026-04-16
**Reviewer**: Zyrkel (service-zyrkel, 4 gesplittete Calls)
**Scope**: vollstaendiger C++20 Core + Tests + Snakemake-Workflow
**Commits reviewed**: `26b1727` (HEAD), `cab0ea8`, `745fbd9`

**Audit (2026-04-17)**: Status markers added per finding. See per-file CODE_REVIEW_*.md for detail.

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
| # | Layer | Datei:Zeile | Problem | Status |
|---|-------|-------------|---------|--------|
| 1 | graph | [lossless_graph.hpp:70-71](../src/graph/lossless_graph.hpp#L70-L71) | `node(NodeId)` ohne Bounds-Check → UB | ✅ Resolved: now throws `std::out_of_range` (commit d815b80 "Consolidate bounds-checks") |
| 2 | graph | [lossless_graph.hpp:83-85](../src/graph/lossless_graph.hpp#L83-L85) | `set_edge_vaf()` ohne Bounds-Check | ✅ Resolved: throws `std::out_of_range` (d815b80) |
| 3 | graph | [reverse_index.hpp:33-37](../src/graph/reverse_index.hpp#L33-L37) | `add()` nach `freeze()` korrumpiert State silently | ✅ Resolved: throws `std::logic_error("add() called after freeze()")` (d815b80) |
| 4 | classify | [cascade.hpp:78-80](../src/classify/cascade.hpp#L78-L80) | Stage 2 setzt `confidence=1.0f` hart, Confidence unkalibriert | ✅ Resolved: FIX 2 uses `min(1, p2/(threshold*1.5))` |
| 5 | workflow | [Snakefile:~roi_catalog](../workflow/Snakefile) | `/tmp/segdup_merged.bed` → Daten nicht auf BeeGFS | ✅ Resolved: uses `{TMP_DIR}` (Snakefile L111,114) derived from `PROJECT_DIR/tmp` |

### MEDIUM — vor v0.2
| # | Layer | Problem | Status |
|---|-------|---------|--------|
| 6 | classify | `Mixed`-Klasse wird nie emittiert — rekursive Dekomposition fehlt komplett | ⚠️ Partial: `BubbleClass::Mixed` wird jetzt emittiert wenn Flank+Depth widersprechen (cascade.hpp FIX 5); rekursive SNP-Cluster-Dekomposition in `mixed_decomposer.{hpp,cpp}` implementiert, aber nicht in `classify_one()` verdrahtet |
| 7 | classify | DepthRatio-Threshold 1.8 biologisch zu niedrig (Duplikation ~2.0x) | ✅ Resolved: FIX 1 setzt `depth_ratio_dup_threshold{2.0f}` |
| 8 | gpu | Per-call `cudaMalloc` in [overlap_kernel.cu:24-34](../src/gpu/overlap_kernel.cu#L24-L34) killt Batch-Perf | ⬜ Open: Kernel wurde auf echtes Sort/Bucket-Enum ausgebaut, aber `cudaMalloc` pro Launch in `launch_overlap_kernel` (L192-195, L209) bleibt |
| 9 | gpu | Keine CUDA-Error-Checks (`cudaGetLastError()` nach Launch) | ✅ Resolved: `CUDA_CHECK` Macro + `cudaGetLastError()` nach jedem Launch (overlap_kernel.cu L227,L258) |
| 10 | gpu | Kein `cudaStream_t`-Parameter — blockiert Multi-GPU-Concurrency auf Hummel-2 H100s | ✅ Resolved: `OverlapKernelLaunchConfig::stream` hinzugefuegt (overlap_kernel.cuh L59) |
| 11 | tests | GPU-Kernel-Logik ungetestet (nur Stub) | ⚠️ Partial: `test_gpu_backend.cpp` existiert, Kernel-Logik teilweise getestet |
| 12 | ci | Keine `.github/workflows/` — ASan/TSan laufen nirgends automatisch | ✅ Resolved: `.github/workflows/ci.yml` mit build-test/asan/tsan-Jobs vorhanden |

### LOW
- ⚠️ Partial: 6 von 12 Features werden nicht in Cascade genutzt — Features 3 (ReadSpanIQR) und 7 (GcContent) sind implementiert im `feature_extractor.cpp`, aber nicht in Cascade-Stages verdrahtet. Features 4,5,6,10 weiterhin offen (fuer LightGBM v0.3).
- ⬜ Open: `CascadeConfig::min_bubble_length_bp` wird jetzt in `classify_one()` geprueft (L86); `min_coverage` + `max_recursion_depth` bleiben Dead-Config fuer non-recursion-pfad.
- ✅ Resolved: Hardcoded Username `bbg6775` in Snakefile — jetzt via `os.environ.get("USER", "bbg6775")` parametrisiert (Snakefile L15); Fallback-String bleibt, aktiver Pfad nutzt $USER.
- ⚠️ Partial: Node-Sequence-Storage — `Node::consensus` string-Member hinzugefuegt (lossless_graph.hpp L30), wird fuer abPOA consensus befuellt.
- ⬜ Open: Adjacency-Index (`edges_from_/edges_to_`) weiterhin nur im Kommentar erwaehnt (lossless_graph.hpp L47), nicht implementiert.

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

**Fehlend (Status):**
- ⚠️ Partial: Mixed-Klasse-Dekomposition — `mixed_decomposer` (SNP-Vector-Clustering) existiert, aber nicht in `classify_one()` verdrahtet
- ⬜ Open: Memory-Pool fuer GPU-Workspace (ersetzt per-call cudaMalloc)
- ⬜ Open: Fat-Binary-Runtime-Device-Selection (`cudaGetDeviceCount`/`cudaSetDevice`)
- ✅ Resolved: Integration-Tests — `test_e2e_assemble.cpp` laeuft End-to-End auf realem FASTA
- ⬜ Open: Fuzz-/Property-Tests (rapidcheck)
- ⚠️ Partial: Bounds-Error-Tests — `EXPECT_THROW` tests existieren in `test_lossless_graph.cpp` und `test_reverse_index.cpp`, keine `EXPECT_DEATH`+ASan-Variante

---

## Phase-0-Readiness Checklist

Vor empirischem Baseline-Run (1 Genom, hifiasm auf Hummel-2) fixen:

- [x] ✅ Resolved: Bounds-Checks in `node()`, `set_edge_vaf()`, `add()` post-freeze
- [x] ✅ Resolved: `/tmp/segdup_merged.bed` → `{PROJECT_DIR}/tmp/` in Snakefile
- [x] ✅ Resolved: `.github/workflows/ci.yml` mit ASan/TSan-Matrix (CPU-only Target)
- [x] ✅ Resolved: Hardcoded `SHARED`-Pfad via `config.yaml` parametrisieren (Snakefile nutzt `$USER`)
- [x] ✅ Resolved: `min_bubble_length_bp` wird in `classify_one()` geprueft

Non-blocker fuer Phase 0 (Paper-Review laut Projekt-Regel nicht Ziel):
- ⚠️ Partial: CUDA-Perf-Fixes — Kernel-Logik implementiert, Memory-Pool fehlt
- ⚠️ Partial: Mixed-Klasse — Emission ja, rekursive Dekomposition noch nicht im Cascade-Pfad
- ✅ Resolved: Integration-Tests (test_e2e_assemble.cpp)

---

## Konkrete Naechste Schritte (empfohlene Reihenfolge)

1. **PR 1**: HIGH-Issues 1-3 fixen + Death-Tests → Graph-Layer memory-safe — ✅ Resolved (d815b80, EXPECT_THROW-Tests)
2. **PR 2**: Snakefile `/tmp`-Bug + `SHARED`-Parametrisierung → Hummel-2-ready — ✅ Resolved
3. **PR 3**: `.github/workflows/ci.yml` mit ASan/TSan → CI aktiv — ✅ Resolved
4. **Phase 0 Baseline-Run** auf Hummel-2 — ✅ Resolved (HG002 SLURM 1823405 durchgelaufen 2026-04-17; 8 HPRC-Genome nachgelaufen)
5. **PR 4** (nach Phase 0): Classification-Kalibrierung + Mixed-Klasse — ⚠️ Partial (Kalibrierung done, Mixed-Recursion nicht verdrahtet)
6. **PR 5** (v0.2): CUDA-Memory-Pool + Stream-API + Multi-GPU-Dispatch — ⚠️ Partial (Stream-API done, Memory-Pool + Device-Selection offen)

## Known-Still-Open Post-Phase-0 (per user, 2026-04-17)

- ⬜ Open: Unitig-Compaction nicht final (graph_compactor wired in b0c76d0, aber noch nicht vollstaendig collapsing)
- ⬜ Open: RAM non-deterministisch (trotz bubble_detector-Fix 71922aa — weitere Quellen vermutet)
- ⬜ Open: CPU-Utilisierung ~6% (Parallelisierung greift nicht durch)
- ⬜ Open: `chrom=NA` in per-contig-Output (linear_mapper-Projektion unvollstaendig in manchen Faellen)
- ⬜ Open: Low-freq-CNV-Klassifikationsroadmap (siehe `docs/cnv_roadmap.md`)
