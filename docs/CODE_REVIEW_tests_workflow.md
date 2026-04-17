# BRANCH Code Review: Tests + Workflow (Teil 4/4)

**Datum**: 2026-04-16
**Scope**: `tests/*.cpp`, `tests/CMakeLists.txt`, `workflow/Snakefile`, Root `CMakeLists.txt`

Status audit: 2026-04-17 — markers added per finding.

## Files

| File | LOC | Was getestet |
|------|-----|--------------|
| test_backend_vtable.cpp | 111 | VTable dispatch, move semantics, destroy lifecycle |
| test_cascade.cpp | 78 | 4-Stage cascade, threshold config, NonSeparable fallback |
| test_delta_read.cpp | 49 | PositionalDelta packing, ReadPath/View basic |
| test_features.cpp | 34 | kNumFeatures=12, FeatureVector size, enum indexing |
| test_gpu_stub.cpp | 21 | MinimizerEntry 16B, stub returns 0 |
| test_lossless_graph.cpp | 59 | add_node/edge, set_copy_count/vaf |
| test_reverse_index.cpp | 65 | Mutable add, freeze(), CSR lookup |
| tests/CMakeLists.txt | 53 | GoogleTest FetchContent, 8 test executables |
| workflow/Snakefile | 117 | align_hifi, roi_catalog, coverage_profile |
| CMakeLists.txt | 74 | C++20, CUDA fat-binary, ASan/TSan/UBSan options |

## Coverage Gaps

| Gap | Severity | Details | Status |
|-----|----------|---------|--------|
| No bounds-error tests | HIGH | `node(invalid_id)` UB not tested | ✅ Resolved: `EXPECT_THROW`-Tests in `test_lossless_graph.cpp` und `test_reverse_index.cpp` |
| No negative/edge-case tests | HIGH | Empty graph edge traversal, freeze()-then-add | ✅ Resolved: Negative-Tests inkl. freeze-then-add-Logic vorhanden |
| No fuzz/property tests | MEDIUM | rapidcheck/hypothesis missing | ⬜ Open: kein rapidcheck in `tests/CMakeLists.txt` |
| Mixed class never tested | MEDIUM | Cascade only tests Branch/Dup/NonSeparable | ✅ Resolved: `test_mixed_decomposer.cpp` + Mixed-Branch in test_cascade |
| GPU kernel logic untested | MEDIUM | Stub-only — real kernel has 0 tests | ⚠️ Partial: `test_gpu_backend.cpp` existiert und testet GpuBackend, Kernel-interne Logik nur indirekt |
| No integration test | LOW | End-to-end BAM→Graph→Classify | ✅ Resolved: `test_e2e_assemble.cpp` (E2E auf realem FASTA), `test_project_smoke.cpp` |

## ASan/TSan Status

✅ **CMake Options present**: `BRANCH_USE_ASAN`, `BRANCH_USE_TSAN`, `BRANCH_USE_UBSAN`
- ✅ Resolved: CI matrix aktiviert ASan/TSan-Jobs in `.github/workflows/ci.yml` (L32 `asan:`, L70 `-DBRANCH_USE_TSAN=ON`)
- ✅ Resolved: CI-Config vorhanden — `.github/workflows/ci.yml` mit build-test/asan/tsan-Jobs
- ⬜ Open: CUDA + ASan-Konflikt — CI baut mit `BRANCH_BUILD_CUDA=OFF` (L25), dadurch gemieden, aber kein dedizierter CUDA-test-Job

## Workflow Issues (Hummel-2)

| Issue | Severity | Details | Status |
|-------|----------|---------|--------|
| `/tmp/` merge files | HIGH | `roi_catalog` writes to `/tmp/segdup_merged.bed` — not BeeGFS, lost on job end | ✅ Resolved: `roi_catalog` schreibt nach `{TMP_DIR}/segdup_merged.bed` (Snakefile L111,114) mit `TMP_DIR = config.get("tmp_dir", f"{PROJECT_DIR}/tmp")` |
| Hardcoded `bbg6775` | MEDIUM | `SHARED` path hardcoded, should use `{USER}` | ✅ Resolved: `USER = os.environ.get("USER", "bbg6775")` (Snakefile L15), `SHARED`/`BEEGFS_BASE` nutzen `{USER}` |
| `mem_mb = 0` | LOW | Works (Hummel ignores), but unclear intent | ⬜ Open: weiterhin `mem_mb = 0` in Regeln, keine Klaerung im Kommentar |
| No `--use-conda` | INFO | Containers via Apptainer — OK, but Snakemake 9 prefers `--software-deployment-method` | ⬜ Open: keine explizite Anpassung im Snakefile, Container-Bindings direkt via `apptainer exec` |
| `set +u` pattern | OK | Correctly handles RRZ init.sh unbound vars | ✅ Resolved (already OK) |
| `TMPDIR=/tmp` | OK | Correct per Hummel-2 docs | ✅ Resolved (already OK) |

## TODOs

1. **[HIGH]** Add bounds-error death tests: `EXPECT_DEATH(g.node(999), "")` with ASan — ⚠️ Partial: `EXPECT_THROW` statt `EXPECT_DEATH`, ASan-Job matched darauf nicht direkt
2. **[HIGH]** Fix `/tmp/` in roi_catalog → use `{PROJECT_DIR}/tmp/` — ✅ Resolved
3. **[HIGH]** Create `.github/workflows/ci.yml` with ASan/TSan matrix — ✅ Resolved
4. **[MEDIUM]** Add `test_cascade_mixed.cpp` for Mixed-class recursive decomposition — ✅ Resolved (`test_mixed_decomposer.cpp`)
5. **[MEDIUM]** Property-based tests with rapidcheck for ReverseIndex invariants — ⬜ Open
6. **[MEDIUM]** Parametrize `SHARED` path in Snakefile via config.yaml — ✅ Resolved (via `{USER}` + `config.get("shared_dir", ...)`)
7. **[LOW]** Add integration test: synthetic BAM → full pipeline — ✅ Resolved (test_e2e_assemble.cpp)

---
*Review 4/4 abgeschlossen.*
