# BRANCH Code Review: Tests + Workflow (Teil 4/4)

**Datum**: 2026-04-16
**Scope**: `tests/*.cpp`, `tests/CMakeLists.txt`, `workflow/Snakefile`, Root `CMakeLists.txt`

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

| Gap | Severity | Details |
|-----|----------|---------|
| No bounds-error tests | HIGH | `node(invalid_id)` UB not tested |
| No negative/edge-case tests | HIGH | Empty graph edge traversal, freeze()-then-add |
| No fuzz/property tests | MEDIUM | rapidcheck/hypothesis missing |
| Mixed class never tested | MEDIUM | Cascade only tests Branch/Dup/NonSeparable |
| GPU kernel logic untested | MEDIUM | Stub-only — real kernel has 0 tests |
| No integration test | LOW | End-to-end BAM→Graph→Classify |

## ASan/TSan Status

✅ **CMake Options present**: `BRANCH_USE_ASAN`, `BRANCH_USE_TSAN`, `BRANCH_USE_UBSAN`
⚠️ **NOT enabled by default** — CI must explicitly pass `-DBRANCH_USE_ASAN=ON`
⚠️ **No CI config found** — `.github/workflows/` missing
⚠️ **CUDA + ASan conflict** — nvcc doesn't support ASan; need separate CPU-only test target

## Workflow Issues (Hummel-2)

| Issue | Severity | Details |
|-------|----------|---------|
| `/tmp/` merge files | HIGH | `roi_catalog` writes to `/tmp/segdup_merged.bed` — not BeeGFS, lost on job end |
| Hardcoded `bbg6775` | MEDIUM | `SHARED` path hardcoded, should use `{USER}` |
| `mem_mb = 0` | LOW | Works (Hummel ignores), but unclear intent |
| No `--use-conda` | INFO | Containers via Apptainer — OK, but Snakemake 9 prefers `--software-deployment-method` |
| `set +u` pattern | OK | Correctly handles RRZ init.sh unbound vars |
| `TMPDIR=/tmp` | OK | Correct per Hummel-2 docs |

## TODOs

1. **[HIGH]** Add bounds-error death tests: `EXPECT_DEATH(g.node(999), "")` with ASan
2. **[HIGH]** Fix `/tmp/` in roi_catalog → use `{PROJECT_DIR}/tmp/`
3. **[HIGH]** Create `.github/workflows/ci.yml` with ASan/TSan matrix
4. **[MEDIUM]** Add `test_cascade_mixed.cpp` for Mixed-class recursive decomposition
5. **[MEDIUM]** Property-based tests with rapidcheck for ReverseIndex invariants
6. **[MEDIUM]** Parametrize `SHARED` path in Snakefile via config.yaml
7. **[LOW]** Add integration test: synthetic BAM → full pipeline

---
*Review 4/4 abgeschlossen.*
