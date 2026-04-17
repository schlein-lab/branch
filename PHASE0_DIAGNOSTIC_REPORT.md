# Phase-0 Critical Findings Diagnostic Report
**BRANCH Pangenome Assembler — SLURM Job 1823405 (Hummel-2)**
**8 HPRC genomes + HG002**

---

## 1. NO UNITIG COLLAPSE

### (a) Root Cause
The `is_chain_start()` lambda in `compact_unitigs()` accesses `unique_in_edge[v]` without guarding against uninitialized values for source nodes (in_degree == 0). The array `unique_in_edge` is only populated when `in_degree == 1` (line 160), but the lookup at line 172 occurs unconditionally after the `in_degree != 1` check returns `true`. This causes **undefined behavior** (out-of-bounds access with index `max(size_t)`), corrupting compactor state.

### (b) File:Line Evidence
| File | Line | Code |
|------|------|------|
| `/home/christian/branch/src/graph/graph_compactor.cpp` | 154 | `std::vector<std::size_t> unique_in_edge(n, std::numeric_limits<std::size_t>::max());` |
| `/home/christian/branch/src/graph/graph_compactor.cpp` | 160 | `if (e.to < n && in_degree[e.to] == 1) { unique_in_edge[e.to] = i; }` — only sets when in_degree==1 |
| `/home/christian/branch/src/graph/graph_compactor.cpp` | 169-175 | `is_chain_start()` lambda accesses `edges[unique_in_edge[v]]` without bounds check |
| `/home/christian/branch/src/graph/graph_compactor.cpp` | 299-309 | Same bug duplicated in `compact_unitigs_with_sequences()` |
| `/home/christian/branch/src/cli/assemble_cmd.cpp` | 241-245 | Invocation after `filter_graph()` — correct order, but compactor fails |

### (c) Proposed Fix
Add early-return guard in `is_chain_start()`:
```cpp
auto is_chain_start = [&](NodeId v) -> bool {
    if (in_degree[v] != 1) return true;  // source or branching → chain start
    // Safe: in_degree[v] == 1, so unique_in_edge[v] was set
    std::size_t ein = unique_in_edge[v];
    ...
};
```
Apply same fix to both overloads (lines ~170 and ~300).

### (d) Effort Estimate
**1.5 hours** — straightforward guard insertion, add test case for source-node chains.

---

## 2. RAM NON-DETERMINISTIC

### (a) Root Cause
Unmitigated `std::unordered_map` iteration in **bubble_detector.cpp** (lines 118-150) iterates over `reachable_within()` results without sorting keys. Bubble detection output order varies by hash distribution, causing downstream vector allocations in different orders. Secondary: `std::async` future join order in cpu_backend.cpp depends on thread scheduling.

### (b) File:Line Evidence
| File | Line | Issue |
|------|------|-------|
| `/home/christian/branch/src/detect/bubble_detector.cpp` | 51-82 | `reachable_within()` returns `unordered_map<NodeId, PathRecord>` |
| `/home/christian/branch/src/detect/bubble_detector.cpp` | 118-150 | Direct iteration: `for (const auto& [exit_node, path_i] : *small)` — **UNMITIGATED** |
| `/home/christian/branch/src/detect/bubble_detector.cpp` | 105 | `unordered_map<EntryExitKey, Bubble, EntryExitHash> bubbles;` — no deterministic iteration |
| `/home/christian/branch/src/classify/feature_extractor.cpp` | 25-34 | `unordered_set<string_view>` iteration in Jaccard similarity |
| `/home/christian/branch/src/backend/cpu_backend.cpp` | 114-143 | `std::async` futures joined in completion order, not submission order |

**Mitigated (not a problem):** `cpu_backend.cpp` lines 166-174 and 220-229 correctly sort bucket keys before iteration.

### (c) Proposed Fix
1. **bubble_detector.cpp**: Extract keys from `reachable_from_i/j` maps, sort, then iterate:
   ```cpp
   std::vector<NodeId> exits;
   for (const auto& kv : *small) exits.push_back(kv.first);
   std::sort(exits.begin(), exits.end());
   for (NodeId exit_node : exits) { ... }
   ```
2. **feature_extractor.cpp**: Return `std::set` or sort k-mers before comparison.
3. **cpu_backend.cpp**: Use indexed futures or thread pool with deterministic task ordering.

### (d) Effort Estimate
**3 hours** — bubble_detector is the critical path; others are secondary.

---

## 3. CPU ~6% UTILIZATION

### (a) Root Cause
Three sequential bottlenecks dominate the pipeline with **zero OpenMP parallelization**:
1. **Transitive reduction** (graph_filter.cpp lines 160-213): O(V·E²) triple-nested loop, entirely single-threaded.
2. **Unitig consensus** (graph_compactor.cpp lines 352-376): Per-unitig `simple_majority_consensus()` in sequential for-loop.
3. **Containment check** (graph_filter.cpp lines 82-133): O(V·E_in·E_out²) with inner `has_out_edge()` linear scan.

### (b) File:Line Evidence
| File | Line | Bottleneck |
|------|------|-----------|
| `/home/christian/branch/src/graph/graph_filter.cpp` | 160-213 | Transitive reduction: `for (NodeId a...) { for (kb...) { for (kc...) { }}}` — **NO OMP** |
| `/home/christian/branch/src/graph/graph_filter.cpp` | 98-103 | `has_out_edge()` does linear scan per check |
| `/home/christian/branch/src/graph/graph_compactor.cpp` | 352-376 | `for (const auto& members : unitigs) { ... consensus ...}` — **SEQUENTIAL** |
| `/home/christian/branch/src/backend/cpu_backend.cpp` | 182-217 | Hash-bucket pair enumeration O(B²) per hash — sequential |
| `/home/christian/branch/src/backend/cpu_backend.cpp` | 138-142 | `f.get()` blocks sequentially on each future |

### (c) Proposed Fix
1. **graph_filter.cpp:160**: Add `#pragma omp parallel for schedule(dynamic)` around node loop, use thread-local `remove_set`, merge at end.
2. **graph_compactor.cpp:352**: Parallelize unitig consensus with task pool or `#pragma omp parallel for`.
3. **graph_filter.cpp:98**: Replace `has_out_edge()` linear scan with `std::unordered_set` for O(1) lookup.
4. **cpu_backend.cpp:182**: Partition bucket_keys across threads.

### (d) Effort Estimate
**8 hours** — transitive reduction parallelization requires careful handling of shared `remove_set`; consensus is straightforward; has_out_edge refactor is medium.

---

## 4. chrom=NA IN BED OUTPUT

### (a) Root Cause
Two code paths produce `chrom=NA`:
1. **No reference provided**: `write_bed()` (line 268) hardcodes `"NA"` when `--ref-linear` is not specified.
2. **Unmapped nodes**: `write_bed_with_refs()` (line 382) falls back to `"NA"` when a node has no consensus or minimap2 produces no alignment.

### (b) File:Line Evidence
| File | Line | Code |
|------|------|------|
| `/home/christian/branch/src/graph/graph_io.cpp` | 268 | `out << "NA" << '\t' ...` — hardcoded in `write_bed()` |
| `/home/christian/branch/src/graph/graph_io.cpp` | 307-310 | `if (node.consensus.empty()) continue;` — skips nodes without consensus |
| `/home/christian/branch/src/graph/graph_io.cpp` | 378-382 | `if (it == best_by_node.end()) { rows.push_back({"NA", ...}); }` — fallback |
| `/home/christian/branch/src/cli/assemble_cmd.cpp` | 262-271 | Chooses `write_bed()` vs `write_bed_with_refs()` based on `--ref-linear` |
| `/home/christian/branch/src/project/linear_mapper.cpp` | 61 | `out.target = target_name;` — source of chrom from minimap2 PAF col 6 |

### (c) Proposed Fix
1. **Ensure consensus is populated**: Check `graph_compactor.cpp` line 374 — if `seqs.empty()`, consensus stays empty. Investigate why nodes lack sequences.
2. **Propagate reference from graph_build**: If building from reference-anchored reads, store contig-of-origin in node metadata during `graph_build()`, then use as fallback in BED output.
3. **Improve logging**: Add warning when nodes fall back to `NA` so user knows to provide `--ref-linear`.

### (d) Effort Estimate
**2 hours** — if root cause is missing `--ref-linear` flag, it's a documentation fix; if consensus is missing due to compactor bug (#1), fixing #1 resolves this.

---

## Summary

| Finding | Root Cause | Critical File | Effort |
|---------|-----------|---------------|--------|
| No unitig collapse | UB in `is_chain_start()` accessing uninitialized array | graph_compactor.cpp:169-175 | 1.5h |
| RAM non-deterministic | Unmitigated unordered_map iteration in bubble_detector | bubble_detector.cpp:118-150 | 3h |
| CPU 6% utilization | Sequential triple-nested loops, no OMP | graph_filter.cpp:160-213 | 8h |
| chrom=NA in BED | Hardcoded fallback + missing consensus | graph_io.cpp:268,382 | 2h |

**Total estimated effort: 14.5 hours**
