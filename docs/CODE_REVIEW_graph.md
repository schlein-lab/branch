# Graph-Layer Review (Part 1/4)

Review Date: 2026-04-16
Status audit: 2026-04-17 — markers added per finding.

## Files
| File | LOC | Purpose |
|------|-----|---------|
| delta_read.hpp | 56 | Read-as-path + delta representation |
| lossless_graph.hpp | 93 | CN-aware nodes, VAF-tagged edges |
| reverse_index.hpp | 97 | Node→Reads CSR index |

## Critical Issues

1. ✅ Resolved — **[lossless_graph.hpp:70-71]** No bounds check on `node(NodeId id)`.
   - Impact: UB on invalid NodeId. Add `assert(id < nodes_.size())` or throw.
   - Why resolved: `node()` now throws `std::out_of_range` (lossless_graph.hpp L77-79, L84-86). Fixed in commit d815b80 "Consolidate bounds-checks and Zyrkel-review integrations".

2. ✅ Resolved — **[lossless_graph.hpp:83-85]** No bounds check on `set_edge_vaf(edge_index)`.
   - Impact: UB on out-of-range index.
   - Why resolved: `set_edge_vaf()` now throws `std::out_of_range` (L104-106), `set_copy_count()` too (L95-97). Commit d815b80.

3. ✅ Resolved — **[reverse_index.hpp:33-37]** `add()` after `freeze()` silently corrupts state.
   - Impact: Violates invariant. Add `assert(!frozen_)`.
   - Why resolved: `add()` throws `std::logic_error("add() called after freeze()")` (reverse_index.hpp L35-37); `reads_for()` similarly guards pre-freeze. Commit d815b80.

4. ⬜ Open — **[reverse_index.hpp:55]** Potential truncation `std::uint32_t` from `size_t`.
   - Impact: >4B entries wraps. Use `std::uint64_t` or static_assert limit.
   - Why open: offsets_ still `std::vector<std::uint32_t>` (L102); no static_assert added. Not yet reached >4B in practice.

## Design Feedback

- **Lossless Property OK**: No data loss paths. Reads stored as path+delta, no collapse.
- **CN-aware Nodes**: `copy_count` per-node is correct design (vs node duplication).
- **VAF-tagged Edges**: Per-edge VAF + confidence implemented correctly.
- **C++20 Idioms**: Good use of designated initializers, `std::span`, `[[nodiscard]]` missing.
- **Missing**: Node consensus sequence storage (only `length_bp`, no actual bases).
- **Missing**: Edge orientation (strand, overlap length) for proper string-graph.
- **Missing**: Adjacency lookup (edges_from/edges_to mentioned in comment but not implemented).

## TODOs

- [x] ✅ Resolved: Add bounds checks to `node()`, `set_edge_vaf()`, `add()` post-freeze (d815b80)
- [x] ✅ Resolved: Add `[[nodiscard]]` to accessor methods (present on `node_count()`, `edge_count()`, `node()`, `nodes()`, `edges()`)
- [ ] ⬜ Open: Implement adjacency index (edges_from_, edges_to_) per header comment
- [x] ⚠️ Partial: Add Node sequence storage or reference to consensus pool — `Node::consensus` string-Member vorhanden (lossless_graph.hpp L30), wird befuellt, aber kein zentraler Consensus-Pool
- [ ] ⬜ Open: Add Edge overlap metadata (overhang, strand) — Edge hat nur `flags`, keine expliziten overhang/strand-Felder
- [ ] ⬜ Open: Consider `std::uint64_t` for CSR offsets if >4B read mappings expected

## Post-Phase-0 additions (not in original review)

graph_compactor (unitig compaction, commit b0c76d0) + graph_filter (transitive reduction + containment removal, commit a8dfb1d) are now wired into the assemble pipeline. Both are Phase-0 extensions of the graph layer beyond the v0.1 scope reviewed here.
