# Graph-Layer Review (Part 1/4)

Review Date: 2026-04-16

## Files
| File | LOC | Purpose |
|------|-----|---------|
| delta_read.hpp | 56 | Read-as-path + delta representation |
| lossless_graph.hpp | 93 | CN-aware nodes, VAF-tagged edges |
| reverse_index.hpp | 97 | Node→Reads CSR index |

## Critical Issues

1. **[lossless_graph.hpp:70-71]** No bounds check on `node(NodeId id)`.
   - Impact: UB on invalid NodeId. Add `assert(id < nodes_.size())` or throw.

2. **[lossless_graph.hpp:83-85]** No bounds check on `set_edge_vaf(edge_index)`.
   - Impact: UB on out-of-range index.

3. **[reverse_index.hpp:33-37]** `add()` after `freeze()` silently corrupts state.
   - Impact: Violates invariant. Add `assert(!frozen_)`.

4. **[reverse_index.hpp:55]** Potential truncation `std::uint32_t` from `size_t`.
   - Impact: >4B entries wraps. Use `std::uint64_t` or static_assert limit.

## Design Feedback

- **Lossless Property OK**: No data loss paths. Reads stored as path+delta, no collapse.
- **CN-aware Nodes**: `copy_count` per-node is correct design (vs node duplication).
- **VAF-tagged Edges**: Per-edge VAF + confidence implemented correctly.
- **C++20 Idioms**: Good use of designated initializers, `std::span`, `[[nodiscard]]` missing.
- **Missing**: Node consensus sequence storage (only `length_bp`, no actual bases).
- **Missing**: Edge orientation (strand, overlap length) for proper string-graph.
- **Missing**: Adjacency lookup (edges_from/edges_to mentioned in comment but not implemented).

## TODOs

- [ ] Add bounds checks to `node()`, `set_edge_vaf()`, `add()` post-freeze
- [ ] Add `[[nodiscard]]` to accessor methods
- [ ] Implement adjacency index (edges_from_, edges_to_) per header comment
- [ ] Add Node sequence storage or reference to consensus pool
- [ ] Add Edge overlap metadata (overhang, strand)
- [ ] Consider `std::uint64_t` for CSR offsets if >4B read mappings expected
