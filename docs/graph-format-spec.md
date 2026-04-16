# BRANCH GFA-1.2 Extension Specification v0.1

## Motivation

Standard GFA-1.2 has no native concept of copy-number-aware nodes, VAF-tagged edges, or branch/duplication classification. BRANCH defines a minimal extension using custom optional tags that downstream tools can ignore without breaking standard parsing.

## Node Custom Tags

| Tag  | Type   | Meaning                                                                 |
|------|--------|-------------------------------------------------------------------------|
| CN:i | int    | Copy-number estimate (integer ≥1)                                       |
| CV:f | float  | CN confidence value in [0,1]                                            |
| CC:i | int    | Raw read coverage count at node                                         |
| HP:Z | string | Hot-node materialized consensus sequence hash (SHA256 first 16 hex) if present |

## Edge Custom Tags

| Tag  | Type   | Meaning                                                        |
|------|--------|----------------------------------------------------------------|
| VF:f | float  | Variant allele frequency in [0,1]                              |
| VC:i | int    | Reads supporting this edge                                     |
| BT:A | char   | Branch type: B=branch, D=duplication, M=mixed, N=non-separable |
| CT:f | float  | Classification confidence in [0,1]                             |
| FS:Z | string | Flank hash pair as "left_sha:right_sha" (16 hex each)          |

## Sidecar Files

Optional companion TSV files with same basename:

- `basename.cn.tsv`: per-node CN distribution (node_id, cn_estimate, ci_low, ci_high, features...)
- `basename.branches.tsv`: per-edge branch classifications (edge_id, class, confidence, features...)
- `basename.classifier.json`: model metadata (type, version, training-set hash, feature list)

## Reader Requirements

- Must gracefully ignore unknown tags (GFA-1.2 compliant behavior)
- Sidecar TSV parsing is optional but recommended for full BRANCH semantics
