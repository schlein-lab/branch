#pragma once

// BRANCH v0.1 — GFA-1.2 output with BRANCH custom tags.
// BRANCH v0.2 — FASTA / BED / PAF sidecar writers.
//
// Writes a LosslessGraph to disk in GFA-1.2 format. Custom optional
// tags follow the specification in docs/graph-format-spec.md:
//   CN:i, CV:f (node copy-number)
//   VF:f, VC:i, BT:A, CT:f (edge branch-type + confidence)
//
// Nodes are written as S-lines, edges as L-lines (no containment or
// walk lines).
//
// Sidecar export formats:
//   write_fasta           — one record per input read (80-col wrap).
//   write_fasta_consensus — one record per Node carrying a non-empty
//                           consensus sequence (post-compaction).
//   write_bed             — one line per node; chrom=NA placeholder.
//   write_bed_with_refs   — BED6 with real chrom/start/end sourced
//                           from linear_mapper (minimap2 asm20).
//   write_paf             — one PAF-12 record per backend overlap pair.

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "graph/lossless_graph.hpp"

namespace branch::graph {

// Write graph in GFA-1.2 format with BRANCH custom tags.
// Sequence content is omitted (S-lines carry "*" for length-only).
// Returns false on write error.
bool write_gfa(const LosslessGraph& graph, std::ostream& out);

// Convenience: write to a file path. Creates or truncates the file.
bool write_gfa(const LosslessGraph& graph, const std::string& path);

// Read a BRANCH GFA-1.2 file. Populates graph; returns false on parse
// error. Custom tags (CN:i, VF:f, etc.) are restored if present.
bool read_gfa(LosslessGraph& graph, std::istream& in);
bool read_gfa(LosslessGraph& graph, const std::string& path);

// ---------------------------------------------------------------------
// FASTA writer (reads).
//
// Emits one FASTA record per read with the sequence wrapped at 80
// columns. Kept for backwards compatibility with callers that want
// the input-reads dump.
//
// Returns false on write error or size mismatch between sequences/names.
bool write_fasta(const LosslessGraph& graph,
                 const std::vector<std::string>& sequences,
                 const std::vector<std::string>& names,
                 std::ostream& out);

bool write_fasta(const LosslessGraph& graph,
                 const std::vector<std::string>& sequences,
                 const std::vector<std::string>& names,
                 const std::string& path);

// ---------------------------------------------------------------------
// FASTA writer (node consensus).
//
// Emits one FASTA record per Node that has a non-empty consensus
// (populated by graph_compactor after unitig compaction). Record name
// is node_<id>. Nodes without consensus are skipped.
//
// Returns false on write error.
bool write_fasta_consensus(const LosslessGraph& graph, std::ostream& out);
bool write_fasta_consensus(const LosslessGraph& graph,
                           const std::string& path);

// ---------------------------------------------------------------------
// BED writer.
//
// The placeholder variant emits one BED line per node with
// chrom=NA, start=0, end=length_bp — kept for backwards compatibility
// where no reference is supplied.
//
// write_bed_with_refs shells out to linear_mapper (minimap2) with the
// per-Node consensus sequences and emits real genomic coordinates. A
// node without a high-MAPQ hit falls back to chrom=NA.
//
// Returns false on write error.
bool write_bed(const LosslessGraph& graph, std::ostream& out);

bool write_bed(const LosslessGraph& graph, const std::string& path);

// One (name, path) pair per linear reference. First hit by MAPQ wins
// when multiple references map a node. refs listed earlier are tried
// first (ties broken by order).
struct BedLinearRef {
    std::string name;
    std::string path;
};

bool write_bed_with_refs(const LosslessGraph& graph,
                         const std::vector<BedLinearRef>& refs,
                         const std::string& path,
                         int threads = 4);

// ---------------------------------------------------------------------
// PAF writer.
//
// Emits one standard PAF-12 record per OverlapPair produced by the
// backend (`read_a`, `read_b` are 0-based indices into `sequences`
// and `names`). Columns:
//   qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq
//
// Behaviour:
//   - strand '+' or '-' is the cluster-majority XOR of per-minimizer
//     strand bits from OverlapPair::strand (0 = same strand / '+',
//     1 = opposite strand / '-'). Minimizer strand is already per-hit
//     via canonical_hash_with_strand, so pair-level strand is accurate.
//   - matches is estimated as overlap_len - diff_count (best-effort;
//     the real aligner computes a precise value later).
//   - mapq is derived from identity = matches/overlap_len via
//     min(60, int(identity * 60)).
//
// Returns false on write error.
bool write_paf(const std::vector<branch::backend::OverlapPair>& pairs,
               const std::vector<std::string>& sequences,
               const std::vector<std::string>& names,
               std::ostream& out);

bool write_paf(const std::vector<branch::backend::OverlapPair>& pairs,
               const std::vector<std::string>& sequences,
               const std::vector<std::string>& names,
               const std::string& path);

}  // namespace branch::graph
