#pragma once

// BRANCH v0.1 — GFA-1.2 output with BRANCH custom tags.
// BRANCH v0.2 — FASTA / BED / PAF sidecar writers.
//
// Writes a LosslessGraph to disk in GFA-1.2 format. Custom optional
// tags follow the specification in docs/graph-format-spec.md:
//   CN:i, CV:f (node copy-number)
//   VF:f, VC:i, BT:A, CT:f (edge branch-type + confidence)
//
// v0.1 writes nodes as S-lines and edges as L-lines (no containment
// or walk lines). Sidecar TSV files are not yet emitted; they arrive
// with the classifier in v0.3.
//
// v0.2 adds three additional export formats for downstream tooling:
//   write_fasta — one record per input read (v0.3 will switch to
//                 per-node consensus once consensus_seq is tracked).
//   write_bed   — one line per node; placeholder NA chrom until v0.3
//                 reference alignment supplies genomic coordinates.
//   write_paf   — one PAF-12 record per overlap pair from the backend.

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
// v0.2 assumptions:
//   - strand '+' or '-' is taken from OverlapPair::strand (0 = +, 1 = -).
//     TODO(v0.3): propagate per-side strand from minimizer_sketcher.
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
