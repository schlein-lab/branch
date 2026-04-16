#pragma once

// BRANCH v0.1 — GFA-1.2 output with BRANCH custom tags.
//
// Writes a LosslessGraph to disk in GFA-1.2 format. Custom optional
// tags follow the specification in docs/graph-format-spec.md:
//   CN:i, CV:f (node copy-number)
//   VF:f, VC:i, BT:A, CT:f (edge branch-type + confidence)
//
// v0.1 writes nodes as S-lines and edges as L-lines (no containment
// or walk lines). Sidecar TSV files are not yet emitted; they arrive
// with the classifier in v0.3.

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>

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

}  // namespace branch::graph
