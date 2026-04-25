#pragma once

// BRANCH — Read database (lossless per-read storage).
//
// Completes the v0.1 design committed to in delta_read.hpp:
// "every read becomes a path through the string-graph plus a small list
//  of positional differences (deltas) from the path consensus".
//
// Until this module landed, a ReadPath was computed once before unitig
// compaction, written to a sidecar TSV, and never updated. After
// compaction the paths referenced dead node IDs; after containment the
// dropped reads simply vanished. Downstream callers had no way to ask
// "which reads traverse this alt?" — the bubble detector was forced to
// approximate VAF by summing edge overlap counts, which loses signal
// at sub-clonal frequencies (see HG008 SV_5 pilot).
//
// ReadDB is the persistent home of every input read, before AND after
// every graph mutation. Three guarantees:
//
//   (1) Every input read has exactly one ReadEntry. Containment never
//       deletes a read; it re-routes the read's path through the
//       absorbing predecessor.
//   (2) Reads that fail to overlap anything land in an explicit
//       Bucket (UnmappedNoOverlap) with their full sequence preserved
//       — could be contamination, could be dark genome, the caller
//       decides.
//   (3) After unitig compaction every path is reindexed to the new
//       unitig IDs with start/end offsets. Round-trip read
//       reconstruction (path consensus + deltas → original sequence)
//       is bit-exact for any read whose deltas were computed.
//
// Persistence uses standard formats:
//   - mapped reads written as GAF (Graph Alignment Format, the canonical
//     "read-as-graph-path" format used by vg, minigraph, GraphAligner).
//     The cs:Z tag carries deltas in the same notation minimap2 uses.
//   - unmapped + filtered reads written as plain FASTQ next to the GAF.
//
// Memory model: ReadDB owns the entries; views/spans are handed out for
// hot loops. Sequence storage for mapped reads is OPTIONAL — once
// deltas are computed against the post-compaction path, the original
// sequence is reconstructible from path + deltas, so the redundant
// in-memory copy can be dropped (huge win for WGS-scale inputs).

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "graph/delta_read.hpp"

namespace branch::graph {

class LosslessGraph;

// Status of a read in the database. Mapped reads have a non-empty
// path; unmapped/filtered reads keep their original sequence so a
// downstream pass can revisit them (e.g. align against a different
// reference, or flag as contamination).
enum class ReadBucket : std::uint8_t {
    Mapped = 0,                  // read has a valid path through the graph
    UnmappedNoOverlap = 1,       // read failed to overlap any other read
    UnmappedLowComplexity = 2,   // (reserved) future low-complexity filter
    Filtered = 3,                // explicitly filtered (e.g. too short)
};

// One step of a read's path through a unitig. start/end are 0-based
// offsets within the unitig consensus, end-exclusive. reverse=true
// means the read traverses the unitig in reverse-complement orientation
// — important once strand-aware overlap lands; safe default false.
struct PathStep {
    NodeId unitig{};
    std::uint32_t start{};
    std::uint32_t end{};
    bool reverse{false};
};

// One read in the database. For mapped reads `path` is non-empty and
// `deltas` describes substitutions / small indels relative to the
// concatenated unitig consensus. For unmapped reads `path` is empty
// and `original_sequence` carries the raw FASTQ bases (so the read
// is not destroyed; the caller can revisit it).
struct ReadEntry {
    ReadId id{};
    std::string name;
    std::uint32_t length_bp{};
    ReadBucket bucket{ReadBucket::UnmappedNoOverlap};
    std::vector<PathStep> path;
    std::vector<DeltaOp> deltas;
    // Held only for unmapped/filtered reads, or for mapped reads
    // until their deltas have been computed (and the sequence is
    // reconstructible). After delta computation a caller may clear
    // this string to free memory.
    std::string original_sequence;
};

// Compact lookup tables produced by graph mutators (containment,
// compaction). Each pair tells the database how to remap a node ID.
//
// remap[old_id] == new_id          -> simple rename
// remap[old_id] == kDropped        -> node was dropped without coverer
// coverer[old_id] != kNoCoverer    -> node was dropped, reroute path
//                                       through the coverer's segment
struct GraphRemap {
    static constexpr NodeId kDropped = std::numeric_limits<NodeId>::max();
    static constexpr NodeId kNoCoverer = std::numeric_limits<NodeId>::max();
    std::vector<NodeId> remap;     // size = old node count
    std::vector<NodeId> coverer;   // size = old node count
};

class ReadDB {
public:
    ReadDB() = default;

    // Append a read in any bucket. Returns the assigned ReadId
    // (the next sequential id). Caller relinquishes ownership of
    // `original_sequence` via std::move.
    ReadId add(std::string name,
               std::uint32_t length_bp,
               ReadBucket bucket,
               std::string original_sequence,
               std::vector<PathStep> path = {},
               std::vector<DeltaOp> deltas = {});

    // Total entries (mapped + unmapped + filtered).
    [[nodiscard]] std::size_t size() const noexcept { return entries_.size(); }

    // Counts per bucket. Cheap, O(N) over the entries vector.
    struct BucketCounts {
        std::size_t mapped{0};
        std::size_t unmapped_no_overlap{0};
        std::size_t unmapped_low_complexity{0};
        std::size_t filtered{0};
    };
    [[nodiscard]] BucketCounts bucket_counts() const noexcept;

    [[nodiscard]] const ReadEntry& get(ReadId id) const { return entries_.at(id); }
    [[nodiscard]] ReadEntry& get_mut(ReadId id) { return entries_.at(id); }
    [[nodiscard]] std::span<const ReadEntry> all() const noexcept { return entries_; }
    [[nodiscard]] std::span<ReadEntry> all_mut() noexcept { return entries_; }

    // Apply a node-id remap produced by graph_filter / graph_compactor.
    // Walks every mapped read's path; nodes mapped to kDropped redirect
    // through the coverer if available, otherwise the path step is
    // erased. A read whose entire path becomes empty stays in the DB
    // with bucket switched to UnmappedNoOverlap so the read is still
    // recoverable.
    void apply_remap(const GraphRemap& remap);

    // Drop original_sequence on mapped reads whose deltas are known to
    // be sufficient to reconstruct the sequence. Saves N * read_len
    // bytes of memory after assemble. No-op for unmapped reads.
    void release_mapped_sequences() noexcept;

    // GAF (Graph Alignment Format) writer + reader. Mapped reads only.
    // Companion FASTQ written to <gaf_path>.unmapped.fastq.gz for
    // unmapped + filtered reads — never silently lost.
    bool write_gaf(const std::string& gaf_path,
                   const LosslessGraph& graph) const;
    bool read_gaf(const std::string& gaf_path,
                  const LosslessGraph& graph);

private:
    std::vector<ReadEntry> entries_;
};

}  // namespace branch::graph
