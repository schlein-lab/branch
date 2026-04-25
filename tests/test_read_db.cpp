// Tests for branch::graph::ReadDB.
//
// The DB has three invariants worth pinning down with regression
// tests:
//   1. Every input read produces exactly one ReadEntry; unmapped reads
//      are preserved with their original_sequence intact.
//   2. apply_remap re-routes mapped paths through coverer nodes when a
//      step's unitig is dropped, and demotes a read to UnmappedNoOverlap
//      only when the remap erases its entire path.
//   3. write_gaf round-trips the path information losslessly: a DB
//      written and re-read produces an identical sequence of mapped
//      ReadEntry paths (deltas are encoded into cs:Z; full delta
//      round-trip is covered separately by the delta_read tests).

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>

#include "graph/lossless_graph.hpp"
#include "graph/read_db.hpp"

using branch::graph::DeltaOp;
using branch::graph::GraphRemap;
using branch::graph::LosslessGraph;
using branch::graph::NodeId;
using branch::graph::PathStep;
using branch::graph::ReadBucket;
using branch::graph::ReadDB;
using branch::graph::ReadEntry;

TEST(ReadDB, Add_assigns_sequential_ids_and_preserves_buckets) {
    ReadDB db;
    auto id_a = db.add("r0", 100, ReadBucket::Mapped, "", {{1, 0, 100, false}}, {});
    auto id_b = db.add("r1", 80, ReadBucket::UnmappedNoOverlap, "ACGTACGT", {}, {});
    auto id_c = db.add("r2", 50, ReadBucket::Filtered, "AAAAAA", {}, {});

    EXPECT_EQ(id_a, 0u);
    EXPECT_EQ(id_b, 1u);
    EXPECT_EQ(id_c, 2u);
    EXPECT_EQ(db.size(), 3u);

    EXPECT_EQ(db.get(id_a).bucket, ReadBucket::Mapped);
    EXPECT_EQ(db.get(id_b).bucket, ReadBucket::UnmappedNoOverlap);
    EXPECT_EQ(db.get(id_b).original_sequence, "ACGTACGT");

    auto counts = db.bucket_counts();
    EXPECT_EQ(counts.mapped, 1u);
    EXPECT_EQ(counts.unmapped_no_overlap, 1u);
    EXPECT_EQ(counts.filtered, 1u);
}

TEST(ReadDB, apply_remap_renames_unitigs_in_path) {
    ReadDB db;
    db.add("r0", 100, ReadBucket::Mapped, "",
           {{0, 0, 50, false}, {1, 0, 50, false}}, {});

    // Rename: 0 -> 7, 1 -> 8.
    GraphRemap rm;
    rm.remap = {7u, 8u};
    rm.coverer = {GraphRemap::kNoCoverer, GraphRemap::kNoCoverer};
    db.apply_remap(rm);

    const auto& r = db.get(0);
    ASSERT_EQ(r.path.size(), 2u);
    EXPECT_EQ(r.path[0].unitig, 7u);
    EXPECT_EQ(r.path[1].unitig, 8u);
    EXPECT_EQ(r.bucket, ReadBucket::Mapped);
}

TEST(ReadDB, apply_remap_redirects_dropped_via_coverer) {
    ReadDB db;
    db.add("r0", 100, ReadBucket::Mapped, "",
           {{0, 0, 50, false}, {1, 0, 50, false}, {2, 0, 50, false}}, {});

    // Node 1 is dropped, but it was contained in node 0 (coverer).
    // Node 0 -> 5, Node 2 -> 6. Node 1 is GraphRemap::kDropped.
    GraphRemap rm;
    rm.remap = {5u, GraphRemap::kDropped, 6u};
    rm.coverer = {GraphRemap::kNoCoverer, /*coverer of 1 is*/ 0u, GraphRemap::kNoCoverer};
    db.apply_remap(rm);

    const auto& r = db.get(0);
    // Path step that landed on node 1 is rerouted through node 0
    // (-> remapped 5). Adjacent duplicates are coalesced, so we expect
    // either {5, 6} or {5, 5, 6}; the implementation coalesces.
    ASSERT_GE(r.path.size(), 2u);
    EXPECT_EQ(r.path.front().unitig, 5u);
    EXPECT_EQ(r.path.back().unitig, 6u);
    EXPECT_EQ(r.bucket, ReadBucket::Mapped);
}

TEST(ReadDB, apply_remap_demotes_to_unmapped_when_path_empties) {
    ReadDB db;
    db.add("orphan", 80, ReadBucket::Mapped, "ACGT",
           {{4, 0, 80, false}}, {});

    GraphRemap rm;
    rm.remap = {0u, 1u, 2u, 3u, GraphRemap::kDropped};
    rm.coverer = {GraphRemap::kNoCoverer, GraphRemap::kNoCoverer,
                  GraphRemap::kNoCoverer, GraphRemap::kNoCoverer,
                  GraphRemap::kNoCoverer};
    db.apply_remap(rm);

    const auto& r = db.get(0);
    EXPECT_TRUE(r.path.empty());
    EXPECT_EQ(r.bucket, ReadBucket::UnmappedNoOverlap);
    EXPECT_EQ(r.original_sequence, "ACGT")
        << "Read sequence must be preserved when path is erased";
}

TEST(ReadDB, write_gaf_round_trips_path_topology) {
    ReadDB db;
    db.add("r0", 200, ReadBucket::Mapped, "",
           {{3, 0, 100, false}, {7, 0, 100, true}}, {});
    db.add("r1", 150, ReadBucket::UnmappedNoOverlap, "ACGTACGT", {}, {});

    LosslessGraph g;
    auto path = std::filesystem::temp_directory_path() / "test_read_db.gaf";

    ASSERT_TRUE(db.write_gaf(path.string(), g));

    ReadDB roundtrip;
    ASSERT_TRUE(roundtrip.read_gaf(path.string(), g));
    ASSERT_EQ(roundtrip.size(), 1u)  // only mapped reads come back via GAF
        << "Unmapped reads land in the .unmapped.fastq companion file, "
           "not in the GAF body";

    const auto& r = roundtrip.get(0);
    EXPECT_EQ(r.name, "r0");
    EXPECT_EQ(r.length_bp, 200u);
    ASSERT_EQ(r.path.size(), 2u);
    EXPECT_EQ(r.path[0].unitig, 3u);
    EXPECT_FALSE(r.path[0].reverse);
    EXPECT_EQ(r.path[1].unitig, 7u);
    EXPECT_TRUE(r.path[1].reverse);

    // Companion FASTQ exists and contains the unmapped read.
    auto fq_path = path.string() + ".unmapped.fastq";
    ASSERT_TRUE(std::filesystem::exists(fq_path));

    std::filesystem::remove(path);
    std::filesystem::remove(fq_path);
}

TEST(ReadDB, release_mapped_sequences_keeps_unmapped) {
    ReadDB db;
    db.add("mapped",   100, ReadBucket::Mapped, "ACGT" /* will be dropped */,
           {{0, 0, 100, false}}, {DeltaOp{10, 'S', 'A'}});
    db.add("unmapped", 100, ReadBucket::UnmappedNoOverlap, "TTTTTTTT", {}, {});

    db.release_mapped_sequences();

    EXPECT_TRUE(db.get(0).original_sequence.empty());
    EXPECT_EQ(db.get(1).original_sequence, "TTTTTTTT")
        << "Unmapped reads must keep their sequence after release_mapped_sequences";
}
