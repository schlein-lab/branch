#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "graph/graph_io.hpp"

using branch::graph::LosslessGraph;
using branch::graph::NodeId;
using branch::graph::read_gfa;
using branch::graph::write_bed;
using branch::graph::write_fasta;
using branch::graph::write_gfa;
using branch::graph::write_paf;

TEST(GraphIOTest, WriteGfa_produces_header_and_S_lines) {
    LosslessGraph g;
    g.add_node(1000, 2);
    g.add_node(500);
    g.add_edge(0, 1, 30);

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));

    std::string gfa = oss.str();
    EXPECT_NE(gfa.find("H\tVN:Z:1.2"), std::string::npos);
    EXPECT_NE(gfa.find("S\t0\t*\tLN:i:1000"), std::string::npos);
    EXPECT_NE(gfa.find("CN:i:2"), std::string::npos);
    EXPECT_NE(gfa.find("L\t0\t+\t1\t+\t0M"), std::string::npos);
    EXPECT_NE(gfa.find("VC:i:30"), std::string::npos);
}

TEST(GraphIOTest, Roundtrip_preserves_graph_structure) {
    LosslessGraph g;
    g.add_node(2000, 3);
    g.add_node(1500);
    g.add_edge(0, 1, 42);
    g.set_copy_count(0, 3, 0.95f);
    g.set_edge_vaf(0, 0.07f, 0.88f);

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));

    LosslessGraph g2;
    std::istringstream iss(oss.str());
    ASSERT_TRUE(read_gfa(g2, iss));

    EXPECT_EQ(g2.node_count(), 2u);
    EXPECT_EQ(g2.edge_count(), 1u);
    EXPECT_EQ(g2.node(0).length_bp, 2000u);
    EXPECT_EQ(g2.node(0).copy_count, 3u);
    EXPECT_NEAR(g2.node(0).copy_count_confidence, 0.95f, 0.01f);
    EXPECT_EQ(g2.edges()[0].from, 0u);
    EXPECT_EQ(g2.edges()[0].to, 1u);
    EXPECT_EQ(g2.edges()[0].read_support, 42u);
    EXPECT_NEAR(g2.edges()[0].vaf, 0.07f, 0.001f);
    EXPECT_NEAR(g2.edges()[0].vaf_confidence, 0.88f, 0.01f);
}

TEST(GraphIOTest, WriteGfa_emits_RC_tag_per_S_line) {
    // Every S-line must carry RC:i:<node.read_support>.
    LosslessGraph g;
    g.add_node(1000);
    g.add_node(500);
    g.add_node(250);
    g.node(0).read_support = 7;
    g.node(1).read_support = 0;
    g.node(2).read_support = 123;

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));
    const std::string gfa = oss.str();

    // Exact format — node_id followed by sequence placeholder then tags.
    EXPECT_NE(gfa.find("S\t0\t*\tLN:i:1000\tRC:i:7\t"), std::string::npos);
    EXPECT_NE(gfa.find("S\t1\t*\tLN:i:500\tRC:i:0\t"),  std::string::npos);
    EXPECT_NE(gfa.find("S\t2\t*\tLN:i:250\tRC:i:123\t"), std::string::npos);
}

TEST(GraphIOTest, Roundtrip_preserves_RC_read_support) {
    // write -> read -> every node's read_support must match.
    LosslessGraph g;
    g.add_node(1000);
    g.add_node(2000);
    g.add_node(500);
    g.node(0).read_support = 42;
    g.node(1).read_support = 0;
    g.node(2).read_support = 1'000'000;

    std::ostringstream oss;
    ASSERT_TRUE(write_gfa(g, oss));

    LosslessGraph g2;
    std::istringstream iss(oss.str());
    ASSERT_TRUE(read_gfa(g2, iss));

    ASSERT_EQ(g2.node_count(), 3u);
    EXPECT_EQ(g2.node(0).read_support, 42u);
    EXPECT_EQ(g2.node(1).read_support, 0u);
    EXPECT_EQ(g2.node(2).read_support, 1'000'000u);
}

TEST(GraphIOTest, ReadGfa_skips_unknown_lines_gracefully) {
    std::string gfa = "H\tVN:Z:1.2\n"
                      "# this is a comment\n"
                      "S\t0\t*\tLN:i:100\n"
                      "W\tunknown\twalk\tline\n"
                      "S\t1\t*\tLN:i:200\n";
    LosslessGraph g;
    std::istringstream iss(gfa);
    ASSERT_TRUE(read_gfa(g, iss));
    EXPECT_EQ(g.node_count(), 2u);
}

// ---------------------------------------------------------------------
// FASTA writer
// ---------------------------------------------------------------------

TEST(GraphIOTest, FASTA_writer_emits_records) {
    LosslessGraph g;
    g.add_node(10);
    g.add_node(4);

    std::vector<std::string> seqs = {
        "ACGTACGTAC",  // 10 bp
        "TTGA",        // 4 bp
    };
    std::vector<std::string> names = {"read0", "read1"};

    std::ostringstream oss;
    ASSERT_TRUE(write_fasta(g, seqs, names, oss));

    const std::string fa = oss.str();
    // Two header lines.
    auto count_char = [&](char c) {
        return static_cast<std::size_t>(std::count(fa.begin(), fa.end(), c));
    };
    EXPECT_EQ(count_char('>'), 2u);
    EXPECT_NE(fa.find(">read0"), std::string::npos);
    EXPECT_NE(fa.find(">read1"), std::string::npos);
    // Sequence content preserved.
    EXPECT_NE(fa.find("ACGTACGTAC"), std::string::npos);
    EXPECT_NE(fa.find("TTGA"), std::string::npos);
}

TEST(GraphIOTest, FASTA_writer_wraps_at_80_columns) {
    LosslessGraph g;
    g.add_node(200);
    // 200 'A' characters → must wrap to exactly 3 lines: 80 + 80 + 40.
    std::vector<std::string> seqs = {std::string(200, 'A')};
    std::vector<std::string> names = {"long"};

    std::ostringstream oss;
    ASSERT_TRUE(write_fasta(g, seqs, names, oss));

    const std::string fa = oss.str();
    // Header + 3 sequence lines = 4 newlines.
    EXPECT_EQ(std::count(fa.begin(), fa.end(), '\n'), 4);

    // No line (other than the header) should exceed 80 chars.
    std::istringstream iss(fa);
    std::string line;
    bool first = true;
    while (std::getline(iss, line)) {
        if (first) { first = false; continue; }
        EXPECT_LE(line.size(), 80u);
    }
}

TEST(GraphIOTest, FASTA_writer_rejects_mismatched_sizes) {
    LosslessGraph g;
    g.add_node(4);
    std::vector<std::string> seqs = {"ACGT", "TTTT"};
    std::vector<std::string> names = {"only_one"};
    std::ostringstream oss;
    EXPECT_FALSE(write_fasta(g, seqs, names, oss));
}

// ---------------------------------------------------------------------
// BED writer
// ---------------------------------------------------------------------

TEST(GraphIOTest, BED_writer_emits_lines) {
    LosslessGraph g;
    g.add_node(1000, 1);
    g.add_node(500, 3);
    g.add_node(250, 2);

    std::ostringstream oss;
    ASSERT_TRUE(write_bed(g, oss));

    const std::string bed = oss.str();
    // 3 newline-terminated lines.
    EXPECT_EQ(std::count(bed.begin(), bed.end(), '\n'), 3);

    // Walk each line and check fields.
    std::istringstream iss(bed);
    std::string line;
    std::vector<std::vector<std::string>> rows;
    while (std::getline(iss, line)) {
        std::vector<std::string> cols;
        std::string col;
        std::istringstream ls(line);
        while (std::getline(ls, col, '\t')) cols.push_back(col);
        rows.push_back(std::move(cols));
    }
    ASSERT_EQ(rows.size(), 3u);
    for (const auto& row : rows) {
        ASSERT_EQ(row.size(), 6u);
        EXPECT_EQ(row[0], "NA");
        EXPECT_EQ(row[1], "0");
        EXPECT_EQ(row[5], ".");
    }
    EXPECT_EQ(rows[0][2], "1000");
    EXPECT_EQ(rows[0][3], "node_0");
    EXPECT_EQ(rows[0][4], "1");
    EXPECT_EQ(rows[1][2], "500");
    EXPECT_EQ(rows[1][3], "node_1");
    EXPECT_EQ(rows[1][4], "3");
    EXPECT_EQ(rows[2][2], "250");
    EXPECT_EQ(rows[2][3], "node_2");
    EXPECT_EQ(rows[2][4], "2");
}

// ---------------------------------------------------------------------
// PAF writer
// ---------------------------------------------------------------------

namespace {

// Minimal split-by-tab parser for PAF-12 records.
std::vector<std::vector<std::string>> parse_paf(const std::string& text) {
    std::vector<std::vector<std::string>> rows;
    std::istringstream iss(text);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.empty()) continue;
        std::vector<std::string> cols;
        std::string col;
        std::istringstream ls(line);
        while (std::getline(ls, col, '\t')) cols.push_back(col);
        rows.push_back(std::move(cols));
    }
    return rows;
}

}  // namespace

TEST(GraphIOTest, PAF_writer_roundtrip) {
    std::vector<std::string> seqs = {
        std::string(1000, 'A'),   // read_a index 0
        std::string(800, 'C'),    // read_b index 1
        std::string(1200, 'G'),   // read index 2
    };
    std::vector<std::string> names = {"readA", "readB", "readC"};

    std::vector<branch::backend::OverlapPair> pairs;
    pairs.push_back(branch::backend::OverlapPair{
        .read_a = 0,
        .read_b = 1,
        .offset_a = 100,
        .offset_b = 0,
        .overlap_len = 500,
        .diff_count = 25,
        .strand = 0,
        ._pad = 0,
    });
    pairs.push_back(branch::backend::OverlapPair{
        .read_a = 2,
        .read_b = 0,
        .offset_a = 0,
        .offset_b = 200,
        .overlap_len = 600,
        .diff_count = 0,
        .strand = 1,
        ._pad = 0,
    });

    std::ostringstream oss;
    ASSERT_TRUE(write_paf(pairs, seqs, names, oss));

    const auto rows = parse_paf(oss.str());
    ASSERT_EQ(rows.size(), 2u);
    for (const auto& r : rows) ASSERT_EQ(r.size(), 12u);

    // Row 0: readA vs readB.
    EXPECT_EQ(rows[0][0], "readA");            // qname
    EXPECT_EQ(rows[0][1], "1000");             // qlen
    EXPECT_EQ(rows[0][2], "100");              // qstart
    EXPECT_EQ(rows[0][3], "600");              // qend = 100 + 500
    EXPECT_EQ(rows[0][4], "+");                // strand
    EXPECT_EQ(rows[0][5], "readB");            // tname
    EXPECT_EQ(rows[0][6], "800");              // tlen
    EXPECT_EQ(rows[0][7], "0");                // tstart
    EXPECT_EQ(rows[0][8], "500");              // tend = 0 + 500
    EXPECT_EQ(rows[0][9], "475");              // matches = 500 - 25
    EXPECT_EQ(rows[0][10], "500");             // alnlen
    // mapq = int(475/500 * 60) = 57.
    EXPECT_EQ(rows[0][11], "57");

    // Row 1: readC vs readA, opposite strand, perfect match.
    EXPECT_EQ(rows[1][0], "readC");
    EXPECT_EQ(rows[1][1], "1200");
    EXPECT_EQ(rows[1][2], "0");
    EXPECT_EQ(rows[1][3], "600");
    EXPECT_EQ(rows[1][4], "-");
    EXPECT_EQ(rows[1][5], "readA");
    EXPECT_EQ(rows[1][6], "1000");
    EXPECT_EQ(rows[1][7], "200");
    EXPECT_EQ(rows[1][8], "800");
    EXPECT_EQ(rows[1][9], "600");              // matches = 600 - 0
    EXPECT_EQ(rows[1][10], "600");             // alnlen
    EXPECT_EQ(rows[1][11], "60");              // mapq clamped to 60
}

TEST(GraphIOTest, PAF_writer_rejects_out_of_range_read_id) {
    std::vector<std::string> seqs = {"ACGT"};
    std::vector<std::string> names = {"only"};
    std::vector<branch::backend::OverlapPair> pairs;
    pairs.push_back(branch::backend::OverlapPair{
        .read_a = 0,
        .read_b = 5,            // out of range
        .offset_a = 0,
        .offset_b = 0,
        .overlap_len = 4,
        .diff_count = 0,
        .strand = 0,
        ._pad = 0,
    });
    std::ostringstream oss;
    EXPECT_FALSE(write_paf(pairs, seqs, names, oss));
}
