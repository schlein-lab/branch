// Unit tests for write_fasta_consensus — verifies node consensus is
// written out as FASTA records, and consensus-less nodes are skipped.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "graph/graph_io.hpp"
#include "graph/lossless_graph.hpp"

namespace {

std::string read_file(const std::string& path) {
    std::ifstream f(path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(WriteFastaConsensus, SkipsNodesWithoutConsensus) {
    branch::graph::LosslessGraph g;
    g.add_node(500);  // no consensus
    g.add_node(300);  // no consensus

    auto tmp = std::filesystem::temp_directory_path() /
        "fasta_consensus_empty_test.fa";
    ASSERT_TRUE(branch::graph::write_fasta_consensus(g, tmp.string()));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_TRUE(content.empty())
        << "Expected empty FASTA when no node has consensus, got: " << content;
}

TEST(WriteFastaConsensus, WritesOneRecordPerNodeWithConsensus) {
    branch::graph::LosslessGraph g;
    auto id0 = g.add_node(6);
    auto id1 = g.add_node(4);
    g.node(id0).consensus = "ACGTAC";
    g.node(id1).consensus = "TTGA";

    auto tmp = std::filesystem::temp_directory_path() /
        "fasta_consensus_two_test.fa";
    ASSERT_TRUE(branch::graph::write_fasta_consensus(g, tmp.string()));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find(">node_0"), std::string::npos);
    EXPECT_NE(content.find("ACGTAC"), std::string::npos);
    EXPECT_NE(content.find(">node_1"), std::string::npos);
    EXPECT_NE(content.find("TTGA"), std::string::npos);
}

TEST(WriteFastaConsensus, WrapsSequenceAt80Columns) {
    branch::graph::LosslessGraph g;
    auto id = g.add_node(200);
    g.node(id).consensus = std::string(200, 'A');

    auto tmp = std::filesystem::temp_directory_path() /
        "fasta_consensus_wrap_test.fa";
    ASSERT_TRUE(branch::graph::write_fasta_consensus(g, tmp.string()));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    // Find first line after header '>node_0\n'
    auto first_nl = content.find('\n');
    ASSERT_NE(first_nl, std::string::npos);
    auto second_nl = content.find('\n', first_nl + 1);
    ASSERT_NE(second_nl, std::string::npos);

    // The first sequence line should be exactly 80 chars.
    std::size_t first_seq_len = second_nl - first_nl - 1;
    EXPECT_EQ(first_seq_len, 80u);
}

TEST(WriteFastaConsensus, SkipsSomeNodesKeepsOrder) {
    branch::graph::LosslessGraph g;
    auto id0 = g.add_node(3);
    /*id1*/ g.add_node(3);          // no consensus, should be skipped
    auto id2 = g.add_node(3);
    g.node(id0).consensus = "AAA";
    g.node(id2).consensus = "CCC";

    auto tmp = std::filesystem::temp_directory_path() /
        "fasta_consensus_gaps_test.fa";
    ASSERT_TRUE(branch::graph::write_fasta_consensus(g, tmp.string()));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    // node_1 must not appear.
    EXPECT_EQ(content.find("node_1"), std::string::npos);
    // node_0 and node_2 both appear, in ID order.
    auto p0 = content.find(">node_0");
    auto p2 = content.find(">node_2");
    ASSERT_NE(p0, std::string::npos);
    ASSERT_NE(p2, std::string::npos);
    EXPECT_LT(p0, p2);
}
