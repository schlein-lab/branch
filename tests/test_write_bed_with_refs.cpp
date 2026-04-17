// Unit test for write_bed_with_refs — verifies chrom=NA fallback when
// refs is empty and when a consensus-less node is encountered.

#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "graph/graph_io.hpp"
#include "graph/lossless_graph.hpp"

namespace {

bool minimap2_available() {
    return std::system("which minimap2 > /dev/null 2>&1") == 0;
}

std::string read_file(const std::string& path) {
    std::ifstream f(path);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

TEST(WriteBedWithRefs, EmptyRefsFallsBackToNA) {
    branch::graph::LosslessGraph g;
    auto id = g.add_node(500);
    (void)id;

    auto tmp = std::filesystem::temp_directory_path() /
        "bed_emptyrefs_test.bed";
    std::vector<branch::graph::BedLinearRef> refs;

    ASSERT_TRUE(branch::graph::write_bed_with_refs(
        g, refs, tmp.string(), 1));

    std::string content = read_file(tmp.string());
    std::filesystem::remove(tmp);

    EXPECT_NE(content.find("NA"), std::string::npos);
    EXPECT_NE(content.find("node_0"), std::string::npos);
}

TEST(WriteBedWithRefs, NodeWithoutConsensusFallsBack) {
    if (!minimap2_available()) {
        GTEST_SKIP() << "minimap2 not in PATH";
    }

    branch::graph::LosslessGraph g;
    g.add_node(500);  // no consensus set

    auto tmp_ref = std::filesystem::temp_directory_path() /
        "bed_test_ref.fa";
    {
        std::ofstream r(tmp_ref);
        r << ">chr1\n"
          << std::string(1000, 'A') << "\n";
    }

    auto tmp_bed = std::filesystem::temp_directory_path() /
        "bed_consensusless_test.bed";
    std::vector<branch::graph::BedLinearRef> refs = {
        {"TestRef", tmp_ref.string()}};

    ASSERT_TRUE(branch::graph::write_bed_with_refs(
        g, refs, tmp_bed.string(), 1));

    std::string content = read_file(tmp_bed.string());
    std::filesystem::remove(tmp_bed);
    std::filesystem::remove(tmp_ref);

    EXPECT_NE(content.find("NA"), std::string::npos);
    EXPECT_NE(content.find("node_0"), std::string::npos);
}
