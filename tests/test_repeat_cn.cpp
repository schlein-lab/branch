#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "analysis/repeat_cn.hpp"

void test_parse_repeat_bed() {
    // Create synthetic BED file
    std::ofstream bed("/tmp/test_repeats.bed");
    bed << "chr1\t100\t200\tLINE/L1\n";
    bed << "chr1\t300\t400\tSINE/Alu\n";
    bed << "chr1\t500\t600\tLINE/L1\n";
    bed << "chr2\t100\t150\tSINE/Alu\n";
    bed.close();

    auto repeats = branch::analysis::parse_repeat_bed("/tmp/test_repeats.bed");
    assert(repeats.size() == 4);
    assert(repeats[0].family == "L1");
    assert(repeats[0].class_ == "LINE");
    assert(repeats[1].family == "Alu");
    assert(repeats[1].class_ == "SINE");
    std::cout << "[PASS] test_parse_repeat_bed\n";
}

void test_compute_family_cn() {
    // Synthetic repeats: 2 families
    std::vector<branch::analysis::RepeatEntry> repeats = {
        {"L1", "LINE", "chr1", 100, 200},  // 100bp
        {"L1", "LINE", "chr1", 300, 400},  // 100bp -> total 200bp L1
        {"Alu", "SINE", "chr1", 500, 600}, // 100bp Alu
    };

    // Synthetic coverage: baseline=30, L1 regions have cov=60 (CN=4), Alu has cov=15 (CN=1)
    std::map<std::string, std::vector<uint32_t>> per_base_cov;
    per_base_cov["chr1"].resize(700, 0);
    for (size_t i = 100; i < 200; ++i) per_base_cov["chr1"][i] = 60;
    for (size_t i = 300; i < 400; ++i) per_base_cov["chr1"][i] = 60;
    for (size_t i = 500; i < 600; ++i) per_base_cov["chr1"][i] = 15;

    double baseline = 30.0;
    auto stats = branch::analysis::compute_family_cn(repeats, per_base_cov, baseline);

    // L1: mean_cov=60, cn=60/30*2=4.0
    assert(stats.count("L1") == 1);
    assert(stats["L1"].total_bp == 200);
    assert(std::abs(stats["L1"].mean_cov - 60.0) < 0.01);
    assert(std::abs(stats["L1"].cn_estimate - 4.0) < 0.01);

    // Alu: mean_cov=15, cn=15/30*2=1.0
    assert(stats.count("Alu") == 1);
    assert(stats["Alu"].total_bp == 100);
    assert(std::abs(stats["Alu"].mean_cov - 15.0) < 0.01);
    assert(std::abs(stats["Alu"].cn_estimate - 1.0) < 0.01);

    std::cout << "[PASS] test_compute_family_cn\n";
}

void test_write_cn_tsv() {
    std::map<std::string, branch::analysis::CnStat> stats;
    stats["L1"] = {200, 60.0, 4.0};
    stats["Alu"] = {100, 15.0, 1.0};

    branch::analysis::write_cn_tsv(stats, "/tmp/test_cn_out.tsv");

    std::ifstream in("/tmp/test_cn_out.tsv");
    std::string line;
    std::getline(in, line); // header
    assert(line.find("family") != std::string::npos);
    int count = 0;
    while (std::getline(in, line)) count++;
    assert(count == 2);
    std::cout << "[PASS] test_write_cn_tsv\n";
}

int main() {
    test_parse_repeat_bed();
    test_compute_family_cn();
    test_write_cn_tsv();
    std::cout << "All repeat_cn tests passed!\n";
    return 0;
}
