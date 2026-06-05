// Unit test for FastqIndex — round-trip sequential and random access.
//
// Builds a small synthetic FASTQ on disk, verifies:
//   1. n_reads matches input
//   2. read_seq(rid) returns the expected sequence
//   3. for_each visits all reads in order with matching sequences
//   4. Same for a gzipped input (via decompress-to-scratch)

#include "wg/phaser/fastq_index.hpp"

#include <zlib.h>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace branch::wg::phaser;

namespace {

const char* kFastqContent =
    "@read0\nACGTACGTAC\n+\n!!!!!!!!!!\n"
    "@read1\nCCCCAAAATTTTGGGG\n+\n????????????????\n"
    "@read2 desc with spaces\nNNNN\n+\n!!!!\n"
    "@read3\nAAAAAAAAAAAAAAAAAAAA\n+\n....................\n";

const std::vector<std::string> kExpectedSeqs = {
    "ACGTACGTAC",
    "CCCCAAAATTTTGGGG",
    "NNNN",
    "AAAAAAAAAAAAAAAAAAAA",
};

// Resolve a writable scratch directory. /tmp is often tiny/full on HPC login
// and compute nodes (our site caps it at 50 MiB), so honour $TMPDIR first.
std::string tmp_path(const std::string& name) {
    if (const char* p = std::getenv("TMPDIR")) {
        if (*p) return std::string(p) + "/" + name;
    }
    return "/tmp/" + name;
}

void write_plain(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    assert(f);
    std::fwrite(kFastqContent, 1, std::strlen(kFastqContent), f);
    std::fclose(f);
}

void write_gz(const std::string& path) {
    gzFile f = gzopen(path.c_str(), "wb");
    assert(f);
    int n = static_cast<int>(std::strlen(kFastqContent));
    int wrote = gzwrite(f, kFastqContent, n);
    assert(wrote == n);
    gzclose(f);
}

void check_roundtrip(FastqIndex& idx) {
    assert(idx.n_reads() == kExpectedSeqs.size());

    // Random access in arbitrary order.
    std::string s;
    idx.read_seq(2, s); assert(s == kExpectedSeqs[2]);
    idx.read_seq(0, s); assert(s == kExpectedSeqs[0]);
    idx.read_seq(3, s); assert(s == kExpectedSeqs[3]);
    idx.read_seq(1, s); assert(s == kExpectedSeqs[1]);

    // Sequential streaming.
    std::vector<std::string> got;
    idx.for_each([&](ReadId rid, const std::string& seq) {
        (void)rid;
        got.push_back(seq);
    });
    assert(got.size() == kExpectedSeqs.size());
    for (std::size_t i = 0; i < got.size(); ++i) {
        assert(got[i] == kExpectedSeqs[i]);
    }

    // ID strings: "@read2 desc with spaces" should yield "read2".
    const auto& ids = idx.read_id_strings();
    assert(ids.size() == kExpectedSeqs.size());
    assert(ids[0] == "read0");
    assert(ids[2] == "read2");
}

void test_plain_fastq() {
    std::string p = tmp_path("test_fastq_index_plain.fq");
    write_plain(p);
    FastqIndex idx;
    idx.build(p);
    check_roundtrip(idx);
    std::remove(p.c_str());
}

void test_gz_fastq() {
    std::string p = tmp_path("test_fastq_index_gz.fq.gz");
    write_gz(p);
    FastqIndex idx;
    idx.build(p);
    check_roundtrip(idx);
    std::remove(p.c_str());
}

}  // namespace

int main() {
    test_plain_fastq();
    std::cout << "plain ok\n";
    test_gz_fastq();
    std::cout << "gz ok\n";
    std::cout << "test_fastq_index: OK\n";
    return 0;
}
