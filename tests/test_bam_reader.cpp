// Tests for branch::io::BamReader.
//
// We build a 3-record BAM at test-time using htslib directly, then verify
// the reader yields the sequences in the expected orientation:
//   record 0 — unmapped read, sequence returned as-is.
//   record 1 — mapped forward-strand read, sequence returned as-is.
//   record 2 — mapped reverse-strand read, sequence returned
//              reverse-complemented (i.e. restored to original
//              sequencing orientation).
//
// The fixture also includes:
//   record 3 — a secondary alignment of record 1. Must be skipped.
//   record 4 — a supplementary alignment of record 1. Must be skipped.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include "io/bam_reader.hpp"

namespace {

// Tiny helper: write one record via bam_set1 + sam_write1.
void write_record(samFile* fp, sam_hdr_t* hdr,
                  const char* qname, std::uint16_t flag,
                  std::int32_t tid, hts_pos_t pos, std::uint8_t mapq,
                  const std::string& seq) {
    bam1_t* b = bam_init1();
    ASSERT_NE(b, nullptr);

    // For mapped reads we need a minimal CIGAR = l_seq M so htslib
    // accepts the record. For unmapped reads we pass no CIGAR.
    std::vector<std::uint32_t> cigar;
    std::size_t n_cigar = 0;
    if ((flag & BAM_FUNMAP) == 0) {
        cigar.push_back(static_cast<std::uint32_t>(seq.size()) << 4 | 0u);  // M
        n_cigar = 1;
    }

    const int rc = bam_set1(
        b,
        std::strlen(qname), qname,
        flag,
        tid, pos, mapq,
        n_cigar, cigar.empty() ? nullptr : cigar.data(),
        /*mtid=*/-1, /*mpos=*/-1, /*isize=*/0,
        seq.size(), seq.c_str(), /*qual=*/nullptr,
        /*l_aux=*/0);
    ASSERT_GE(rc, 0);

    const int wrc = sam_write1(fp, hdr, b);
    ASSERT_GE(wrc, 0);
    bam_destroy1(b);
}

// Build a 5-record BAM at `path`. Returns the three canonical sequences
// (unmapped, fwd, rev) in the orientation we expect to read them back.
struct Fixture {
    std::string unmapped_seq;
    std::string fwd_seq;
    std::string rev_recorded_seq;  // as stored in the BAM (reference-forward)
    std::string rev_expected_seq;  // what BamReader should yield (revcomp)
};

Fixture build_fixture(const std::string& path) {
    Fixture f;
    f.unmapped_seq   = "ACGTACGTACGTACGTAAAA";     // 20 bp
    f.fwd_seq        = "GGGGCCCCAAAATTTTGGGG";     // 20 bp
    f.rev_recorded_seq = "ACACACACGTGTGTGTACGT";   // what's in the BAM
    // Expected = reverse-complement of rev_recorded_seq.
    // ACACACACGTGTGTGTACGT -> revcomp -> ACGTACACACACGTGTGTGT
    f.rev_expected_seq = "ACGTACACACACGTGTGTGT";

    samFile* fp = sam_open(path.c_str(), "wb");
    EXPECT_NE(fp, nullptr);
    sam_hdr_t* hdr = sam_hdr_init();
    EXPECT_NE(hdr, nullptr);
    EXPECT_EQ(sam_hdr_add_line(hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL), 0);
    EXPECT_EQ(sam_hdr_add_line(hdr, "SQ", "SN", "chr1", "LN", "1000", NULL), 0);
    EXPECT_EQ(sam_hdr_write(fp, hdr), 0);

    // Record 0: unmapped read. BAM stores seq as-provided.
    write_record(fp, hdr, "read_unmapped", BAM_FUNMAP,
                 /*tid=*/-1, /*pos=*/-1, /*mapq=*/0, f.unmapped_seq);

    // Record 1: mapped, forward strand. Stored as-is.
    write_record(fp, hdr, "read_fwd", /*flag=*/0,
                 /*tid=*/0, /*pos=*/100, /*mapq=*/60, f.fwd_seq);

    // Record 2: mapped, reverse strand. BAM stores reference-forward seq;
    // the reader must revcomp it back to original orientation.
    write_record(fp, hdr, "read_rev", BAM_FREVERSE,
                 /*tid=*/0, /*pos=*/200, /*mapq=*/60, f.rev_recorded_seq);

    // Record 3: secondary alignment of read_fwd. Must be skipped.
    write_record(fp, hdr, "read_fwd", BAM_FSECONDARY,
                 /*tid=*/0, /*pos=*/300, /*mapq=*/0, f.fwd_seq);

    // Record 4: supplementary alignment of read_fwd. Must be skipped.
    write_record(fp, hdr, "read_fwd", BAM_FSUPPLEMENTARY,
                 /*tid=*/0, /*pos=*/400, /*mapq=*/0, f.fwd_seq);

    sam_hdr_destroy(hdr);
    EXPECT_EQ(sam_close(fp), 0);
    return f;
}

std::string make_tmp_bam_path() {
    namespace fs = std::filesystem;
    fs::path p = fs::temp_directory_path() /
                 ("branch_bam_reader_test_" +
                  std::to_string(::getpid()) + ".bam");
    return p.string();
}

}  // namespace

TEST(BamReader, YieldsThreeReads_SkipsSecondaryAndSupplementary_RevCompsReverseStrand) {
    const std::string path = make_tmp_bam_path();
    std::filesystem::remove(path);
    const Fixture fx = build_fixture(path);

    branch::io::BamReader r;
    ASSERT_TRUE(r.open(path));

    std::string name, seq;

    ASSERT_TRUE(r.next(name, seq));
    EXPECT_EQ(name, "read_unmapped");
    EXPECT_EQ(seq, fx.unmapped_seq);

    ASSERT_TRUE(r.next(name, seq));
    EXPECT_EQ(name, "read_fwd");
    EXPECT_EQ(seq, fx.fwd_seq);

    ASSERT_TRUE(r.next(name, seq));
    EXPECT_EQ(name, "read_rev");
    EXPECT_EQ(seq, fx.rev_expected_seq)
        << "reverse-strand read must be reverse-complemented";

    // No further records: secondary + supplementary filtered out.
    EXPECT_FALSE(r.next(name, seq));
    EXPECT_EQ(r.reads_read(), 3u);

    r.close();
    std::filesystem::remove(path);
}

TEST(BamReader, OpenReturnsFalseForMissingFile) {
    branch::io::BamReader r;
    EXPECT_FALSE(r.open("/nonexistent/path/does_not_exist.bam"));

    // Safe to call next() on a never-opened reader.
    std::string name, seq;
    EXPECT_FALSE(r.next(name, seq));
}
