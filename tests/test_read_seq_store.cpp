// Tests for branch::io::ReadSeqStoreWriter / ReadSeqStoreReader.
//
// The two properties worth the most scrutiny: open() must stay cheap
// regardless of how much sequence data is behind it (verified via
// /proc/self/status VmRSS, not just record counts), and fetch_window's
// clip-vs-throw split at the sequence boundary has to land exactly
// where documented in read_seq_store.hpp, not just "somewhere near it".

#include "io/read_seq_store.hpp"

#include <gtest/gtest.h>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using branch::io::ReadSeqStoreReader;
using branch::io::ReadSeqStoreWriter;

namespace {

// Returns VmRSS in kB from /proc/self/status, or -1 if unavailable
// (e.g. non-Linux) so the caller can skip rather than fail outright.
long vm_rss_kb() {
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.compare(0, 6, "VmRSS:") == 0) {
            std::istringstream iss(line.substr(6));
            long kb = -1;
            iss >> kb;
            return kb;
        }
    }
    return -1;
}

}  // namespace

// --- roundtrip -----------------------------------------------------------

TEST(ReadSeqStoreTest, roundtrip_basic_acgt_reads) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_basic.rsq";
    const std::vector<std::string> seqs = {
        "ACGT", "TTTTTTTT", "A", "ACGTACGTACGTACGTACGT", "GATTACA",
    };

    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    for (const auto& s : seqs) w.add(s);
    ASSERT_TRUE(w.close());
    EXPECT_EQ(w.reads_written(), static_cast<std::uint32_t>(seqs.size()));

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    ASSERT_EQ(r.read_count(), static_cast<std::uint32_t>(seqs.size()));
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        const auto id = static_cast<std::uint32_t>(i);
        EXPECT_EQ(r.sequence_length(id), static_cast<std::uint32_t>(seqs[i].size()));
        EXPECT_EQ(r.fetch(id), seqs[i]);
    }

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, roundtrip_n_containing_reads) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_n.rsq";
    const std::vector<std::string> seqs = {
        "ACGTNACGT",  // single N in the middle
        "NNNNN",      // all exceptions
        "N",          // one-base exception read
        "ACGTN",      // trailing N
        "NACGT",      // leading N
    };

    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    for (const auto& s : seqs) w.add(s);
    ASSERT_TRUE(w.close());

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    for (std::size_t i = 0; i < seqs.size(); ++i) {
        const auto id = static_cast<std::uint32_t>(i);
        EXPECT_EQ(r.fetch(id), seqs[i]) << "mismatch at read " << i;
    }

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, roundtrip_empty_sequence) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_empty.rsq";

    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.add("");            // read 0: empty
    w.add("ACGT");         // read 1: non-empty, brackets the empty one
    w.add("");            // read 2: empty again
    ASSERT_TRUE(w.close());
    EXPECT_EQ(w.reads_written(), 3u);

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_EQ(r.sequence_length(0), 0u);
    EXPECT_EQ(r.fetch(0), "");
    EXPECT_EQ(r.fetch(1), "ACGT");
    EXPECT_EQ(r.sequence_length(2), 0u);
    EXPECT_EQ(r.fetch(2), "");

    std::filesystem::remove(path);
}

// --- fetch_window edge cases ---------------------------------------------
//
// Documented behavior (read_seq_store.hpp): start > length throws
// std::out_of_range; start <= length with start+len past the end clips
// to whatever remains instead of throwing.

TEST(ReadSeqStoreTest, fetch_window_normal_range) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_window_normal.rsq";
    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.add("ACGTACGTACGT");  // length 12
    ASSERT_TRUE(w.close());

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_EQ(r.fetch_window(0, 0, 4), "ACGT");
    EXPECT_EQ(r.fetch_window(0, 4, 4), "ACGT");
    EXPECT_EQ(r.fetch_window(0, 2, 5), "GTACG");
    EXPECT_EQ(r.fetch_window(0, 0, 12), "ACGTACGTACGT");

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, fetch_window_clips_when_len_runs_past_end) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_window_clip.rsq";
    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.add("ACGTACGT");  // length 8
    ASSERT_TRUE(w.close());

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_EQ(r.fetch_window(0, 6, 100), "GT");   // clipped from 100 to 2
    EXPECT_EQ(r.fetch_window(0, 8, 5), "");        // start == length: empty, no throw
    EXPECT_EQ(r.fetch_window(0, 0, 8), "ACGTACGT");  // exact length: no clip needed

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, fetch_window_throws_when_start_past_end) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_window_throw.rsq";
    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.add("ACGT");  // length 4
    ASSERT_TRUE(w.close());

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_THROW(r.fetch_window(0, 5, 1), std::out_of_range);
    EXPECT_THROW(r.fetch_window(0, 1000, 0), std::out_of_range);

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, fetch_throws_on_invalid_read_id) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_bad_id.rsq";
    ReadSeqStoreWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.add("ACGT");
    ASSERT_TRUE(w.close());

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_THROW(r.fetch(1), std::out_of_range);
    EXPECT_THROW(r.fetch_window(1, 0, 1), std::out_of_range);
    EXPECT_THROW((void)r.sequence_length(1), std::out_of_range);

    std::filesystem::remove(path);
}

// --- scale: 100k reads -----------------------------------------------------

TEST(ReadSeqStoreTest, roundtrip_100k_random_access_sample) {
    constexpr std::uint32_t kCount = 100000;
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_100k_sample.rsq";

    // Deterministic per-read content: read r's sequence is length
    // (20 + r % 30) and base i is chosen by (r + i) % 4, with every
    // 13th read carrying a single N near its midpoint so the sample
    // below also touches the exception path at scale.
    auto make_seq = [](std::uint32_t r) {
        static constexpr char kBases[4] = {'A', 'C', 'G', 'T'};
        const std::uint32_t len = 20 + (r % 30);
        std::string s(len, 'A');
        for (std::uint32_t i = 0; i < len; ++i) s[i] = kBases[(r + i) % 4];
        if (r % 13 == 0) s[len / 2] = 'N';
        return s;
    };

    {
        ReadSeqStoreWriter w;
        ASSERT_TRUE(w.open(path.string()));
        for (std::uint32_t r = 0; r < kCount; ++r) w.add(make_seq(r));
        ASSERT_TRUE(w.close());
        ASSERT_EQ(w.reads_written(), kCount);
    }

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    ASSERT_EQ(r.read_count(), kCount);

    // Prime stride: touches a scattered, non-periodic sample without
    // fetching all 100k reads.
    std::uint32_t checked = 0;
    for (std::uint32_t id = 0; id < kCount; id += 997) {
        const std::string expected = make_seq(id);
        EXPECT_EQ(r.sequence_length(id), static_cast<std::uint32_t>(expected.size()));
        EXPECT_EQ(r.fetch(id), expected) << "mismatch at read_id " << id;
        if (expected.size() >= 6) {
            EXPECT_EQ(r.fetch_window(id, 2, 4), expected.substr(2, 4))
                << "window mismatch at read_id " << id;
        }
        ++checked;
    }
    EXPECT_GT(checked, 0u);
    // First and last reads explicitly, as boundary checks.
    EXPECT_EQ(r.fetch(0), make_seq(0));
    EXPECT_EQ(r.fetch(kCount - 1), make_seq(kCount - 1));

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, index_ram_check_100k_reads) {
    constexpr std::uint32_t kCount = 100000;
    constexpr std::uint32_t kReadLen = 1000;
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_ram_check.rsq";

    std::string seq(kReadLen, 'A');
    for (std::uint32_t i = 0; i < kReadLen; ++i) {
        static constexpr char kBases[4] = {'A', 'C', 'G', 'T'};
        seq[i] = kBases[i % 4];
    }

    {
        ReadSeqStoreWriter w;
        ASSERT_TRUE(w.open(path.string()));
        for (std::uint32_t r = 0; r < kCount; ++r) w.add(seq);
        ASSERT_TRUE(w.close());
    }

    const long rss_before = vm_rss_kb();
    if (rss_before < 0) GTEST_SKIP() << "VmRSS not available from /proc/self/status here";

    ReadSeqStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    ASSERT_EQ(r.read_count(), kCount);

    const long rss_after = vm_rss_kb();
    ASSERT_GE(rss_after, 0);
    const long grew_kb = rss_after - rss_before;

    // Theoretical minimum is 16B/read (sizeof(IndexEntry) once loaded,
    // per the u64-alignment padding documented in read_seq_store.hpp) --
    // ~1.6MB for 100k reads. The real regression this guards against is
    // open() accidentally loading the packed sequence data too, which
    // for this fixture would be ~24MB (kCount * kReadLen / 4 bytes) --
    // so comparing against a quarter of that leaves generous allocator
    // slack while still catching that regression clearly.
    const std::uint64_t packed_data_bytes = (static_cast<std::uint64_t>(kCount) * kReadLen) / 4;
    const auto packed_data_kb = static_cast<long>(packed_data_bytes / 1024);
    EXPECT_LT(grew_kb, packed_data_kb / 4)
        << "VmRSS grew " << grew_kb << " kB opening a store whose full packed sequence data is "
        << packed_data_kb << " kB -- open() should load only the index, not sequence data";

    std::filesystem::remove(path);
}

// --- error handling --------------------------------------------------------

TEST(ReadSeqStoreTest, reader_open_rejects_missing_file) {
    ReadSeqStoreReader r;
    EXPECT_FALSE(r.open("/nonexistent/path/that/should/not/exist.rsq"));
    EXPECT_FALSE(r.last_error().empty());
}

TEST(ReadSeqStoreTest, reader_open_rejects_bad_magic) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_bad_magic.rsq";
    {
        // Deliberately >> kFooterSize (37B) so this exercises the bad-magic
        // path specifically, not the too-small-to-contain-a-footer path --
        // a fixed repeat count keeps that margin obvious rather than
        // resting on a hand-counted string literal.
        std::ofstream f(path, std::ios::binary);
        f << std::string(64, 'X');
    }

    ReadSeqStoreReader r;
    EXPECT_FALSE(r.open(path.string()));
    EXPECT_NE(r.last_error().find("magic"), std::string::npos) << r.last_error();

    std::filesystem::remove(path);
}

TEST(ReadSeqStoreTest, reader_open_rejects_truncated_file) {
    auto path = std::filesystem::temp_directory_path() / "read_seq_store_truncated.rsq";
    {
        ReadSeqStoreWriter w;
        ASSERT_TRUE(w.open(path.string()));
        w.add("ACGTACGTACGT");
        w.add("TTTTGGGGCCCC");
        ASSERT_TRUE(w.close());
    }
    // Chop off the last few bytes so the footer itself is incomplete.
    const auto full_size = std::filesystem::file_size(path);
    ASSERT_GT(full_size, 10u);
    std::filesystem::resize_file(path, full_size - 10);

    ReadSeqStoreReader r;
    EXPECT_FALSE(r.open(path.string()));
    EXPECT_FALSE(r.last_error().empty());

    std::filesystem::remove(path);
}
