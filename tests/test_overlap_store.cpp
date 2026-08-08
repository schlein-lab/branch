// Tests for branch::io::OverlapSpillWriter / merge_spills / OverlapStoreReader.
//
// merge_spills is the part worth the most scrutiny: it has two entirely
// different code paths (a single in-RAM sort vs. an external chunked
// k-way merge) that must agree byte-for-byte, and a dedup rule (keep
// the larger overlap_len per (read_a, read_b)) that has to survive both
// a single-chunk sort and duplicates split across chunk boundaries.

#include "io/overlap_store.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using branch::io::merge_spills;
using branch::io::OverlapPair;
using branch::io::OverlapSpillWriter;
using branch::io::OverlapStoreReader;

namespace {

OverlapPair make_overlap(std::uint32_t read_a, std::uint32_t read_b, std::uint32_t offset_a,
                          std::uint32_t offset_b, std::uint32_t overlap_len,
                          std::int16_t diff_count = 0, std::uint8_t strand = 0) {
    OverlapPair p{};
    p.read_a = read_a;
    p.read_b = read_b;
    p.offset_a = offset_a;
    p.offset_b = offset_b;
    p.overlap_len = overlap_len;
    p.diff_count = diff_count;
    p.strand = strand;
    return p;
}

bool overlap_eq(const OverlapPair& a, const OverlapPair& b) {
    return a.read_a == b.read_a && a.read_b == b.read_b && a.offset_a == b.offset_a &&
           a.offset_b == b.offset_b && a.overlap_len == b.overlap_len &&
           a.diff_count == b.diff_count && a.strand == b.strand;
}

std::string read_whole_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

}  // namespace

// --- OverlapSpillWriter / OverlapStoreReader roundtrip ---------------------

TEST(OverlapStoreTest, roundtrip_zero_records) {
    auto path = std::filesystem::temp_directory_path() / "overlap_store_roundtrip_zero.brov1";

    OverlapSpillWriter w;
    ASSERT_TRUE(w.open(path.string()));
    ASSERT_TRUE(w.close());
    EXPECT_EQ(w.records_written(), 0u);

    OverlapStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_EQ(r.record_count(), 0u);
    OverlapPair rec;
    EXPECT_FALSE(r.next(rec));
    EXPECT_TRUE(r.last_error().empty());

    std::filesystem::remove(path);
}

TEST(OverlapStoreTest, roundtrip_one_record) {
    auto path = std::filesystem::temp_directory_path() / "overlap_store_roundtrip_one.brov1";
    const OverlapPair original = make_overlap(42, 99, 1000, 2000, 350, -7, 1);

    OverlapSpillWriter w;
    ASSERT_TRUE(w.open(path.string()));
    w.append(original);
    ASSERT_TRUE(w.close());
    EXPECT_EQ(w.records_written(), 1u);

    OverlapStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    EXPECT_EQ(r.record_count(), 1u);
    OverlapPair got;
    ASSERT_TRUE(r.next(got));
    EXPECT_TRUE(overlap_eq(got, original));
    EXPECT_FALSE(r.next(got));

    std::filesystem::remove(path);
}

TEST(OverlapStoreTest, roundtrip_one_million_records) {
    constexpr std::uint32_t kCount = 1'000'000;
    auto path = std::filesystem::temp_directory_path() / "overlap_store_roundtrip_1m.brov1";

    auto record_for = [](std::uint32_t i) {
        return make_overlap(i, i + 1, i * 3, i * 5, 100 + (i % 400),
                             static_cast<std::int16_t>(i % 50),
                             static_cast<std::uint8_t>(i % 2));
    };

    OverlapSpillWriter w;
    ASSERT_TRUE(w.open(path.string()));
    for (std::uint32_t i = 0; i < kCount; ++i) w.append(record_for(i));
    ASSERT_TRUE(w.close());
    ASSERT_EQ(w.records_written(), kCount);

    OverlapStoreReader r;
    ASSERT_TRUE(r.open(path.string())) << r.last_error();
    ASSERT_EQ(r.record_count(), kCount);

    std::uint32_t n = 0;
    std::uint32_t first_mismatch = kCount;  // sentinel: "none found"
    OverlapPair rec;
    while (r.next(rec)) {
        if (!overlap_eq(rec, record_for(n)) && first_mismatch == kCount) first_mismatch = n;
        ++n;
    }
    EXPECT_EQ(first_mismatch, kCount) << "first mismatching record index";
    EXPECT_EQ(n, kCount);
    EXPECT_TRUE(r.last_error().empty());

    std::filesystem::remove(path);
}

// --- merge_spills ------------------------------------------------------

TEST(OverlapStoreTest, merge_three_spills_with_duplicates_in_ram) {
    const auto dir = std::filesystem::temp_directory_path();
    const std::vector<std::string> spill_paths = {
        (dir / "overlap_store_merge_spill_a.brov1").string(),
        (dir / "overlap_store_merge_spill_b.brov1").string(),
        (dir / "overlap_store_merge_spill_c.brov1").string(),
    };
    const auto out_path = (dir / "overlap_store_merge_out.brov1").string();
    const auto tmp_dir = (dir / "overlap_store_merge_tmp").string();

    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[0]));
        w.append(make_overlap(1, 2, 0, 0, 50));   // (1,2) dup, smaller overlap_len
        w.append(make_overlap(3, 4, 10, 20, 80));
        ASSERT_TRUE(w.close());
    }
    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[1]));
        w.append(make_overlap(1, 2, 5, 5, 90));   // (1,2) dup, larger -> must win
        w.append(make_overlap(5, 6, 0, 0, 30));
        ASSERT_TRUE(w.close());
    }
    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[2]));
        w.append(make_overlap(2, 9, 0, 0, 15));
        ASSERT_TRUE(w.close());
    }

    ASSERT_NO_THROW(merge_spills(spill_paths, out_path, tmp_dir, std::uint64_t{1} << 30));

    OverlapStoreReader r;
    ASSERT_TRUE(r.open(out_path)) << r.last_error();
    // (1,2) (2,9) (3,4) (5,6): the (1,2) duplicate collapses into one.
    EXPECT_EQ(r.record_count(), 4u);

    std::vector<OverlapPair> got;
    OverlapPair rec;
    while (r.next(rec)) got.push_back(rec);
    ASSERT_EQ(got.size(), 4u);

    for (std::size_t i = 1; i < got.size(); ++i) {
        const auto& p = got[i - 1];
        const auto& c = got[i];
        EXPECT_TRUE(std::tie(p.read_a, p.read_b) < std::tie(c.read_a, c.read_b))
            << "records out of order at index " << i;
    }

    auto it = std::find_if(got.begin(), got.end(), [](const OverlapPair& p) {
        return p.read_a == 1 && p.read_b == 2;
    });
    ASSERT_NE(it, got.end());
    EXPECT_EQ(it->overlap_len, 90u);
    EXPECT_EQ(it->offset_a, 5u);

    std::filesystem::remove(spill_paths[0]);
    std::filesystem::remove(spill_paths[1]);
    std::filesystem::remove(spill_paths[2]);
    std::filesystem::remove(out_path);
    std::filesystem::remove_all(tmp_dir);
}

TEST(OverlapStoreTest, merge_forces_external_sort_under_tiny_ram_budget) {
    const auto dir = std::filesystem::temp_directory_path();
    const auto spill_path = (dir / "overlap_store_merge_big_spill.brov1").string();
    const auto out_path = (dir / "overlap_store_merge_big_out.brov1").string();
    const auto tmp_dir = (dir / "overlap_store_merge_big_tmp").string();

    constexpr std::uint32_t kCount = 200;
    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_path));
        // The winning record and its duplicate are kCount apart in the
        // stream, so a small chunk_capacity guarantees they land in
        // different sorted runs -- this exercises dedup across chunk
        // boundaries during the k-way merge, not just within one chunk.
        for (std::uint32_t i = 0; i < kCount; ++i) w.append(make_overlap(i, i + 1, 0, 0, 200));
        for (std::uint32_t i = 0; i < kCount; ++i) w.append(make_overlap(i, i + 1, 1, 1, 50));
        ASSERT_TRUE(w.close());
    }

    // 8 records/chunk over 400 total records forces ~50 sorted runs.
    const std::uint64_t tiny_ram_budget = 8 * sizeof(OverlapPair);
    ASSERT_NO_THROW(merge_spills({spill_path}, out_path, tmp_dir, tiny_ram_budget));

    OverlapStoreReader r;
    ASSERT_TRUE(r.open(out_path)) << r.last_error();
    ASSERT_EQ(r.record_count(), kCount);

    std::uint32_t n = 0;
    OverlapPair prev{};
    bool have_prev = false;
    OverlapPair rec;
    while (r.next(rec)) {
        EXPECT_EQ(rec.overlap_len, 200u) << "dedup should have kept the larger overlap_len";
        if (have_prev) {
            EXPECT_TRUE(std::tie(prev.read_a, prev.read_b) < std::tie(rec.read_a, rec.read_b));
        }
        prev = rec;
        have_prev = true;
        ++n;
    }
    EXPECT_EQ(n, kCount);

    if (std::filesystem::exists(tmp_dir)) {
        EXPECT_TRUE(std::filesystem::is_empty(tmp_dir)) << "temp run files must not survive close";
    }

    std::filesystem::remove(spill_path);
    std::filesystem::remove(out_path);
    std::filesystem::remove_all(tmp_dir);
}

TEST(OverlapStoreTest, merge_twice_is_byte_identical_in_ram_path) {
    const auto dir = std::filesystem::temp_directory_path();
    const std::vector<std::string> spill_paths = {
        (dir / "overlap_store_det_spill_a.brov1").string(),
        (dir / "overlap_store_det_spill_b.brov1").string(),
    };
    const auto out1 = (dir / "overlap_store_det_out1.brov1").string();
    const auto out2 = (dir / "overlap_store_det_out2.brov1").string();
    const auto tmp_dir = (dir / "overlap_store_det_tmp").string();

    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[0]));
        for (std::uint32_t i = 0; i < 300; ++i) {
            w.append(make_overlap(i % 37, (i * 13) % 41, i, i * 2, 100 + (i % 20)));
        }
        ASSERT_TRUE(w.close());
    }
    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[1]));
        for (std::uint32_t i = 0; i < 300; ++i) {
            w.append(make_overlap((i * 3) % 37, (i * 17) % 41, i + 1, i * 2 + 1, 100 + (i % 25)));
        }
        ASSERT_TRUE(w.close());
    }

    ASSERT_NO_THROW(merge_spills(spill_paths, out1, tmp_dir, std::uint64_t{1} << 30));
    ASSERT_NO_THROW(merge_spills(spill_paths, out2, tmp_dir, std::uint64_t{1} << 30));

    const std::string bytes1 = read_whole_file(out1);
    const std::string bytes2 = read_whole_file(out2);
    EXPECT_FALSE(bytes1.empty());
    EXPECT_EQ(bytes1, bytes2);

    std::filesystem::remove(spill_paths[0]);
    std::filesystem::remove(spill_paths[1]);
    std::filesystem::remove(out1);
    std::filesystem::remove(out2);
    std::filesystem::remove_all(tmp_dir);
}

TEST(OverlapStoreTest, merge_twice_is_byte_identical_external_path) {
    const auto dir = std::filesystem::temp_directory_path();
    const std::vector<std::string> spill_paths = {
        (dir / "overlap_store_det_ext_spill_a.brov1").string(),
        (dir / "overlap_store_det_ext_spill_b.brov1").string(),
    };
    const auto out1 = (dir / "overlap_store_det_ext_out1.brov1").string();
    const auto out2 = (dir / "overlap_store_det_ext_out2.brov1").string();
    const auto tmp_dir = (dir / "overlap_store_det_ext_tmp").string();

    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[0]));
        for (std::uint32_t i = 0; i < 300; ++i) {
            w.append(make_overlap(i % 37, (i * 13) % 41, i, i * 2, 100 + (i % 20)));
        }
        ASSERT_TRUE(w.close());
    }
    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(spill_paths[1]));
        for (std::uint32_t i = 0; i < 300; ++i) {
            w.append(make_overlap((i * 3) % 37, (i * 17) % 41, i + 1, i * 2 + 1, 100 + (i % 25)));
        }
        ASSERT_TRUE(w.close());
    }

    const std::uint64_t tiny_budget = 16 * sizeof(OverlapPair);
    ASSERT_NO_THROW(merge_spills(spill_paths, out1, tmp_dir, tiny_budget));
    ASSERT_NO_THROW(merge_spills(spill_paths, out2, tmp_dir, tiny_budget));

    const std::string bytes1 = read_whole_file(out1);
    const std::string bytes2 = read_whole_file(out2);
    EXPECT_FALSE(bytes1.empty());
    EXPECT_EQ(bytes1, bytes2);

    std::filesystem::remove(spill_paths[0]);
    std::filesystem::remove(spill_paths[1]);
    std::filesystem::remove(out1);
    std::filesystem::remove(out2);
    std::filesystem::remove_all(tmp_dir);
}

// --- error handling: corrupt magic / truncated files ------------------

TEST(OverlapStoreTest, reader_open_rejects_missing_file) {
    OverlapStoreReader r;
    EXPECT_FALSE(r.open("/nonexistent/path/that/should/not/exist.brov1"));
    EXPECT_FALSE(r.last_error().empty());
}

TEST(OverlapStoreTest, reader_open_rejects_bad_magic) {
    auto path = std::filesystem::temp_directory_path() / "overlap_store_bad_magic.brov1";
    {
        std::ofstream f(path, std::ios::binary);
        f << "NOTAVALIDHEADERATALL";
    }

    OverlapStoreReader r;
    EXPECT_FALSE(r.open(path.string()));
    EXPECT_NE(r.last_error().find("magic"), std::string::npos) << r.last_error();

    std::filesystem::remove(path);
}

TEST(OverlapStoreTest, reader_open_rejects_truncated_file) {
    auto path = std::filesystem::temp_directory_path() / "overlap_store_truncated.brov1";

    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(path.string()));
        for (std::uint32_t i = 0; i < 10; ++i) w.append(make_overlap(i, i + 1, 0, 0, 10));
        ASSERT_TRUE(w.close());
    }
    // Header (13B) plus a fraction of one 24B record -- header claims
    // 10 records, file holds none of them.
    std::filesystem::resize_file(path, 20);

    OverlapStoreReader r;
    EXPECT_FALSE(r.open(path.string()));
    EXPECT_NE(r.last_error().find("truncated"), std::string::npos) << r.last_error();

    std::filesystem::remove(path);
}

TEST(OverlapStoreTest, merge_throws_on_corrupt_spill) {
    const auto dir = std::filesystem::temp_directory_path();
    const auto good_path = (dir / "overlap_store_merge_good.brov1").string();
    const auto bad_path = (dir / "overlap_store_merge_bad.brov1").string();
    const auto out_path = (dir / "overlap_store_merge_should_not_exist.brov1").string();
    const auto tmp_dir = (dir / "overlap_store_merge_corrupt_tmp").string();

    {
        OverlapSpillWriter w;
        ASSERT_TRUE(w.open(good_path));
        w.append(make_overlap(1, 2, 0, 0, 10));
        ASSERT_TRUE(w.close());
    }
    {
        std::ofstream f(bad_path, std::ios::binary);
        f << "GARBAGE";
    }

    EXPECT_THROW(merge_spills({good_path, bad_path}, out_path, tmp_dir), std::runtime_error);

    std::filesystem::remove(good_path);
    std::filesystem::remove(bad_path);
    std::filesystem::remove(out_path);
    std::filesystem::remove_all(tmp_dir);
}
