// Tests for branch::common::memory helpers.
//
// parse_memory_size is the critical thing — it feeds setrlimit and a
// wrong multiplier would either refuse every valid budget or cap the
// process at a byte, so we cover every suffix form we accept.

#include <gtest/gtest.h>

#include <string>

#include "common/memory.hpp"

using branch::common::format_bytes;
using branch::common::parse_memory_size;
using branch::common::set_memory_budget;

TEST(Memory, parse_plain_bytes) {
    EXPECT_EQ(parse_memory_size("0").value(), 0u);
    EXPECT_EQ(parse_memory_size("1").value(), 1u);
    EXPECT_EQ(parse_memory_size("1024").value(), 1024u);
    EXPECT_EQ(parse_memory_size("1024B").value(), 1024u);
}

TEST(Memory, parse_kmgt_suffixes_binary) {
    EXPECT_EQ(parse_memory_size("1K").value(),   1024u);
    EXPECT_EQ(parse_memory_size("1KB").value(),  1024u);
    EXPECT_EQ(parse_memory_size("1KiB").value(), 1024u);
    EXPECT_EQ(parse_memory_size("1M").value(),   1024u * 1024u);
    EXPECT_EQ(parse_memory_size("1MB").value(),  1024u * 1024u);
    EXPECT_EQ(parse_memory_size("1MiB").value(), 1024u * 1024u);
    EXPECT_EQ(parse_memory_size("1G").value(),   1024u * 1024u * 1024u);
    EXPECT_EQ(parse_memory_size("1T").value(),   1024ull * 1024 * 1024 * 1024);
}

TEST(Memory, parse_case_insensitive_suffix) {
    EXPECT_EQ(parse_memory_size("8g").value(),   8ull * 1024 * 1024 * 1024);
    EXPECT_EQ(parse_memory_size("8G").value(),   8ull * 1024 * 1024 * 1024);
    EXPECT_EQ(parse_memory_size("8gb").value(),  8ull * 1024 * 1024 * 1024);
    EXPECT_EQ(parse_memory_size("8GIB").value(), 8ull * 1024 * 1024 * 1024);
}

TEST(Memory, parse_decimal_multiplier) {
    // 1.5 GiB = 1.5 * 2^30 = 1610612736 bytes (truncated toward zero).
    EXPECT_EQ(parse_memory_size("1.5G").value(), 1610612736u);
    // Locale-robustness: comma as decimal point.
    EXPECT_EQ(parse_memory_size("1,5G").value(), 1610612736u);
}

TEST(Memory, parse_strips_whitespace) {
    EXPECT_EQ(parse_memory_size("  8G  ").value(),
              8ull * 1024 * 1024 * 1024);
}

TEST(Memory, parse_rejects_malformed) {
    EXPECT_FALSE(parse_memory_size("").has_value());
    EXPECT_FALSE(parse_memory_size("G").has_value());
    EXPECT_FALSE(parse_memory_size("abc").has_value());
    EXPECT_FALSE(parse_memory_size("8X").has_value());
    EXPECT_FALSE(parse_memory_size("-1G").has_value());
}

TEST(Memory, format_bytes_picks_best_unit) {
    EXPECT_EQ(format_bytes(0), "0.0 B");
    EXPECT_EQ(format_bytes(1), "1.0 B");
    EXPECT_EQ(format_bytes(1024), "1.0 KiB");
    EXPECT_EQ(format_bytes(1024 * 1024), "1.0 MiB");
    EXPECT_EQ(format_bytes(1ull * 1024 * 1024 * 1024), "1.0 GiB");
    EXPECT_EQ(format_bytes(1536 * 1024 * 1024ull), "1.5 GiB");
}

TEST(Memory, set_budget_is_noop_on_zero) {
    // Zero bytes = no cap. Must succeed without touching limits.
    EXPECT_TRUE(set_memory_budget(0));
}

TEST(Memory, set_budget_clamps_to_hard_limit_safely) {
    // Setting a very large soft limit (e.g. 2 EiB) gets clamped to the
    // inherited hard limit rather than failing. The call returns true
    // and we don't assert any specific limit change — this test is
    // chiefly a crash-safety check.
    EXPECT_TRUE(set_memory_budget(1ull << 62));
}
