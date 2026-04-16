#include "consensus/majority_voter.hpp"

#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using branch::consensus::build_consensus;
using branch::consensus::ConsensusConfig;
using branch::consensus::ConsensusResult;

namespace {

std::vector<std::string_view>
views(const std::vector<std::string>& owners) {
    std::vector<std::string_view> out;
    out.reserve(owners.size());
    for (const auto& s : owners) {
        out.emplace_back(s);
    }
    return out;
}

}  // namespace

TEST(MajorityVoter, EmptyInputReturnsEmptyResult) {
    const std::vector<std::string_view> empty;
    const ConsensusResult r = build_consensus(empty);
    EXPECT_TRUE(r.sequence.empty());
    EXPECT_TRUE(r.depths.empty());
    EXPECT_TRUE(r.majority_fraction.empty());
}

TEST(MajorityVoter, ThreeIdenticalSequencesProducesInput) {
    const std::vector<std::string> owners{
        "ACGTACGT",
        "ACGTACGT",
        "ACGTACGT",
    };
    const auto ins = views(owners);
    const ConsensusResult r = build_consensus(ins);
    EXPECT_EQ(r.sequence, "ACGTACGT");
    ASSERT_EQ(r.depths.size(), 8u);
    for (auto d : r.depths) {
        EXPECT_EQ(d, 3u);
    }
    for (auto f : r.majority_fraction) {
        EXPECT_FLOAT_EQ(f, 1.0f);
    }
}

TEST(MajorityVoter, MajorityOverridesOneDifferingPosition) {
    // Column 3 has two 'T' and one 'A' → majority 'T'.
    const std::vector<std::string> owners{
        "ACGTACGT",
        "ACGAACGT",  // differs at col 3
        "ACGTACGT",
    };
    const auto ins = views(owners);
    const ConsensusResult r = build_consensus(ins);
    EXPECT_EQ(r.sequence, "ACGTACGT");
    ASSERT_EQ(r.majority_fraction.size(), 8u);
    // Col 3: 2/3 T.
    EXPECT_NEAR(r.majority_fraction[3], 2.0f / 3.0f, 1e-6f);
    EXPECT_EQ(r.depths[3], 3u);
}

TEST(MajorityVoter, TieAtColumnEmitsAmbiguousChar) {
    // Two sequences, one 'A' one 'C' → tie at col 0.
    // min_coverage=2 satisfied, min_majority_fraction=0.5 satisfied
    // (0.5 >= 0.5), but tie must still resolve to 'N'.
    const std::vector<std::string> owners{
        "A",
        "C",
    };
    const auto ins = views(owners);
    ConsensusConfig cfg;
    cfg.min_coverage = 2;
    cfg.min_majority_fraction = 0.5f;
    const ConsensusResult r = build_consensus(ins, cfg);
    ASSERT_EQ(r.sequence.size(), 1u);
    EXPECT_EQ(r.sequence[0], cfg.ambiguous_char);
    EXPECT_EQ(r.depths[0], 2u);
    EXPECT_FLOAT_EQ(r.majority_fraction[0], 0.5f);
}

TEST(MajorityVoter, MismatchedLengthsThrow) {
    const std::vector<std::string> owners{
        "ACGT",
        "ACG",  // shorter
    };
    const auto ins = views(owners);
    EXPECT_THROW(build_consensus(ins), std::invalid_argument);
}

TEST(MajorityVoter, GapsAreIgnoredForCoverage) {
    // Col 0: A, A, '-' → coverage 2, majority 'A' at 100%.
    // Col 1: C, '-', '-' → coverage 1, below min_coverage=2 → 'N'.
    const std::vector<std::string> owners{
        "AC",
        "A-",
        "--",
    };
    const auto ins = views(owners);
    ConsensusConfig cfg;
    cfg.min_coverage = 2;
    const ConsensusResult r = build_consensus(ins, cfg);
    ASSERT_EQ(r.sequence.size(), 2u);
    EXPECT_EQ(r.sequence[0], 'A');
    EXPECT_EQ(r.sequence[1], cfg.ambiguous_char);
    EXPECT_EQ(r.depths[0], 2u);
    EXPECT_EQ(r.depths[1], 1u);
}
