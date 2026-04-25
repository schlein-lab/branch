// Tests for branch::phase::snv_phaser. Coarse-but-correct: a tiny
// hand-built ReadDB lets us exercise the informative-SNV detector
// and the greedy block-builder without depending on real assembly.

#include <gtest/gtest.h>

#include "graph/read_db.hpp"
#include "phase/snv_phaser.hpp"

using branch::graph::DeltaOp;
using branch::graph::PathStep;
using branch::graph::ReadBucket;
using branch::graph::ReadDB;
using branch::phase::build_phase_blocks;
using branch::phase::find_informative_snvs;
using branch::phase::SnvPhaserConfig;

namespace {

// Helper: add a mapped read with a single path step on `unitig`
// covering positions [0, length), and a list of SNV deltas at
// (read-local pos, alt_base).
void add_read_with_snvs(ReadDB& db,
                         const char* name,
                         branch::graph::NodeId unitig,
                         std::uint32_t length,
                         std::vector<std::pair<std::uint32_t, char>> snvs) {
    std::vector<DeltaOp> deltas;
    deltas.reserve(snvs.size());
    for (auto [pos, base] : snvs) {
        deltas.push_back(DeltaOp{pos, 'S', base});
    }
    db.add(name, length, ReadBucket::Mapped, "",
           {{unitig, 0, length, false}}, std::move(deltas));
}

}  // namespace

TEST(SnvPhaser, Singleton_HiFi_errors_are_NOT_informative) {
    // Two reads carry an alt at the same position — but only ONE has
    // it. With min_alt_reads=2 the singleton is dismissed as a HiFi
    // sequencing artefact (the user-validated assumption: artefacts
    // are typically n=1 and never form a real branch).
    ReadDB db;
    for (int i = 0; i < 5; ++i) add_read_with_snvs(db, "ref_only", 0, 1000, {});
    add_read_with_snvs(db, "singleton_error", 0, 1000, {{500, 'A'}});

    auto snvs = find_informative_snvs(db, SnvPhaserConfig{.min_alt_reads = 2});
    EXPECT_TRUE(snvs.empty())
        << "n=1 alt at a position must NOT be flagged as informative";
}

TEST(SnvPhaser, Real_SNV_with_n_geq_2_alt_reads_is_informative) {
    ReadDB db;
    for (int i = 0; i < 5; ++i) add_read_with_snvs(db, "ref", 0, 1000, {});
    for (int i = 0; i < 3; ++i) add_read_with_snvs(db, "alt", 0, 1000, {{500, 'A'}});

    auto snvs = find_informative_snvs(db, SnvPhaserConfig{.min_alt_reads = 2});
    ASSERT_EQ(snvs.size(), 1u);
    EXPECT_EQ(snvs[0].unitig, 0u);
    EXPECT_EQ(snvs[0].pos, 500u);
    EXPECT_EQ(snvs[0].alt_base, 'A');
    EXPECT_EQ(snvs[0].alt_reads.size(), 3u);
}

TEST(SnvPhaser, Two_linked_SNVs_form_a_single_phase_block) {
    // 5 ref-ref reads, 3 alt-alt reads. Both SNVs at positions 200
    // and 800 of unitig 0. min_link_reads=2, the 3 alt-alt reads
    // satisfy it.
    ReadDB db;
    for (int i = 0; i < 5; ++i) {
        add_read_with_snvs(db, "haploA", 0, 1000, {});
    }
    for (int i = 0; i < 3; ++i) {
        add_read_with_snvs(db, "haploB", 0, 1000, {{200, 'A'}, {800, 'C'}});
    }

    auto snvs = find_informative_snvs(db);
    ASSERT_EQ(snvs.size(), 2u);

    auto blocks = build_phase_blocks(snvs);
    ASSERT_EQ(blocks.size(), 1u);
    EXPECT_EQ(blocks[0].snvs.size(), 2u);
    EXPECT_FALSE(blocks[0].read_haplotype.empty());
}

TEST(SnvPhaser, SNVs_on_different_unitigs_split_into_separate_blocks) {
    // Definitionally unlinkable: two SNVs sit on different unitigs, so
    // no path step can span both. The phaser must split them.
    ReadDB db;
    for (int i = 0; i < 5; ++i) add_read_with_snvs(db, "u0_ref", 0, 1000, {});
    for (int i = 0; i < 3; ++i) add_read_with_snvs(db, "u0_alt", 0, 1000, {{200, 'A'}});
    for (int i = 0; i < 5; ++i) add_read_with_snvs(db, "u1_ref", 1, 1000, {});
    for (int i = 0; i < 3; ++i) add_read_with_snvs(db, "u1_alt", 1, 1000, {{300, 'C'}});

    auto snvs = find_informative_snvs(db);
    ASSERT_EQ(snvs.size(), 2u);
    auto blocks = build_phase_blocks(snvs, SnvPhaserConfig{.min_link_reads = 2});
    EXPECT_EQ(blocks.size(), 2u);
}

// Documented limitation: when alt-carriers at SNV1 and alt-carriers at
// SNV2 are disjoint sets but the bulk of ref-only reads spans both
// positions, the greedy phaser links the SNVs into the same block
// through ref-ref evidence. Recovering trans phasing (alt at SNV1 on
// haplotype 1, alt at SNV2 on haplotype 2 — never co-occurring on one
// read) requires bipartite-min-cut-style logic. Tracked as a follow-up
// improvement; this test pins the current behaviour so any future
// fix shows up as an intentional change.
TEST(SnvPhaser, KnownLimit_disjoint_alt_sets_still_link_via_ref_reads) {
    ReadDB db;
    for (int i = 0; i < 4; ++i) add_read_with_snvs(db, "ref", 0, 1000, {});
    for (int i = 0; i < 3; ++i) {
        add_read_with_snvs(db, std::string("snv1_carrier_" + std::to_string(i)).c_str(),
                            0, 1000, {{200, 'A'}});
    }
    for (int i = 0; i < 3; ++i) {
        add_read_with_snvs(db, std::string("snv2_carrier_" + std::to_string(i)).c_str(),
                            0, 1000, {{800, 'C'}});
    }
    auto snvs = find_informative_snvs(db);
    ASSERT_EQ(snvs.size(), 2u);
    auto blocks = build_phase_blocks(snvs);
    EXPECT_EQ(blocks.size(), 1u);
}
