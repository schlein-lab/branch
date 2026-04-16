// BRANCH v0.2 — Strand-handling regression tests.
//
// Guards the v0.2 fix for the canonical-hash strand-collapse bug: every
// overlap MUST carry a strand bit so downstream graph building can emit
// correct "+/+" vs "+/-" GFA L-lines and classifier logic can reason
// about dovetail geometry.
//
// The tests also pin the struct sizes (OverlapPair=24, MinimizerHit=16)
// because the GPU minimizer-entry layout depends on those values and
// any regression silently breaks device-side uploads.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "graph/kmer_sketch.hpp"

using branch::backend::make_cpu_backend;
using branch::backend::OverlapPair;
using branch::backend::ReadBatch;
using branch::graph::MinimizerHit;

namespace {

// Produce a deterministic DNA string of the given length.
std::string make_seq(std::size_t len, unsigned seed) {
    const char bases[] = "ACGT";
    std::string s;
    s.reserve(len);
    for (std::size_t i = 0; i < len; ++i) {
        seed = seed * 1103515245u + 12345u;
        s.push_back(bases[(seed >> 16) & 3u]);
    }
    return s;
}

// Reverse-complement a DNA string.
std::string revcomp(const std::string& s) {
    std::string r;
    r.resize(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        char c = s[s.size() - 1 - i];
        switch (c) {
            case 'A': r[i] = 'T'; break;
            case 'C': r[i] = 'G'; break;
            case 'G': r[i] = 'C'; break;
            case 'T': r[i] = 'A'; break;
            default:  r[i] = c;   break;
        }
    }
    return r;
}

}  // namespace

// ---------------------------------------------------------------------
// Layout invariants
// ---------------------------------------------------------------------

TEST(StrandHandling, OverlapPair_packs_to_24_bytes) {
    // The struct picked up a 1-byte strand field in v0.2. Alignment must
    // still keep the whole thing at 24 bytes so SIMD batches and GPU
    // zero-copy transfers don't regress.
    EXPECT_EQ(sizeof(OverlapPair), 24u);
}

TEST(StrandHandling, MinimizerHit_packs_to_16_bytes) {
    // MinimizerHit layout is shared with branch::gpu::MinimizerEntry for
    // host→device upload. A 1-byte strand field was added by shrinking
    // `pos` from 32 to 24 bits; total size must stay 16.
    EXPECT_EQ(sizeof(MinimizerHit), 16u);
}

// ---------------------------------------------------------------------
// End-to-end overlap strand detection
// ---------------------------------------------------------------------

TEST(StrandHandling, Same_strand_overlap_emits_strand_zero) {
    // Two reads sharing a 400 bp identical substring on the same strand
    // should produce >= 1 overlap with strand == 0.
    std::string shared = make_seq(400, 42);

    std::string read1 = make_seq(50, 100) + shared + make_seq(50, 200);
    std::string read2 = make_seq(60, 300) + shared + make_seq(40, 400);

    ReadBatch batch;
    batch.reads.push_back({std::string_view(read1), 0});
    batch.reads.push_back({std::string_view(read2), 1});

    auto b = make_cpu_backend();
    std::array<OverlapPair, 16> out{};
    std::size_t count = 0;
    b.compute_overlaps(&batch, out, &count);

    ASSERT_GE(count, 1u)
        << "Expected at least one overlap for 400bp same-strand region";

    bool found_same = false;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& p = out[i];
        const bool pair_is_01 =
            (p.read_a == 0u && p.read_b == 1u) ||
            (p.read_a == 1u && p.read_b == 0u);
        if (pair_is_01 && p.strand == 0u) {
            found_same = true;
            break;
        }
    }
    EXPECT_TRUE(found_same)
        << "Expected an overlap between reads 0 and 1 with strand == 0";
}

TEST(StrandHandling, Opposite_strand_overlap_emits_strand_one) {
    // Read 2 embeds the REVERSE COMPLEMENT of read 1's shared region.
    // Same minimizer k-mers (canonical hash collapses forward/rc), but
    // the sketcher's strand-bit disagrees between the two reads → the
    // backend must mark the OverlapPair with strand == 1.
    std::string shared = make_seq(400, 7);
    std::string shared_rc = revcomp(shared);

    std::string read1 = make_seq(50, 111) + shared     + make_seq(50, 222);
    std::string read2 = make_seq(60, 333) + shared_rc  + make_seq(40, 444);

    ReadBatch batch;
    batch.reads.push_back({std::string_view(read1), 0});
    batch.reads.push_back({std::string_view(read2), 1});

    auto b = make_cpu_backend();
    std::array<OverlapPair, 16> out{};
    std::size_t count = 0;
    b.compute_overlaps(&batch, out, &count);

    ASSERT_GE(count, 1u)
        << "Expected at least one overlap for 400bp opposite-strand region";

    bool found_opp = false;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& p = out[i];
        const bool pair_is_01 =
            (p.read_a == 0u && p.read_b == 1u) ||
            (p.read_a == 1u && p.read_b == 0u);
        if (pair_is_01 && p.strand == 1u) {
            found_opp = true;
            break;
        }
    }
    EXPECT_TRUE(found_opp)
        << "Expected an overlap between reads 0 and 1 with strand == 1";
}
