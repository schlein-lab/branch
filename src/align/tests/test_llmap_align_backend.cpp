// Unit tests for the LLmap alignment backend (the minimap2/ksw2 replacement).
//
// Standalone (no GoogleTest) so it runs under ctest with only BUILD_TESTING,
// no network FetchContent. Exercises both seams:
//   1. align_pair / pairwise_edit_distance  (former ksw2)
//   2. ReferenceMapper                       (former minimap2 mm_idx + mm_map)

#include "align/llmap_align_backend.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace bk = branch::align::llmap_backend;

static int g_failures = 0;

#define CHECK(cond, msg)                                                      \
    do {                                                                      \
        if (!(cond)) {                                                        \
            std::fprintf(stderr, "FAIL: %s (%s:%d)\n", msg, __FILE__,         \
                         __LINE__);                                           \
            ++g_failures;                                                     \
        }                                                                     \
    } while (0)

// Deterministic, k-mer-diverse sequence (LCG over the 4 bases).
static std::string make_seq(std::size_t n, std::uint64_t seed) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    std::string s;
    s.reserve(n);
    std::uint64_t x = seed | 1ULL;
    for (std::size_t i = 0; i < n; ++i) {
        x = x * 6364136223846793005ULL + 1442695040888963407ULL;
        s.push_back(B[(x >> 33) & 3]);
    }
    return s;
}

static void test_pairwise() {
    const std::string a = make_seq(300, 42);

    // Identical → zero edits, full CIGAR M.
    std::string cig;
    CHECK(bk::pairwise_edit_distance(a, a, {}, &cig) == 0, "identical edits==0");

    // One substitution → exactly one edit.
    std::string b = a;
    b[150] = (b[150] == 'A') ? 'C' : 'A';
    CHECK(bk::pairwise_edit_distance(a, b) == 1, "1 substitution edits==1");

    // A two-base deletion in the query → two indel edits.
    std::string c = a.substr(0, 100) + a.substr(102);
    CHECK(bk::pairwise_edit_distance(a, c) == 2, "2bp deletion edits==2");

    // Full alignment record: identity high, cigar non-empty.
    auto aln = bk::align_pair(a, b);
    CHECK(aln.has_value(), "align_pair returns a result");
    if (aln) {
        CHECK(!aln->cigar.empty(), "cigar non-empty");
        CHECK(aln->identity > 0.99f, "identity > 0.99 for 1-mismatch/300");
        CHECK(aln->num_mismatches == 1, "exactly 1 mismatch counted");
    }
}

static void test_reference_mapper() {
    // Two reference contigs; map a near-identical substring of contig 1.
    const std::string ref0 = make_seq(800, 7);
    const std::string ref1 = make_seq(800, 99);
    std::vector<std::string> names = {"contig0", "contig1"};
    std::vector<std::string> seqs  = {ref0, ref1};

    bk::ReferenceMapper mapper(names, seqs, "map-hifi");
    CHECK(mapper.valid(), "mapper built from sequences");

    // Query: a 300 bp window of contig1 with two substitutions.
    std::string q = ref1.substr(250, 300);
    q[40]  = (q[40]  == 'A') ? 'G' : 'A';
    q[200] = (q[200] == 'T') ? 'C' : 'T';

    auto hit = mapper.best(q);
    CHECK(hit.has_value(), "query maps somewhere");
    if (hit) {
        CHECK(hit->ref_name == "contig1", "maps to the correct contig");
        CHECK(hit->identity > 0.95f, "high identity to source contig");
        // Mapped near the true origin (allow chain/extension slack).
        CHECK(hit->ref_start >= 200 && hit->ref_start <= 300,
              "ref_start near true window origin (~250)");
    }

    // A random sequence unrelated to either contig should not map (or map
    // far below the identity floor → no primary hit).
    std::string noise = make_seq(300, 123456);
    auto miss = mapper.best(noise);
    CHECK(!miss.has_value(), "unrelated sequence does not map");
}

int main() {
    test_pairwise();
    test_reference_mapper();
    if (g_failures == 0) {
        std::printf("llmap_align_backend: all checks passed\n");
        return 0;
    }
    std::fprintf(stderr, "llmap_align_backend: %d checks FAILED\n", g_failures);
    return 1;
}
