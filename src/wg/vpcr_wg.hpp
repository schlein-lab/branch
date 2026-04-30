// BRANCH v0.5 — Phase 3: two-stage vPCR (cheap whole-genome screen +
// precise per-region amplicon counter).
//
// Stage 1 (in stage1_screen.{hpp,cpp}) flags candidate regions over the
// whole haplotype using the Phase 0 k-mer counter; only ~10⁴-10⁵
// windows survive instead of the ~10⁹ amplicons that naive design over
// 3 Gbp would emit.
//
// Stage 2 here:
//   - design 1-N primer pairs PER candidate region (low-entropy windows
//     paired at PCR-able distance), tagged with the region's flag class
//   - count via a k-mer-seed inverted index across all primers,
//     OpenMP-parallel over reads on CPU and a CUDA kernel on GPU.
//
// Output: every amplicon's read_count + CN estimate + Wilson 95 % CI,
// labelled with `region_class` so downstream consumers can distinguish
// elevated / depleted / repetitive / paralog loci without re-running.

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"
#include "wg/stage1_screen.hpp"

namespace branch::wg {

struct VpcrAmplicon {
    std::string amplicon_name;
    std::string fwd_primer;
    std::string rev_primer;
    std::uint32_t hap_idx = 0;
    std::uint32_t start_bp = 0;       // amplicon start in haplotype coords
    std::uint32_t expected_size_bp = 0;
    std::uint8_t  region_flags = 0;   // copied from CandidateRegion.flags
    std::uint32_t read_count = 0;
    double cn_estimate = 0.0;
    double ci_low = 0.0;
    double ci_high = 0.0;
};

class VpcrAutoDesigner {
public:
    explicit VpcrAutoDesigner(const ::branch::graph::TechProfile& profile);

    /// Stage 2 design: pick low-entropy primer windows inside each
    /// candidate region and pair them at PCR-able distance.
    /// `n_per_region` caps how many amplicons to emit per region.
    void design(const std::vector<std::string>& hap_seqs,
                const std::vector<CandidateRegion>& regions,
                std::vector<VpcrAmplicon>& out,
                std::uint32_t n_per_region = 3) const;

    /// Indexed-scan counter (CPU OpenMP path; GPU dispatch in
    /// wg_stage2_kernel.cu when CUDA is available).
    void count(std::vector<VpcrAmplicon>& amplicons,
               const std::vector<std::pair<std::string, std::string>>& reads,
               int expected_single_copy_count) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
