// BRANCH v0.5 — Phase 3: auto-primer vPCR design + counting.

#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "graph/kmer_sketch.hpp"

namespace branch::wg {

struct VpcrAmplicon {
    std::string amplicon_name;
    std::string fwd_primer;
    std::string rev_primer;
    std::uint32_t expected_size_bp;
    std::uint32_t read_count = 0;
    double cn_estimate = 0.0;
    double ci_low = 0.0;
    double ci_high = 0.0;
};

class VpcrAutoDesigner {
public:
    explicit VpcrAutoDesigner(const ::branch::graph::TechProfile& profile);

    /// Slide a window over haplotype sequences, find low-entropy primer
    /// candidates, pair them at PCR-able distance.
    void design(const std::vector<std::string>& hap_seqs,
                std::vector<VpcrAmplicon>& out) const;

    /// For each amplicon: scan reads, count those carrying both primers
    /// at the expected distance + orientation. Updates read_count + cn.
    void count(std::vector<VpcrAmplicon>& amplicons,
               const std::vector<std::pair<std::string, std::string>>& reads,
               int expected_single_copy_count) const;

private:
    ::branch::graph::TechProfile profile_;
};

}  // namespace branch::wg
