#pragma once
#include <string>
#include <string_view>
#include <vector>

namespace branch::vpcr {

struct PrimerPair {
    std::string name;
    std::string fwd_seq;
    std::string rev_seq;
};

struct AmpliconHit {
    std::string read_id;
    size_t start;
    size_t end;
    std::string amplicon_seq;
};

/// Returns the reverse complement of a DNA sequence
std::string reverse_complement(std::string_view seq);

/// Scans a read for amplicons bracketed by the primer pair
/// Returns all hits where fwd primer is found followed by rev primer (or its revcomp)
std::vector<AmpliconHit> scan_amplicons(
    const PrimerPair& primers,
    std::string_view read_id,
    std::string_view read_seq,
    int max_mismatches = 2
);

} // namespace branch::vpcr
