#include "vpcr.hpp"
#include <algorithm>
#include <cctype>

namespace branch::vpcr {

std::string reverse_complement(std::string_view seq) {
    std::string result;
    result.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
        char c = static_cast<char>(std::toupper(static_cast<unsigned char>(*it)));
        switch (c) {
            case 'A': result += 'T'; break;
            case 'T': result += 'A'; break;
            case 'C': result += 'G'; break;
            case 'G': result += 'C'; break;
            default:  result += 'N'; break;
        }
    }
    return result;
}

namespace {

// Simple mismatch count between pattern and text at position
int count_mismatches(std::string_view text, size_t pos, std::string_view pattern) {
    if (pos + pattern.size() > text.size()) return static_cast<int>(pattern.size()) + 1;
    int mm = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        if (std::toupper(static_cast<unsigned char>(text[pos + i])) != 
            std::toupper(static_cast<unsigned char>(pattern[i]))) {
            ++mm;
        }
    }
    return mm;
}

// Find all positions where pattern matches with <= max_mismatches
std::vector<size_t> find_matches(std::string_view seq, std::string_view pattern, int max_mm) {
    std::vector<size_t> hits;
    if (pattern.empty() || pattern.size() > seq.size()) return hits;
    for (size_t i = 0; i + pattern.size() <= seq.size(); ++i) {
        if (count_mismatches(seq, i, pattern) <= max_mm) {
            hits.push_back(i);
        }
    }
    return hits;
}

} // anonymous namespace

std::vector<AmpliconHit> scan_amplicons(
    const PrimerPair& primers,
    std::string_view read_id,
    std::string_view read_seq,
    int max_mismatches
) {
    std::vector<AmpliconHit> results;
    
    // Forward primer matches
    auto fwd_hits = find_matches(read_seq, primers.fwd_seq, max_mismatches);
    
    // Reverse primer: we look for its reverse complement in the read
    std::string rev_rc = reverse_complement(primers.rev_seq);
    auto rev_hits = find_matches(read_seq, rev_rc, max_mismatches);
    
    // For each fwd hit, find rev hits that come after it
    for (size_t fwd_pos : fwd_hits) {
        size_t amplicon_start = fwd_pos;
        for (size_t rev_pos : rev_hits) {
            // rev_pos is where rev_rc starts; amplicon ends at rev_pos + rev_rc.size()
            if (rev_pos > fwd_pos + primers.fwd_seq.size()) {
                size_t amplicon_end = rev_pos + rev_rc.size();
                AmpliconHit hit;
                hit.read_id = std::string(read_id);
                hit.start = amplicon_start;
                hit.end = amplicon_end;
                hit.amplicon_seq = std::string(read_seq.substr(amplicon_start, amplicon_end - amplicon_start));
                results.push_back(std::move(hit));
            }
        }
    }
    
    return results;
}

} // namespace branch::vpcr
