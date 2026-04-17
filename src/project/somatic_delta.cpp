// BRANCH v0.4.3 — Somatic delta computation implementation.
//
// Uses ksw2 (via minimap2 FetchContent) for edit distance and CIGAR generation.

#include "somatic_delta.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

// ksw2 from minimap2
extern "C" {
#include <ksw2.h>
}

// kseq for FASTA parsing
#include <zlib.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
extern "C" {
#include <kseq.h>
}
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

namespace branch::project {

namespace {

// Convert ACGT to 0-3 encoding for ksw2
uint8_t base_to_num(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;  // N or other
    }
}

// Convert sequence to numeric encoding
std::vector<uint8_t> encode_sequence(const std::string& seq) {
    std::vector<uint8_t> encoded(seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        encoded[i] = base_to_num(seq[i]);
    }
    return encoded;
}

// Build scoring matrix for ksw2
void build_scoring_matrix(int8_t* mat, int match, int mismatch) {
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            if (i == 4 || j == 4) {
                mat[i * 5 + j] = 0;  // N matches nothing
            } else if (i == j) {
                mat[i * 5 + j] = static_cast<int8_t>(match);
            } else {
                mat[i * 5 + j] = static_cast<int8_t>(mismatch);
            }
        }
    }
}

// Parse CIGAR from ksw2 result and compute edit distance
// query and target are needed to count mismatches within M operations
std::pair<std::string, int> parse_cigar_and_edits(
    const uint32_t* cigar, int n_cigar,
    const std::string& query, const std::string& target) {
    
    std::ostringstream cigar_ss;
    int edit_distance = 0;
    size_t q_pos = 0;  // position in query
    size_t t_pos = 0;  // position in target
    
    for (int i = 0; i < n_cigar; ++i) {
        uint32_t len = cigar[i] >> 4;
        uint32_t op = cigar[i] & 0xf;
        
        char op_char;
        switch (op) {
            case 0:  // M: match/mismatch - need to count actual mismatches
                op_char = 'M';
                for (uint32_t j = 0; j < len && q_pos < query.size() && t_pos < target.size(); ++j) {
                    if (query[q_pos] != target[t_pos]) {
                        ++edit_distance;
                    }
                    ++q_pos;
                    ++t_pos;
                }
                break;
            case 1:  // I: insertion in query
                op_char = 'I';
                edit_distance += static_cast<int>(len);
                q_pos += len;
                break;
            case 2:  // D: deletion from query
                op_char = 'D';
                edit_distance += static_cast<int>(len);
                t_pos += len;
                break;
            case 7:  // =: sequence match
                op_char = '=';
                q_pos += len;
                t_pos += len;
                break;
            case 8:  // X: mismatch
                op_char = 'X';
                edit_distance += static_cast<int>(len);
                q_pos += len;
                t_pos += len;
                break;
            default:
                op_char = '?';
                break;
        }
        cigar_ss << len << op_char;
    }
    
    return {cigar_ss.str(), edit_distance};
}

// Fallback: count differences directly when alignment fails
int count_differences_direct(const std::string& query, const std::string& target) {
    int diff = 0;
    size_t min_len = std::min(query.size(), target.size());
    
    for (size_t i = 0; i < min_len; ++i) {
        if (query[i] != target[i]) {
            ++diff;
        }
    }
    
    // Length difference counts as edits
    diff += static_cast<int>(std::abs(static_cast<long>(query.size()) - static_cast<long>(target.size())));
    
    return diff;
}

// Read sequences from FASTA file
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string& path) {
    std::vector<std::pair<std::string, std::string>> sequences;
    
    gzFile fp = gzopen(path.c_str(), "r");
    if (!fp) return sequences;
    
    kseq_t* seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        sequences.emplace_back(seq->name.s, seq->seq.s);
    }
    kseq_destroy(seq);
    gzclose(fp);
    
    return sequences;
}

// Check for shell-unsafe characters
bool is_safe_path(const std::string& path) {
    for (char c : path) {
        if (c == '\0' || c == '\'' || c == '"' || c == '$' || c == '`' || c == '\\') {
            return false;
        }
    }
    return true;
}

}  // namespace

int ksw2_edit_distance(
    const std::string& query,
    const std::string& target,
    const SomaticDeltaOptions& opts,
    std::string* cigar_out) {
    
    if (query.empty() || target.empty()) {
        if (cigar_out) *cigar_out = "";
        return static_cast<int>(std::max(query.size(), target.size()));
    }
    
    // Fast path: identical sequences
    if (query == target) {
        if (cigar_out) *cigar_out = std::to_string(query.size()) + "M";
        return 0;
    }
    
    // Encode sequences
    auto q_enc = encode_sequence(query);
    auto t_enc = encode_sequence(target);
    
    // Build scoring matrix
    int8_t mat[25];
    build_scoring_matrix(mat, opts.match, opts.mismatch);
    
    // Allocate result structure
    ksw_extz_t ez;
    std::memset(&ez, 0, sizeof(ez));
    
    // Run ksw2 extension alignment with CIGAR
    // Use flag 0 for standard alignment (not just extension)
    ksw_extz2_sse(
        nullptr,  // km allocator (nullptr = malloc)
        static_cast<int>(q_enc.size()),
        q_enc.data(),
        static_cast<int>(t_enc.size()),
        t_enc.data(),
        5,  // alphabet size (A,C,G,T,N)
        mat,
        static_cast<int8_t>(opts.gap_open),
        static_cast<int8_t>(opts.gap_extend),
        -1,  // band width (-1 = full)
        -1,  // zdrop (-1 = disabled)
        0,   // end_bonus
        0,   // flag: 0 for standard alignment with CIGAR
        &ez
    );
    
    int edits = 0;
    std::string cigar;
    
    // Check if we got a valid CIGAR
    if (ez.n_cigar > 0 && ez.cigar != nullptr) {
        auto [parsed_cigar, parsed_edits] = parse_cigar_and_edits(ez.cigar, ez.n_cigar, query, target);
        cigar = parsed_cigar;
        edits = parsed_edits;
    } else {
        // Fallback: direct comparison
        edits = count_differences_direct(query, target);
        cigar = std::to_string(std::min(query.size(), target.size())) + "M";
        if (query.size() > target.size()) {
            cigar += std::to_string(query.size() - target.size()) + "I";
        } else if (target.size() > query.size()) {
            cigar += std::to_string(target.size() - query.size()) + "D";
        }
    }
    
    if (cigar_out) {
        *cigar_out = cigar;
    }
    
    // Free CIGAR memory
    if (ez.cigar) free(ez.cigar);
    
    return edits;
}

std::vector<SomaticDelta> compute_somatic_deltas(
    const std::string& branches_fasta,
    const std::vector<std::pair<std::string, std::string>>& ref_paths,
    const SomaticDeltaOptions& opts,
    std::string* err_out) {
    
    std::vector<SomaticDelta> results;
    
    // Early return for empty refs
    if (ref_paths.empty()) {
        return results;
    }
    
    // Validate paths
    if (!is_safe_path(branches_fasta)) {
        if (err_out) *err_out = "branches_fasta contains unsafe characters";
        return results;
    }
    
    for (const auto& [ref_name, ref_path] : ref_paths) {
        if (!is_safe_path(ref_path)) {
            if (err_out) *err_out = "ref_path contains unsafe characters: " + ref_name;
            return results;
        }
        if (!is_safe_path(ref_name)) {
            if (err_out) *err_out = "ref_name contains unsafe characters";
            return results;
        }
    }
    
    // Read branch sequences
    auto branches = read_fasta(branches_fasta);
    if (branches.empty()) {
        if (err_out) *err_out = "Failed to read branches FASTA or empty: " + branches_fasta;
        return results;
    }
    
    // Process each reference
    for (const auto& [ref_name, ref_path] : ref_paths) {
        auto refs = read_fasta(ref_path);
        if (refs.empty()) {
            continue;  // Skip missing references
        }
        
        // For simplicity, use first sequence in ref FASTA
        const std::string& ref_seq = refs[0].second;
        
        // Align each branch against this reference
        for (const auto& [branch_id, branch_seq] : branches) {
            SomaticDelta delta;
            delta.branch_id = branch_id;
            delta.ref_name = ref_name;
            
            std::string cigar;
            delta.edit_distance = ksw2_edit_distance(branch_seq, ref_seq, opts, &cigar);
            delta.cigar = cigar;
            
            delta.aligned_query_len = static_cast<std::int64_t>(branch_seq.size());
            delta.aligned_ref_len = static_cast<std::int64_t>(ref_seq.size());
            
            // Compute identity
            std::int64_t max_len = std::max(delta.aligned_query_len, delta.aligned_ref_len);
            if (max_len > 0) {
                delta.identity = 1.0 - static_cast<double>(delta.edit_distance) / static_cast<double>(max_len);
            }
            
            results.push_back(std::move(delta));
        }
    }
    
    // Sort by (ref_name, branch_id)
    std::sort(results.begin(), results.end(),
              [](const SomaticDelta& a, const SomaticDelta& b) {
                  if (a.ref_name != b.ref_name) return a.ref_name < b.ref_name;
                  return a.branch_id < b.branch_id;
              });
    
    return results;
}

}  // namespace branch::project
