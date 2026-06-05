// BRANCH v0.4.3 — Somatic delta computation implementation.
//
// Edit distance + CIGAR via the LLmap classical alignment engine
// (branch::align::llmap_backend). No minimap2/ksw2 dependency.

#include "somatic_delta.hpp"
#include "../align/llmap_align_backend.hpp"

#include <algorithm>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>

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

int aligned_edit_distance(
    const std::string& query,
    const std::string& target,
    const SomaticDeltaOptions& opts,
    std::string* cigar_out) {

    // Map BRANCH's score-style options onto the LLmap gap-affine penalty
    // model. opts.mismatch is stored as a negative score; the backend wants
    // a positive penalty. opts.match (a positive score) has no penalty
    // analogue and is intentionally dropped — edit distance is derived from
    // alignment op-counts, not the scoring scheme.
    branch::align::llmap_backend::PairwiseOptions popts;
    popts.mismatch_penalty = std::abs(opts.mismatch);
    popts.gap_open         = opts.gap_open;
    popts.gap_extend       = opts.gap_extend;
    popts.global           = true;  // end-to-end ⇒ true edit distance

    return branch::align::llmap_backend::pairwise_edit_distance(
        query, target, popts, cigar_out);
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

            // Skip branches that are too long for the WFA2/Gotoh DP.
            if (opts.max_branch_len_bp > 0
                && static_cast<int>(branch_seq.size()) > opts.max_branch_len_bp) {
                delta.edit_distance = -1;  // sentinel: skipped
                delta.cigar = std::string("SKIP_LEN>") +
                              std::to_string(opts.max_branch_len_bp);
                delta.aligned_query_len = static_cast<std::int64_t>(branch_seq.size());
                delta.aligned_ref_len = static_cast<std::int64_t>(ref_seq.size());
                delta.identity = 0.0;
                results.push_back(std::move(delta));
                continue;
            }

            std::string cigar;
            delta.edit_distance = aligned_edit_distance(branch_seq, ref_seq, opts, &cigar);
            delta.cigar = cigar;
            delta.aligned_query_len = static_cast<std::int64_t>(branch_seq.size());
            delta.aligned_ref_len = static_cast<std::int64_t>(ref_seq.size());

            // Identity over the aligned region only.
            //
            // ksw2 in NW mode gives a global alignment whose CIGAR contains
            // leading + trailing D-runs when the query is shorter than the
            // target (i.e. branch fits inside a much bigger ref window).
            // Counting those Ds against identity collapses the metric to
            // ~0 % regardless of how good the actual hit was. So we strip
            // leading + trailing D from the CIGAR before computing identity,
            // matching the convention used by minimap2 (NM / aligned_len).
            // Edges-of-the-alignment behaviour (S/H clip ops) is also
            // skipped so the metric stays meaningful for partial hits.
            std::int64_t aln_q = 0, aln_t = 0;
            std::int64_t lead_d = 0, trail_d = 0;
            std::vector<std::pair<int,char>> ops_parsed;
            {
                std::size_t i = 0;
                while (i < cigar.size()) {
                    int n = 0;
                    while (i < cigar.size() && std::isdigit(static_cast<unsigned char>(cigar[i]))) {
                        n = n * 10 + (cigar[i] - '0'); ++i;
                    }
                    if (i >= cigar.size()) break;
                    char op = cigar[i++];
                    ops_parsed.emplace_back(n, op);
                }
            }
            std::size_t lo = 0, hi = ops_parsed.size();
            while (lo < hi && ops_parsed[lo].second == 'D') {
                lead_d += ops_parsed[lo].first; ++lo;
            }
            while (hi > lo && ops_parsed[hi-1].second == 'D') {
                trail_d += ops_parsed[hi-1].first; --hi;
            }
            for (std::size_t i = lo; i < hi; ++i) {
                int n = ops_parsed[i].first; char op = ops_parsed[i].second;
                if (op == 'M' || op == '=' || op == 'X') {
                    aln_q += n; aln_t += n;
                } else if (op == 'I') {
                    aln_q += n;
                } else if (op == 'D') {
                    aln_t += n;
                }
            }
            std::int64_t edits_in_region =
                static_cast<std::int64_t>(delta.edit_distance) - lead_d - trail_d;
            if (edits_in_region < 0) edits_in_region = 0;
            std::int64_t denom = std::max(aln_q, aln_t);
            if (denom > 0) {
                delta.identity = 1.0 - static_cast<double>(edits_in_region)
                                       / static_cast<double>(denom);
                if (delta.identity < 0.0) delta.identity = 0.0;
            } else {
                delta.identity = 0.0;
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
