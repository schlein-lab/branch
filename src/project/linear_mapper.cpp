// BRANCH v0.4 — Linear reference mapper implementation.
//
// Maps every branch consensus contig against each linear reference
// (CHM13, GRCh38, ...) via the in-process LLmap classical engine
// (branch::align::llmap_backend::ReferenceMapper). No minimap2 subprocess,
// no PAF round-trip.

#include "linear_mapper.hpp"
#include "../align/llmap_align_backend.hpp"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <zlib.h>

namespace branch::project {

namespace {

// Minimal gz-aware FASTA reader → (name, sequence) records. Name is the
// header up to the first whitespace.
std::vector<std::pair<std::string, std::string>>
read_fasta_named(const std::string& path) {
    std::vector<std::pair<std::string, std::string>> out;
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) return out;

    constexpr int kBuf = 1 << 16;
    std::vector<char> buf(kBuf);
    std::string line, name, seq;
    auto flush = [&]() {
        if (!name.empty() || !seq.empty()) out.emplace_back(name, seq);
        name.clear();
        seq.clear();
    };
    auto take_header = [&](const std::string& l) {
        flush();
        std::size_t sp = l.find_first_of(" \t");
        name = l.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
    };

    int n;
    while ((n = gzread(fp, buf.data(), kBuf)) > 0) {
        for (int i = 0; i < n; ++i) {
            char c = buf[i];
            if (c == '\n') {
                if (!line.empty()) {
                    if (line[0] == '>') take_header(line);
                    else seq += line;
                }
                line.clear();
            } else if (c != '\r') {
                line += c;
            }
        }
    }
    if (!line.empty()) {
        if (line[0] == '>') take_header(line);
        else seq += line;
    }
    flush();
    gzclose(fp);
    return out;
}

}  // namespace

std::vector<LinearMapping> map_branches_linear(
    const std::string& fasta_path,
    const std::vector<LinearRef>& refs,
    const LinearMapOptions& opts,
    std::string* err_out) {

    std::vector<LinearMapping> result;

    // Load the branch contigs once; they are queried against every reference.
    auto branches = read_fasta_named(fasta_path);
    if (branches.empty()) {
        if (err_out) *err_out = "no branch sequences read from: " + fasta_path;
        return result;
    }

    for (const auto& ref : refs) {
        // One minimizer index per reference. The preset name is passed
        // through; ReferenceMapper maps non-ONT presets (asm5/asm20/map-hifi)
        // onto its high-identity defaults, appropriate for same-species
        // contig-vs-reference projection.
        branch::align::llmap_backend::ReferenceMapper mapper(ref.path, opts.preset);
        if (!mapper.valid()) {
            if (err_out) *err_out = "failed to index reference: " + ref.name;
            continue;  // try remaining refs
        }

        for (const auto& [branch_id, seq] : branches) {
            for (const auto& hit : mapper.all(seq)) {
                if (static_cast<int>(hit.mapq) < opts.min_mapq) continue;
                LinearMapping m;
                m.branch_id    = branch_id;
                m.ref_name     = ref.name;
                m.target       = hit.ref_name;
                m.target_start = hit.ref_start;
                m.target_end   = hit.ref_end;
                m.mapq         = static_cast<int>(hit.mapq);
                m.strand       = hit.is_forward ? '+' : '-';
                m.query_len    = static_cast<std::int64_t>(seq.size());
                m.query_start  = hit.query_start;
                m.query_end    = hit.query_end;
                result.push_back(std::move(m));
            }
        }
    }

    // Sort by (ref_name, branch_id, query_start).
    std::sort(result.begin(), result.end(),
              [](const LinearMapping& a, const LinearMapping& b) {
                  if (a.ref_name != b.ref_name) return a.ref_name < b.ref_name;
                  if (a.branch_id != b.branch_id) return a.branch_id < b.branch_id;
                  return a.query_start < b.query_start;
              });

    return result;
}

}  // namespace branch::project
