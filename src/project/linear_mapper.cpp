// BRANCH v0.4 — Linear reference mapper implementation.
//
// Shell-out to minimap2 per (branch FASTA x ref) pair, parse PAF output,
// return merged mappings sorted by (ref_name, query_start).

#include "linear_mapper.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>

namespace branch::project {

namespace {

// Check if a path contains unsafe shell characters.
// Returns true if safe, false if dangerous.
bool is_safe_shell_path(const std::string& path) {
    for (char c : path) {
        if (c == '\0' || c == '\'' || c == '"' || c == '$' || c == '`' || c == '\\') {
            return false;
        }
    }
    return true;
}

// Parse a single PAF line into a LinearMapping.
// PAF format (tab-separated):
//   0: query_name, 1: query_len, 2: qstart, 3: qend, 4: strand (+/-),
//   5: target_name, 6: target_len, 7: tstart, 8: tend, 9: n_matches,
//   10: alignment_block_len, 11: mapq, 12+: optional tags
// Returns true on success, false on parse error.
bool parse_paf_line(const std::string& line, const std::string& ref_name,
                    LinearMapping& out) {
    if (line.empty() || line[0] == '#') {
        return false;
    }

    std::istringstream iss(line);
    std::string query_name, strand_str, target_name;
    std::int64_t query_len, qstart, qend, target_len, tstart, tend;
    std::int64_t n_matches, block_len;
    int mapq;

    // Read all 12 mandatory PAF columns
    if (!(iss >> query_name >> query_len >> qstart >> qend >> strand_str
              >> target_name >> target_len >> tstart >> tend
              >> n_matches >> block_len >> mapq)) {
        return false;
    }

    // Validate strand
    if (strand_str.empty() || (strand_str[0] != '+' && strand_str[0] != '-')) {
        return false;
    }

    out.branch_id = query_name;
    out.ref_name = ref_name;
    out.target = target_name;
    out.target_start = tstart;
    out.target_end = tend;
    out.mapq = mapq;
    out.strand = strand_str[0];
    out.query_len = query_len;
    out.query_start = qstart;
    out.query_end = qend;

    return true;
}

// Build the minimap2 command line.
std::string build_minimap2_cmd(const std::string& minimap2_path,
                                const std::string& preset,
                                int threads,
                                const std::string& ref_path,
                                const std::string& fasta_path) {
    std::ostringstream cmd;
    cmd << minimap2_path
        << " -x " << preset
        << " -t " << threads
        << " '" << ref_path << "'"
        << " '" << fasta_path << "'"
        << " 2>/dev/null";
    return cmd.str();
}

}  // namespace

std::vector<LinearMapping> map_branches_linear(
    const std::string& fasta_path,
    const std::vector<LinearRef>& refs,
    const LinearMapOptions& opts,
    std::string* err_out) {

    std::vector<LinearMapping> result;

    // Validate paths for shell safety
    if (!is_safe_shell_path(fasta_path)) {
        if (err_out) {
            *err_out = "fasta_path contains unsafe shell characters";
        }
        return result;
    }

    if (!is_safe_shell_path(opts.minimap2_path)) {
        if (err_out) {
            *err_out = "minimap2_path contains unsafe shell characters";
        }
        return result;
    }

    for (const auto& ref : refs) {
        if (!is_safe_shell_path(ref.path)) {
            if (err_out) {
                *err_out = "ref path '" + ref.name + "' contains unsafe shell characters";
            }
            return result;
        }
        if (!is_safe_shell_path(ref.name)) {
            if (err_out) {
                *err_out = "ref name '" + ref.name + "' contains unsafe shell characters";
            }
            return result;
        }
    }

    // Map against each reference
    for (const auto& ref : refs) {
        std::string cmd = build_minimap2_cmd(
            opts.minimap2_path, opts.preset, opts.threads,
            ref.path, fasta_path);

        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            if (err_out) {
                *err_out = "minimap2 not found or failed to execute: " + opts.minimap2_path;
            }
            return {};  // Return empty on failure
        }

        // Read and parse PAF output line by line
        char buffer[8192];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            // Remove trailing newline
            std::size_t len = std::strlen(buffer);
            if (len > 0 && buffer[len - 1] == '\n') {
                buffer[len - 1] = '\0';
            }

            LinearMapping mapping;
            if (parse_paf_line(buffer, ref.name, mapping)) {
                if (mapping.mapq >= opts.min_mapq) {
                    result.push_back(std::move(mapping));
                }
            }
        }

        int status = pclose(pipe);
        if (status != 0) {
            // minimap2 exited non-zero — could be missing ref, but we continue
            // to try other refs. Only error out if ALL refs fail.
        }
    }

    // Sort by (ref_name, branch_id, query_start)
    std::sort(result.begin(), result.end(),
              [](const LinearMapping& a, const LinearMapping& b) {
                  if (a.ref_name != b.ref_name) {
                      return a.ref_name < b.ref_name;
                  }
                  if (a.branch_id != b.branch_id) {
                      return a.branch_id < b.branch_id;
                  }
                  return a.query_start < b.query_start;
              });

    return result;
}

}  // namespace branch::project
