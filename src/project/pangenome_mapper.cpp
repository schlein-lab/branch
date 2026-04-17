// BRANCH v0.4 — Pangenome mapper implementation.
//
// Shell-out to GraphAligner or minigraph per (branch FASTA x ref) pair,
// parse GAF output, return merged mappings sorted by (ref_name, branch_id, path_start).

#include "pangenome_mapper.hpp"

#include <algorithm>
#include <cstdint>
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

// Parse a single GAF line into a PangenomeMapping.
// GAF format (tab-separated):
//   0: query_name, 1: query_len, 2: qstart, 3: qend, 4: strand (+/-),
//   5: path_name (string starting with > or <), 6: path_len, 7: path_start,
//   8: path_end, 9: n_matches, 10: alignment_block_len, 11: mapq, 12+: optional tags
// Returns true on success, false on parse error.
bool parse_gaf_line(const std::string& line, const std::string& ref_name,
                    PangenomeMapping& out) {
    if (line.empty() || line[0] == '#') {
        return false;
    }

    std::istringstream iss(line);
    std::string query_name, strand_str, path_name;
    std::int64_t query_len, qstart, qend, path_len, path_start, path_end;
    std::int64_t n_matches, block_len;
    int mapq;

    // Read all 12 mandatory GAF columns
    if (!(iss >> query_name >> query_len >> qstart >> qend >> strand_str
              >> path_name >> path_len >> path_start >> path_end
              >> n_matches >> block_len >> mapq)) {
        return false;
    }

    // Validate strand
    if (strand_str.empty() || (strand_str[0] != '+' && strand_str[0] != '-')) {
        return false;
    }

    // Compute identity from matches / block_len
    double identity = 0.0;
    if (block_len > 0) {
        identity = static_cast<double>(n_matches) / static_cast<double>(block_len);
    }

    out.branch_id = query_name;
    out.ref_name = ref_name;
    out.path_names = path_name;
    out.path_start = path_start;
    out.path_end = path_end;
    out.mapq = mapq;
    out.identity = identity;
    out.query_len = query_len;

    return true;
}

// Build the GraphAligner command line.
std::string build_graphaligner_cmd(const std::string& tool_path,
                                    int threads,
                                    const std::string& graph_path,
                                    const std::string& fasta_path) {
    std::ostringstream cmd;
    cmd << tool_path
        << " -g '" << graph_path << "'"
        << " -f '" << fasta_path << "'"
        << " -t " << threads
        << " -a /dev/stdout"
        << " -x vg"
        << " 2>/dev/null";
    return cmd.str();
}

// Build the minigraph command line.
std::string build_minigraph_cmd(const std::string& tool_path,
                                 int threads,
                                 const std::string& graph_path,
                                 const std::string& fasta_path) {
    std::ostringstream cmd;
    cmd << tool_path
        << " -t " << threads
        << " -c"
        << " '" << graph_path << "'"
        << " '" << fasta_path << "'"
        << " 2>/dev/null";
    return cmd.str();
}

}  // namespace

std::vector<PangenomeMapping> map_branches_pangenome(
    const std::string& fasta_path,
    const std::vector<PangenomeRef>& refs,
    const PangenomeMapOptions& opts,
    std::string* err_out) {

    std::vector<PangenomeMapping> result;

    // Early return for empty refs
    if (refs.empty()) {
        return result;
    }

    // Validate paths for shell safety
    if (!is_safe_shell_path(fasta_path)) {
        if (err_out) {
            *err_out = "fasta_path contains unsafe shell characters";
        }
        return result;
    }

    if (!is_safe_shell_path(opts.tool_path)) {
        if (err_out) {
            *err_out = "tool_path contains unsafe shell characters";
        }
        return result;
    }

    for (const auto& ref : refs) {
        if (!is_safe_shell_path(ref.gbz_path)) {
            if (err_out) {
                *err_out = "ref gbz_path '" + ref.name + "' contains unsafe shell characters";
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
        std::string cmd;
        if (opts.tool == "minigraph") {
            cmd = build_minigraph_cmd(opts.tool_path, opts.threads,
                                       ref.gbz_path, fasta_path);
        } else {
            // Default: GraphAligner
            cmd = build_graphaligner_cmd(opts.tool_path, opts.threads,
                                          ref.gbz_path, fasta_path);
        }

        FILE* pipe = popen(cmd.c_str(), "r");
        if (!pipe) {
            if (err_out) {
                *err_out = opts.tool + " not found or failed to execute: " + opts.tool_path;
            }
            return {};  // Return empty on failure
        }

        // Read and parse GAF output line by line
        char buffer[16384];  // GAF paths can be long
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            // Remove trailing newline
            std::size_t len = std::strlen(buffer);
            if (len > 0 && buffer[len - 1] == '\n') {
                buffer[len - 1] = '\0';
            }

            PangenomeMapping mapping;
            if (parse_gaf_line(buffer, ref.name, mapping)) {
                if (mapping.mapq >= opts.min_mapq) {
                    result.push_back(std::move(mapping));
                }
            }
        }

        int status = pclose(pipe);
        if (status != 0) {
            // Tool exited non-zero — could be missing ref, but we continue
            // to try other refs. Only error out if ALL refs fail.
        }
    }

    // Sort by (ref_name, branch_id, path_start)
    std::sort(result.begin(), result.end(),
              [](const PangenomeMapping& a, const PangenomeMapping& b) {
                  if (a.ref_name != b.ref_name) {
                      return a.ref_name < b.ref_name;
                  }
                  if (a.branch_id != b.branch_id) {
                      return a.branch_id < b.branch_id;
                  }
                  return a.path_start < b.path_start;
              });

    return result;
}

}  // namespace branch::project
