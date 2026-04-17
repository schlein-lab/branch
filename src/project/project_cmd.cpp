// BRANCH v0.4 — `branch project` subcommand.
//
// Linear-layer implementation (v0.4.1):
// - Maps branch consensus FASTA against linear references (CHM13, GRCh38)
// - Shell-out to minimap2 with asm20 preset
// - Generates JSON branch report with per-branch mappings
//
// Future layers (v0.4.2, v0.4.3):
// - Pangenome (HPRC v1.1 via minigraph/GraphAligner)
// - Somatic delta (edit distance to closest HPRC path)

#include "project_cmd.hpp"
#include "project/linear_mapper.hpp"
#include "project/project_report.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace branch::cli {

namespace {

void print_project_usage(std::ostream& os) {
    os << "branch project — linear reference projection (v0.4.1)\n"
          "\nUsage:\n"
          "  branch project --fasta <path.fasta> --ref-linear <name=path>\n"
          "                 --out-prefix <path> [--threads <N>] [--min-mapq <N>]\n"
          "\nRequired:\n"
          "  --fasta <path>         Branch consensus FASTA (from branch assemble)\n"
          "  --ref-linear <spec>    Linear reference (name=path, repeatable)\n"
          "                         Example: --ref-linear CHM13=/refs/chm13.mmi\n"
          "  --out-prefix <path>    Output prefix for result files\n"
          "\nOptional:\n"
          "  --threads <N>          Thread count for minimap2 (default: 4)\n"
          "  --min-mapq <N>         Minimum MAPQ for annotation (default: 0)\n"
          "                         Branches with no mapping >= min-mapq are marked unannotated.\n"
          "  --help                 Show this help\n"
          "\nOutput files:\n"
          "  <prefix>.branch-report.json  Per-branch analysis with linear mappings\n"
          "\nExample:\n"
          "  branch project --fasta branches.fa \\\n"
          "                 --ref-linear CHM13=/refs/chm13.mmi \\\n"
          "                 --ref-linear GRCh38=/refs/grch38.mmi \\\n"
          "                 --out-prefix output/sample1 \\\n"
          "                 --threads 8 --min-mapq 20\n"
          "\nNote: Pangenome layer and somatic delta will be added in v0.4.2/v0.4.3.\n";
}

/// Parsed arguments for branch project.
struct ProjectArgs {
    std::string fasta;
    std::vector<std::pair<std::string, std::string>> ref_linear;  // name, path
    std::string out_prefix;
    int threads{4};
    int min_mapq{0};
    bool ok{false};
    std::string err;
};

/// Parse name=path format for --ref-linear.
std::pair<std::string, std::string> parse_ref_spec(const std::string& spec) {
    auto eq = spec.find('=');
    if (eq == std::string::npos) {
        // No name given, derive from filename
        auto slash = spec.rfind('/');
        std::string name = (slash != std::string::npos) ? spec.substr(slash + 1) : spec;
        auto dot = name.rfind('.');
        if (dot != std::string::npos) name = name.substr(0, dot);
        return {name, spec};
    }
    return {spec.substr(0, eq), spec.substr(eq + 1)};
}

ProjectArgs parse_project_args(int argc, char** argv) {
    ProjectArgs a;

    for (int i = 1; i < argc; ++i) {
        std::string_view k = argv[i];

        auto needs_val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                a.err = std::string("missing value for ") + name;
                return nullptr;
            }
            return argv[++i];
        };

        if (k == "--fasta") {
            auto v = needs_val("--fasta");
            if (!v) return a;
            a.fasta = v;
        } else if (k == "--ref-linear") {
            auto v = needs_val("--ref-linear");
            if (!v) return a;
            a.ref_linear.push_back(parse_ref_spec(v));
        } else if (k == "--out-prefix") {
            auto v = needs_val("--out-prefix");
            if (!v) return a;
            a.out_prefix = v;
        } else if (k == "--threads") {
            auto v = needs_val("--threads");
            if (!v) return a;
            a.threads = std::stoi(v);
        } else if (k == "--min-mapq") {
            auto v = needs_val("--min-mapq");
            if (!v) return a;
            a.min_mapq = std::stoi(v);
        } else if (k == "--help" || k == "-h") {
            a.err = "HELP";
            return a;
        } else {
            a.err = std::string("unknown argument: ") + std::string(k);
            return a;
        }
    }

    // Validation
    if (a.fasta.empty()) {
        a.err = "--fasta is required";
        return a;
    }
    if (a.ref_linear.empty()) {
        a.err = "at least one --ref-linear is required";
        return a;
    }
    if (a.out_prefix.empty()) {
        a.err = "--out-prefix is required";
        return a;
    }

    a.ok = true;
    return a;
}

/// Read FASTA file and extract branch IDs and lengths.
/// Returns: map of branch_id -> length_bp
std::map<std::string, std::int64_t> read_fasta_entries(
    const std::string& fasta_path, std::string* err_out) {

    std::map<std::string, std::int64_t> entries;

    std::ifstream fs(fasta_path);
    if (!fs) {
        if (err_out) {
            *err_out = "cannot open fasta file: " + fasta_path;
        }
        return entries;
    }

    std::string line;
    std::string current_id;
    std::int64_t current_len = 0;

    while (std::getline(fs, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // Save previous entry
            if (!current_id.empty()) {
                entries[current_id] = current_len;
            }
            // Parse new header: >branch_id [description...]
            std::size_t space_pos = line.find_first_of(" \t", 1);
            if (space_pos != std::string::npos) {
                current_id = line.substr(1, space_pos - 1);
            } else {
                current_id = line.substr(1);
            }
            current_len = 0;
        } else {
            // Sequence line
            current_len += static_cast<std::int64_t>(line.size());
        }
    }

    // Save last entry
    if (!current_id.empty()) {
        entries[current_id] = current_len;
    }

    return entries;
}

/// Derive sample name from out_prefix (basename without extension).
std::string derive_sample_name(const std::string& out_prefix) {
    auto slash = out_prefix.rfind('/');
    std::string name = (slash != std::string::npos)
        ? out_prefix.substr(slash + 1)
        : out_prefix;
    // Remove common extensions if present
    for (const char* ext : {".branch-report", ".json", ".txt"}) {
        auto pos = name.rfind(ext);
        if (pos != std::string::npos && pos + std::strlen(ext) == name.size()) {
            name = name.substr(0, pos);
        }
    }
    return name.empty() ? "unknown" : name;
}

}  // namespace

int run_project(int argc, char** argv) {
    ProjectArgs args = parse_project_args(argc, argv);

    if (!args.ok) {
        if (args.err == "HELP") {
            print_project_usage(std::cout);
            return 0;
        }
        std::cerr << "branch project: " << args.err << "\n\n";
        print_project_usage(std::cerr);
        return 2;
    }

    // 1. Read FASTA to get branch IDs and lengths
    std::string fasta_err;
    auto fasta_entries = read_fasta_entries(args.fasta, &fasta_err);
    if (fasta_entries.empty()) {
        std::cerr << "branch project: " << fasta_err << "\n";
        return 3;
    }

    std::cerr << "[branch project] loaded " << fasta_entries.size()
              << " branches from " << args.fasta << "\n";

    // 2. Build LinearRef vector and run mapping
    std::vector<project::LinearRef> refs;
    for (const auto& [name, path] : args.ref_linear) {
        refs.push_back({name, path});
    }

    project::LinearMapOptions opts;
    opts.threads = args.threads;
    opts.min_mapq = 0;  // Get all mappings, filter later for unannotated flag

    std::string map_err;
    auto mappings = project::map_branches_linear(args.fasta, refs, opts, &map_err);

    if (!map_err.empty()) {
        std::cerr << "branch project: mapping error: " << map_err << "\n";
        return 3;
    }

    std::cerr << "[branch project] got " << mappings.size()
              << " mappings across " << refs.size() << " references\n";

    // 3. Group mappings by branch_id
    std::map<std::string, std::vector<project::LinearMapping>> mappings_by_branch;
    for (auto& m : mappings) {
        mappings_by_branch[m.branch_id].push_back(std::move(m));
    }

    // 4. Build BranchEntry for each branch
    std::vector<project::BranchEntry> branches;
    branches.reserve(fasta_entries.size());

    std::size_t unannotated_count = 0;
    for (const auto& [branch_id, length] : fasta_entries) {
        project::BranchEntry entry;
        entry.branch_id = branch_id;
        entry.length_bp = length;
        entry.vaf = -1.0;  // Unknown
        entry.coverage = -1.0;  // Unknown

        // Get mappings for this branch
        auto it = mappings_by_branch.find(branch_id);
        if (it != mappings_by_branch.end()) {
            entry.linear_mappings = std::move(it->second);
        }

        // Check if unannotated: no mapping with mapq >= min_mapq
        bool has_good_mapping = false;
        for (const auto& m : entry.linear_mappings) {
            if (m.mapq >= args.min_mapq) {
                has_good_mapping = true;
                break;
            }
        }
        entry.unannotated = !has_good_mapping;
        if (entry.unannotated) {
            ++unannotated_count;
        }

        branches.push_back(std::move(entry));
    }

    // 5. Write JSON report
    std::string report_path = args.out_prefix + ".branch-report.json";
    std::string sample_name = derive_sample_name(args.out_prefix);
    std::string report_err;

    if (!project::write_branch_report_json(report_path, sample_name, branches, &report_err)) {
        std::cerr << "branch project: " << report_err << "\n";
        return 4;
    }

    std::cerr << "[branch project] mapped " << fasta_entries.size()
              << " branches across " << refs.size() << " refs, "
              << unannotated_count << " unannotated, report=" << report_path << "\n";

    return 0;
}

}  // namespace branch::cli
