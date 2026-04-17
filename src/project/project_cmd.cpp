// BRANCH v0.4 — `branch project` subcommand.
//
// Three-layer reference projection:
// 1. Linear (CHM13 + GRCh38 via minimap2)
// 2. Pangenome (HPRC v1.1 via minigraph/GraphAligner)
// 3. Somatic delta (edit distance to closest HPRC path)
//
// This is a stub implementation. Full functionality in v0.4.

#include "project_cmd.hpp"

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

namespace branch::cli {

namespace {

void print_project_usage(std::ostream& os) {
    os << "branch project — three-layer reference projection (v0.4 in development)\n"
          "\nUsage:\n"
          "  branch project --gfa <path.gfa> --fasta <path.fasta>\n"
          "                 --ref-linear <name=path> [--ref-linear <name=path>]\n"
          "                 --ref-pangenome <path.gbz>\n"
          "                 --out-prefix <path>\n"
          "                 [--mapper linear=minimap2,pangenome=minigraph]\n"
          "                 [--threads <N>] [--emit-bam] [--min-mapq <N>]\n"
          "\nRequired:\n"
          "  --gfa <path>           GFA from branch assemble (nodes = branches)\n"
          "  --fasta <path>         Branch consensus FASTA\n"
          "  --ref-linear <spec>    Linear reference (name=path, repeatable)\n"
          "                         Example: --ref-linear CHM13=/refs/chm13.mmi\n"
          "  --ref-pangenome <path> Pangenome GBZ file (repeatable)\n"
          "  --out-prefix <path>    Output prefix for all result files\n"
          "\nOptional:\n"
          "  --mapper <spec>        Mapper selection (default: linear=minimap2,pangenome=minigraph)\n"
          "                         Alternatives: pangenome=graphaligner\n"
          "  --threads <N>          Thread count (default: available cores)\n"
          "  --emit-bam             Also emit BAM (default: off)\n"
          "  --min-mapq <N>         Minimum MAPQ for annotation (default: 20)\n"
          "  --help                 Show this help\n"
          "\nOutput files:\n"
          "  <prefix>.linear.paf       PAF alignments to linear references\n"
          "  <prefix>.pangenome.gaf    GAF alignments to pangenome\n"
          "  <prefix>.somatic.vcf      Somatic variant calls\n"
          "  <prefix>.branch-report.json  Detailed per-branch analysis\n"
          "\nSee docs/branch-project-design.md for full specification.\n";
}

/// Parsed arguments for branch project.
struct ProjectArgs {
    std::string gfa;
    std::string fasta;
    std::vector<std::pair<std::string, std::string>> ref_linear;  // name, path
    std::vector<std::string> ref_pangenome;
    std::string out_prefix;
    std::string mapper_linear{"minimap2"};
    std::string mapper_pangenome{"minigraph"};
    int threads{0};  // 0 = auto-detect
    int min_mapq{20};
    bool emit_bam{false};
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

/// Parse --mapper spec (e.g., "linear=minimap2,pangenome=graphaligner").
void parse_mapper_spec(const std::string& spec, ProjectArgs& args) {
    std::string::size_type pos = 0;
    while (pos < spec.size()) {
        auto comma = spec.find(',', pos);
        std::string part = (comma != std::string::npos)
            ? spec.substr(pos, comma - pos)
            : spec.substr(pos);

        auto eq = part.find('=');
        if (eq != std::string::npos) {
            std::string key = part.substr(0, eq);
            std::string val = part.substr(eq + 1);
            if (key == "linear") args.mapper_linear = val;
            else if (key == "pangenome") args.mapper_pangenome = val;
        }

        if (comma == std::string::npos) break;
        pos = comma + 1;
    }
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

        if (k == "--gfa") {
            auto v = needs_val("--gfa");
            if (!v) return a;
            a.gfa = v;
        } else if (k == "--fasta") {
            auto v = needs_val("--fasta");
            if (!v) return a;
            a.fasta = v;
        } else if (k == "--ref-linear") {
            auto v = needs_val("--ref-linear");
            if (!v) return a;
            a.ref_linear.push_back(parse_ref_spec(v));
        } else if (k == "--ref-pangenome") {
            auto v = needs_val("--ref-pangenome");
            if (!v) return a;
            a.ref_pangenome.push_back(v);
        } else if (k == "--out-prefix") {
            auto v = needs_val("--out-prefix");
            if (!v) return a;
            a.out_prefix = v;
        } else if (k == "--mapper") {
            auto v = needs_val("--mapper");
            if (!v) return a;
            parse_mapper_spec(v, a);
        } else if (k == "--threads") {
            auto v = needs_val("--threads");
            if (!v) return a;
            a.threads = std::stoi(v);
        } else if (k == "--min-mapq") {
            auto v = needs_val("--min-mapq");
            if (!v) return a;
            a.min_mapq = std::stoi(v);
        } else if (k == "--emit-bam") {
            a.emit_bam = true;
        } else if (k == "--help" || k == "-h") {
            a.err = "HELP";
            return a;
        } else {
            a.err = std::string("unknown argument: ") + std::string(k);
            return a;
        }
    }

    // Validation
    if (a.gfa.empty()) {
        a.err = "--gfa is required";
        return a;
    }
    if (a.fasta.empty()) {
        a.err = "--fasta is required";
        return a;
    }
    if (a.ref_linear.empty()) {
        a.err = "at least one --ref-linear is required";
        return a;
    }
    if (a.ref_pangenome.empty()) {
        a.err = "at least one --ref-pangenome is required";
        return a;
    }
    if (a.out_prefix.empty()) {
        a.err = "--out-prefix is required";
        return a;
    }

    a.ok = true;
    return a;
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

    // Stub: print configuration and exit
    std::cerr << "[branch project] Configuration:\n"
              << "  GFA:           " << args.gfa << "\n"
              << "  FASTA:         " << args.fasta << "\n"
              << "  Linear refs:   " << args.ref_linear.size() << "\n";
    for (const auto& [name, path] : args.ref_linear) {
        std::cerr << "    " << name << " -> " << path << "\n";
    }
    std::cerr << "  Pangenome refs: " << args.ref_pangenome.size() << "\n";
    for (const auto& p : args.ref_pangenome) {
        std::cerr << "    " << p << "\n";
    }
    std::cerr << "  Out prefix:    " << args.out_prefix << "\n"
              << "  Mapper linear: " << args.mapper_linear << "\n"
              << "  Mapper pangenome: " << args.mapper_pangenome << "\n"
              << "  Threads:       " << (args.threads > 0 ? std::to_string(args.threads) : "auto") << "\n"
              << "  Min MAPQ:      " << args.min_mapq << "\n"
              << "  Emit BAM:      " << (args.emit_bam ? "yes" : "no") << "\n";

    std::cerr << "\nbranch project: not yet implemented (v0.4 in development)\n";
    return 2;
}

}  // namespace branch::cli
