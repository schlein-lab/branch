// BRANCH v0.3 — `branch analyze` subcommand.
//
// Two operating modes share one command:
//
//   1. Mosdepth mode (legacy): --regions <regions.bed> [--baseline-...]
//      emits per-target CN estimates + cluster summaries to stdout.
//
//   2. Graph mode (P2.1):     --graph <path.gfa> [--out-bed <path>]
//      loads a BRANCH GFA written by `branch assemble`, detects bubbles,
//      runs coverage-conservation on the loaded graph, classifies each
//      bubble via the hierarchical disambiguator, and writes a per-bubble
//      BED (chrom/start/end/name/confidence). Either mode may be combined
//      with the other on the same invocation.

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "analysis/annotation_bed.hpp"
#include "analysis/coverage_cn.hpp"
#include "analysis/coverage_conservation.hpp"
#include "analysis/repeat_cn.hpp"
#include "common/memory.hpp"
#include "graph/read_db.hpp"
#include "classify/feature_extractor.hpp"
#include "classify/features.hpp"
#include "classify/hierarchical_disambiguator.hpp"
#include "detect/bubble_detector.hpp"
#include "graph/graph_io.hpp"
#include "graph/lossless_graph.hpp"
#include "io/mosdepth_reader.hpp"

namespace branch::cli {

namespace {

void print_analyze_usage(std::ostream& os) {
    os << "branch analyze — CN inference + bubble classification\n"
          "\nUsage (graph mode, P2.1):\n"
          "  branch analyze --graph <path.gfa> [--out-bed <path>]\n"
          "\nUsage (mosdepth mode):\n"
          "  branch analyze --regions <regions.bed>\n"
          "                 [--baseline-prefix <prefix>] [--baseline-cov <float>]\n"
          "                 [--cluster <name:prefix,prefix,...>]\n"
          "                 [--repeat-bed <path>]\n"
          "\nGraph mode consumes the BRANCH GFA written by `branch assemble`,\n"
          "runs bubble detection + coverage-conservation + hierarchical\n"
          "disambiguation, and emits one BED row per detected bubble with its\n"
          "classifier confidence in column 5.\n"
          "\nMosdepth mode: one of --baseline-cov (explicit haploid baseline) or\n"
          "--baseline-prefix (median of regions whose name begins with prefix)\n"
          "must be provided.\n"
          "\nResource caps:\n"
          "  --max-memory <size>  cap process virtual memory (e.g. 8G, 16GiB,\n"
          "                       500M). OOM becomes a clean std::bad_alloc exit\n"
          "                       (code 9) instead of SLURM / kernel SIGKILL.\n"
          "\nQuantitative VAF (lossless ReadDB):\n"
          "  --reads <path.gaf>   GAF file emitted by `branch assemble` describing\n"
          "                       per-read paths through the graph. When given,\n"
          "                       per-bubble alt support is computed by counting\n"
          "                       reads whose path traverses each alt — quantitative\n"
          "                       VAF rather than the legacy edge-overlap-sum proxy.\n";
}

// Simple flag parser: captures --key value pairs.
struct Args {
    std::string regions;
    std::string baseline_prefix;
    float baseline_cov = 0.0f;
    std::vector<std::pair<std::string, std::vector<std::string>>> clusters;
    std::string repeat_bed;
    std::string graph_path;
    std::string out_bed;
    std::string max_memory;  // e.g. "8G"; applied via setrlimit(RLIMIT_AS)
    // Optional path to a GAF file (produced by `branch assemble`) so
    // bubble VAFs can be computed from real read-path intersection
    // instead of the legacy edge-overlap-sum proxy.
    std::string reads_gaf;
    bool ok = false;
    std::string err;
};

std::vector<std::string> split_csv(const std::string& s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') {
            if (!cur.empty()) out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string_view k = argv[i];
        auto needs_val = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                a.err = std::string("missing value for ") + name;
                return nullptr;
            }
            return argv[++i];
        };
        if (k == "--regions") {
            auto v = needs_val("--regions");
            if (!v) return a;
            a.regions = v;
        } else if (k == "--baseline-prefix") {
            auto v = needs_val("--baseline-prefix");
            if (!v) return a;
            a.baseline_prefix = v;
        } else if (k == "--baseline-cov") {
            auto v = needs_val("--baseline-cov");
            if (!v) return a;
            a.baseline_cov = std::strtof(v, nullptr);
        } else if (k == "--cluster") {
            auto v = needs_val("--cluster");
            if (!v) return a;
            std::string s = v;
            auto colon = s.find(':');
            if (colon == std::string::npos) {
                a.err = "--cluster requires NAME:PREFIX1,PREFIX2,...";
                return a;
            }
            a.clusters.emplace_back(s.substr(0, colon), split_csv(s.substr(colon + 1)));
        } else if (k == "--repeat-bed") {
            auto v = needs_val("--repeat-bed");
            if (!v) return a;
            a.repeat_bed = v;
        } else if (k == "--graph") {
            auto v = needs_val("--graph");
            if (!v) return a;
            a.graph_path = v;
        } else if (k == "--out-bed") {
            auto v = needs_val("--out-bed");
            if (!v) return a;
            a.out_bed = v;
        } else if (k == "--max-memory") {
            auto v = needs_val("--max-memory");
            if (!v) return a;
            a.max_memory = v;
        } else if (k == "--reads") {
            auto v = needs_val("--reads");
            if (!v) return a;
            a.reads_gaf = v;
        } else if (k == "--help" || k == "-h") {
            a.err = "HELP";
            return a;
        } else {
            a.err = std::string("unknown arg: ") + std::string(k);
            return a;
        }
    }
    // Exactly one of --regions (mosdepth mode) or --graph (graph mode)
    // must be provided. Both may be combined in one invocation: the
    // mosdepth pass runs first, then the graph pass emits the per-bubble
    // BED.
    if (a.regions.empty() && a.graph_path.empty()) {
        a.err = "--regions (mosdepth mode) or --graph (graph mode) is required";
        return a;
    }
    if (!a.regions.empty()) {
        if (a.baseline_cov <= 0.0f && a.baseline_prefix.empty()) {
            a.err = "one of --baseline-cov or --baseline-prefix is required when --regions is set";
            return a;
        }
    }
    a.ok = true;
    return a;
}

// Convert a detect::Bubble into a classify::BubbleCandidate ready for
// feature extraction. Read-support is partitioned across the first two
// alts (branch vs. alt); extra alts on n-ary bubbles are folded into
// `read_support_alt` so the disambiguator still sees the full depth.
// `bubble_length_bp` = sum of intermediate-node length_bp on the longer
// alt path plus entry+exit flanking node lengths.
branch::classify::BubbleCandidate to_candidate(
    std::uint32_t bubble_id,
    const branch::detect::Bubble& b,
    const branch::graph::LosslessGraph& graph) {
    branch::classify::BubbleCandidate c{};
    c.bubble_id = bubble_id;
    c.chrom_id = 0;
    c.start = 0;
    c.end = 0;
    c.entry_node = b.entry;
    c.exit_node = b.exit;

    std::uint32_t support_branch = 0;
    std::uint32_t support_alt = 0;
    std::uint32_t longest_alt_len = 0;
    for (std::size_t i = 0; i < b.alts.size(); ++i) {
        const auto& alt = b.alts[i];
        if (i == 0) {
            support_branch = alt.total_read_support;
        } else {
            support_alt += alt.total_read_support;
        }
        std::uint32_t alt_len = 0;
        for (auto nid : alt.nodes) {
            if (nid < graph.node_count()) {
                alt_len += graph.node(nid).length_bp;
            }
        }
        longest_alt_len = std::max(longest_alt_len, alt_len);
    }
    c.read_support_branch = support_branch;
    c.read_support_alt = support_alt;

    std::uint32_t entry_len = 0, exit_len = 0;
    if (b.entry < graph.node_count()) entry_len = graph.node(b.entry).length_bp;
    if (b.exit  < graph.node_count()) exit_len  = graph.node(b.exit).length_bp;
    c.bubble_length_bp = entry_len + longest_alt_len + exit_len;
    return c;
}

}  // namespace

int run_analyze(int argc, char** argv) {
    Args a = parse_args(argc, argv);
    if (!a.ok) {
        if (a.err == "HELP") {
            print_analyze_usage(std::cout);
            return 0;
        }
        std::cerr << "branch analyze: " << a.err << "\n\n";
        print_analyze_usage(std::cerr);
        return 2;
    }

    if (!a.max_memory.empty()) {
        auto bytes = branch::common::parse_memory_size(a.max_memory);
        if (!bytes) {
            std::cerr << "branch analyze: cannot parse --max-memory '"
                      << a.max_memory << "'\n";
            return 2;
        }
        if (!branch::common::set_memory_budget(*bytes)) {
            std::cerr << "branch analyze: setrlimit(RLIMIT_AS) failed; "
                         "continuing without budget\n";
        } else {
            std::cerr << "[branch analyze] memory budget = "
                      << branch::common::format_bytes(*bytes) << "\n";
        }
    }

    // --- Mosdepth mode ----------------------------------------------------
    std::vector<branch::analysis::CNEstimate> est;  // shared with graph mode
    if (!a.regions.empty()) {
        auto regions = branch::io::read_mosdepth_regions(a.regions);
        if (regions.empty()) {
            std::cerr << "branch analyze: failed to read regions from " << a.regions << "\n";
            return 3;
        }

        float baseline = a.baseline_cov;
        if (baseline <= 0.0f) {
            std::vector<branch::io::RegionCoverage> baseline_regions;
            for (const auto& r : regions) {
                if (r.name.size() >= a.baseline_prefix.size() &&
                    r.name.compare(0, a.baseline_prefix.size(), a.baseline_prefix) == 0) {
                    baseline_regions.push_back(r);
                }
            }
            if (baseline_regions.empty()) {
                std::cerr << "branch analyze: no baseline regions matched prefix '"
                          << a.baseline_prefix << "'\n";
                return 4;
            }
            baseline = branch::analysis::median_baseline_coverage(baseline_regions);
        }

        std::cout << "# branch analyze\n";
        std::cout << "# baseline_haploid_coverage=" << baseline << "\n";
        std::cout << "# regions=" << regions.size() << "\n";
        std::cout << "chrom\tstart\tend\tname\tmean_cov\trelative_cn\n";
        est = branch::analysis::estimate_cn(regions, baseline);
        for (const auto& e : est) {
            std::cout << e.chrom << '\t' << e.start << '\t' << e.end << '\t'
                      << e.name << '\t' << e.mean_coverage << '\t' << e.relative_cn << '\n';
        }

        for (const auto& [cluster_name, prefixes] : a.clusters) {
            auto s = branch::analysis::summarise_cluster(est, cluster_name, prefixes);
            std::cout << "# cluster " << s.cluster_name
                      << " members=" << s.member_count
                      << " total_rel_cn=" << s.total_relative_cn
                      << " mean_per_member=" << s.mean_per_member << "\n";
        }

        if (!a.repeat_bed.empty()) {
            auto repeats = branch::analysis::parse_repeat_bed(a.repeat_bed);
            std::map<std::string, std::vector<uint32_t>> per_base_cov;
            for (const auto& r : regions) {
                auto& vec = per_base_cov[r.chrom];
                if (vec.size() < r.end) vec.resize(r.end, 0);
                uint32_t cov_val = static_cast<uint32_t>(r.mean_coverage + 0.5);
                for (size_t i = r.start; i < r.end && i < vec.size(); ++i) {
                    vec[i] = cov_val;
                }
            }
            auto cn_stats = branch::analysis::compute_family_cn(repeats, per_base_cov, baseline);
            std::string out_path = a.regions + ".repeat_cn.tsv";
            branch::analysis::write_cn_tsv(cn_stats, out_path);
            std::cerr << "# repeat_cn written to " << out_path << "\n";
        }
    }

    // --- Graph mode (P2.1) ------------------------------------------------
    //
    // Load the GFA emitted by `branch assemble`, run bubble detection on
    // the loaded graph, feed the real graph into the coverage-conservation
    // solver, and classify each bubble via the hierarchical disambiguator.
    if (!a.graph_path.empty()) {
        branch::graph::LosslessGraph graph;
        if (!branch::graph::read_gfa(graph, a.graph_path)) {
            std::cerr << "branch analyze: failed to read graph from "
                      << a.graph_path << "\n";
            return 5;
        }
        std::cerr << "[branch analyze] graph loaded nodes="
                  << graph.node_count() << " edges=" << graph.edge_count()
                  << "\n";

        // Structural bubble enumeration. Empty on linear assemblies.
        // If --reads <gaf> was passed, additionally count real reads
        // per alt by intersecting bubble paths with the ReadDB; this
        // populates AltPath::reads_traversing for downstream
        // quantitative VAF.
        std::vector<branch::detect::Bubble> bubbles;
        bool reads_gaf_loaded = false;
        branch::graph::ReadDB read_db;
        if (!a.reads_gaf.empty()) {
            if (!read_db.read_gaf(a.reads_gaf, graph)) {
                std::cerr << "branch analyze: cannot read GAF " << a.reads_gaf
                          << "; falling back to topological VAF\n";
            } else {
                reads_gaf_loaded = true;
                std::cerr << "[branch analyze] loaded " << read_db.size()
                          << " read paths from " << a.reads_gaf << "\n";
            }
        }
        if (reads_gaf_loaded) {
            bubbles = branch::detect::detect_bubbles_with_reads(graph, read_db);
        } else {
            bubbles = branch::detect::detect_bubbles(graph);
        }
        std::cerr << "[branch analyze] bubbles=" << bubbles.size() << "\n";

        // Run the real solver on the loaded graph.
        auto report = branch::analysis::run_conservation(graph, bubbles, est);
        std::cout << "# conservation iterations=" << report.iterations_used
                  << " converged=" << (report.converged ? 1 : 0)
                  << " node_violations=" << report.per_node_violations.size()
                  << " bubble_violations=" << report.per_bubble_violations.size()
                  << " residual=" << report.final_global_residual << "\n";

        // Classify every bubble via the hierarchical disambiguator and
        // stage one BedEntry per bubble. The per-bubble SNP matrix used
        // by mixed_decomposer is not yet wired through assemble->analyze
        // (P2.2 work), so we pass an empty vector — disambiguate_hierarchical
        // degrades to the single-level disambiguate() call in that case.
        //
        // Each bubble's classification + BedEntry build is independent of
        // every other bubble's, so we pre-size the output vector and fan
        // the loop out with OpenMP — writes to distinct indices never
        // race. Order of output matches input bubbles[] so determinism
        // is preserved without an explicit sort.
        std::vector<branch::analysis::BedEntry> bed_entries(bubbles.size());
        const std::vector<std::vector<char>> empty_snp_matrix;
        const branch::classify::DisambiguatorConfig dcfg{};

        #pragma omp parallel for schedule(dynamic, 64)
        for (std::size_t i = 0; i < bubbles.size(); ++i) {
            const auto& b = bubbles[i];
            auto cand = to_candidate(static_cast<std::uint32_t>(i), b, graph);
            cand.features = branch::classify::extract_features(cand, graph);

            auto r = branch::classify::disambiguate_hierarchical(
                cand.features, empty_snp_matrix, dcfg);

            auto& entry = bed_entries[i];
            entry.chrom = "graph";
            entry.name = "bubble_" + std::to_string(i) + "_n" +
                         std::to_string(b.entry) + "_" + std::to_string(b.exit);
            entry.start = 0;
            entry.end = cand.bubble_length_bp;
            entry.confidence = r.confidence;
            entry.alt_read_supports.reserve(b.alts.size());
            for (const auto& alt : b.alts) {
                // When --reads was passed, alt.reads_traversing was
                // populated by detect_bubbles_with_reads with the
                // exact set of reads whose path traverses this alt.
                // That is the quantity we want for VAF; fall back to
                // the legacy edge-sum proxy when no GAF was provided.
                const std::uint32_t support = !alt.reads_traversing.empty()
                    ? static_cast<std::uint32_t>(alt.reads_traversing.size())
                    : alt.total_read_support;
                entry.alt_read_supports.push_back(support);
            }
        }
        std::cout << "# bubbles_classified=" << bed_entries.size() << "\n";

        // Diagnostic: surface min/max/mean confidence so callers can gate
        // without re-parsing the BED. NaN means "no bubbles detected".
        if (!bed_entries.empty()) {
            float min_c = 1.0f, max_c = 0.0f, sum_c = 0.0f;
            for (const auto& e : bed_entries) {
                const float c = e.confidence.value_or(std::nanf(""));
                if (std::isnan(c)) continue;
                min_c = std::min(min_c, c);
                max_c = std::max(max_c, c);
                sum_c += c;
            }
            std::cout << "# confidence_min=" << min_c
                      << " max=" << max_c
                      << " mean=" << (sum_c / static_cast<float>(bed_entries.size()))
                      << "\n";
        }

        if (!a.out_bed.empty()) {
            std::ofstream ofs(a.out_bed);
            if (!ofs) {
                std::cerr << "branch analyze: cannot open " << a.out_bed
                          << "\n";
                return 6;
            }
            for (const auto& e : bed_entries) {
                branch::analysis::write_bed_entry(ofs, e);
            }
            if (!ofs) {
                std::cerr << "branch analyze: write error on " << a.out_bed
                          << "\n";
                return 6;
            }
            std::cerr << "[branch analyze] wrote BED " << a.out_bed
                      << " (" << bed_entries.size() << " bubbles)\n";
        }
    }

    return 0;
}

}  // namespace branch::cli
