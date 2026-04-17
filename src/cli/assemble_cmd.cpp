// BRANCH v0.2 — `branch assemble` subcommand.
//
// Reads a FASTQ file, sketches every read with the minimizer
// sketcher, runs the CPU backend's all-vs-all overlap computation,
// builds a LosslessGraph from the overlaps, and writes it out as
// BRANCH-extended GFA-1.2.
//
// This is the first end-to-end path through the BRANCH pipeline.
// Limitations: no repeat resolution, no classifier
// pass — pure "reads → graph". v0.3 layers those on top.

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "graph/delta_read.hpp"
#include "graph/graph_builder.hpp"
#include "graph/graph_io.hpp"
#include "graph/graph_filter.hpp"
#include "graph/graph_compactor.hpp"
#include "graph/lossless_graph.hpp"
#include "io/bam_reader.hpp"
#include "io/fastq_reader.hpp"

namespace branch::cli {

namespace {

void print_assemble_usage(std::ostream& os) {
    os << "branch assemble — reads -> minimizer overlap -> lossless graph\n"
          "\nUsage:\n"
          "  branch assemble (--fastq <path.fastq> | --bam <path.bam>) --out <path.gfa>\n"
          "                  [--max-reads <N>] [--max-overlaps <N>]\n"
          "                  [--fasta <path>] [--bed <path>] [--paf <path>]\n"
          "                  [--paths <path>] [--reference <path>]\n"
          "\nInput (exactly one of):\n"
          "  --fastq <path>  Plain-text FASTQ. For .gz, pipe via `zcat | branch ... --fastq /dev/stdin`.\n"
          "  --bam   <path>  BAM/CRAM/SAM via htslib. Secondary + supplementary alignments are\n"
          "                  skipped; reverse-strand reads are reverse-complemented back to the\n"
          "                  read's original sequencing orientation.\n"
          "\nOptional output sidecars:\n"
          "  --fasta <path>  Write input reads as FASTA (one record per read, 80-col wrap).\n"
          "  --fasta-consensus <path>  Write per-node consensus FASTA (post-compaction).\n"
          "  --bed   <path>  Write per-node BED (chrom=NA by default; see --ref-linear).\n"
          "  --ref-linear name=path  Linear reference for BED chrom/start/end (repeatable).\n"
          "  --paf   <path>  Write backend overlap pairs as PAF-12 (pre-graph-build).\n"
          "  --paths <path>  Write per-read graph paths as TSV: read_name\\tnode_ids\\tn_deltas.\n"
          "\nv0.2 notes:\n"
          "  - Unitig compaction enabled. No repeat resolution.\n";
}

struct Args {
    std::string fastq;
    std::string bam;
    std::string out;
    std::size_t max_reads{0};       // 0 = no cap
    std::size_t max_overlaps{10'000'000};
    std::string fasta_path;
    std::string fasta_consensus_path;
    std::string bed_path;
    std::string paf_path;
    std::string paths_path;
    std::string reference_path;  // Reference FASTA for alignment
    std::vector<std::pair<std::string, std::string>> ref_linear;  // name,path pairs
    int threads{4};
    bool ok{false};
    std::string err;
};

Args parse(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string_view k = argv[i];
        auto needs = [&](const char* name) -> const char* {
            if (i + 1 >= argc) { a.err = std::string("missing value for ") + name; return nullptr; }
            return argv[++i];
        };
        if (k == "--fastq")         { auto v = needs("--fastq");         if (!v) return a; a.fastq = v; }
        else if (k == "--bam")      { auto v = needs("--bam");           if (!v) return a; a.bam = v; }
        else if (k == "--out")      { auto v = needs("--out");           if (!v) return a; a.out = v; }
        else if (k == "--max-reads"){ auto v = needs("--max-reads");     if (!v) return a; a.max_reads = std::strtoull(v, nullptr, 10); }
        else if (k == "--max-overlaps"){ auto v = needs("--max-overlaps"); if (!v) return a; a.max_overlaps = std::strtoull(v, nullptr, 10); }
        else if (k == "--fasta")    { auto v = needs("--fasta");         if (!v) return a; a.fasta_path = v; }
        else if (k == "--fasta-consensus") { auto v = needs("--fasta-consensus"); if (!v) return a; a.fasta_consensus_path = v; }
        else if (k == "--bed")      { auto v = needs("--bed");           if (!v) return a; a.bed_path = v; }
        else if (k == "--paf")      { auto v = needs("--paf");           if (!v) return a; a.paf_path = v; }
        else if (k == "--paths")    { auto v = needs("--paths");         if (!v) return a; a.paths_path = v; }
        else if (k == "--reference"){ auto v = needs("--reference");     if (!v) return a; a.reference_path = v; }
        else if (k == "--ref-linear"){ auto v = needs("--ref-linear");   if (!v) return a;
            std::string spec(v);
            auto eq = spec.find('=');
            if (eq == std::string::npos) {
                auto slash = spec.rfind('/');
                std::string name = (slash != std::string::npos) ? spec.substr(slash + 1) : spec;
                auto dot = name.rfind('.');
                if (dot != std::string::npos) name = name.substr(0, dot);
                a.ref_linear.emplace_back(name, spec);
            } else {
                a.ref_linear.emplace_back(spec.substr(0, eq), spec.substr(eq + 1));
            }
        }
        else if (k == "--threads" || k == "-t") { auto v = needs("--threads"); if (!v) return a; a.threads = std::atoi(v); }
        else if (k == "--help" || k == "-h") { a.err = "HELP"; return a; }
        else { a.err = std::string("unknown arg: ") + std::string(k); return a; }
    }
    const bool has_fastq = !a.fastq.empty();
    const bool has_bam = !a.bam.empty();
    if (has_fastq == has_bam) {
        a.err = "exactly one of --fastq or --bam is required";
        return a;
    }
    if (a.out.empty()) {
        a.err = "--out is required";
        return a;
    }
    a.ok = true;
    return a;
}

struct ReadStore {
    std::vector<std::string> names;
    std::vector<std::string> sequences;
};

}  // namespace

int run_assemble(int argc, char** argv) {
    Args a = parse(argc, argv);
    if (!a.ok) {
        if (a.err == "HELP") { print_assemble_usage(std::cout); return 0; }
        std::cerr << "branch assemble: " << a.err << "\n\n";
        print_assemble_usage(std::cerr);
        return 2;
    }

    auto t0 = std::chrono::steady_clock::now();

    // 1. Read input (FASTQ or BAM) into memory. string_views into this
    // backing store are handed to the backend.
    ReadStore store;
    if (!a.fastq.empty()) {
        branch::io::FastqReader reader(a.fastq);
        if (!reader.ok()) {
            std::cerr << "branch assemble: cannot open " << a.fastq << "\n";
            return 3;
        }
        branch::io::FastqRecord rec;
        while (reader.next_record(rec)) {
            store.names.push_back(std::move(rec.name));
            store.sequences.push_back(std::move(rec.sequence));
            if (a.max_reads > 0 && store.names.size() >= a.max_reads) break;
        }
    } else {
        branch::io::BamReader reader;
        if (!reader.open(a.bam)) {
            std::cerr << "branch assemble: cannot open " << a.bam << "\n";
            return 3;
        }
        std::string name, seq;
        while (reader.next(name, seq)) {
            store.names.push_back(std::move(name));
            store.sequences.push_back(std::move(seq));
            name.clear();
            seq.clear();
            if (a.max_reads > 0 && store.names.size() >= a.max_reads) break;
        }
    }
    std::cerr << "[branch assemble] reads=" << store.names.size() << "\n";
    if (store.names.empty()) {
        std::cerr << "branch assemble: no reads found\n";
        return 4;
    }

    // 2. Build ReadBatch for the CPU backend.
    branch::backend::ReadBatch batch;
    batch.reads.reserve(store.names.size());
    for (std::size_t i = 0; i < store.names.size(); ++i) {
        batch.reads.push_back(branch::backend::ReadBatch::Read{
            .seq = store.sequences[i],
            .id = static_cast<branch::graph::ReadId>(i),
        });
    }

    // 3. Compute overlaps.
    branch::backend::set_cpu_overlap_threads(
        static_cast<unsigned int>(a.threads > 0 ? a.threads : 1));
    auto bk = branch::backend::make_cpu_backend();
    std::vector<branch::backend::OverlapPair> overlaps(a.max_overlaps);
    std::size_t n_overlaps = 0;
    bk.compute_overlaps(&batch, overlaps, &n_overlaps);
    overlaps.resize(n_overlaps);
    std::cerr << "[branch assemble] overlaps=" << n_overlaps << "\n";

    // 3b. Optional FASTA / PAF sidecars (before graph-build so the
    // PAF reflects the raw backend overlaps, not post-filter edges).
    if (!a.fasta_path.empty()) {
        if (!branch::graph::write_fasta(branch::graph::LosslessGraph{},
                                        store.sequences, store.names,
                                        a.fasta_path)) {
            std::cerr << "branch assemble: failed to write FASTA " << a.fasta_path << "\n";
            return 6;
        }
        std::cerr << "[branch assemble] wrote FASTA " << a.fasta_path
                  << " (" << store.sequences.size() << " records)\n";
    }
    if (!a.paf_path.empty()) {
        if (!branch::graph::write_paf(overlaps, store.sequences,
                                      store.names, a.paf_path)) {
            std::cerr << "branch assemble: failed to write PAF " << a.paf_path << "\n";
            return 7;
        }
        std::cerr << "[branch assemble] wrote PAF " << a.paf_path
                  << " (" << overlaps.size() << " pairs)\n";
    }

    // 4. Build graph.
    std::vector<branch::graph::ReadMeta> metas;
    metas.reserve(store.names.size());
    for (std::size_t i = 0; i < store.names.size(); ++i) {
        metas.push_back(branch::graph::ReadMeta{
            .read_id = static_cast<branch::graph::ReadId>(i),
            .length_bp = static_cast<std::uint32_t>(store.sequences[i].size()),
        });
    }
    auto build = branch::graph::build_graph(metas, overlaps);
    std::cerr << "[branch assemble] raw_nodes=" << build.graph.node_count()
              << " raw_edges=" << build.graph.edge_count() << "\n";
    // 4a-filter. Apply graph filtering (containment + transitive reduction).
    auto filter_stats = branch::graph::filter_graph(build.graph);
    std::cerr << "[branch assemble] filtered: edges_before=" << filter_stats.edges_before
              << " edges_after=" << filter_stats.edges_after
              << " contained_dropped=" << filter_stats.nodes_dropped_contained
              << " transitive_removed=" << filter_stats.transitive_edges_removed << "\n";

    // 4b. Compact unitigs (collapse linear chains).
    auto compaction = branch::graph::compact_unitigs_with_sequences(
        build.graph, store.sequences);
    std::cerr << "[branch assemble] compacted_nodes=" << compaction.compacted.node_count()
              << " compacted_edges=" << compaction.compacted.edge_count() << "\n";

    // Use compacted graph for output
    const auto& final_graph = compaction.compacted;

    // 5. Write GFA.
    if (!branch::graph::write_gfa(final_graph, a.out)) {
        std::cerr << "branch assemble: failed to write " << a.out << "\n";
        return 5;
    }

    // 5b. Optional BED sidecar (post-graph-build, per-node). With --ref-linear
    // the nodes' consensus sequences are mapped to the supplied references and
    // real chrom/start/end are emitted. Without --ref-linear the chrom=NA
    // placeholder row is written (backwards compat).
    if (!a.bed_path.empty()) {
        bool bed_ok = false;
        if (!a.ref_linear.empty()) {
            std::vector<branch::graph::BedLinearRef> refs;
            refs.reserve(a.ref_linear.size());
            for (const auto& [name, path] : a.ref_linear) {
                refs.push_back({name, path});
            }
            bed_ok = branch::graph::write_bed_with_refs(
                final_graph, refs, a.bed_path, a.threads);
        } else {
            bed_ok = branch::graph::write_bed(final_graph, a.bed_path);
        }
        if (!bed_ok) {
            std::cerr << "branch assemble: failed to write BED " << a.bed_path << "\n";
            return 8;
        }
        std::cerr << "[branch assemble] wrote BED " << a.bed_path
                  << " (" << final_graph.node_count() << " nodes"
                  << (a.ref_linear.empty() ? "" : ", ref-linear mapped")
                  << ")\n";
    }

    // 5b-bis. Optional consensus-FASTA sidecar (one record per Node
    // with a non-empty consensus). Emitted post-compaction so the
    // consensus reflects unitig-collapsed sequences.
    if (!a.fasta_consensus_path.empty()) {
        if (!branch::graph::write_fasta_consensus(final_graph,
                                                   a.fasta_consensus_path)) {
            std::cerr << "branch assemble: failed to write consensus FASTA "
                      << a.fasta_consensus_path << "\n";
            return 9;
        }
        std::size_t consensus_nodes = 0;
        for (const auto& n : final_graph.nodes()) {
            if (!n.consensus.empty()) ++consensus_nodes;
        }
        std::cerr << "[branch assemble] wrote consensus FASTA "
                  << a.fasta_consensus_path
                  << " (" << consensus_nodes << " nodes)\n";
    }

    // 5c. Populate ReadPaths from overlap pairs (v0.2 best-effort; deltas
    // stay empty until Node carries consensus). Emit TSV only if requested,
    // but always compute so we can print a coverage summary.
    {
        std::vector<std::uint32_t> read_lengths(store.sequences.size());
        for (std::size_t i = 0; i < store.sequences.size(); ++i) {
            read_lengths[i] = static_cast<std::uint32_t>(store.sequences[i].size());
        }
        std::vector<branch::graph::ReadPath> paths;
        branch::graph::populate_read_paths(
            store.sequences.size(), build.read_to_node, read_lengths,
            overlaps, paths);

        std::size_t mapped = 0;
        for (const auto& p : paths) {
            if (!p.path.empty()) ++mapped;
        }
        const double pct = paths.empty() ? 0.0 :
            (100.0 * static_cast<double>(mapped) /
             static_cast<double>(paths.size()));
        std::cerr << "[branch assemble] read_paths=" << paths.size()
                  << " mapped=" << mapped
                  << " (" << pct << "% of input reads mapped to >=1 node)\n";

        if (!a.paths_path.empty()) {
            std::ofstream os(a.paths_path);
            if (!os) {
                std::cerr << "branch assemble: failed to open " << a.paths_path << "\n";
                return 9;
            }
            os << "read_name\tnode_ids\tn_deltas\n";
            for (const auto& p : paths) {
                const std::string& name = (p.read_id < store.names.size())
                                              ? store.names[p.read_id]
                                              : std::string{"<unknown>"};
                os << name << '\t';
                for (std::size_t i = 0; i < p.path.size(); ++i) {
                    if (i) os << ',';
                    os << p.path[i];
                }
                os << '\t' << p.deltas.size() << '\n';
            }
            if (!os) {
                std::cerr << "branch assemble: write error on " << a.paths_path << "\n";
                return 9;
            }
            std::cerr << "[branch assemble] wrote paths TSV " << a.paths_path
                      << " (" << paths.size() << " rows)\n";
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    const double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cerr << "[branch assemble] wrote " << a.out
              << " in " << secs << "s\n";
    return 0;
}

}  // namespace branch::cli
