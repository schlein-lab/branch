// BRANCH v0.2 — `branch assemble` subcommand.
//
// Reads a FASTQ file, sketches every read with the minimizer
// sketcher, runs the CPU backend's all-vs-all overlap computation,
// builds a LosslessGraph from the overlaps, and writes it out as
// BRANCH-extended GFA-1.2.
//
// This is the first end-to-end path through the BRANCH pipeline.
// Limitations: no compaction, no repeat resolution, no classifier
// pass — pure "reads → graph". v0.3 layers those on top.

#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "backend/backend_vtable.hpp"
#include "backend/cpu_backend.hpp"
#include "graph/graph_builder.hpp"
#include "graph/graph_io.hpp"
#include "graph/lossless_graph.hpp"
#include "io/fastq_reader.hpp"

namespace branch::cli {

namespace {

void print_assemble_usage(std::ostream& os) {
    os << "branch assemble — reads -> minimizer overlap -> lossless graph\n"
          "\nUsage:\n"
          "  branch assemble --fastq <path.fastq> --out <path.gfa>\n"
          "                  [--max-reads <N>] [--max-overlaps <N>]\n"
          "\nv0.2 notes:\n"
          "  - Accepts plain-text FASTQ only. Use `zcat x.fastq.gz | branch assemble --fastq /dev/stdin ...`\n"
          "  - No compaction, no repeat resolution. Raw read-level graph.\n";
}

struct Args {
    std::string fastq;
    std::string out;
    std::size_t max_reads{0};       // 0 = no cap
    std::size_t max_overlaps{10'000'000};
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
        else if (k == "--out")      { auto v = needs("--out");           if (!v) return a; a.out = v; }
        else if (k == "--max-reads"){ auto v = needs("--max-reads");     if (!v) return a; a.max_reads = std::strtoull(v, nullptr, 10); }
        else if (k == "--max-overlaps"){ auto v = needs("--max-overlaps"); if (!v) return a; a.max_overlaps = std::strtoull(v, nullptr, 10); }
        else if (k == "--help" || k == "-h") { a.err = "HELP"; return a; }
        else { a.err = std::string("unknown arg: ") + std::string(k); return a; }
    }
    if (a.fastq.empty() || a.out.empty()) {
        a.err = "--fastq and --out are required";
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

    // 1. Read FASTQ into memory (string_views into this backing store
    // are handed to the backend).
    branch::io::FastqReader reader(a.fastq);
    if (!reader.ok()) {
        std::cerr << "branch assemble: cannot open " << a.fastq << "\n";
        return 3;
    }
    ReadStore store;
    branch::io::FastqRecord rec;
    while (reader.next_record(rec)) {
        store.names.push_back(std::move(rec.name));
        store.sequences.push_back(std::move(rec.sequence));
        if (a.max_reads > 0 && store.names.size() >= a.max_reads) break;
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
    auto bk = branch::backend::make_cpu_backend();
    std::vector<branch::backend::OverlapPair> overlaps(a.max_overlaps);
    std::size_t n_overlaps = 0;
    bk.compute_overlaps(&batch, overlaps, &n_overlaps);
    overlaps.resize(n_overlaps);
    std::cerr << "[branch assemble] overlaps=" << n_overlaps << "\n";

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
    std::cerr << "[branch assemble] nodes=" << build.graph.node_count()
              << " edges=" << build.graph.edge_count() << "\n";

    // 5. Write GFA.
    if (!branch::graph::write_gfa(build.graph, a.out)) {
        std::cerr << "branch assemble: failed to write " << a.out << "\n";
        return 5;
    }
    auto t1 = std::chrono::steady_clock::now();
    const double secs = std::chrono::duration<double>(t1 - t0).count();
    std::cerr << "[branch assemble] wrote " << a.out
              << " in " << secs << "s\n";
    return 0;
}

}  // namespace branch::cli
