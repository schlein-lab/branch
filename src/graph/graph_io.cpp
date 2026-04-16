// BRANCH v0.1 — GFA-1.2 I/O implementation.
// BRANCH v0.2 — FASTA / BED / PAF sidecar writers.

#include "graph/graph_io.hpp"

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace branch::graph {

bool write_gfa(const LosslessGraph& graph, std::ostream& out) {
    out << "H\tVN:Z:1.2\n";

    for (const auto& node : graph.nodes()) {
        out << "S\t" << node.id
            << "\t*"  // placeholder sequence
            << "\tLN:i:" << node.length_bp
            << "\tRC:i:" << node.read_support
            << "\tCN:i:" << node.copy_count
            << "\tCV:f:" << std::fixed << std::setprecision(3) << node.copy_count_confidence;
        if (!node.consensus.empty()) {
            out << "\tCS:Z:" << node.consensus;
        }
        out << '\n';
    }

    for (std::size_t i = 0; i < graph.edges().size(); ++i) {
        const auto& edge = graph.edges()[i];
        // GFA L-line: L <from> <from_orient> <to> <to_orient> <cigar>
        //
        // v0.2 strand-aware output: bit 4 (0x10) of edge.flags marks an
        // opposite-strand overlap. When set we emit "+/-" so the target
        // node enters on its reverse-complement side; otherwise "+/+".
        const bool opp_strand = (edge.flags & 0x10) != 0;
        const char to_orient = opp_strand ? '-' : '+';
        out << "L\t" << edge.from << "\t+\t" << edge.to << '\t' << to_orient
            << "\t0M"
            << "\tVC:i:" << edge.read_support
            << "\tVF:f:" << std::fixed << std::setprecision(4) << edge.vaf
            << "\tCT:f:" << std::fixed << std::setprecision(4) << edge.vaf_confidence;

        // BT:A — branch-type character
        if (edge.flags & 0x01) {
            out << "\tBT:A:B";  // Branch
        } else if (edge.flags & 0x02) {
            out << "\tBT:A:D";  // Duplication
        } else if (edge.flags & 0x04) {
            out << "\tBT:A:M";  // Mixed
        } else if (edge.flags & 0x08) {
            out << "\tBT:A:N";  // NonSeparable
        }

        out << '\n';
    }

    return out.good();
}

bool write_gfa(const LosslessGraph& graph, const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    return write_gfa(graph, ofs);
}

namespace {

// Parse "TAG:TYPE:VALUE" optional fields from GFA lines.
bool parse_tag_int(std::string_view token, std::string_view tag, std::uint32_t& out) {
    if (token.size() > tag.size() && token.substr(0, tag.size()) == tag) {
        auto val = token.substr(tag.size());
        auto r = std::from_chars(val.data(), val.data() + val.size(), out);
        return r.ec == std::errc{};
    }
    return false;
}

bool parse_tag_float(std::string_view token, std::string_view tag, float& out) {
    if (token.size() > tag.size() && token.substr(0, tag.size()) == tag) {
        auto val = token.substr(tag.size());
        // std::from_chars for float is not always available pre-C++26;
        // fall back to strtof for portability.
        char* end = nullptr;
        std::string s(val);
        out = std::strtof(s.c_str(), &end);
        return end != s.c_str();
    }
    return false;
}

bool parse_tag_string(std::string_view token, std::string_view tag, std::string& out) {
    if (token.size() > tag.size() && token.substr(0, tag.size()) == tag) {
        out = std::string(token.substr(tag.size()));
        return true;
    }
    return false;
}

}  // namespace

bool read_gfa(LosslessGraph& graph, std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        char record_type{};
        iss >> record_type;

        if (record_type == 'S') {
            std::uint32_t id{};
            std::string seq;
            iss >> id >> seq;

            std::uint32_t len = 0, rc = 0, cn = 1;
            float cv = 1.0f;
            std::string cs;
            std::string token;
            while (iss >> token) {
                parse_tag_int(token, "LN:i:", len);
                parse_tag_int(token, "RC:i:", rc);
                parse_tag_int(token, "CN:i:", cn);
                parse_tag_float(token, "CV:f:", cv);
                parse_tag_string(token, "CS:Z:", cs);
            }

            // Grow node vector to fit id (nodes must arrive in order or
            // sparse IDs are auto-filled with empties).
            while (graph.node_count() <= id) {
                graph.add_node(0);
            }
            auto& n = graph.node(static_cast<NodeId>(id));
            n.length_bp = len;
            n.read_support = rc;
            n.copy_count = cn;
            n.copy_count_confidence = cv;
            n.consensus = std::move(cs);

        } else if (record_type == 'L') {
            std::uint32_t from{}, to{};
            std::string orient_from, orient_to, cigar;
            iss >> from >> orient_from >> to >> orient_to >> cigar;

            std::uint32_t vc = 0;
            float vf = 1.0f, ct = 0.0f;
            std::uint16_t flags = 0;
            std::string token;
            while (iss >> token) {
                parse_tag_int(token, "VC:i:", vc);
                parse_tag_float(token, "VF:f:", vf);
                parse_tag_float(token, "CT:f:", ct);
                if (token.size() == 6 && token.substr(0, 5) == "BT:A:") {
                    char c = token[5];
                    if (c == 'B') flags |= 0x01;
                    else if (c == 'D') flags |= 0x02;
                    else if (c == 'M') flags |= 0x04;
                    else if (c == 'N') flags |= 0x08;
                }
            }
            graph.add_edge(static_cast<NodeId>(from), static_cast<NodeId>(to), vc);
            auto& e = graph.edges().back();  // just added
            const_cast<Edge&>(e).vaf = vf;
            const_cast<Edge&>(e).vaf_confidence = ct;
            const_cast<Edge&>(e).flags = flags;
        }
        // H-lines and unknown record types are silently skipped.
    }
    return !in.bad();
}

bool read_gfa(LosslessGraph& graph, const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    return read_gfa(graph, ifs);
}

// =====================================================================
// FASTA writer
// =====================================================================
//
// TODO(v0.3): write node consensus instead of raw reads once the
// compactor stores a per-Node consensus sequence.

bool write_fasta(const LosslessGraph& /*graph*/,
                 const std::vector<std::string>& sequences,
                 const std::vector<std::string>& names,
                 std::ostream& out) {
    if (sequences.size() != names.size()) return false;

    constexpr std::size_t kWrap = 80;
    for (std::size_t i = 0; i < sequences.size(); ++i) {
        out << '>' << names[i] << '\n';
        const auto& seq = sequences[i];
        for (std::size_t p = 0; p < seq.size(); p += kWrap) {
            const std::size_t len = std::min(kWrap, seq.size() - p);
            out.write(seq.data() + p, static_cast<std::streamsize>(len));
            out << '\n';
        }
    }
    return out.good();
}

bool write_fasta(const LosslessGraph& graph,
                 const std::vector<std::string>& sequences,
                 const std::vector<std::string>& names,
                 const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    return write_fasta(graph, sequences, names, ofs);
}

// =====================================================================
// BED writer
// =====================================================================
//
// TODO(v0.3): use real chrom/start/end once reference alignment exists.
// BED6 layout: chrom  start  end  name  score  strand

bool write_bed(const LosslessGraph& graph, std::ostream& out) {
    for (const auto& node : graph.nodes()) {
        out << "NA" << '\t'
            << 0 << '\t'
            << node.length_bp << '\t'
            << "node_" << node.id << '\t'
            << node.copy_count << '\t'
            << '.' << '\n';
    }
    return out.good();
}

bool write_bed(const LosslessGraph& graph, const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    return write_bed(graph, ofs);
}

// =====================================================================
// PAF writer
// =====================================================================
//
// PAF-12 columns (minimap2 standard):
//   qname qlen qstart qend strand tname tlen tstart tend matches alnlen mapq
//
// Our OverlapPair carries (read_a, read_b, offset_a, offset_b,
// overlap_len, diff_count, strand). Map them directly:
//   qstart = offset_a, qend = offset_a + overlap_len
//   tstart = offset_b, tend = offset_b + overlap_len
//   matches = max(0, overlap_len - diff_count)
//   alnlen  = overlap_len
//   mapq    = min(60, int((matches / overlap_len) * 60))
//
// TODO(v0.3): propagate per-side strand from minimizer_sketcher so
// the '+'/'-' column is accurate at sub-read granularity.

bool write_paf(const std::vector<branch::backend::OverlapPair>& pairs,
               const std::vector<std::string>& sequences,
               const std::vector<std::string>& names,
               std::ostream& out) {
    if (sequences.size() != names.size()) return false;

    for (const auto& p : pairs) {
        if (p.read_a >= sequences.size() || p.read_b >= sequences.size()) {
            return false;
        }
        const auto qlen = static_cast<std::uint64_t>(sequences[p.read_a].size());
        const auto tlen = static_cast<std::uint64_t>(sequences[p.read_b].size());

        std::uint64_t qstart = p.offset_a;
        std::uint64_t tstart = p.offset_b;
        std::uint64_t aln    = p.overlap_len;

        // Clamp against read lengths in case the backend emits an
        // overlap that extends past the end of either sequence. The
        // real aligner enforces this invariant; we defend here so the
        // PAF stays consistent for downstream parsers.
        if (qstart > qlen) qstart = qlen;
        if (tstart > tlen) tstart = tlen;
        if (qstart + aln > qlen) aln = qlen - qstart;
        if (tstart + aln > tlen) aln = tlen - tstart;

        const std::uint64_t qend = qstart + aln;
        const std::uint64_t tend = tstart + aln;

        const std::int64_t diff = p.diff_count;
        std::int64_t matches = static_cast<std::int64_t>(aln) - diff;
        if (matches < 0) matches = 0;

        int mapq = 0;
        if (aln > 0) {
            const double identity = static_cast<double>(matches) /
                                    static_cast<double>(aln);
            mapq = static_cast<int>(identity * 60.0);
            if (mapq > 60) mapq = 60;
            if (mapq < 0) mapq = 0;
        }

        const char strand_c = (p.strand == 0) ? '+' : '-';

        out << names[p.read_a] << '\t'
            << qlen   << '\t'
            << qstart << '\t'
            << qend   << '\t'
            << strand_c << '\t'
            << names[p.read_b] << '\t'
            << tlen   << '\t'
            << tstart << '\t'
            << tend   << '\t'
            << matches << '\t'
            << aln    << '\t'
            << mapq   << '\n';
    }
    return out.good();
}

bool write_paf(const std::vector<branch::backend::OverlapPair>& pairs,
               const std::vector<std::string>& sequences,
               const std::vector<std::string>& names,
               const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs) return false;
    return write_paf(pairs, sequences, names, ofs);
}

}  // namespace branch::graph
