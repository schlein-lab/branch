// BRANCH v0.1 — GFA-1.2 I/O implementation.

#include "graph/graph_io.hpp"

#include <charconv>
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
            << "\tCV:f:" << std::fixed << std::setprecision(3) << node.copy_count_confidence
            << '\n';
    }

    for (std::size_t i = 0; i < graph.edges().size(); ++i) {
        const auto& edge = graph.edges()[i];
        // GFA L-line: L <from> <from_orient> <to> <to_orient> <cigar>
        out << "L\t" << edge.from << "\t+\t" << edge.to << "\t+\t0M"
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
            std::string token;
            while (iss >> token) {
                parse_tag_int(token, "LN:i:", len);
                parse_tag_int(token, "RC:i:", rc);
                parse_tag_int(token, "CN:i:", cn);
                parse_tag_float(token, "CV:f:", cv);
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

}  // namespace branch::graph
