#include "gfa_writer.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>

namespace branch::wg::phaser {

namespace {

const char* kind_str(NodeKind k) {
    switch (k) {
        case NodeKind::SHARED:       return "shared";
        case NodeKind::BUBBLE_H1:    return "bubble_h1";
        case NodeKind::BUBBLE_H2:    return "bubble_h2";
        case NodeKind::BRANCH:       return "branch";
        case NodeKind::UNCLASSIFIED: return "unc";
    }
    return "?";
}

std::string node_name(const UnitigNode& u) {
    std::ostringstream s;
    s << kind_str(u.kind) << "_" << u.id;
    return s.str();
}

// Forward-direction walk from a starting node, following next_nodes
// edges until we hit a fork or a node already on the path. Visits each
// node at most once. The walker stops at any node that doesn't match the
// allowed kind set — that's how h1_walk skips BUBBLE_H2 nodes and vice
// versa.
std::vector<NodeId> walk_path(
    const std::vector<UnitigNode>& unitigs,
    NodeId start,
    bool allow_h1, bool allow_h2,
    bool allow_branch)
{
    std::vector<NodeId> path;
    std::vector<bool> seen(unitigs.size(), false);
    NodeId cur = start;
    while (cur != ~0u && cur < unitigs.size() && !seen[cur]) {
        const auto& u = unitigs[cur];
        // Filter by allowed kinds.
        bool ok = false;
        switch (u.kind) {
            case NodeKind::SHARED:    ok = true; break;
            case NodeKind::BUBBLE_H1: ok = allow_h1; break;
            case NodeKind::BUBBLE_H2: ok = allow_h2; break;
            case NodeKind::BRANCH:    ok = allow_branch; break;
            default:                  ok = true; break;
        }
        if (!ok) break;
        path.push_back(cur);
        seen[cur] = true;
        // Choose next node: prefer SHARED; else preferred kind; else stop.
        NodeId next = ~0u;
        for (NodeId n : u.next_nodes) {
            if (n >= unitigs.size() || seen[n]) continue;
            const auto& nu = unitigs[n];
            if (nu.kind == NodeKind::SHARED) { next = n; break; }
        }
        if (next == ~0u) {
            for (NodeId n : u.next_nodes) {
                if (n >= unitigs.size() || seen[n]) continue;
                const auto& nu = unitigs[n];
                if ((allow_h1 && nu.kind == NodeKind::BUBBLE_H1) ||
                    (allow_h2 && nu.kind == NodeKind::BUBBLE_H2)) {
                    next = n; break;
                }
            }
        }
        cur = next;
    }
    return path;
}

void write_fasta_from_path(const std::string& path_out,
                           const std::string& contig_name,
                           const std::vector<NodeId>& path,
                           const std::vector<UnitigNode>& unitigs)
{
    std::ofstream f(path_out);
    if (!f) {
        std::cerr << "[phaser/gfa] cannot open " << path_out << "\n";
        return;
    }
    f << '>' << contig_name << "\n";
    std::string buf;
    buf.reserve(64'000);
    for (NodeId nid : path) {
        if (nid >= unitigs.size()) continue;
        for (char c : unitigs[nid].seq) {
            buf.push_back(c);
            if (buf.size() >= 60) { f << buf << "\n"; buf.clear(); }
        }
    }
    if (!buf.empty()) f << buf << "\n";
}

}  // namespace

void GfaWriter::write(
    const std::vector<UnitigNode>& unitigs,
    const std::vector<Bubble>& bubbles,
    const std::vector<ReadAssignment>& assignments,
    const std::vector<std::string>& read_id_strings,
    const GfaWriterOpts& opts)
{
    (void)bubbles;

    // ----- GFA -----
    if (!opts.out_gfa_path.empty()) {
        std::ofstream g(opts.out_gfa_path);
        if (!g) {
            std::cerr << "[phaser/gfa] cannot open " << opts.out_gfa_path << "\n";
        } else {
            g << "H\tVN:Z:1.0\n";
            for (const auto& u : unitigs) {
                if (u.kind == NodeKind::UNCLASSIFIED) continue;
                g << "S\t" << node_name(u) << '\t' << u.seq
                  << "\tLN:i:" << u.seq.size()
                  << "\tKD:Z:" << kind_str(u.kind) << "\n";
            }
            for (const auto& u : unitigs) {
                for (NodeId n : u.next_nodes) {
                    if (n >= unitigs.size()) continue;
                    if (unitigs[n].kind == NodeKind::UNCLASSIFIED) continue;
                    g << "L\t" << node_name(u) << "\t+\t"
                      << node_name(unitigs[n]) << "\t+\t0M\n";
                }
            }
            // Find a starting SHARED or BUBBLE_H1 node with no incoming
            // SHARED predecessor → walk start.
            auto find_start = [&](bool h1) {
                for (const auto& u : unitigs) {
                    if (u.kind == NodeKind::UNCLASSIFIED) continue;
                    if (u.kind == NodeKind::BRANCH) continue;
                    if (h1 && u.kind == NodeKind::BUBBLE_H2) continue;
                    if (!h1 && u.kind == NodeKind::BUBBLE_H1) continue;
                    bool has_main_pred = false;
                    for (NodeId p : u.prev_nodes) {
                        if (p >= unitigs.size()) continue;
                        const auto& pu = unitigs[p];
                        if (pu.kind == NodeKind::SHARED) { has_main_pred = true; break; }
                        if (h1 && pu.kind == NodeKind::BUBBLE_H1) { has_main_pred = true; break; }
                        if (!h1 && pu.kind == NodeKind::BUBBLE_H2) { has_main_pred = true; break; }
                    }
                    if (!has_main_pred) return u.id;
                }
                return ~0u;
            };
            auto write_path = [&](const char* name, NodeId start, bool h1) {
                if (start == ~0u) return;
                auto path = walk_path(unitigs, start, h1, !h1, false);
                if (path.empty()) return;
                g << "P\t" << name << '\t';
                for (std::size_t i = 0; i < path.size(); ++i) {
                    if (i) g << ',';
                    g << node_name(unitigs[path[i]]) << '+';
                }
                g << '\t';
                for (std::size_t i = 0; i + 1 < path.size(); ++i) {
                    if (i) g << ',';
                    g << "0M";
                }
                g << "\n";
            };
            write_path("hap1", find_start(true),  true);
            write_path("hap2", find_start(false), false);
            std::cerr << "[phaser/gfa] wrote " << opts.out_gfa_path << "\n";
        }
    }

    // ----- FASTA outputs (optional, for callers without gfatools) -----
    auto find_start_for_fa = [&](bool h1) {
        for (const auto& u : unitigs) {
            if (u.kind == NodeKind::UNCLASSIFIED) continue;
            if (u.kind == NodeKind::BRANCH) continue;
            if (h1 && u.kind == NodeKind::BUBBLE_H2) continue;
            if (!h1 && u.kind == NodeKind::BUBBLE_H1) continue;
            bool has_main_pred = false;
            for (NodeId p : u.prev_nodes) {
                if (p >= unitigs.size()) continue;
                const auto& pu = unitigs[p];
                if (pu.kind == NodeKind::SHARED) { has_main_pred = true; break; }
                if (h1 && pu.kind == NodeKind::BUBBLE_H1) { has_main_pred = true; break; }
                if (!h1 && pu.kind == NodeKind::BUBBLE_H2) { has_main_pred = true; break; }
            }
            if (!has_main_pred) return u.id;
        }
        return ~0u;
    };
    if (!opts.out_h1_fa_path.empty()) {
        NodeId s = find_start_for_fa(true);
        if (s != ~0u) {
            auto p = walk_path(unitigs, s, true, false, false);
            write_fasta_from_path(opts.out_h1_fa_path, "hap1", p, unitigs);
            std::cerr << "[phaser/gfa] wrote " << opts.out_h1_fa_path
                      << " (" << p.size() << " unitigs)\n";
        }
    }
    if (!opts.out_h2_fa_path.empty()) {
        NodeId s = find_start_for_fa(false);
        if (s != ~0u) {
            auto p = walk_path(unitigs, s, false, true, false);
            write_fasta_from_path(opts.out_h2_fa_path, "hap2", p, unitigs);
            std::cerr << "[phaser/gfa] wrote " << opts.out_h2_fa_path
                      << " (" << p.size() << " unitigs)\n";
        }
    }

    // ----- Read assignment TSV -----
    if (!opts.out_assignments_tsv_path.empty()) {
        std::ofstream t(opts.out_assignments_tsv_path);
        if (!t) {
            std::cerr << "[phaser/gfa] cannot open "
                      << opts.out_assignments_tsv_path << "\n";
        } else {
            t << "read_id\ttag\tconfidence\thome_node\tbranch_idx\n";
            for (const auto& a : assignments) {
                const std::string& nm =
                    (a.read_id < read_id_strings.size())
                        ? read_id_strings[a.read_id]
                        : std::to_string(a.read_id);
                t << nm << '\t' << read_tag_str(a.tag)
                  << '\t' << a.confidence
                  << '\t' << a.home_node
                  << '\t' << a.branch_idx << "\n";
            }
            std::cerr << "[phaser/gfa] wrote "
                      << opts.out_assignments_tsv_path << "\n";
        }
    }
}

}  // namespace branch::wg::phaser
