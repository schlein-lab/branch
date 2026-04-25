// BRANCH — ReadDB implementation.

#include "graph/read_db.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "graph/lossless_graph.hpp"

namespace branch::graph {

ReadId ReadDB::add(std::string name,
                   std::uint32_t length_bp,
                   ReadBucket bucket,
                   std::string original_sequence,
                   std::vector<PathStep> path,
                   std::vector<DeltaOp> deltas) {
    const ReadId id = static_cast<ReadId>(entries_.size());
    ReadEntry e;
    e.id = id;
    e.name = std::move(name);
    e.length_bp = length_bp;
    e.bucket = bucket;
    e.path = std::move(path);
    e.deltas = std::move(deltas);
    e.original_sequence = std::move(original_sequence);
    entries_.push_back(std::move(e));
    return id;
}

ReadDB::BucketCounts ReadDB::bucket_counts() const noexcept {
    BucketCounts c{};
    for (const auto& e : entries_) {
        switch (e.bucket) {
            case ReadBucket::Mapped:                 ++c.mapped; break;
            case ReadBucket::UnmappedNoOverlap:      ++c.unmapped_no_overlap; break;
            case ReadBucket::UnmappedLowComplexity:  ++c.unmapped_low_complexity; break;
            case ReadBucket::Filtered:               ++c.filtered; break;
        }
    }
    return c;
}

void ReadDB::apply_remap(const GraphRemap& remap) {
    const std::size_t N = remap.remap.size();
    auto resolve = [&](NodeId old) -> std::pair<NodeId, bool> {
        // Returns (new_id, valid). valid=false when the node was
        // dropped without a coverer; caller drops the path step.
        if (old >= N) return {old, false};
        const NodeId mapped = remap.remap[old];
        if (mapped != GraphRemap::kDropped) return {mapped, true};
        if (old < remap.coverer.size()
            && remap.coverer[old] != GraphRemap::kNoCoverer) {
            const NodeId via = remap.coverer[old];
            if (via < N) {
                const NodeId via_new = remap.remap[via];
                if (via_new != GraphRemap::kDropped) return {via_new, true};
            }
        }
        return {old, false};
    };

    for (auto& e : entries_) {
        if (e.bucket != ReadBucket::Mapped) continue;
        std::vector<PathStep> remapped;
        remapped.reserve(e.path.size());
        for (const auto& step : e.path) {
            auto [nu, ok] = resolve(step.unitig);
            if (!ok) continue;
            // Coalesce consecutive duplicates that arise after
            // multiple path steps land on the same coverer.
            if (!remapped.empty() && remapped.back().unitig == nu &&
                remapped.back().reverse == step.reverse) {
                remapped.back().end = std::max(remapped.back().end, step.end);
                continue;
            }
            PathStep s = step;
            s.unitig = nu;
            remapped.push_back(s);
        }
        e.path = std::move(remapped);
        if (e.path.empty()) {
            // The graph mutator removed every node this read traversed
            // (rare — usually means the read was wholly contained AND
            // its coverer was also subsequently dropped without a
            // grandcoverer). Demote to unmapped rather than delete:
            // the read is still in the DB and still recoverable from
            // its original_sequence.
            e.bucket = ReadBucket::UnmappedNoOverlap;
        }
    }
}

void ReadDB::release_mapped_sequences() noexcept {
    for (auto& e : entries_) {
        if (e.bucket == ReadBucket::Mapped && !e.deltas.empty()) {
            e.original_sequence.clear();
            e.original_sequence.shrink_to_fit();
        }
    }
}

namespace {

// GAF cs:Z encoding for a list of DeltaOp entries. Format follows
// minimap2's short-form convention: ":N" for N matches, "*xy" for a
// substitution from x to y, "+xx" for insertion, "-xx" for deletion.
std::string encode_cs_tag(std::span<const DeltaOp> deltas,
                          std::uint32_t total_len) {
    std::ostringstream os;
    os << "cs:Z:";
    std::uint32_t pos = 0;
    for (const auto& d : deltas) {
        if (d.pos > pos) os << ':' << (d.pos - pos);
        switch (d.op) {
            case 'S': os << '*' << char(std::tolower(d.base)) << char(std::tolower(d.base)); break;
            case 'I': os << '+' << char(std::tolower(d.base)); break;
            case 'D': os << '-' << char(std::tolower(d.base)); break;
            default:  break;
        }
        pos = d.pos + (d.op == 'D' ? 0 : 1);
    }
    if (total_len > pos) os << ':' << (total_len - pos);
    return os.str();
}

}  // namespace

bool ReadDB::write_gaf(const std::string& gaf_path,
                       const LosslessGraph& graph) const {
    std::ofstream gaf(gaf_path);
    if (!gaf) return false;

    // Companion FASTQ for unmapped + filtered reads.
    const std::string fastq_path = gaf_path + ".unmapped.fastq";
    std::ofstream fq;
    bool fq_open = false;

    for (const auto& e : entries_) {
        if (e.bucket == ReadBucket::Mapped) {
            // Compute path total length from unitig spans.
            std::uint32_t path_total = 0;
            for (const auto& s : e.path) path_total += (s.end - s.start);
            // GAF columns:
            // 1=name 2=length 3=qstart 4=qend 5=strand
            // 6=path 7=path_len 8=path_start 9=path_end
            // 10=matches 11=block 12=mapq 13+=tags
            std::ostringstream path_str;
            for (const auto& s : e.path) {
                path_str << (s.reverse ? '<' : '>') << s.unitig;
                (void)graph;  // currently unused, reserved for future tag enrichment
            }
            const std::uint32_t matches = path_total > static_cast<std::uint32_t>(e.deltas.size())
                ? path_total - static_cast<std::uint32_t>(e.deltas.size())
                : 0;
            gaf << e.name << '\t'
                << e.length_bp << '\t'
                << 0 << '\t' << e.length_bp << '\t'
                << '+' << '\t'
                << path_str.str() << '\t'
                << path_total << '\t' << 0 << '\t' << path_total << '\t'
                << matches << '\t' << path_total << '\t' << 60 << '\t'
                << encode_cs_tag(e.deltas, e.length_bp)
                << '\n';
        } else {
            if (!fq_open) {
                fq.open(fastq_path);
                if (!fq) return false;
                fq_open = true;
            }
            const char bucket_tag = e.bucket == ReadBucket::UnmappedNoOverlap ? 'U'
                                  : e.bucket == ReadBucket::UnmappedLowComplexity ? 'L'
                                  : 'F';
            fq << '@' << e.name << " bucket=" << bucket_tag
               << " length=" << e.length_bp << '\n'
               << e.original_sequence << '\n'
               << '+' << '\n'
               << std::string(e.original_sequence.size(), '!') << '\n';
        }
    }
    return static_cast<bool>(gaf);
}

bool ReadDB::read_gaf(const std::string& gaf_path, const LosslessGraph& graph) {
    (void)graph;
    std::ifstream gaf(gaf_path);
    if (!gaf) return false;
    entries_.clear();

    std::string line;
    while (std::getline(gaf, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string name, qlen_s, qstart_s, qend_s, strand, path_s,
                    path_len_s, ps, pe, matches_s, block_s, mapq_s;
        if (!std::getline(iss, name, '\t')) continue;
        if (!std::getline(iss, qlen_s, '\t')) continue;
        if (!std::getline(iss, qstart_s, '\t')) continue;
        if (!std::getline(iss, qend_s, '\t')) continue;
        if (!std::getline(iss, strand, '\t')) continue;
        if (!std::getline(iss, path_s, '\t')) continue;
        if (!std::getline(iss, path_len_s, '\t')) continue;
        // Remaining fields not required for round-trip of the path.

        std::vector<PathStep> path;
        std::size_t i = 0;
        while (i < path_s.size()) {
            char dir = path_s[i++];
            const bool reverse = (dir == '<');
            std::string num;
            while (i < path_s.size() && std::isdigit(static_cast<unsigned char>(path_s[i]))) {
                num.push_back(path_s[i++]);
            }
            if (num.empty()) break;
            PathStep st;
            st.unitig = static_cast<NodeId>(std::stoul(num));
            st.start = 0;
            st.end = 0;  // span left at 0; caller computes from graph if needed
            st.reverse = reverse;
            path.push_back(st);
        }

        ReadEntry e;
        e.id = static_cast<ReadId>(entries_.size());
        e.name = name;
        e.length_bp = static_cast<std::uint32_t>(std::stoul(qlen_s));
        e.bucket = ReadBucket::Mapped;
        e.path = std::move(path);
        entries_.push_back(std::move(e));
    }
    return true;
}

}  // namespace branch::graph
