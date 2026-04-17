// BRANCH v0.4 — JSON report writer implementation.
//
// Hand-rolled JSON writer (no external dependencies).
// Produces branch-report.json with sample metadata and per-branch entries.

#include "project_report.hpp"

#include <cstdio>
#include <fstream>
#include <sstream>

namespace branch::project {

namespace {

// Escape a string for JSON output.
// Handles: backslash, double-quote, control characters.
std::string json_escape(const std::string& s) {
    std::ostringstream out;
    for (unsigned char c : s) {
        switch (c) {
            case '"':  out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (c < 0x20) {
                    // Control character: use \uXXXX
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", c);
                    out << buf;
                } else {
                    out << c;
                }
                break;
        }
    }
    return out.str();
}

// Write a JSON string value (quoted).
void write_json_string(std::ostream& os, const std::string& s) {
    os << '"' << json_escape(s) << '"';
}

// Write a LinearMapping as JSON object.
void write_mapping_json(std::ostream& os, const LinearMapping& m, int indent) {
    std::string pad(static_cast<std::size_t>(indent), ' ');
    std::string pad2(static_cast<std::size_t>(indent + 2), ' ');
    os << pad << "{\n";
    os << pad2 << "\"ref\": "; write_json_string(os, m.ref_name); os << ",\n";
    os << pad2 << "\"target\": "; write_json_string(os, m.target); os << ",\n";
    os << pad2 << "\"start\": " << m.target_start << ",\n";
    os << pad2 << "\"end\": " << m.target_end << ",\n";
    os << pad2 << "\"mapq\": " << m.mapq << ",\n";
    os << pad2 << "\"strand\": \"" << m.strand << "\",\n";
    os << pad2 << "\"query_len\": " << m.query_len << ",\n";
    os << pad2 << "\"query_range\": [" << m.query_start << ", " << m.query_end << "]\n";
    os << pad << "}";
}

// Write a BranchEntry as JSON object.
void write_branch_json(std::ostream& os, const BranchEntry& b, int indent) {
    std::string pad(static_cast<std::size_t>(indent), ' ');
    std::string pad2(static_cast<std::size_t>(indent + 2), ' ');

    os << pad << "{\n";
    os << pad2 << "\"branch_id\": "; write_json_string(os, b.branch_id); os << ",\n";
    os << pad2 << "\"length_bp\": " << b.length_bp << ",\n";

    // VAF: -1 means unknown, write as null
    os << pad2 << "\"vaf\": ";
    if (b.vaf < 0) {
        os << "null";
    } else {
        os << b.vaf;
    }
    os << ",\n";

    // Coverage: -1 means unknown, write as null
    os << pad2 << "\"coverage\": ";
    if (b.coverage < 0) {
        os << "null";
    } else {
        os << b.coverage;
    }
    os << ",\n";

    os << pad2 << "\"unannotated\": " << (b.unannotated ? "true" : "false") << ",\n";

    // Linear mappings array
    os << pad2 << "\"linear_mappings\": [";
    if (b.linear_mappings.empty()) {
        os << "]";
    } else {
        os << "\n";
        for (std::size_t i = 0; i < b.linear_mappings.size(); ++i) {
            write_mapping_json(os, b.linear_mappings[i], indent + 4);
            if (i + 1 < b.linear_mappings.size()) {
                os << ",";
            }
            os << "\n";
        }
        os << pad2 << "]";
    }
    os << "\n";
    os << pad << "}";
}

}  // namespace

bool write_branch_report_json(
    const std::string& out_path,
    const std::string& sample_name,
    const std::vector<BranchEntry>& branches,
    std::string* err_out) {

    std::ofstream os(out_path);
    if (!os) {
        if (err_out) {
            *err_out = "failed to open output file: " + out_path;
        }
        return false;
    }

    // Count unannotated
    std::size_t unannotated_count = 0;
    for (const auto& b : branches) {
        if (b.unannotated) {
            ++unannotated_count;
        }
    }

    os << "{\n";
    os << "  \"sample\": "; write_json_string(os, sample_name); os << ",\n";
    os << "  \"version\": \"0.4.1\",\n";
    os << "  \"branch_count\": " << branches.size() << ",\n";
    os << "  \"unannotated_count\": " << unannotated_count << ",\n";
    os << "  \"branches\": [";

    if (branches.empty()) {
        os << "]\n";
    } else {
        os << "\n";
        for (std::size_t i = 0; i < branches.size(); ++i) {
            write_branch_json(os, branches[i], 4);
            if (i + 1 < branches.size()) {
                os << ",";
            }
            os << "\n";
        }
        os << "  ]\n";
    }

    os << "}\n";

    if (!os) {
        if (err_out) {
            *err_out = "write error on output file: " + out_path;
        }
        return false;
    }

    return true;
}

}  // namespace branch::project
