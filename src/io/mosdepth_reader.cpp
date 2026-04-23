#include "io/mosdepth_reader.hpp"

#include <charconv>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

#include <zlib.h>

namespace branch::io {

namespace {

bool parse_u32(std::string_view sv, std::uint32_t& out) {
    auto r = std::from_chars(sv.data(), sv.data() + sv.size(), out);
    return r.ec == std::errc{};
}

bool ends_with_gz(const std::string& path) {
    if (path.size() < 3) return false;
    return path.compare(path.size() - 3, 3, ".gz") == 0;
}

// Read all lines from a gzip-compressed file. Returns empty vector on failure.
std::vector<std::string> read_all_lines_gz(const std::string& path) {
    std::vector<std::string> lines;
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz) return lines;
    constexpr std::size_t BUF_SIZE = 1 << 16;
    std::string buffer;
    buffer.resize(BUF_SIZE);
    std::string current;
    int n = 0;
    while ((n = gzread(gz, buffer.data(), static_cast<unsigned>(BUF_SIZE))) > 0) {
        for (int i = 0; i < n; ++i) {
            char c = buffer[static_cast<std::size_t>(i)];
            if (c == '\n') {
                lines.push_back(std::move(current));
                current.clear();
            } else if (c != '\r') {
                current.push_back(c);
            }
        }
    }
    if (!current.empty()) lines.push_back(std::move(current));
    gzclose(gz);
    return lines;
}

// Read all lines from a plain-text file.
std::vector<std::string> read_all_lines_plain(const std::string& path) {
    std::vector<std::string> lines;
    std::ifstream ifs(path);
    if (!ifs) return lines;
    std::string line;
    while (std::getline(ifs, line)) {
        // Strip trailing \r (Windows line endings)
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(std::move(line));
    }
    return lines;
}

bool parse_region_line(const std::string& line, RegionCoverage& rc) {
    if (line.empty() || line[0] == '#') return false;
    std::vector<std::string> fields;
    fields.reserve(5);
    std::string tok;
    std::istringstream iss(line);
    while (std::getline(iss, tok, '\t')) fields.push_back(tok);

    if (fields.size() < 4) return false;
    rc.chrom = fields[0];
    if (!parse_u32(fields[1], rc.start)) return false;
    if (!parse_u32(fields[2], rc.end)) return false;

    if (fields.size() >= 5) {
        rc.name = fields[3];
        rc.mean_coverage = std::strtof(fields[4].c_str(), nullptr);
    } else {
        rc.name.clear();
        rc.mean_coverage = std::strtof(fields[3].c_str(), nullptr);
    }
    return true;
}

}  // namespace

std::vector<RegionCoverage> read_mosdepth_regions(const std::string& path) {
    std::vector<RegionCoverage> out;
    const auto lines = ends_with_gz(path) ? read_all_lines_gz(path)
                                          : read_all_lines_plain(path);
    out.reserve(lines.size());
    for (const auto& line : lines) {
        RegionCoverage rc;
        if (parse_region_line(line, rc)) out.push_back(std::move(rc));
    }
    return out;
}

}  // namespace branch::io
