#include "io/mosdepth_reader.hpp"

#include <charconv>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>

namespace branch::io {

namespace {

bool parse_u32(std::string_view sv, std::uint32_t& out) {
    auto r = std::from_chars(sv.data(), sv.data() + sv.size(), out);
    return r.ec == std::errc{};
}

}  // namespace

std::vector<RegionCoverage> read_mosdepth_regions(const std::string& path) {
    std::vector<RegionCoverage> out;
    std::ifstream ifs(path);
    if (!ifs) return out;

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        RegionCoverage rc;
        // Up to 5 tab-separated fields: chrom start end name mean_cov
        // If only 4 fields (no name), mean_cov shifts to field 4.
        std::vector<std::string> fields;
        fields.reserve(5);
        std::string tok;
        std::istringstream iss(line);
        while (std::getline(iss, tok, '\t')) fields.push_back(tok);

        if (fields.size() < 4) continue;
        rc.chrom = fields[0];
        if (!parse_u32(fields[1], rc.start)) continue;
        if (!parse_u32(fields[2], rc.end)) continue;

        if (fields.size() >= 5) {
            rc.name = fields[3];
            rc.mean_coverage = std::strtof(fields[4].c_str(), nullptr);
        } else {
            rc.name.clear();
            rc.mean_coverage = std::strtof(fields[3].c_str(), nullptr);
        }
        out.push_back(std::move(rc));
    }
    return out;
}

}  // namespace branch::io
