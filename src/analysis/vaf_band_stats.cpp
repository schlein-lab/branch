// BRANCH P1.3 — VAF-band benchmark statistics (implementation).

#include "analysis/vaf_band_stats.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

namespace branch::analysis {

namespace {

// Split a line by `delim` into fields, no quote handling (TSV/GFA are
// tab-separated without embedded tabs).
std::vector<std::string> split(std::string_view line, char delim) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == delim) {
            out.push_back(std::move(cur));
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(std::move(cur));
    return out;
}

// Parse a "TAG:TYPE:VALUE" GFA optional field. Returns the VALUE part
// as a string_view, or empty if the tag doesn't match.
bool parse_tag_uint(std::string_view token, std::string_view tag,
                    std::uint32_t& out) {
    if (token.size() > tag.size() &&
        token.substr(0, tag.size()) == tag) {
        // strtoull tolerates trailing non-digits; fine for our tags.
        std::string val(token.substr(tag.size()));
        char* endp = nullptr;
        const auto v = std::strtoull(val.c_str(), &endp, 10);
        if (endp == val.c_str()) return false;
        out = static_cast<std::uint32_t>(v);
        return true;
    }
    return false;
}

double median_of(std::vector<double> v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2 == 1) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

}  // namespace

// ---------------------------------------------------------------------
// assign_band
// ---------------------------------------------------------------------

std::size_t assign_band(double true_vaf) noexcept {
    if (!(true_vaf > 0.0) || true_vaf > 1.0) return kVafBands.size();
    for (std::size_t i = 0; i < kVafBands.size(); ++i) {
        const auto& b = kVafBands[i];
        if (true_vaf > b.lo && true_vaf <= b.hi) return i;
    }
    return kVafBands.size();
}

// ---------------------------------------------------------------------
// BandStats derived metrics
// ---------------------------------------------------------------------

double BandStats::recall() const noexcept {
    if (n_truth == 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(n_true_positive) / static_cast<double>(n_truth);
}

double BandStats::precision() const noexcept {
    const std::size_t denom = n_true_positive + n_false_positive;
    if (denom == 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(n_true_positive) / static_cast<double>(denom);
}

double BandStats::median_abs_vaf_residual() const noexcept {
    return median_of(vaf_abs_residuals);
}

double BandStats::copy_count_accuracy() const noexcept {
    if (n_true_positive == 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(n_copy_count_correct) /
           static_cast<double>(n_true_positive);
}

// ---------------------------------------------------------------------
// compute_band_stats
// ---------------------------------------------------------------------

BandBenchmarkReport compute_band_stats(
    const std::vector<TruthAllele>& truth,
    const std::vector<CalledBubble>& calls,
    const std::string& sample_name) {

    BandBenchmarkReport report;
    report.sample_name = sample_name;
    report.n_truth_total = truth.size();
    report.n_called_total = calls.size();

    for (std::size_t i = 0; i < kVafBands.size(); ++i) {
        report.bands[i].band_name = kVafBands[i].name;
        report.bands[i].lo = kVafBands[i].lo;
        report.bands[i].hi = kVafBands[i].hi;
    }

    // Band assignment for each truth allele.
    std::vector<std::size_t> truth_band(truth.size(), kVafBands.size());
    for (std::size_t i = 0; i < truth.size(); ++i) {
        truth_band[i] = assign_band(truth[i].vaf_target);
        if (truth_band[i] < kVafBands.size()) {
            ++report.bands[truth_band[i]].n_truth;
        }
    }

    // One-to-one greedy matching by |est_vaf - truth.vaf_target|.
    // Scan truth in order; each truth claims the nearest unclaimed
    // call within kVafMatchTolerance.
    std::vector<char> call_claimed(calls.size(), 0);

    for (std::size_t ti = 0; ti < truth.size(); ++ti) {
        if (truth_band[ti] >= kVafBands.size()) continue;  // OOB truth VAF
        double best_dist = kVafMatchTolerance + 1e-12;
        std::size_t best_ci = static_cast<std::size_t>(-1);
        for (std::size_t ci = 0; ci < calls.size(); ++ci) {
            if (call_claimed[ci]) continue;
            const double d = std::fabs(calls[ci].est_vaf - truth[ti].vaf_target);
            if (d < best_dist) {
                best_dist = d;
                best_ci = ci;
            }
        }
        if (best_ci == static_cast<std::size_t>(-1)) continue;  // missed
        auto& band = report.bands[truth_band[ti]];
        call_claimed[best_ci] = 1;
        ++band.n_true_positive;
        ++band.n_called;
        band.vaf_abs_residuals.push_back(best_dist);
        if (calls[best_ci].est_copy_count == truth[ti].copy_count) {
            ++band.n_copy_count_correct;
        }
    }

    // Any call left unclaimed is a false positive attributed to the
    // band its est_vaf falls into.
    for (std::size_t ci = 0; ci < calls.size(); ++ci) {
        if (call_claimed[ci]) continue;
        const std::size_t b = assign_band(calls[ci].est_vaf);
        if (b < kVafBands.size()) {
            ++report.bands[b].n_false_positive;
        }
        // Calls with out-of-range est_vaf (0.0 or > 1.0) are dropped
        // from both TP and FP counts — they indicate degenerate nodes
        // (e.g. read_support=0) and shouldn't affect any band.
    }

    return report;
}

// ---------------------------------------------------------------------
// parse_fixture_meta_tsv
// ---------------------------------------------------------------------

std::vector<TruthAllele> parse_fixture_meta_tsv(const std::string& path,
                                                std::string* err) {
    std::vector<TruthAllele> out;
    std::ifstream ifs(path);
    if (!ifs) {
        if (err) *err = "cannot open meta TSV: " + path;
        return out;
    }
    std::string line;
    std::size_t lineno = 0;
    while (std::getline(ifs, line)) {
        ++lineno;
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        auto fields = split(line, '\t');
        if (fields.size() < 5) {
            if (err) {
                *err = "meta TSV line " + std::to_string(lineno) +
                       " has <5 fields";
            }
            return {};
        }
        TruthAllele t;
        t.allele_id = fields[0];
        t.copy_count = std::atoi(fields[1].c_str());
        t.vaf_target = std::atof(fields[2].c_str());
        t.n_reads = static_cast<std::uint64_t>(
            std::strtoull(fields[3].c_str(), nullptr, 10));
        t.cassette_len_bp = static_cast<std::uint64_t>(
            std::strtoull(fields[4].c_str(), nullptr, 10));
        out.push_back(std::move(t));
    }
    return out;
}

// ---------------------------------------------------------------------
// parse_gfa_nodes
// ---------------------------------------------------------------------

std::vector<CalledBubble> parse_gfa_nodes(const std::string& path,
                                          std::string* err) {
    std::vector<CalledBubble> out;
    std::ifstream ifs(path);
    if (!ifs) {
        if (err) *err = "cannot open GFA: " + path;
        return out;
    }
    std::string line;
    std::uint64_t rc_sum = 0;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] != 'S') continue;
        auto fields = split(line, '\t');
        // S-line: "S", id, seq_or_*, LN:i:..., RC:i:..., CN:i:..., CV:f:...
        if (fields.size() < 3) continue;
        CalledBubble b;
        b.node_id = static_cast<std::uint32_t>(
            std::strtoull(fields[1].c_str(), nullptr, 10));
        std::uint32_t len = 0, rc = 0, cn = 1;
        for (std::size_t i = 3; i < fields.size(); ++i) {
            std::string_view tok(fields[i]);
            parse_tag_uint(tok, "LN:i:", len);
            parse_tag_uint(tok, "RC:i:", rc);
            parse_tag_uint(tok, "CN:i:", cn);
        }
        b.length_bp = len;
        b.read_support = rc;
        b.est_copy_count = static_cast<int>(cn);
        rc_sum += rc;
        out.push_back(b);
    }
    if (rc_sum > 0) {
        for (auto& b : out) {
            b.est_vaf = static_cast<double>(b.read_support) /
                        static_cast<double>(rc_sum);
        }
    }
    return out;
}

// ---------------------------------------------------------------------
// aggregate_nodes_by_length
// ---------------------------------------------------------------------

std::vector<CalledBubble> aggregate_nodes_by_length(
    const std::vector<CalledBubble>& nodes,
    std::uint32_t bucket_size_bp) {

    if (nodes.empty() || bucket_size_bp == 0) return {};

    struct Bucket {
        std::uint32_t first_node_id = 0;
        std::uint32_t bucket_idx = 0;
        std::uint64_t sum_len = 0;
        std::uint64_t sum_rc  = 0;
        std::size_t count = 0;
        // Tally est_copy_count mode.
        std::vector<int> cc_tally;  // small vector; duplicates ok for mode
    };

    std::vector<Bucket> buckets;
    auto find_or_make = [&](std::uint32_t idx,
                            std::uint32_t node_id) -> Bucket& {
        for (auto& b : buckets) {
            if (b.bucket_idx == idx) return b;
        }
        buckets.push_back({node_id, idx, 0, 0, 0, {}});
        return buckets.back();
    };

    for (const auto& n : nodes) {
        const std::uint32_t idx =
            (n.length_bp + bucket_size_bp / 2) / bucket_size_bp;
        auto& b = find_or_make(idx, n.node_id);
        if (n.node_id < b.first_node_id) b.first_node_id = n.node_id;
        b.sum_len += n.length_bp;
        b.sum_rc  += n.read_support;
        ++b.count;
        b.cc_tally.push_back(n.est_copy_count);
    }

    std::vector<CalledBubble> out;
    out.reserve(buckets.size());
    for (const auto& b : buckets) {
        // Mode of cc_tally (small set; O(n^2) is fine for a handful).
        int mode = 1;
        std::size_t best_ct = 0;
        for (std::size_t i = 0; i < b.cc_tally.size(); ++i) {
            std::size_t ct = 0;
            for (std::size_t j = 0; j < b.cc_tally.size(); ++j) {
                if (b.cc_tally[j] == b.cc_tally[i]) ++ct;
            }
            if (ct > best_ct) { best_ct = ct; mode = b.cc_tally[i]; }
        }
        CalledBubble cb;
        cb.node_id = b.first_node_id;
        cb.est_copy_count = mode;
        cb.length_bp = static_cast<std::uint32_t>(
            b.sum_len / (b.count ? b.count : 1));
        cb.read_support = static_cast<std::uint32_t>(b.sum_rc);
        cb.est_vaf = static_cast<double>(b.count) /
                     static_cast<double>(nodes.size());
        out.push_back(cb);
    }
    return out;
}

// ---------------------------------------------------------------------
// write_report_json
// ---------------------------------------------------------------------

namespace {

void write_number_or_null(std::ostream& os, double v) {
    if (std::isnan(v) || std::isinf(v)) {
        os << "null";
    } else {
        // Fixed 6-digit precision keeps the output diffable.
        os << std::fixed << std::setprecision(6) << v;
        os.unsetf(std::ios::fixed);
    }
}

void write_json_string(std::ostream& os, const std::string& s) {
    os << '"';
    for (char c : s) {
        switch (c) {
            case '"':  os << "\\\""; break;
            case '\\': os << "\\\\"; break;
            case '\n': os << "\\n"; break;
            case '\r': os << "\\r"; break;
            case '\t': os << "\\t"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x",
                                  static_cast<unsigned char>(c));
                    os << buf;
                } else {
                    os << c;
                }
        }
    }
    os << '"';
}

}  // namespace

bool write_report_json(const BandBenchmarkReport& report,
                       const std::string& path,
                       std::string* err) {
    std::ofstream os(path);
    if (!os) {
        if (err) *err = "cannot open " + path + " for writing";
        return false;
    }
    os << "{\n";
    os << "  \"schema_version\": ";
    write_json_string(os, report.schema_version);
    os << ",\n";
    os << "  \"sample\": ";
    write_json_string(os, report.sample_name);
    os << ",\n";
    os << "  \"n_truth\": " << report.n_truth_total << ",\n";
    os << "  \"n_called\": " << report.n_called_total << ",\n";
    os << "  \"bands\": [\n";
    for (std::size_t i = 0; i < report.bands.size(); ++i) {
        const auto& b = report.bands[i];
        os << "    {\n";
        os << "      \"name\": "; write_json_string(os, b.band_name);
        os << ",\n";
        os << "      \"lo\": " << b.lo << ",\n";
        os << "      \"hi\": " << b.hi << ",\n";
        os << "      \"n_truth\": " << b.n_truth << ",\n";
        os << "      \"n_called\": " << b.n_called << ",\n";
        os << "      \"tp\": " << b.n_true_positive << ",\n";
        os << "      \"fp\": " << b.n_false_positive << ",\n";
        os << "      \"recall\": ";
        write_number_or_null(os, b.recall());
        os << ",\n";
        os << "      \"precision\": ";
        write_number_or_null(os, b.precision());
        os << ",\n";
        os << "      \"median_abs_vaf_residual\": ";
        write_number_or_null(os, b.median_abs_vaf_residual());
        os << ",\n";
        os << "      \"copy_count_accuracy\": ";
        write_number_or_null(os, b.copy_count_accuracy());
        os << "\n";
        os << "    }";
        if (i + 1 < report.bands.size()) os << ",";
        os << "\n";
    }
    os << "  ]\n";
    os << "}\n";
    if (!os) {
        if (err) *err = "write error on " + path;
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------
// write_report_table
// ---------------------------------------------------------------------

void write_report_table(const BandBenchmarkReport& report, std::ostream& os) {
    auto fmt = [](double v) -> std::string {
        if (std::isnan(v)) return std::string("   NaN");
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%6.3f", v);
        return std::string(buf);
    };
    os << "VAF-band benchmark — sample=" << report.sample_name
       << "  n_truth=" << report.n_truth_total
       << "  n_called=" << report.n_called_total << "\n";
    os << "schema=" << report.schema_version << "\n";
    os << "------------------------------------------------------------------"
          "----------\n";
    os << "band        n_truth n_call   tp   fp  recall   prec   |dVAF|   "
          "CN_acc\n";
    os << "------------------------------------------------------------------"
          "----------\n";
    for (const auto& b : report.bands) {
        char line[256];
        std::snprintf(line, sizeof(line),
                      "%-10s %7zu %6zu %4zu %4zu %s %s %s %s\n",
                      b.band_name.c_str(),
                      b.n_truth, b.n_called,
                      b.n_true_positive, b.n_false_positive,
                      fmt(b.recall()).c_str(),
                      fmt(b.precision()).c_str(),
                      fmt(b.median_abs_vaf_residual()).c_str(),
                      fmt(b.copy_count_accuracy()).c_str());
        os << line;
    }
    os << "------------------------------------------------------------------"
          "----------\n";
}

}  // namespace branch::analysis
