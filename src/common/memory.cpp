// BRANCH — Memory-budget + peak-tracking implementation.

#include "common/memory.hpp"

#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>

#include <sys/resource.h>
#include <unistd.h>

namespace branch::common {

namespace {

constexpr std::size_t kKiB = 1024ULL;
constexpr std::size_t kMiB = kKiB * 1024ULL;
constexpr std::size_t kGiB = kMiB * 1024ULL;
constexpr std::size_t kTiB = kGiB * 1024ULL;

std::atomic<bool> g_reporter_installed{false};

void peak_rss_atexit() noexcept {
    const std::size_t rss = peak_rss_bytes();
    if (rss == 0) return;
    std::fprintf(stderr, "[branch] peak RSS = %s\n",
                 format_bytes(rss).c_str());
}

}  // namespace

std::optional<std::size_t>
parse_memory_size(std::string_view s) noexcept {
    // Trim leading / trailing whitespace.
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))  s.remove_suffix(1);
    if (s.empty()) return std::nullopt;

    // Split numeric prefix from suffix letters.
    std::size_t split = 0;
    while (split < s.size()) {
        const char c = s[split];
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.' || c == ',') {
            ++split;
        } else {
            break;
        }
    }
    if (split == 0) return std::nullopt;

    std::string num_str(s.substr(0, split));
    // Allow comma as a decimal separator for locale-robustness.
    for (auto& ch : num_str) if (ch == ',') ch = '.';
    char* end = nullptr;
    const double val = std::strtod(num_str.c_str(), &end);
    if (end == num_str.c_str() || val < 0.0) return std::nullopt;

    std::string_view suf = s.substr(split);
    // Lowercase + strip 'B' / 'iB' suffix variants so "G", "GB",
    // "GiB", "gb" all mean the same. IEC and SI are treated as
    // base-1024 here — the caller cares about "fits in RAM", not
    // about packaging-marketing precision.
    std::size_t mult = 1;
    if (suf.empty()) {
        mult = 1;
    } else {
        char c = static_cast<char>(std::tolower(static_cast<unsigned char>(suf.front())));
        switch (c) {
            case 'k': mult = kKiB; break;
            case 'm': mult = kMiB; break;
            case 'g': mult = kGiB; break;
            case 't': mult = kTiB; break;
            case 'b': mult = 1; break;
            default: return std::nullopt;
        }
    }

    const double bytes_d = val * static_cast<double>(mult);
    if (!std::isfinite(bytes_d) || bytes_d < 0.0) return std::nullopt;
    // Cap at SIZE_MAX - 1.
    constexpr double kCap = 18446744073709550000.0;  // slightly below 2^64
    if (bytes_d >= kCap) return std::nullopt;
    return static_cast<std::size_t>(bytes_d);
}

bool set_memory_budget(std::size_t bytes) noexcept {
    if (bytes == 0) return true;  // no-op

    struct rlimit cur{};
    if (getrlimit(RLIMIT_AS, &cur) != 0) return false;

    rlim_t want = (bytes > static_cast<std::size_t>(cur.rlim_max))
                      ? cur.rlim_max
                      : static_cast<rlim_t>(bytes);

    struct rlimit set{};
    set.rlim_cur = want;
    set.rlim_max = cur.rlim_max;  // never raise hard limit
    return setrlimit(RLIMIT_AS, &set) == 0;
}

std::size_t peak_rss_bytes() noexcept {
    struct rusage ru{};
    if (getrusage(RUSAGE_SELF, &ru) != 0) return 0;
    // Linux reports ru_maxrss in KiB; macOS in bytes. We live on
    // Linux here; ru_maxrss is long, can be up to ~2 PiB worth of
    // KiB before overflow — beyond any realistic use.
    return static_cast<std::size_t>(ru.ru_maxrss) * 1024ULL;
}

std::string format_bytes(std::size_t bytes) {
    const double b = static_cast<double>(bytes);
    double v; const char* unit;
    if (b >= static_cast<double>(kTiB))      { v = b / kTiB; unit = "TiB"; }
    else if (b >= static_cast<double>(kGiB)) { v = b / kGiB; unit = "GiB"; }
    else if (b >= static_cast<double>(kMiB)) { v = b / kMiB; unit = "MiB"; }
    else if (b >= static_cast<double>(kKiB)) { v = b / kKiB; unit = "KiB"; }
    else                                     { v = b;        unit = "B"; }

    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << v << ' ' << unit;
    return os.str();
}

void install_peak_rss_reporter() noexcept {
    bool expected = false;
    if (g_reporter_installed.compare_exchange_strong(expected, true)) {
        std::atexit(peak_rss_atexit);
    }
}

}  // namespace branch::common
