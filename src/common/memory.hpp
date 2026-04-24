#pragma once

// BRANCH — Memory-budget + peak-tracking utilities.
//
// Three layers of defence against OOM:
//
//  (1) parse_memory_size("8G" | "500M" | "1.5T" | "16GiB" | "16GB")
//      -> std::optional<std::size_t>. Case-insensitive suffix; both
//      decimal (SI) and binary (IEC) forms accepted and treated
//      identically as base-1024. Empty / malformed returns nullopt.
//
//  (2) set_memory_budget(bytes). Calls setrlimit(RLIMIT_AS, bytes)
//      so any allocation that would push the process past the cap
//      fails at the OS layer. glibc then throws std::bad_alloc,
//      which callers catch at the CLI top level and turn into a
//      graceful exit with the partial output flushed. Idempotent:
//      calling twice re-sets the cap. Never lowers an existing cap
//      inherited from the environment (e.g. SLURM cgroup) — if the
//      caller-specified bytes exceed the existing hard limit, we
//      clamp to the hard limit so the call can't fail due to
//      insufficient privilege.
//
//  (3) peak_rss_bytes() + format_bytes(). Zero-cost telemetry:
//      getrusage(RUSAGE_SELF).ru_maxrss reports the high-water mark
//      of resident memory (Linux reports in KiB). Pair with
//      format_bytes() for human-readable stderr summaries at process
//      exit.
//
// The helpers here do not depend on any other BRANCH module; they
// are safe to pull into main.cpp and every CLI subcommand.

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace branch::common {

// Parse a human-friendly byte-size string. Accepted suffixes:
//   (none)    -> bytes
//   k / K / kB / KiB / KB
//   m / M / mB / MiB / MB
//   g / G / gB / GiB / GB
//   t / T / tB / TiB / TB
// Decimal numbers are accepted ("1.5G" -> 1.5 * 2^30 bytes, rounded down).
// Returns std::nullopt on empty input, invalid number, or bad suffix.
[[nodiscard]] std::optional<std::size_t>
parse_memory_size(std::string_view s) noexcept;

// Apply `bytes` as an RLIMIT_AS cap for the current process. Returns
// true on success. On false, the caller should surface the error
// (getrlimit(2) failure, or caller asked to raise the cap above the
// hard limit without CAP_SYS_RESOURCE).
//
// Pass 0 to mean "no cap" — the function is a no-op in that case.
bool set_memory_budget(std::size_t bytes) noexcept;

// High-water resident set size in bytes. 0 on unsupported platforms.
[[nodiscard]] std::size_t peak_rss_bytes() noexcept;

// Pretty-print a byte count with a KiB/MiB/GiB/TiB suffix, 1 decimal
// place ("8.0 GiB", "517.4 MiB").
[[nodiscard]] std::string format_bytes(std::size_t bytes);

// Install an atexit handler that writes a single summary line to
// stderr: "[branch] peak RSS = X.Y GiB". Idempotent.
void install_peak_rss_reporter() noexcept;

}  // namespace branch::common
