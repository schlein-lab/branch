#pragma once

// BRANCH — Backend factory.
//
// Hides the "GPU or CPU?" decision behind a single make_backend(mode)
// entry point. Call sites (CLI commands, tests) specify intent, not
// physical hardware:
//
//   BackendMode::Auto  — try GPU, fall back to CPU on any failure.
//                        The common case.
//   BackendMode::Cpu   — force CPU, even on GPU hardware. Used by
//                        deterministic tests and CPU-vs-GPU regression
//                        comparisons.
//   BackendMode::Gpu   — require GPU. Returns an empty Backend if GPU
//                        is unavailable so callers can surface the
//                        error explicitly instead of silently running
//                        on CPU.
//
// Never throws. Empty Backend returns mean "requested mode not
// available" — callers decide whether to retry with a different mode
// or abort.

#include <string_view>

#include "backend/backend_vtable.hpp"

namespace branch::backend {

enum class BackendMode : std::uint8_t {
    Auto = 0,
    Cpu = 1,
    Gpu = 2,
};

// Parse a CLI string into a BackendMode. Recognised tokens: "auto",
// "cpu", "gpu". Empty / unknown input falls back to Auto and writes
// a note to *err so the CLI can surface it.
[[nodiscard]] BackendMode parse_backend_mode(std::string_view s,
                                              std::string* err = nullptr) noexcept;

// Produce a Backend honouring the requested mode. Logs the chosen
// backend name on stderr so benchmark runs are self-documenting.
[[nodiscard]] Backend make_backend(BackendMode mode) noexcept;

}  // namespace branch::backend
