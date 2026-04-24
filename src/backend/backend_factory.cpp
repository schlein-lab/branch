// BRANCH — Backend factory implementation.

#include "backend/backend_factory.hpp"

#include <cstdio>

#include "backend/cpu_backend.hpp"
#include "backend/gpu_backend.hpp"

namespace branch::backend {

BackendMode parse_backend_mode(std::string_view s, std::string* err) noexcept {
    if (s == "auto" || s.empty()) return BackendMode::Auto;
    if (s == "cpu") return BackendMode::Cpu;
    if (s == "gpu") return BackendMode::Gpu;
    if (err) {
        *err = "unrecognised backend mode '";
        err->append(s);
        err->append("'; valid: auto, cpu, gpu. Falling back to auto.");
    }
    return BackendMode::Auto;
}

Backend make_backend(BackendMode mode) noexcept {
    switch (mode) {
        case BackendMode::Cpu: {
            Backend b = make_cpu_backend();
            std::fprintf(stderr, "[branch backend] selected cpu-reference (forced)\n");
            return b;
        }
        case BackendMode::Gpu: {
            Backend b = make_gpu_backend();
            if (b.empty()) {
                std::fprintf(stderr,
                             "[branch backend] gpu requested but not available "
                             "(compiled_in=%d). Returning empty backend; caller "
                             "should decide whether to fail or retry with auto.\n",
                             gpu_backend_compiled_in() ? 1 : 0);
            } else {
                std::fprintf(stderr, "[branch backend] selected %s (forced)\n",
                             b.name());
            }
            return b;
        }
        case BackendMode::Auto:
        default: {
            Backend gpu = make_gpu_backend();
            if (!gpu.empty()) {
                std::fprintf(stderr, "[branch backend] selected %s (auto, GPU available)\n",
                             gpu.name());
                return gpu;
            }
            Backend cpu = make_cpu_backend();
            std::fprintf(stderr,
                         "[branch backend] selected cpu-reference (auto, GPU %s)\n",
                         gpu_backend_compiled_in()
                             ? "compiled but no device"
                             : "not compiled in");
            return cpu;
        }
    }
}

}  // namespace branch::backend
