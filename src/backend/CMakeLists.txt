# branch_backend: header-only for v0.1.
# Actual CPU/GPU backend implementations arrive later (v0.1+ CUDA kernel).
add_library(branch_backend INTERFACE)
target_include_directories(branch_backend INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_backend INTERFACE cxx_std_20)
