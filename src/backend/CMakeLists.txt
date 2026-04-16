# branch_backend: Type-Erased VTable headers + CPU reference implementation.
# v0.1 ships a CPU-only backend. GPU kernels live in branch_gpu.
add_library(branch_backend STATIC
    cpu_backend.cpp
)
target_include_directories(branch_backend PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_backend PUBLIC cxx_std_20)
target_link_libraries(branch_backend PUBLIC branch_graph branch_classify)
