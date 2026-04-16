# branch_classify: header-only for v0.1; classifier implementation arrives v0.3.
add_library(branch_classify INTERFACE)
target_include_directories(branch_classify INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_classify INTERFACE cxx_std_20)
