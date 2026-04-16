# branch_graph: header-only library for v0.1.
# Will gain .cpp files once implementation starts.
add_library(branch_graph INTERFACE)
target_include_directories(branch_graph INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_graph INTERFACE cxx_std_20)
