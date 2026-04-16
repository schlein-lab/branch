# branch_graph: header-only types + GFA I/O implementation.
add_library(branch_graph STATIC
    graph_io.cpp
)
target_include_directories(branch_graph PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_graph PUBLIC cxx_std_20)
