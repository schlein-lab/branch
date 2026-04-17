# branch_graph: data types + GFA I/O + minimizer sketcher + builder + compactor.
add_library(branch_graph STATIC
    graph_io.cpp
    minimizer_sketcher.cpp
    graph_builder.cpp
    graph_compactor.cpp
    graph_filter.cpp
    delta_read.cpp
)
target_include_directories(branch_graph PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_graph PUBLIC cxx_std_20)
target_link_libraries(branch_graph PUBLIC branch_project)
