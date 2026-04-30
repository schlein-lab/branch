# BRANCH v0.5 — whole-genome pipeline modules.
add_library(branch_wg STATIC
    bloom_filter.cpp
    kmer_counter.cpp
    kmer_hist.cpp
)
target_include_directories(branch_wg PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(branch_wg PUBLIC branch_graph branch_common)
target_compile_features(branch_wg PUBLIC cxx_std_20)
