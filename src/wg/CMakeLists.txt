# BRANCH v0.5 — whole-genome pipeline modules.
add_library(branch_wg STATIC
    bloom_filter.cpp
    kmer_counter.cpp
    kmer_hist.cpp
    haplotype_router.cpp
    master_tiler.cpp
    master_curator.cpp
    branch_attacher.cpp
    stage1_screen.cpp
    vpcr_wg.cpp
    orphan_clusterer.cpp
    vpf_writer.cpp
)
target_include_directories(branch_wg PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(branch_wg PUBLIC branch_graph branch_common branch_vpcr)

# zstd is needed by vpf_writer.cpp. Reuse the system find done at top-level.
target_compile_features(branch_wg PUBLIC cxx_std_20)
