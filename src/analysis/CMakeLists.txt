add_library(branch_analysis STATIC
    coverage_cn.cpp
    coverage_conservation.cpp
    repeat_cn.cpp
    annotation_bed.cpp
    vaf_band_stats.cpp
)
target_include_directories(branch_analysis PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_analysis PUBLIC cxx_std_20)
target_link_libraries(branch_analysis PUBLIC branch_io branch_graph branch_detect)
