# branch_classify: cascade header + v0.2 feature extractor.
#
# v0.1 was INTERFACE-only (cascade is header-inlineable). v0.2 adds
# feature_extractor.cpp which touches graph::LosslessGraph, so the
# target is now a STATIC library that links branch_graph transitively.
# The LightGBM classifier implementation still arrives with v0.3.
add_library(branch_classify STATIC
    feature_extractor.cpp
)
target_include_directories(branch_classify PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_classify PUBLIC cxx_std_20)
target_link_libraries(branch_classify PUBLIC branch_graph)
