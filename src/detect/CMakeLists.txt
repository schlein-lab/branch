add_library(branch_detect STATIC
    bubble_detector.cpp
)
target_include_directories(branch_detect PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_detect PUBLIC cxx_std_20)
target_link_libraries(branch_detect PUBLIC branch_graph)
