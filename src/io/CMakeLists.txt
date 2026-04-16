add_library(branch_io STATIC
    mosdepth_reader.cpp
)
target_include_directories(branch_io PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_io PUBLIC cxx_std_20)
