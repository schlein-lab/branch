# branch_common: process-wide utilities (memory budget, peak-RSS
# tracking, byte-size parsing). Header-light module that every CLI
# subcommand can depend on without dragging in the heavier graph /
# backend libraries.

add_library(branch_common STATIC
    memory.cpp
)
target_include_directories(branch_common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_common PUBLIC cxx_std_20)
