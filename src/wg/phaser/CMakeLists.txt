# BRANCH v0.6 phaser — modular replacement for monolithic Phase 1.1 routing.
#
# Each translation unit is independently testable. No sub-module reaches
# into another's internals; they communicate through the typed structures
# in types.hpp.

add_library(branch_phaser STATIC
    overlap_graph.cpp
    bubble_finder.cpp
    phase_engine.cpp
    refiner.cpp
    balance_qc.cpp
    gfa_writer.cpp
    phaser.cpp
)
target_include_directories(branch_phaser PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_features(branch_phaser PUBLIC cxx_std_20)

# OpenMP is optional; overlap_graph parallelises sketching when present.
find_package(OpenMP QUIET)
if (OpenMP_CXX_FOUND)
    target_link_libraries(branch_phaser PUBLIC OpenMP::OpenMP_CXX)
endif()

# Tests — compiled only when BUILD_TESTING is on.
if (BUILD_TESTING)
    enable_testing()
    add_executable(test_phaser_toy tests/test_phaser_toy.cpp)
    target_link_libraries(test_phaser_toy PRIVATE branch_phaser)
    add_test(NAME phaser_toy COMMAND test_phaser_toy)

    add_executable(test_balance_qc tests/test_balance_qc.cpp)
    target_link_libraries(test_balance_qc PRIVATE branch_phaser)
    add_test(NAME balance_qc COMMAND test_balance_qc)
endif()
