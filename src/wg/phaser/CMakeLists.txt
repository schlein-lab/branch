# BRANCH v0.6 phaser — modular replacement for monolithic Phase 1.1 routing.
#
# Each translation unit is independently testable. No sub-module reaches
# into another's internals; they communicate through the typed structures
# in types.hpp.

add_library(branch_phaser STATIC
    fastq_index.cpp
    overlap_graph.cpp
    bubble_finder.cpp
    bubble_classifier.cpp
    bubble_voter.cpp
    phase_engine.cpp
    refiner.cpp
    balance_qc.cpp
    gfa_writer.cpp
    phaser.cpp
    # v0.8 dual-phaser scaffold
    anchored/anchored_phaser.cpp
    reconcile/reconciler.cpp
    reconcile/mlp.cpp
)
target_include_directories(branch_phaser PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_compile_features(branch_phaser PUBLIC cxx_std_20)

# zlib for FastqIndex gz handling.
find_package(ZLIB REQUIRED)
target_link_libraries(branch_phaser PUBLIC ZLIB::ZLIB)

# anchored_phaser maps reads via the LLmap backend (no minimap2 subprocess).
target_link_libraries(branch_phaser PUBLIC branch_llmap_align)

# OpenMP is optional; overlap_graph parallelises sketching when present.
find_package(OpenMP QUIET)
if (OpenMP_CXX_FOUND)
    target_link_libraries(branch_phaser PUBLIC OpenMP::OpenMP_CXX)
endif()

# Tests — compiled only when BUILD_TESTING is on.
if (BUILD_TESTING)
    enable_testing()
    # These tests assert() their expectations. Keep assertions live even in
    # Release builds (which otherwise define NDEBUG and turn assert() into a
    # no-op, leaving the checked variables unused under -Werror).
    add_compile_options(-UNDEBUG)
    add_executable(test_phaser_toy tests/test_phaser_toy.cpp)
    target_link_libraries(test_phaser_toy PRIVATE branch_phaser)
    add_test(NAME phaser_toy COMMAND test_phaser_toy)

    add_executable(test_balance_qc tests/test_balance_qc.cpp)
    target_link_libraries(test_balance_qc PRIVATE branch_phaser)
    add_test(NAME balance_qc COMMAND test_balance_qc)

    add_executable(test_fastq_index tests/test_fastq_index.cpp)
    target_link_libraries(test_fastq_index PRIVATE branch_phaser)
    add_test(NAME fastq_index COMMAND test_fastq_index)

    add_executable(test_transitive_reduction tests/test_transitive_reduction.cpp)
    target_link_libraries(test_transitive_reduction PRIVATE branch_phaser)
    add_test(NAME transitive_reduction COMMAND test_transitive_reduction)

    add_executable(test_reconciler tests/test_reconciler.cpp)
    target_link_libraries(test_reconciler PRIVATE branch_phaser)
    add_test(NAME reconciler COMMAND test_reconciler)

    add_executable(test_bubble_classifier tests/test_bubble_classifier.cpp)
    target_link_libraries(test_bubble_classifier PRIVATE branch_phaser)
    add_test(NAME bubble_classifier COMMAND test_bubble_classifier)

    add_executable(test_bubble_voter tests/test_bubble_voter.cpp)
    target_link_libraries(test_bubble_voter PRIVATE branch_phaser)
    add_test(NAME bubble_voter COMMAND test_bubble_voter)

    add_executable(test_mlp tests/test_mlp.cpp)
    target_link_libraries(test_mlp PRIVATE branch_phaser)
    add_test(NAME mlp COMMAND test_mlp)
endif()
