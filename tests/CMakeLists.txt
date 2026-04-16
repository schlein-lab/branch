# GoogleTest via FetchContent (self-contained, no system dep).
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0
    SYSTEM
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

add_executable(test_delta_read test_delta_read.cpp)
target_link_libraries(test_delta_read PRIVATE branch_graph branch_backend GTest::gtest_main)
target_compile_features(test_delta_read PRIVATE cxx_std_20)
gtest_discover_tests(test_delta_read)

add_executable(test_backend_vtable test_backend_vtable.cpp)
target_link_libraries(test_backend_vtable PRIVATE branch_backend GTest::gtest_main)
target_compile_features(test_backend_vtable PRIVATE cxx_std_20)
gtest_discover_tests(test_backend_vtable)

add_executable(test_features test_features.cpp)
target_link_libraries(test_features PRIVATE branch_classify GTest::gtest_main)
target_compile_features(test_features PRIVATE cxx_std_20)
gtest_discover_tests(test_features)

add_executable(test_lossless_graph test_lossless_graph.cpp)
target_link_libraries(test_lossless_graph PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_lossless_graph PRIVATE cxx_std_20)
gtest_discover_tests(test_lossless_graph)

add_executable(test_reverse_index test_reverse_index.cpp)
target_link_libraries(test_reverse_index PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_reverse_index PRIVATE cxx_std_20)
gtest_discover_tests(test_reverse_index)

add_executable(test_cascade test_cascade.cpp)
target_link_libraries(test_cascade PRIVATE branch_classify branch_backend GTest::gtest_main)
target_compile_features(test_cascade PRIVATE cxx_std_20)
gtest_discover_tests(test_cascade)

add_executable(test_gpu_stub test_gpu_stub.cpp)
target_link_libraries(test_gpu_stub PRIVATE branch_gpu GTest::gtest_main)
target_compile_features(test_gpu_stub PRIVATE cxx_std_20)
gtest_discover_tests(test_gpu_stub)

add_executable(test_kmer_sketch test_kmer_sketch.cpp)
target_link_libraries(test_kmer_sketch PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_kmer_sketch PRIVATE cxx_std_20)
gtest_discover_tests(test_kmer_sketch)

add_executable(test_minimizer_sketcher test_minimizer_sketcher.cpp)
target_link_libraries(test_minimizer_sketcher PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_minimizer_sketcher PRIVATE cxx_std_20)
gtest_discover_tests(test_minimizer_sketcher)

add_executable(test_graph_builder test_graph_builder.cpp)
target_link_libraries(test_graph_builder PRIVATE branch_graph branch_backend GTest::gtest_main)
target_compile_features(test_graph_builder PRIVATE cxx_std_20)
gtest_discover_tests(test_graph_builder)

add_executable(test_bubble_detector test_bubble_detector.cpp)
target_link_libraries(test_bubble_detector PRIVATE branch_detect GTest::gtest_main)
target_compile_features(test_bubble_detector PRIVATE cxx_std_20)
gtest_discover_tests(test_bubble_detector)

add_executable(test_cpu_backend test_cpu_backend.cpp)
target_link_libraries(test_cpu_backend PRIVATE branch_backend GTest::gtest_main)
target_compile_features(test_cpu_backend PRIVATE cxx_std_20)
gtest_discover_tests(test_cpu_backend)

add_executable(test_graph_io test_graph_io.cpp)
target_link_libraries(test_graph_io PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_graph_io PRIVATE cxx_std_20)
gtest_discover_tests(test_graph_io)

add_executable(test_graph_compactor test_graph_compactor.cpp)
target_link_libraries(test_graph_compactor PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_graph_compactor PRIVATE cxx_std_20)
gtest_discover_tests(test_graph_compactor)

add_executable(test_majority_voter test_majority_voter.cpp)
target_link_libraries(test_majority_voter PRIVATE branch_consensus GTest::gtest_main)
target_compile_features(test_majority_voter PRIVATE cxx_std_20)
gtest_discover_tests(test_majority_voter)

add_executable(test_graph_filter test_graph_filter.cpp)
target_link_libraries(test_graph_filter PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_graph_filter PRIVATE cxx_std_20)
gtest_discover_tests(test_graph_filter)

add_executable(test_strand_handling test_strand_handling.cpp)
target_link_libraries(test_strand_handling PRIVATE branch_backend branch_graph GTest::gtest_main)
target_compile_features(test_strand_handling PRIVATE cxx_std_20)
gtest_discover_tests(test_strand_handling)

# --- E2E integration test -------------------------------------------------
# Spawns the real `branch` binary on a synthetic FASTQ and validates the
# emitted GFA. Gated on the `branch` target being defined so a standalone
# tests build (unlikely, but possible) still configures.
if(TARGET branch)
    add_executable(generate_toy_fastq fixtures/generate_toy_fastq.cpp)
    target_compile_features(generate_toy_fastq PRIVATE cxx_std_20)

    add_executable(test_e2e_assemble test_e2e_assemble.cpp)
    target_link_libraries(test_e2e_assemble PRIVATE GTest::gtest GTest::gtest_main)
    target_compile_features(test_e2e_assemble PRIVATE cxx_std_20)
    target_compile_definitions(test_e2e_assemble PRIVATE
        BRANCH_BINARY="$<TARGET_FILE:branch>"
        GEN_FASTQ_BINARY="$<TARGET_FILE:generate_toy_fastq>")
    add_dependencies(test_e2e_assemble branch generate_toy_fastq)
    add_test(NAME E2E_Assemble COMMAND test_e2e_assemble)
endif()
