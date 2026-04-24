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

add_executable(test_backend_factory test_backend_factory.cpp)
target_link_libraries(test_backend_factory PRIVATE branch_backend GTest::gtest_main)
target_compile_features(test_backend_factory PRIVATE cxx_std_20)
gtest_discover_tests(test_backend_factory)

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

add_executable(test_gpu_backend test_gpu_backend.cpp)
target_link_libraries(test_gpu_backend PRIVATE branch_gpu GTest::gtest_main)
target_compile_features(test_gpu_backend PRIVATE cxx_std_20)
gtest_discover_tests(test_gpu_backend)

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

add_executable(test_coverage_conservation test_coverage_conservation.cpp)
target_link_libraries(test_coverage_conservation PRIVATE branch_analysis branch_detect branch_graph GTest::gtest_main)
target_compile_features(test_coverage_conservation PRIVATE cxx_std_20)
gtest_discover_tests(test_coverage_conservation)

add_executable(test_cpu_backend test_cpu_backend.cpp)
target_link_libraries(test_cpu_backend PRIVATE branch_backend GTest::gtest_main)
target_compile_features(test_cpu_backend PRIVATE cxx_std_20)
gtest_discover_tests(test_cpu_backend)

add_executable(test_graph_io test_graph_io.cpp)
target_link_libraries(test_graph_io PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_graph_io PRIVATE cxx_std_20)
gtest_discover_tests(test_graph_io)

add_executable(test_write_bed_with_refs test_write_bed_with_refs.cpp)
target_link_libraries(test_write_bed_with_refs PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_write_bed_with_refs PRIVATE cxx_std_20)
gtest_discover_tests(test_write_bed_with_refs)

add_executable(test_write_fasta_consensus test_write_fasta_consensus.cpp)
target_link_libraries(test_write_fasta_consensus PRIVATE branch_graph GTest::gtest_main)
target_compile_features(test_write_fasta_consensus PRIVATE cxx_std_20)
gtest_discover_tests(test_write_fasta_consensus)

add_executable(test_cpu_thread_determinism test_cpu_thread_determinism.cpp)
target_link_libraries(test_cpu_thread_determinism PRIVATE branch_backend branch_graph GTest::gtest_main)
target_compile_features(test_cpu_thread_determinism PRIVATE cxx_std_20)
gtest_discover_tests(test_cpu_thread_determinism)

add_executable(test_project_report test_project_report.cpp)
target_link_libraries(test_project_report PRIVATE branch_project GTest::gtest_main)
target_compile_features(test_project_report PRIVATE cxx_std_20)
target_include_directories(test_project_report PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_project_report)

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

add_executable(test_bam_reader test_bam_reader.cpp)
target_link_libraries(test_bam_reader PRIVATE branch_io branch_htslib GTest::gtest_main)
target_compile_features(test_bam_reader PRIVATE cxx_std_20)
# Silence -Werror=conversion from vendored htslib headers included here.
target_include_directories(test_bam_reader SYSTEM PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party/htslib)
gtest_discover_tests(test_bam_reader)

# --- E2E integration test -------------------------------------------------
# Spawns the real `branch` binary on a synthetic FASTQ and validates the
# emitted GFA. Gated on the `branch` target being defined so a standalone
# tests build (unlikely, but possible) still configures.
if(TARGET branch)
    add_executable(generate_toy_fastq fixtures/generate_toy_fastq.cpp)
    target_compile_features(generate_toy_fastq PRIVATE cxx_std_20)

    add_executable(generate_branch_fixture fixtures/generate_branch_fixture.cpp)
    target_compile_features(generate_branch_fixture PRIVATE cxx_std_20)

    add_executable(generate_branch_complex_fixture fixtures/generate_branch_complex_fixture.cpp)
    target_compile_features(generate_branch_complex_fixture PRIVATE cxx_std_20)

    add_executable(test_e2e_assemble test_e2e_assemble.cpp)
    target_link_libraries(test_e2e_assemble PRIVATE GTest::gtest GTest::gtest_main)
    target_compile_features(test_e2e_assemble PRIVATE cxx_std_20)
    target_compile_definitions(test_e2e_assemble PRIVATE
        BRANCH_BINARY="$<TARGET_FILE:branch>"
        GEN_FASTQ_BINARY="$<TARGET_FILE:generate_toy_fastq>"
        GEN_BRANCH_FIXTURE_BINARY="$<TARGET_FILE:generate_branch_fixture>")
    add_dependencies(test_e2e_assemble branch generate_toy_fastq generate_branch_fixture)
    add_test(NAME E2E_Assemble COMMAND test_e2e_assemble)

    add_executable(test_project_smoke test_project_smoke.cpp)
    target_link_libraries(test_project_smoke PRIVATE GTest::gtest GTest::gtest_main)
    target_compile_features(test_project_smoke PRIVATE cxx_std_20)
    target_compile_definitions(test_project_smoke PRIVATE
        BRANCH_BINARY="$<TARGET_FILE:branch>")
    add_dependencies(test_project_smoke branch)
    add_test(NAME ProjectSmoke COMMAND test_project_smoke)

    # P2.1 assemble->analyze wiring. Exercises the GFA handoff:
    # `branch assemble` writes a GFA, then `branch analyze --graph`
    # loads it, detects bubbles, runs conservation, classifies via the
    # hierarchical disambiguator, and writes per-bubble BED with
    # confidences. Separate test because bubble detection on real
    # assemblies is slow (bubble_detector complexity is quadratic in
    # out-degree); this test uses a hand-crafted minimal bubble graph.
    add_executable(test_e2e_analyze_wiring test_e2e_analyze_wiring.cpp)
    target_link_libraries(test_e2e_analyze_wiring
        PRIVATE branch_graph GTest::gtest GTest::gtest_main)
    target_compile_features(test_e2e_analyze_wiring PRIVATE cxx_std_20)
    target_include_directories(test_e2e_analyze_wiring
        PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_compile_definitions(test_e2e_analyze_wiring PRIVATE
        BRANCH_BINARY="$<TARGET_FILE:branch>"
        GEN_FASTQ_BINARY="$<TARGET_FILE:generate_toy_fastq>")
    add_dependencies(test_e2e_analyze_wiring branch generate_toy_fastq)
    add_test(NAME E2E_AnalyzeWiring COMMAND test_e2e_analyze_wiring)

    # P1.3 VAF-band benchmark. Generates the 4-way bubble fixture,
    # runs branch assemble, parses the GFA, and reports per-band
    # retrieval metrics with a permissive low-band recall gate.
    add_executable(test_vaf_band_benchmark test_vaf_band_benchmark.cpp)
    target_link_libraries(test_vaf_band_benchmark
        PRIVATE branch_analysis GTest::gtest GTest::gtest_main)
    target_compile_features(test_vaf_band_benchmark PRIVATE cxx_std_20)
    target_include_directories(test_vaf_band_benchmark
        PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_compile_definitions(test_vaf_band_benchmark PRIVATE
        BRANCH_BINARY="$<TARGET_FILE:branch>"
        GEN_BRANCH_FIXTURE_BINARY="$<TARGET_FILE:generate_branch_fixture>")
    add_dependencies(test_vaf_band_benchmark branch generate_branch_fixture)
    add_test(NAME VafBandBenchmark COMMAND test_vaf_band_benchmark)
endif()

add_executable(test_reference_aligner test_reference_aligner.cpp)
target_link_libraries(test_reference_aligner PRIVATE GTest::gtest_main)
target_compile_features(test_reference_aligner PRIVATE cxx_std_20)
gtest_discover_tests(test_reference_aligner)

add_executable(test_linear_mapper test_linear_mapper.cpp)
target_link_libraries(test_linear_mapper PRIVATE branch_project GTest::gtest_main)
target_compile_features(test_linear_mapper PRIVATE cxx_std_20)
target_include_directories(test_linear_mapper PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_linear_mapper)

add_executable(test_pangenome_mapper test_pangenome_mapper.cpp)
target_link_libraries(test_pangenome_mapper PRIVATE branch_project GTest::gtest_main)
target_compile_features(test_pangenome_mapper PRIVATE cxx_std_20)
target_include_directories(test_pangenome_mapper PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_pangenome_mapper)


add_executable(test_vpcr test_vpcr.cpp)
target_link_libraries(test_vpcr PRIVATE branch_vpcr GTest::gtest_main)
target_compile_features(test_vpcr PRIVATE cxx_std_20)
target_include_directories(test_vpcr PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_vpcr)

add_executable(test_repeat_cn test_repeat_cn.cpp)
target_link_libraries(test_repeat_cn PRIVATE branch_analysis)
add_test(NAME test_repeat_cn COMMAND test_repeat_cn)

add_executable(test_annotation_bed test_annotation_bed.cpp)
target_link_libraries(test_annotation_bed PRIVATE branch_analysis GTest::gtest_main)
target_compile_features(test_annotation_bed PRIVATE cxx_std_20)
target_include_directories(test_annotation_bed PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_annotation_bed)

add_executable(test_mixed_decomposer test_mixed_decomposer.cpp)
target_link_libraries(test_mixed_decomposer PRIVATE branch_classify GTest::gtest_main)
target_compile_features(test_mixed_decomposer PRIVATE cxx_std_20)
target_include_directories(test_mixed_decomposer PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_mixed_decomposer)

add_executable(test_hierarchical_disambiguator test_hierarchical_disambiguator.cpp)
target_link_libraries(test_hierarchical_disambiguator PRIVATE branch_classify branch_backend GTest::gtest_main)
target_compile_features(test_hierarchical_disambiguator PRIVATE cxx_std_20)
target_include_directories(test_hierarchical_disambiguator PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_hierarchical_disambiguator)

add_executable(test_somatic_delta test_somatic_delta.cpp)
target_link_libraries(test_somatic_delta PRIVATE branch_project GTest::gtest_main)
target_compile_features(test_somatic_delta PRIVATE cxx_std_20)
target_include_directories(test_somatic_delta PRIVATE ${CMAKE_SOURCE_DIR}/src)
gtest_discover_tests(test_somatic_delta)
