# BRANCH — Alignment Module

# LLmap alignment backend — the single seam to the LLmap engine. This is the
# ONLY target compiled as C++23 (LLmap requires it) and the only one that
# includes LLmap headers / links its static libs. Everything else stays C++20
# and consumes the plain C++20 header llmap_align_backend.hpp.
add_library(branch_llmap_align STATIC
    llmap_align_backend.cpp
)
target_include_directories(branch_llmap_align PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(branch_llmap_align PUBLIC cxx_std_23)
target_link_libraries(branch_llmap_align
    PUBLIC z              # zlib for the FASTA reader
    PRIVATE branch_llmap  # LLmap classical static libs (C++23)
)

# Reference aligner — query-vs-reference mapping built on the LLmap backend.
add_library(branch_align STATIC
    reference_aligner.cpp
)
target_include_directories(branch_align PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(branch_align PUBLIC cxx_std_20)
target_link_libraries(branch_align PUBLIC branch_llmap_align)

# Tests — standalone (no GoogleTest), built when BUILD_TESTING is on.
if (BUILD_TESTING)
    enable_testing()
    add_executable(test_llmap_align_backend tests/test_llmap_align_backend.cpp)
    target_include_directories(test_llmap_align_backend PRIVATE ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(test_llmap_align_backend PRIVATE branch_llmap_align)
    add_test(NAME llmap_align_backend COMMAND test_llmap_align_backend)
endif()
