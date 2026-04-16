# BRANCH — Alignment Module

add_library(branch_align STATIC
    reference_aligner.cpp
)

target_include_directories(branch_align PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link minimap2 headers (via interface library from root CMakeLists.txt)
target_link_libraries(branch_align PUBLIC
    branch_minimap2
    z  # zlib for kseq.h
)

# Suppress warnings in minimap2 C code compiled into our wrapper
target_compile_options(branch_align PRIVATE
    -Wno-sign-compare
    -Wno-unused-parameter
)
