# BRANCH v0.4.3 — Project module (three-layer reference projection)
#
# Linear layer (v0.4.1), Pangenome (v0.4.2), Somatic delta (v0.4.3).

find_package(ZLIB REQUIRED)

add_library(branch_project STATIC
    project_cmd.cpp
    linear_mapper.cpp
    pangenome_mapper.cpp
    somatic_delta.cpp
    project_report.cpp
)

target_include_directories(branch_project PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(branch_project PUBLIC cxx_std_20)

# Link minimap2 for ksw2 alignment (somatic_delta) and zlib for gzip support
target_link_libraries(branch_project PUBLIC branch_minimap2 ZLIB::ZLIB)
