# BRANCH v0.4 — Project module (three-layer reference projection)
#
# Linear layer implemented in v0.4.1. Pangenome and somatic delta in v0.4.2+.

add_library(branch_project STATIC
    project_cmd.cpp
    linear_mapper.cpp
    pangenome_mapper.cpp
    project_report.cpp
)

target_include_directories(branch_project PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(branch_project PUBLIC cxx_std_20)
