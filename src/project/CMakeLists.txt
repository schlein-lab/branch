# BRANCH v0.4 — Project module (three-layer reference projection)
#
# Stub library for v0.3. Full implementation in v0.4.

add_library(branch_project STATIC
    project_cmd.cpp
    # TODO v0.4: Add implementation files
    # linear_mapper.cpp
    # pangenome_mapper.cpp
    # somatic_delta.cpp
    # project_report.cpp
)

target_include_directories(branch_project PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_compile_features(branch_project PRIVATE cxx_std_20)

# TODO v0.4: Link dependencies
# target_link_libraries(branch_project PRIVATE branch_graph branch_io)
