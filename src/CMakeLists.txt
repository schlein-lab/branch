# BRANCH source tree
# Sub-modules are added individually as they grow; v0.1 starts with graph/ + backend/.
# project/ is added before graph/ because branch_graph now depends on
# branch_project (write_bed_with_refs uses linear_mapper).
add_subdirectory(common)
add_subdirectory(project)
add_subdirectory(graph)
add_subdirectory(classify)
add_subdirectory(backend)
add_subdirectory(gpu)
add_subdirectory(io)
add_subdirectory(detect)
add_subdirectory(analysis)
add_subdirectory(consensus)
add_subdirectory(align)
add_subdirectory(cli)
add_subdirectory(vpcr)
