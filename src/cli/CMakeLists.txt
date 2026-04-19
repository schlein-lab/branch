add_executable(branch
    main.cpp
    analyze_cmd.cpp
    assemble_cmd.cpp
)
target_link_libraries(branch PRIVATE branch_analysis branch_detect branch_io branch_backend branch_graph branch_project)
target_compile_features(branch PRIVATE cxx_std_20)
