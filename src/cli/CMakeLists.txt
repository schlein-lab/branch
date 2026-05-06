add_executable(branch
    main.cpp
    wg_cmd.cpp
    analyze_cmd.cpp
    assemble_cmd.cpp
)
target_link_libraries(branch PRIVATE branch_analysis branch_classify branch_detect branch_io branch_backend branch_graph branch_project branch_common branch_phase branch_wg branch_phaser)
target_compile_features(branch PRIVATE cxx_std_20)
