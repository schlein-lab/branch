add_executable(branch
    main.cpp
    analyze_cmd.cpp
)
target_link_libraries(branch PRIVATE branch_analysis branch_io)
target_compile_features(branch PRIVATE cxx_std_20)
