add_library(branch_vpcr STATIC vpcr.cpp)
target_include_directories(branch_vpcr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(branch_vpcr PUBLIC cxx_std_20)
