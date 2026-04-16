add_library(branch_consensus STATIC
    majority_voter.cpp
)
target_include_directories(branch_consensus PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_consensus PUBLIC cxx_std_20)
