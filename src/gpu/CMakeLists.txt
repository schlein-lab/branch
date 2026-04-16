# branch_gpu: overlap-kernel skeleton. Compiles differently depending on
# BRANCH_BUILD_CUDA: a CUDA .cu file when ON, a CPU stub .cpp when OFF.
# Either way the header branch/gpu/overlap_kernel.cuh is usable.

if(BRANCH_BUILD_CUDA)
    add_library(branch_gpu STATIC overlap_kernel.cu)
    set_target_properties(branch_gpu PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )
else()
    add_library(branch_gpu STATIC overlap_kernel_stub.cpp)
endif()

target_include_directories(branch_gpu PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_gpu PUBLIC cxx_std_20)
target_link_libraries(branch_gpu PUBLIC branch_backend)
