# branch_gpu: GPU kernels for the BRANCH pipeline.
#
# Compiles differently depending on the resolved BRANCH_CUDA_ENABLED
# flag (set by the root CMakeLists from BRANCH_BUILD_CUDA=AUTO|ON|OFF):
# CUDA .cu files when ON, CPU stub .cpp files when OFF. Either way the
# headers (overlap_kernel.cuh, wg_kernels.cuh) stay usable.

if(BRANCH_CUDA_ENABLED)
    add_library(branch_gpu STATIC
        overlap_kernel.cu
        wg_kernels.cu
    )
    set_target_properties(branch_gpu PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
    )
else()
    add_library(branch_gpu STATIC
        overlap_kernel_stub.cpp
        wg_kernels_stub.cpp
    )
endif()

target_include_directories(branch_gpu PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_gpu PUBLIC cxx_std_20)
# branch_gpu does not link other branch libs; it only consumes their
# headers (e.g. wg/stage1_screen.hpp for Stage1Window struct layout).
# branch_backend already links branch_gpu, and branch_wg links into the
# CLI; the new wg_kernels are reached via that chain.
