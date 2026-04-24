# branch_gpu: overlap-kernel skeleton. Compiles differently depending on
# the resolved BRANCH_CUDA_ENABLED flag (set by the root CMakeLists from
# BRANCH_BUILD_CUDA=AUTO|ON|OFF): a CUDA .cu file when ON, a CPU stub
# .cpp when OFF. Either way the header branch/gpu/overlap_kernel.cuh
# stays usable.

if(BRANCH_CUDA_ENABLED)
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
# branch_gpu only uses branch_backend's backend_vtable.hpp for the
# OverlapPair struct layout — header-only, no linked symbols — so we
# do not add branch_backend as a link dependency. This is deliberate:
# branch_backend links branch_gpu for launch_overlap_kernel, so the
# opposite link direction would create a cycle.
