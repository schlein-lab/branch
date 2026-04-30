# branch_gpu: GPU kernels for the BRANCH pipeline.
#
# Compiles differently depending on the resolved BRANCH_CUDA_ENABLED
# flag (set by the root CMakeLists from BRANCH_BUILD_CUDA=AUTO|ON|OFF):
# CUDA .cu files when ON, CPU stub .cpp files when OFF. Either way the
# headers (overlap_kernel.cuh, wg_kernels.cuh) stay usable.
#
# overlap_kernel.cu (the v0.2 minimizer all-vs-all kernel) is NOT
# rebuilt under CUDA 12.9 — its CUB::DeviceRadixSort over a custom
# MinimizerEntry type needs a CUB-2.x decomposer that the v0.2 code
# doesn't provide. The CPU stub `overlap_kernel_stub.cpp` is linked
# in instead so callers (gpu_backend.cpp) still resolve the symbol;
# they do not invoke the kernel because BRANCH_CUDA_ENABLED on the
# new wg_kernels path is what matters for v0.5 GPU acceleration.
# A future commit can revive the v0.2 kernel with a proper decomposer.

if(BRANCH_CUDA_ENABLED)
    add_library(branch_gpu STATIC
        overlap_kernel_stub.cpp
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
