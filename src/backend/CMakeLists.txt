# branch_backend: Type-Erased VTable + concrete backends.
#
# Always includes:
#   cpu_backend.cpp        — CPU reference implementation.
#   backend_factory.cpp    — make_backend(mode) selection logic.
#
# GPU implementation is picked between two mutually-exclusive TUs
# based on BRANCH_CUDA_ENABLED (resolved by the root CMakeLists from
# BRANCH_BUILD_CUDA=AUTO|ON|OFF):
#   gpu_backend.cpp        — real CUDA impl, links libcudart.
#   gpu_backend_stub.cpp   — empty make_gpu_backend() returning {}.
# Callers always get a linkable make_gpu_backend() symbol; emptiness
# of the returned Backend tells them GPU is unavailable.

set(branch_backend_sources
    cpu_backend.cpp
    backend_factory.cpp
)

if(BRANCH_CUDA_ENABLED)
    list(APPEND branch_backend_sources gpu_backend.cpp)
else()
    list(APPEND branch_backend_sources gpu_backend_stub.cpp)
endif()

add_library(branch_backend STATIC ${branch_backend_sources})
target_include_directories(branch_backend PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_backend PUBLIC cxx_std_20)
target_link_libraries(branch_backend PUBLIC branch_graph branch_classify)

if(BRANCH_CUDA_ENABLED)
    # gpu_backend.cpp needs CUDA runtime headers + libcudart. branch_gpu
    # already carries the kernel translation unit; linking it here also
    # pulls in the driver runtime transitively.
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(branch_backend PUBLIC branch_gpu CUDA::cudart)
    target_compile_definitions(branch_backend PUBLIC BRANCH_HAVE_CUDA=1)
else()
    # Still need the overlap_kernel header path for gpu_backend_stub.cpp
    # (empty, but includes gpu_backend.hpp which has no CUDA include).
    # branch_gpu carries the stub launch_overlap_kernel symbol.
    target_link_libraries(branch_backend PUBLIC branch_gpu)
endif()
