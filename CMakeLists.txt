# BRANCH - Breakpoint-Resolved Assembly of Non-diploid Copy-number Heterogeneity
# Somatic-mosaicism-aware genome assembler for PacBio HiFi bulk-blood
# https://github.com/schlein-lab/branch

cmake_minimum_required(VERSION 3.28)
project(branch VERSION 0.1.0 LANGUAGES C CXX)

# ------------------------------------------------------------------------------
# C++ Standard
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# Build Options
# ------------------------------------------------------------------------------
# BRANCH_BUILD_CUDA: AUTO = probe for a CUDA toolchain and enable the
# GPU build when found, silently fall back to CPU-only otherwise.
# Explicit ON / OFF overrides the probe. Default AUTO so the same build
# recipe works on GPU login nodes (enables GPU path) and plain CPU
# workstations / CI (stays CPU-only) without any flag flipping.
set(BRANCH_BUILD_CUDA "AUTO" CACHE STRING
    "Build CUDA backends: AUTO (detect), ON (required), OFF (force-disable)")
set_property(CACHE BRANCH_BUILD_CUDA PROPERTY STRINGS AUTO ON OFF)

option(BRANCH_BUILD_TESTS "Build test suite" ON)
option(BRANCH_USE_ASAN "Enable AddressSanitizer" OFF)
option(BRANCH_USE_TSAN "Enable ThreadSanitizer" OFF)
option(BRANCH_USE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

# ------------------------------------------------------------------------------
# CUDA Configuration
# ------------------------------------------------------------------------------
# Resolve BRANCH_BUILD_CUDA=AUTO via CUDAToolkit probe. We do the probe
# in a separate block so the ON / OFF paths below see a boolean value.
set(BRANCH_CUDA_ENABLED OFF)
if(BRANCH_BUILD_CUDA STREQUAL "AUTO")
    # find_package in QUIET mode so missing toolchain is silent on CPU boxes.
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        message(STATUS "BRANCH_BUILD_CUDA=AUTO: CUDA toolchain found (${CMAKE_CUDA_COMPILER}), enabling GPU build")
        set(BRANCH_CUDA_ENABLED ON)
    else()
        message(STATUS "BRANCH_BUILD_CUDA=AUTO: no CUDA toolchain found, GPU build disabled (CPU-only binary)")
    endif()
elseif(BRANCH_BUILD_CUDA)
    # Explicit ON — require CUDA and fail loudly if missing.
    message(STATUS "BRANCH_BUILD_CUDA=ON: GPU build required")
    set(BRANCH_CUDA_ENABLED ON)
else()
    message(STATUS "BRANCH_BUILD_CUDA=OFF: GPU build force-disabled")
endif()

if(BRANCH_CUDA_ENABLED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    # Fat binary for A100 (sm_80) and H100 (sm_90)
    set(CMAKE_CUDA_ARCHITECTURES 80 90)
endif()

# ------------------------------------------------------------------------------
# Compiler Warnings (CXX only — nvcc chokes on many host-compiler flags)
# ------------------------------------------------------------------------------
add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wextra>
    $<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>
    $<$<COMPILE_LANGUAGE:CXX>:-Werror>
    $<$<COMPILE_LANGUAGE:CXX>:-Wconversion>
    $<$<COMPILE_LANGUAGE:CXX>:-Wshadow>
)

# ------------------------------------------------------------------------------
# Sanitizers
# ------------------------------------------------------------------------------
if(BRANCH_USE_ASAN)
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
endif()

if(BRANCH_USE_TSAN)
    add_compile_options(-fsanitize=thread)
    add_link_options(-fsanitize=thread)
endif()

if(BRANCH_USE_UBSAN)
    add_compile_options(-fsanitize=undefined)
    add_link_options(-fsanitize=undefined)
endif()

# ------------------------------------------------------------------------------
# HTSlib (vendored under third_party/htslib — built by its own Makefile)
# ------------------------------------------------------------------------------
# BAM/CRAM/VCF I/O. We deliberately vendor rather than depend on a system
# package because Ubuntu ships libhts3 without matching dev headers in some
# configurations. `third_party/htslib/libhts.a` is the pre-built static
# library (the bundled Makefile is invoked once at setup time); we only
# consume its headers + .a here.
set(BRANCH_HTSLIB_DIR "${CMAKE_SOURCE_DIR}/third_party/htslib")
if(NOT EXISTS "${BRANCH_HTSLIB_DIR}/libhts.a")
    message(FATAL_ERROR
        "Vendored htslib not built. Run:\n"
        "  make -C ${BRANCH_HTSLIB_DIR} -j4 libhts.a")
endif()
add_library(branch_htslib INTERFACE)
target_include_directories(branch_htslib INTERFACE "${BRANCH_HTSLIB_DIR}")
# htslib transitive deps. htscodecs and the rANS/RLE codecs are compiled
# into libhts.a from third_party/htslib's bundled htscodecs sources (see
# `htscodecs_bundled.mk` and the configure output `linking
# htscodecs_bundled.mk to htscodecs.mk`). libdeflate is autodetected at
# htslib configure time; probe for it here and only link if present so the
# build works on hosts without libdeflate (e.g. HPC login nodes).
target_link_libraries(branch_htslib INTERFACE
    "${BRANCH_HTSLIB_DIR}/libhts.a"
    z m pthread)
find_library(BRANCH_LIBDEFLATE NAMES deflate)
if(BRANCH_LIBDEFLATE)
    target_link_libraries(branch_htslib INTERFACE ${BRANCH_LIBDEFLATE})
endif()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# OpenMP — required for parallel hotpaths in graph_filter / graph_compactor.
# Without this the transitive-reduction + unitig-consensus passes sit at
# ~6 % CPU on 32-thread nodes (one core used, the rest idle).
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------------------------
# minimap2 (FetchContent for ksw2 alignment)
# ------------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(
    minimap2
    GIT_REPOSITORY https://github.com/lh3/minimap2.git
    GIT_TAG        v2.28
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(minimap2)

# Create static library for ksw2 alignment (compiled from minimap2 C sources)
# We need ksw2_extz2_sse.c and kalloc.c for the ksw_extz2_sse function.
add_library(branch_minimap2 STATIC
    "${minimap2_SOURCE_DIR}/ksw2_extz2_sse.c"
    "${minimap2_SOURCE_DIR}/kalloc.c"
)
target_include_directories(branch_minimap2 SYSTEM PUBLIC "${minimap2_SOURCE_DIR}")
# Disable warnings for vendored C code
target_compile_options(branch_minimap2 PRIVATE -w)

# Subdirectories
# ------------------------------------------------------------------------------
add_subdirectory(src)

if(BRANCH_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
