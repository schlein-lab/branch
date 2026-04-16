# BRANCH - Breakpoint-Resolved Assembly of Non-diploid Copy-number Heterogeneity
# Somatic-mosaicism-aware genome assembler for PacBio HiFi bulk-blood
# https://github.com/schlein-lab/branch

cmake_minimum_required(VERSION 3.28)
project(branch VERSION 0.1.0 LANGUAGES CXX)

# ------------------------------------------------------------------------------
# C++ Standard
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------------------------------------------------------------------
# Build Options
# ------------------------------------------------------------------------------
option(BRANCH_BUILD_CUDA "Build CUDA backends" ON)
option(BRANCH_BUILD_TESTS "Build test suite" ON)
option(BRANCH_USE_ASAN "Enable AddressSanitizer" OFF)
option(BRANCH_USE_TSAN "Enable ThreadSanitizer" OFF)
option(BRANCH_USE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

# ------------------------------------------------------------------------------
# CUDA Configuration (optional)
# ------------------------------------------------------------------------------
if(BRANCH_BUILD_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 20)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    # Fat binary for A100 (sm_80) and H100 (sm_90)
    set(CMAKE_CUDA_ARCHITECTURES 80 90)
endif()

# ------------------------------------------------------------------------------
# Compiler Warnings
# ------------------------------------------------------------------------------
add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -Werror
    -Wconversion
    -Wshadow
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
# Subdirectories
# ------------------------------------------------------------------------------
add_subdirectory(src)

if(BRANCH_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
