# branch_phase: SNV-aware read phasing layered on top of the ReadDB.
# Greedy single-pass phaser (HapCUT2-lite) — coarse haplotype assignment
# from per-read deltas, used by bubble-VAF analysis to split per-alt
# read counts into haplotype-aware fractions.

add_library(branch_phase STATIC
    snv_phaser.cpp
)
target_include_directories(branch_phase PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
target_compile_features(branch_phase PUBLIC cxx_std_20)
target_link_libraries(branch_phase PUBLIC branch_graph)
