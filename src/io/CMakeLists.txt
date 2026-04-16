add_library(branch_io STATIC
    mosdepth_reader.cpp
    fastq_reader.cpp
    bam_reader.cpp
)
target_include_directories(branch_io PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
# htslib headers are vendored in third_party/htslib and trip -Werror=conversion
# (kstring.h, hts.h, sam.h contain legitimate narrowing conversions inside the
# upstream implementation). Treat them as SYSTEM so their warnings are silenced
# without loosening our own strict flags.
target_include_directories(branch_io SYSTEM PRIVATE
    ${CMAKE_SOURCE_DIR}/third_party/htslib
)
target_compile_features(branch_io PUBLIC cxx_std_20)
# BAM reader uses htslib. Everything else in branch_io is header-only
# w.r.t. third-party deps, but a PUBLIC link is fine — consumers that
# only include fastq_reader.hpp won't pay any runtime cost.
target_link_libraries(branch_io PUBLIC branch_htslib)
