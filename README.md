# BRANCH

**Breakpoint-Resolved Assembly of Non-diploid Copy-number Heterogeneity**

## What

Somatic-mosaicism-aware genome assembler for human PacBio HiFi bulk-blood data with vPCR-based CNV quantification.

BRANCH builds a lossless Delta-Graph representation where reads are stored as path + delta-list, not as sequences. This enables efficient classification of parallel graph paths as Branch vs. Duplication vs. Mixed using flanking identity, depth-sum, read-span, and reference priors.

## Status

**v0.1** — In development, alpha.

## Out of Scope

- Plants
- Tumors
- Cancer subclones

## Dependencies

- C++20 compiler (GCC 13+ or Clang 17+)
- CMake 3.28+
- CUDA 12+ (optional, for GPU acceleration)
- htslib
- ksw2
- abPOA

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

## License

MIT — see [LICENSE](LICENSE)
