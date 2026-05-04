// BRANCH v0.5 — GPU-kernel CPU stub.
//
// Compiled when BRANCH_CUDA_ENABLED=OFF. Every launcher returns false
// so callers fall back to the CPU path. Symbols match wg_kernels.cu so
// CPU-only builds link cleanly.

#include "gpu/wg_kernels.cuh"

#include "wg/stage1_screen.hpp"
#include "wg/vpcr_wg.hpp"
#include "wg/master_curator.hpp"

namespace branch::gpu::wg_kernels {

bool is_gpu_available() noexcept { return false; }

bool launch_stage1(
    std::uint32_t,
    const std::string&,
    const std::vector<std::uint64_t>&,
    const std::vector<std::uint32_t>&,
    std::size_t,
    std::size_t,
    double, double, double,
    std::vector<branch::wg::Stage1Window>&) noexcept { return false; }

bool launch_stage2(
    std::vector<branch::wg::VpcrAmplicon>&,
    const std::vector<std::pair<std::string, std::string>>&,
    int, int, int, int) noexcept { return false; }

bool launch_phase15_curation(
    std::uint32_t, const std::string&,
    const std::vector<std::pair<std::string, std::string>>&,
    std::size_t, std::size_t, int, double, int,
    std::vector<branch::wg::CurationEvent>&,
    std::vector<branch::wg::BranchCandidate>&) noexcept { return false; }

bool launch_phase0_count(
    const std::vector<std::pair<std::string, std::string>>&,
    std::size_t,
    std::vector<std::uint64_t>&,
    std::vector<std::uint32_t>&) noexcept { return false; }

bool launch_phase11_het_pairs(
    const std::vector<std::uint64_t>&,
    const std::vector<std::uint32_t>&,
    std::size_t, int, double, double, double, double,
    std::vector<std::uint64_t>&,
    std::vector<std::uint64_t>&) noexcept { return false; }

bool launch_phase2_direct_attach(
    const std::vector<std::string>&,
    const std::vector<std::pair<std::string, std::string>>&,
    std::size_t, double,
    std::vector<std::int32_t>&) noexcept { return false; }

}  // namespace branch::gpu::wg_kernels
