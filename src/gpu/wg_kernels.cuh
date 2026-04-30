// BRANCH v0.5 — GPU kernel header for the whole-genome pipeline.
//
// One header, four kernels, host-callable launcher functions. Each
// launcher is implemented twice:
//   - in wg_kernels.cu       when BRANCH_CUDA_ENABLED=ON
//   - in wg_kernels_stub.cpp when BRANCH_CUDA_ENABLED=OFF
//
// Callers (the orchestrator + Stage1Screener / VpcrAutoDesigner) use
// `is_gpu_available()` to decide whether to dispatch GPU or fall back
// to the CPU OpenMP path. The CPU path is always present.
//
// The four kernels — in priority order:
//   1. wg_stage1     — per-window k-mer-coverage stats (Phase 3 Stage 1)
//   2. wg_stage2     — per-read seed-indexed primer scan (Phase 3 Stage 2)
//   3. wg_phase15    — per-position pile-up curation (Phase 1.5)
//   4. wg_phase0     — bloom-filtered k-mer counter (Phase 0)
//
// Memory model:
//   Each launcher allocates its own device buffers, copies inputs in,
//   launches the kernel, copies results out, frees device memory. No
//   cross-kernel device-resident state — simpler bookkeeping at the
//   cost of a few extra HtoD/DtoH transfers. A future "persistent GPU
//   state" mode could keep the read pool + counter resident across
//   phases.

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace branch::wg {
struct Stage1Window;
struct VpcrAmplicon;
struct CurationEvent;
struct BranchCandidate;
}  // namespace branch::wg

namespace branch::gpu::wg_kernels {

/// True when the binary was built with BRANCH_CUDA_ENABLED=ON AND
/// at least one CUDA-capable device is visible at run time.
[[nodiscard]] bool is_gpu_available() noexcept;

/// Phase 3 Stage 1: window-stats screen.
///   - hap_seq:      uppercase ATCG sequence (host)
///   - counter_keys: sorted ascending vector of canonical_kmer_bits
///   - counter_cnts: parallel vector of counts (same length as keys)
///   - k:            kmer length (from TechProfile)
///   - w:            minimizer window length (from TechProfile)
///   - hap_idx:      stamped into output records
///   - z_threshold / rep_entropy / min_density_z: same defaults as CPU
/// Returns true on success; false → caller falls back to CPU.
[[nodiscard]] bool launch_stage1(
    std::uint32_t hap_idx,
    const std::string& hap_seq,
    const std::vector<std::uint64_t>& counter_keys,
    const std::vector<std::uint32_t>& counter_cnts,
    std::size_t k,
    std::size_t w,
    double z_threshold,
    double rep_entropy,
    double min_density_z,
    std::vector<branch::wg::Stage1Window>& out_windows) noexcept;

/// Phase 3 Stage 2: per-read seed-indexed primer scan.
///   - amplicons:  in/out — read_count + cn populated
///   - reads:      flat (id, seq) pairs (id only used for diagnostics)
///   - max_mm:     max primer-mismatch tolerance (= 2)
///   - amp_min/max: PCR-able distance gate
///   - expected_single_copy_count: same semantics as CPU count()
[[nodiscard]] bool launch_stage2(
    std::vector<branch::wg::VpcrAmplicon>& amplicons,
    const std::vector<std::pair<std::string, std::string>>& reads,
    int max_mm,
    int amp_min,
    int amp_max,
    int expected_single_copy_count) noexcept;

/// Phase 1.5: per-position pile-up curation. (Stretch goal — not yet
/// implemented, returns false unconditionally to force CPU fallback.)
[[nodiscard]] bool launch_phase15_curation(
    std::uint32_t /*hap_idx*/,
    const std::string& /*hap_seq*/,
    const std::vector<std::pair<std::string, std::string>>& /*reads*/,
    std::size_t /*k*/,
    std::size_t /*w*/,
    int /*min_coverage_branch*/,
    double /*min_base_agreement*/,
    int /*default_q*/,
    std::vector<branch::wg::CurationEvent>& /*out_events*/,
    std::vector<branch::wg::BranchCandidate>& /*out_branches*/) noexcept;

/// Phase 0: bloom-filtered k-mer counter. (Stretch goal — not yet
/// implemented, returns false to force CPU fallback.)
[[nodiscard]] bool launch_phase0_count(
    const std::vector<std::pair<std::string, std::string>>& /*reads*/,
    std::size_t /*k*/,
    std::vector<std::uint64_t>& /*out_keys*/,
    std::vector<std::uint32_t>& /*out_counts*/) noexcept;

}  // namespace branch::gpu::wg_kernels
