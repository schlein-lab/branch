#include "anchored_phaser.hpp"

#include "../../../align/llmap_align_backend.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

namespace branch::wg::phaser {

namespace {

// Matched-base weight for a hit: identity × aligned query span. A robust,
// always-non-negative proxy for alignment quality used to compare the two
// haplotypes' fits for one read (replaces minimap2's AS:i score).
int hit_weight(const branch::align::llmap_backend::RefHit& h) {
    const int aligned = h.query_end - h.query_start;
    if (aligned <= 0) return 0;
    return static_cast<int>(h.identity * static_cast<float>(aligned));
}

}  // namespace

void AnchoredPhaser::run(
    FastqIndex& reads,
    std::vector<ReadAssignment>& out_assignments,
    const AnchoredPhaserOpts& opts)
{
    auto t0 = std::chrono::steady_clock::now();
    out_assignments.clear();
    out_assignments.resize(reads.n_reads());

    // Initialize all assignments as UNTAGGED (no hits yet).
    for (std::size_t i = 0; i < reads.n_reads(); ++i) {
        out_assignments[i].read_id = static_cast<ReadId>(i);
        out_assignments[i].tag = ReadTag::UNTAGGED;
        out_assignments[i].confidence = 0.0f;
        out_assignments[i].home_node = 0;
    }

    if (opts.ref_hap1_fa_path.empty() || opts.ref_hap2_fa_path.empty()) {
        std::cerr << "[anchored] no haplotype reference FASTAs given; "
                     "emitting UNTAGGED for all reads\n";
        return;
    }

    // -- Build one minimizer index per haplotype (held for the whole pass).
    branch::align::llmap_backend::ReferenceMapper mapper_h1(
        opts.ref_hap1_fa_path, opts.preset);
    branch::align::llmap_backend::ReferenceMapper mapper_h2(
        opts.ref_hap2_fa_path, opts.preset);
    if (!mapper_h1.valid() || !mapper_h2.valid()) {
        std::cerr << "[anchored] failed to index one or both haplotype "
                     "references; emitting UNTAGGED for all reads\n";
        return;
    }

    // Per-read matched-base weight against each haplotype.
    std::vector<int> weight_h1(reads.n_reads(), 0);
    std::vector<int> weight_h2(reads.n_reads(), 0);

    // -- Stream reads in bounded batches; map each batch in parallel.
    std::vector<ReadId>      batch_rids;
    std::vector<std::string> batch_seqs;
    batch_rids.reserve(opts.batch_size);
    batch_seqs.reserve(opts.batch_size);
    std::size_t n_mapped = 0;

    auto flush_batch = [&]() {
        if (batch_seqs.empty()) return;
        auto hits1 = mapper_h1.best_batch(batch_seqs, opts.n_threads);
        auto hits2 = mapper_h2.best_batch(batch_seqs, opts.n_threads);
        for (std::size_t i = 0; i < batch_rids.size(); ++i) {
            const ReadId rid = batch_rids[i];
            if (hits1[i]) {
                const int w = hit_weight(*hits1[i]);
                if (w >= opts.min_align_score) weight_h1[rid] = w;
            }
            if (hits2[i]) {
                const int w = hit_weight(*hits2[i]);
                if (w >= opts.min_align_score) weight_h2[rid] = w;
            }
        }
        n_mapped += batch_seqs.size();
        batch_rids.clear();
        batch_seqs.clear();
    };

    reads.for_each([&](ReadId rid, const std::string& seq) {
        batch_rids.push_back(rid);
        batch_seqs.push_back(seq);  // copy: for_each reuses its buffer
        if (batch_seqs.size() >= opts.batch_size) flush_batch();
    });
    flush_batch();  // trailing partial batch

    // -- Per-read decision.
    for (std::size_t rid = 0; rid < reads.n_reads(); ++rid) {
        auto& a = out_assignments[rid];
        const int s1 = weight_h1[rid];
        const int s2 = weight_h2[rid];
        const int total = s1 + s2;
        if (total == 0) {
            // No hits anywhere → genuinely unmappable (novel SV, contam, etc.)
            a.tag = ReadTag::UNCERTAIN;
            a.confidence = 0.0f;
            continue;
        }
        const double frac_h1 = static_cast<double>(s1) / static_cast<double>(total);
        const double frac_h2 = static_cast<double>(s2) / static_cast<double>(total);
        if (std::abs(frac_h1 - frac_h2) < opts.hap_call_margin) {
            // Both haplotypes equally well aligned → SHARED.
            a.tag = ReadTag::SHARED;
            a.confidence = static_cast<float>(1.0 - std::abs(frac_h1 - frac_h2));
        } else if (frac_h1 > frac_h2) {
            a.tag = ReadTag::H1_ONLY;
            a.confidence = static_cast<float>(frac_h1);
        } else {
            a.tag = ReadTag::H2_ONLY;
            a.confidence = static_cast<float>(frac_h2);
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[anchored] mapped " << n_mapped << " reads vs 2 haplotypes in "
              << std::chrono::duration<double>(t1 - t0).count() << "s\n";
}

}  // namespace branch::wg::phaser
