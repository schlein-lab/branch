#include "anchored_phaser.hpp"

#include "../../../align/llmap_align_backend.hpp"

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace branch::wg::phaser {

void AnchoredPhaser::run(
    FastqIndex& reads,
    std::vector<ReadAssignment>& out_assignments,
    const AnchoredPhaserOpts& opts)
{
    auto t0 = std::chrono::steady_clock::now();
    out_assignments.clear();
    out_assignments.resize(reads.n_reads());

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

    branch::align::llmap_backend::ReferenceMapper mapper_h1(
        opts.ref_hap1_fa_path, opts.preset);
    branch::align::llmap_backend::ReferenceMapper mapper_h2(
        opts.ref_hap2_fa_path, opts.preset);
    if (!mapper_h1.valid() || !mapper_h2.valid()) {
        std::cerr << "[anchored] failed to index one or both haplotype "
                     "references; emitting UNTAGGED for all reads\n";
        return;
    }

    // Per-read mismatch counts to each haplotype (-1 = did not map).
    constexpr int kUnmapped = -1;
    std::vector<int> mm_h1(reads.n_reads(), kUnmapped);
    std::vector<int> mm_h2(reads.n_reads(), kUnmapped);

    std::vector<ReadId>      batch_rids;
    std::vector<std::string> batch_seqs;
    batch_rids.reserve(opts.batch_size);
    batch_seqs.reserve(opts.batch_size);
    std::size_t n_mapped = 0;

    auto record = [&](const std::optional<branch::align::llmap_backend::RefHit>& hit,
                      std::vector<int>& mm, ReadId rid) {
        if (hit && hit->aligned_len >= opts.min_aligned_bp) {
            mm[rid] = hit->mismatches;
        }
    };

    auto flush_batch = [&]() {
        if (batch_seqs.empty()) return;
        auto hits1 = mapper_h1.best_batch(batch_seqs, opts.n_threads);
        auto hits2 = mapper_h2.best_batch(batch_seqs, opts.n_threads);
        for (std::size_t i = 0; i < batch_rids.size(); ++i) {
            record(hits1[i], mm_h1, batch_rids[i]);
            record(hits2[i], mm_h2, batch_rids[i]);
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
    flush_batch();

    // -- Per-read decision: SNV-aware discrimination by mismatch difference.
    //
    // A read from its true haplotype matches that hap's allele at every
    // differing (het) site and mismatches the other hap there; sequencing
    // errors hit both haps about equally and cancel. So the SIGN of
    // (mm_other - mm_true) isolates the het-site signal, even when the haps
    // are >99.9% identical and the total alignment weights are nearly equal.
    for (std::size_t rid = 0; rid < reads.n_reads(); ++rid) {
        auto& a = out_assignments[rid];
        const int m1 = mm_h1[rid];
        const int m2 = mm_h2[rid];
        const bool map1 = (m1 != kUnmapped);
        const bool map2 = (m2 != kUnmapped);

        if (!map1 && !map2) {
            a.tag = ReadTag::UNCERTAIN;     // maps to neither hap
            a.confidence = 0.0f;
            continue;
        }
        if (map1 != map2) {
            // Aligns to exactly one hap → unambiguous assignment.
            a.tag = map1 ? ReadTag::H1_ONLY : ReadTag::H2_ONLY;
            a.confidence = 1.0f;
            continue;
        }

        const int d = m2 - m1;            // >0 ⇒ fewer mismatches to hap1
        const int ad = std::abs(d);
        if (ad < opts.min_discriminating_sites) {
            // No discriminating mismatch difference → homozygous span.
            a.tag = ReadTag::SHARED;
            a.confidence = 1.0f;
        } else {
            a.tag = (d > 0) ? ReadTag::H1_ONLY : ReadTag::H2_ONLY;
            // Confidence grows with the net discriminating margin (bounded,
            // monotone): 1 site → 0.5, 2 → 0.67, 4 → 0.8, → 1.0.
            a.confidence = static_cast<float>(ad) / static_cast<float>(ad + 1);
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cerr << "[anchored] mapped " << n_mapped
              << " reads vs 2 haplotypes (SNV-aware) in "
              << std::chrono::duration<double>(t1 - t0).count() << "s\n";
}

}  // namespace branch::wg::phaser
