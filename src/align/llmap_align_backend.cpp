// BRANCH — LLmap alignment backend (implementation)
//
// Compiled as C++23 (LLmap requires cxx_std_23) and the ONLY translation
// unit in BRANCH that includes LLmap headers / links its static libs.
// See llmap_align_backend.hpp for the rationale behind this single seam.

#include "llmap_align_backend.hpp"

#include <algorithm>
#include <cstdio>

// LLmap classical alignment engine (gap-affine WFA2/Gotoh + seed-chain-extend).
#include "classical/wfa2_aligner.h"
#include "classical/minimizer_index.h"
#include "classical/classical_pipeline.h"
#include "core/thread_pool.h"

#include <zlib.h>

namespace branch::align::llmap_backend {

namespace {

// --- preset → minimizer/identity knobs -------------------------------------
struct PresetParams {
    std::uint8_t k;
    std::uint8_t w;
    float        min_identity;
};

PresetParams resolve_preset(const std::string& preset) {
    if (preset == "map-ont") {
        return {15, 10, 0.70f};   // noisier ONT R10
    }
    // "asm5", "map-hifi", and anything unrecognized → high-identity defaults.
    return {19, 19, 0.85f};
}

// --- minimal gz-aware FASTA reader -----------------------------------------
// We read FASTA here (rather than pull in another LLmap lib) so the backend
// links only llmap_core/annot/classical. Handles plain or gzipped input.
bool read_fasta(const std::string& path,
                std::vector<std::string>& names,
                std::vector<std::string>& seqs) {
    gzFile fp = gzopen(path.c_str(), "rb");
    if (!fp) return false;

    constexpr int kBuf = 1 << 16;
    std::vector<char> buf(kBuf);
    std::string cur_name, cur_seq;
    auto flush = [&]() {
        if (!cur_name.empty() || !cur_seq.empty()) {
            names.push_back(cur_name);
            seqs.push_back(cur_seq);
        }
        cur_name.clear();
        cur_seq.clear();
    };

    std::string line;
    int n;
    while ((n = gzread(fp, buf.data(), kBuf)) > 0) {
        for (int i = 0; i < n; ++i) {
            char c = buf[i];
            if (c == '\n') {
                if (!line.empty()) {
                    if (line[0] == '>') {
                        flush();
                        // Header up to first whitespace is the name.
                        std::size_t sp = line.find_first_of(" \t");
                        cur_name = line.substr(1, sp == std::string::npos
                                                      ? std::string::npos
                                                      : sp - 1);
                    } else {
                        cur_seq += line;
                    }
                }
                line.clear();
            } else if (c != '\r') {
                line += c;
            }
        }
    }
    // Trailing line with no newline.
    if (!line.empty()) {
        if (line[0] == '>') {
            flush();
            std::size_t sp = line.find_first_of(" \t");
            cur_name = line.substr(1, sp == std::string::npos
                                          ? std::string::npos
                                          : sp - 1);
        } else {
            cur_seq += line;
        }
    }
    flush();
    gzclose(fp);
    return !seqs.empty();
}

}  // namespace

// ===========================================================================
// Pairwise alignment
// ===========================================================================

std::optional<PairwiseAlignment>
align_pair(std::string_view query,
           std::string_view target,
           const PairwiseOptions& opts) {
    if (query.empty() || target.empty()) return std::nullopt;
    if (query.size() > opts.max_query_len || target.size() > opts.max_ref_len) {
        return std::nullopt;
    }

    llmap::classical::WFA2Config cfg;
    cfg.match_score        = 0;
    cfg.mismatch_penalty   = opts.mismatch_penalty;
    cfg.gap_open           = opts.gap_open;
    cfg.gap_extend         = opts.gap_extend;
    cfg.end_to_end         = opts.global;
    cfg.max_pattern_length = opts.max_query_len;
    cfg.max_text_length    = opts.max_ref_len;

    llmap::classical::WFA2Aligner aln{cfg};
    auto res = aln.Align(query, target);
    if (!res) return std::nullopt;

    PairwiseAlignment out;
    out.cigar          = res->CigarString();
    out.score          = res->score;
    out.identity       = res->identity;
    out.query_start    = res->query_start;
    out.query_end      = res->query_end;
    out.ref_start      = res->ref_start;
    out.ref_end        = res->ref_end;
    out.num_matches    = res->num_matches;
    out.num_mismatches = res->num_mismatches;
    out.num_insertions = res->num_insertions;
    out.num_deletions  = res->num_deletions;
    // Edit distance = substitutions + indel bases. Independent of the
    // scoring scheme; derived purely from the alignment's op counts.
    out.edit_distance = static_cast<int>(res->num_mismatches +
                                         res->num_insertions +
                                         res->num_deletions);
    return out;
}

int pairwise_edit_distance(std::string_view query,
                           std::string_view target,
                           const PairwiseOptions& opts,
                           std::string* cigar_out) {
    if (query.empty() || target.empty()) {
        if (cigar_out) *cigar_out = "";
        return static_cast<int>(std::max(query.size(), target.size()));
    }
    if (query == target) {
        if (cigar_out) *cigar_out = std::to_string(query.size()) + "M";
        return 0;
    }

    auto aln = align_pair(query, target, opts);
    if (aln) {
        if (cigar_out) *cigar_out = aln->cigar;
        return aln->edit_distance;
    }

    // Fallback when LLmap aborts (e.g. length guard): positional difference
    // count + length delta. Never throws, so callers always get a number.
    int diff = 0;
    const std::size_t m = std::min(query.size(), target.size());
    for (std::size_t i = 0; i < m; ++i) {
        if (query[i] != target[i]) ++diff;
    }
    diff += static_cast<int>(query.size() > target.size()
                                 ? query.size() - target.size()
                                 : target.size() - query.size());
    if (cigar_out) *cigar_out = "";
    return diff;
}

// ===========================================================================
// Reference mapping
// ===========================================================================

struct ReferenceMapper::Impl {
    llmap::classical::ClassicalPipeline pipeline;
    bool ok = false;

    Impl(const std::vector<std::string>& names,
         const std::vector<std::string>& seqs,
         const std::string& preset)
        : pipeline(make_config(preset)) {
        if (names.size() != seqs.size() || seqs.empty()) return;

        llmap::classical::MinimizerIndex::Builder builder{
            make_config(preset).minimizer_config};
        for (std::size_t i = 0; i < seqs.size(); ++i) {
            builder.AddSequence(names[i], seqs[i]);
        }
        pipeline.SetIndex(builder.Build());
        // Owning copy: the pipeline needs the reference bases for WFA2
        // extension of each chain.
        pipeline.SetReferenceSequences(seqs);
        ok = pipeline.HasIndex();
    }

    static llmap::classical::ClassicalPipelineConfig
    make_config(const std::string& preset) {
        const PresetParams p = resolve_preset(preset);
        llmap::classical::ClassicalPipelineConfig cfg;
        cfg.minimizer_config.k = p.k;
        cfg.minimizer_config.w = p.w;
        cfg.min_identity       = p.min_identity;
        cfg.report_secondary   = false;  // primary-only, minimap2-style
        cfg.num_threads        = 1;      // caller parallelizes across reads
        return cfg;
    }

    static RefHit to_hit(const llmap::classical::ClassicalAlignment& a) {
        RefHit h;
        h.ref_name    = a.ref_name;
        h.ref_id      = a.ref_id;
        h.ref_start   = a.ref_start;
        h.ref_end     = a.ref_end;
        h.query_start = a.query_start;
        h.query_end   = a.query_end;
        h.is_forward  = a.is_forward;
        h.is_primary  = a.is_primary;
        h.cigar       = a.CigarString();
        h.score       = a.score;
        h.identity    = a.identity;
        h.mapq        = a.mapq;
        // Mismatch count within the aligned query span. identity is the
        // match fraction, so (1 - identity) × span ≈ substitutions + indels.
        // This is the SNV-aware discriminator used by the anchored phaser:
        // the *difference* in mismatches between two near-identical haplotype
        // references isolates exactly the bases at their differing sites.
        h.aligned_len = a.query_end - a.query_start;
        if (h.aligned_len < 0) h.aligned_len = 0;
        h.mismatches  = static_cast<int>(
            (1.0f - a.identity) * static_cast<float>(h.aligned_len) + 0.5f);
        return h;
    }
};

ReferenceMapper::ReferenceMapper(const std::vector<std::string>& ref_names,
                                 const std::vector<std::string>& ref_seqs,
                                 const std::string& preset)
    : impl_(std::make_unique<Impl>(ref_names, ref_seqs, preset)) {}

ReferenceMapper::ReferenceMapper(const std::string& ref_fasta_path,
                                 const std::string& preset) {
    std::vector<std::string> names, seqs;
    if (!read_fasta(ref_fasta_path, names, seqs)) {
        // Construct an invalid mapper rather than throw; callers check valid().
        impl_ = std::make_unique<Impl>(std::vector<std::string>{},
                                       std::vector<std::string>{}, preset);
        return;
    }
    impl_ = std::make_unique<Impl>(names, seqs, preset);
}

ReferenceMapper::~ReferenceMapper() = default;
ReferenceMapper::ReferenceMapper(ReferenceMapper&&) noexcept = default;
ReferenceMapper& ReferenceMapper::operator=(ReferenceMapper&&) noexcept = default;

bool ReferenceMapper::valid() const noexcept {
    return impl_ && impl_->ok;
}

std::optional<RefHit> ReferenceMapper::best(std::string_view query) const {
    if (!valid() || query.empty()) return std::nullopt;
    auto rar = impl_->pipeline.AlignRead("q", query);
    if (!rar.HasAlignment()) return std::nullopt;
    return Impl::to_hit(rar.alignments.front());
}

std::vector<RefHit> ReferenceMapper::all(std::string_view query) const {
    std::vector<RefHit> out;
    if (!valid() || query.empty()) return out;
    auto rar = impl_->pipeline.AlignRead("q", query);
    out.reserve(rar.alignments.size());
    for (const auto& a : rar.alignments) out.push_back(Impl::to_hit(a));
    return out;
}

std::vector<std::optional<RefHit>>
ReferenceMapper::best_batch(const std::vector<std::string>& seqs,
                            int n_threads) const {
    std::vector<std::optional<RefHit>> out(seqs.size());
    if (!valid() || seqs.empty()) return out;

    // ClassicalPipeline::AlignReadsParallel wants parallel name/seq spans.
    // Names are unused downstream here, so cheap sequential indices suffice.
    std::vector<std::string> names(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); ++i) names[i] = std::to_string(i);

    llmap::core::ThreadPool pool(n_threads > 0
                                     ? static_cast<std::size_t>(n_threads)
                                     : 1);
    auto rars = impl_->pipeline.AlignReadsParallel(names, seqs, pool);
    for (std::size_t i = 0; i < rars.size() && i < out.size(); ++i) {
        if (rars[i].HasAlignment()) {
            out[i] = Impl::to_hit(rars[i].alignments.front());
        }
    }
    return out;
}

}  // namespace branch::align::llmap_backend
