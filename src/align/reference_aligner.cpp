// BRANCH — Reference Aligner Implementation
// Query-vs-reference alignment via the LLmap classical engine
// (seed → chain → WFA2 extension). Replaces the former minimap2 mm_map.

#include "reference_aligner.hpp"
#include "llmap_align_backend.hpp"

#include <stdexcept>
#include <utility>

namespace branch::align {

struct ReferenceAligner::Impl {
    // "asm5"-equivalent preset: high-identity assembly/consensus vs
    // reference (k=19, w=19, identity floor 0.85). The former code used
    // minimap2's asm5 for exactly this regime.
    explicit Impl(const std::string& ref_path)
        : mapper(ref_path, "asm5") {}

    llmap_backend::ReferenceMapper mapper;
};

ReferenceAligner::ReferenceAligner(const std::string& ref_path)
    : impl_(std::make_unique<Impl>(ref_path)) {
    if (!impl_->mapper.valid()) {
        throw std::runtime_error("Failed to build reference index from: " +
                                 ref_path);
    }
}

ReferenceAligner::~ReferenceAligner() = default;
ReferenceAligner::ReferenceAligner(ReferenceAligner&&) noexcept = default;
ReferenceAligner& ReferenceAligner::operator=(ReferenceAligner&&) noexcept = default;

bool ReferenceAligner::is_valid() const noexcept {
    return impl_ && impl_->mapper.valid();
}

std::optional<AlignmentResult> ReferenceAligner::align(std::string_view query) const {
    if (!is_valid() || query.empty()) {
        return std::nullopt;
    }

    auto hit = impl_->mapper.best(query);
    if (!hit) {
        return std::nullopt;
    }

    AlignmentResult ar;
    ar.score       = hit->score;
    ar.identity    = hit->identity;
    ar.cigar       = std::move(hit->cigar);
    ar.ref_start   = hit->ref_start;
    ar.ref_end     = hit->ref_end;
    ar.query_start = hit->query_start;
    ar.query_end   = hit->query_end;
    return ar;
}

float ReferenceAligner::align_identity(std::string_view query) const {
    auto result = align(query);
    return result ? result->identity : 0.0f;
}

}  // namespace branch::align
