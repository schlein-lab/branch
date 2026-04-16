// BRANCH — Reference Aligner Implementation
// Uses minimap2/ksw2 for sequence alignment.

#include "reference_aligner.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>

// minimap2 headers (C interface)
extern "C" {
#include <minimap.h>
}
#include <zlib.h>

// kseq.h generates code with -Wconversion warnings — suppress them
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
extern "C" {
#include <kseq.h>
}
KSEQ_INIT(gzFile, gzread)
#pragma GCC diagnostic pop

namespace branch::align {

struct ReferenceAligner::Impl {
    mm_idx_t* index{nullptr};
    mm_mapopt_t map_opt;
    mm_idxopt_t idx_opt;
    std::string ref_path;
    
    Impl() {
        mm_set_opt(nullptr, &idx_opt, &map_opt);
        // Use asm5 preset for high-identity alignments (like reference vs consensus)
        mm_set_opt("asm5", &idx_opt, &map_opt);
        map_opt.flag |= MM_F_CIGAR;  // Generate CIGAR
    }
    
    ~Impl() {
        if (index) {
            mm_idx_destroy(index);
        }
    }
};

ReferenceAligner::ReferenceAligner(const std::string& ref_path) 
    : impl_(std::make_unique<Impl>()) {
    
    impl_->ref_path = ref_path;
    
    // Build index from reference FASTA
    mm_idx_reader_t* reader = mm_idx_reader_open(ref_path.c_str(), &impl_->idx_opt, nullptr);
    if (!reader) {
        throw std::runtime_error("Failed to open reference: " + ref_path);
    }
    
    impl_->index = mm_idx_reader_read(reader, 1);  // Read first index
    mm_idx_reader_close(reader);
    
    if (!impl_->index) {
        throw std::runtime_error("Failed to build index from: " + ref_path);
    }
    
    // Map options post-index
    mm_mapopt_update(&impl_->map_opt, impl_->index);
}

ReferenceAligner::~ReferenceAligner() = default;

ReferenceAligner::ReferenceAligner(ReferenceAligner&&) noexcept = default;
ReferenceAligner& ReferenceAligner::operator=(ReferenceAligner&&) noexcept = default;

bool ReferenceAligner::is_valid() const noexcept {
    return impl_ && impl_->index != nullptr;
}

std::optional<AlignmentResult> ReferenceAligner::align(std::string_view query) const {
    if (!is_valid() || query.empty()) {
        return std::nullopt;
    }
    
    // Prepare query
    mm_tbuf_t* tbuf = mm_tbuf_init();
    if (!tbuf) {
        return std::nullopt;
    }
    
    int n_reg = 0;
    mm_reg1_t* regs = mm_map(
        impl_->index,
        static_cast<int>(query.size()),
        query.data(),
        &n_reg,
        tbuf,
        &impl_->map_opt,
        nullptr  // name (not needed)
    );
    
    std::optional<AlignmentResult> result;
    
    if (n_reg > 0 && regs) {
        // Take best (first) alignment
        const mm_reg1_t& r = regs[0];
        
        AlignmentResult ar;
        ar.score = r.score;
        ar.ref_start = r.rs;
        ar.ref_end = r.re;
        ar.query_start = r.qs;
        ar.query_end = r.qe;
        
        // Calculate identity from CIGAR
        int matches = 0;
        int alignment_len = 0;
        
        if (r.p) {
            std::ostringstream cigar_ss;
            for (uint32_t i = 0; i < r.p->n_cigar; ++i) {
                uint32_t op = r.p->cigar[i] & 0xf;
                uint32_t len = r.p->cigar[i] >> 4;
                
                char op_char = "MIDNSHP=X"[op];
                cigar_ss << len << op_char;
                
                // Count matches and alignment length
                switch (op) {
                    case 0:  // M (match/mismatch)
                    case 7:  // = (match)
                        matches += static_cast<int>(len);
                        alignment_len += static_cast<int>(len);
                        break;
                    case 8:  // X (mismatch)
                        alignment_len += static_cast<int>(len);
                        break;
                    case 1:  // I (insertion)
                    case 2:  // D (deletion)
                        alignment_len += static_cast<int>(len);
                        break;
                    default:
                        break;
                }
            }
            ar.cigar = cigar_ss.str();
            
            // For M operations, we need the actual match count from blen/mlen
            // minimap2 stores block length and match length
            if (r.blen > 0) {
                ar.identity = static_cast<float>(r.mlen) / static_cast<float>(r.blen);
            } else if (alignment_len > 0) {
                ar.identity = static_cast<float>(matches) / static_cast<float>(alignment_len);
            }
        } else {
            // No CIGAR, estimate from block info
            if (r.blen > 0) {
                ar.identity = static_cast<float>(r.mlen) / static_cast<float>(r.blen);
            }
        }
        
        result = ar;
        
        // Free CIGAR memory
        for (int i = 0; i < n_reg; ++i) {
            if (regs[i].p) free(regs[i].p);
        }
    }
    
    free(regs);
    mm_tbuf_destroy(tbuf);
    
    return result;
}

float ReferenceAligner::align_identity(std::string_view query) const {
    auto result = align(query);
    return result ? result->identity : 0.0f;
}

}  // namespace branch::align
