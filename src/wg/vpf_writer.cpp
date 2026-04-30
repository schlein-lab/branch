// BRANCH v0.5 — extended-GFA writer implementation.

#include "wg/vpf_writer.hpp"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>


namespace branch::wg {

VpfWriter::VpfWriter(const std::string& out_prefix) : prefix_(out_prefix) {}

namespace {


void write_fasta_record(std::ostream& f, const std::string& name,
                        const std::string& seq) {
    f << '>' << name << '\n';
    for (std::size_t i = 0; i < seq.size(); i += 60) {
        f << seq.substr(i, 60) << '\n';
    }
}

}  // namespace

void VpfWriter::write(const WgOutputs& out) {
    // 1. <prefix>.fa — FASTA with all sequences (haplotypes + branches + orphans).
    {
        std::string path = prefix_ + ".fa";
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot open " + path);
        for (std::size_t i = 0; i < out.haplotype_seqs.size(); ++i) {
            write_fasta_record(f, "hap" + std::to_string(i + 1),
                               out.haplotype_seqs[i]);
        }
        for (const auto& b : out.branches) {
            write_fasta_record(f, b.branch_name, b.seq);
        }
        for (const auto& c : out.orphan_clusters) {
            if (!c.consensus_hint.empty()) {
                write_fasta_record(
                    f, "orphan_oc" + std::to_string(c.cluster_id),
                    c.consensus_hint);
            }
        }
    }

    // 2. <prefix>.bed — BED with branch + amplicon positions.
    {
        std::string path = prefix_ + ".bed";
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot open " + path);
        for (const auto& b : out.branches) {
            f << "hap" << (b.anchor_hap + 1) << '\t'
              << b.anchor_pos << '\t'
              << (b.anchor_pos + static_cast<std::uint32_t>(b.seq.size())) << '\t'
              << b.branch_name << '\t'
              << b.confidence << "\t.\n";
        }
        for (const auto& a : out.vpcr_amplicons) {
            f << "hap" << (a.hap_idx + 1) << '\t'
              << a.start_bp << '\t'
              << (a.start_bp + a.expected_size_bp) << '\t'
              << a.amplicon_name << '\t'
              << a.cn_estimate << "\t.\t"
              << a.read_count << '\t'
              << a.ci_low << '\t' << a.ci_high << '\t'
              << static_cast<int>(a.region_flags) << '\n';
        }
    }

    // 3. <prefix>.gfa.zst — extended GFA with all the v0.5 metadata.
    {
        std::ostringstream ss;
        // Header
        ss << "##VPF v5.0\n";
        ss << "##assembler " << out.meta.assembler_version << '\n';
        ss << "##sample " << out.meta.sample << '\n';
        ss << "##read_tech " << out.meta.read_tech << '\n';
        ss << "##haplotype_count " << out.haplotype_seqs.size() << '\n';
        ss << "##input_reads " << out.meta.input_reads << '\n';
        ss << "##input_bases " << out.meta.input_bases << '\n';
        ss << "##detected_coverage " << out.meta.detected_coverage << '\n';
        ss << "##detected_het_coverage " << out.meta.detected_het_coverage << '\n';

        // #vph lines per haplotype
        for (std::size_t h = 0; h < out.haplotype_tiles.size(); ++h) {
            const auto& tiles = out.haplotype_tiles[h];
            ss << "#vph hap=hap" << (h + 1)
               << " length=" << out.haplotype_seqs[h].size()
               << " master_tiles=" << tiles.size() << '\n';
            for (std::size_t t = 0; t < tiles.size(); ++t) {
                const auto& tt = tiles[t];
                ss << "#vph hap=hap" << (h + 1)
                   << " tile=" << t
                   << " start=" << tt.start_bp
                   << " len=" << tt.length_bp
                   << " master=read_" << tt.master_read_idx
                   << " qmean=" << (tt.q_mean_x10 / 10.0) << '\n';
            }
        }

        // #cur lines (curation log)
        for (const auto& c : out.curation_events) {
            ss << "#cur hap=hap" << (c.hap_idx + 1)
               << " pos=" << c.pos_in_hap
               << " from=" << c.original_base
               << " to=" << c.curated_base
               << " n_support=" << c.n_supporting_reads << '\n';
        }

        // #vpcr lines
        for (const auto& a : out.vpcr_amplicons) {
            ss << "#vpcr amplicon=" << a.amplicon_name
               << " fwd=" << a.fwd_primer
               << " rev=" << a.rev_primer
               << " count=" << a.read_count
               << " cn=" << a.cn_estimate
               << " ci=[" << a.ci_low << ',' << a.ci_high << "]\n";
        }

        // #disp lines
        for (const auto& [id, info] : out.read_disposition) {
            ss << "#disp read=" << id << ' ' << info << '\n';
        }

        // #orphan lines
        for (const auto& c : out.orphan_clusters) {
            ss << "#orphan cluster=oc_" << c.cluster_id
               << " size=" << c.member_read_ids.size()
               << " ref=" << c.ref_attempt_chrom
               << " best_score=" << c.ref_attempt_score << '\n';
        }

        // #log lines
        for (const auto& [phase, secs] : out.phase_log) {
            ss << "#log phase=" << phase << " wall_seconds=" << secs << '\n';
        }

        // Standard GFA records
        for (std::size_t i = 0; i < out.haplotype_seqs.size(); ++i) {
            ss << "S\thap" << (i + 1) << '\t' << out.haplotype_seqs[i]
               << "\tLN:i:" << out.haplotype_seqs[i].size() << '\n';
        }
        for (const auto& b : out.branches) {
            ss << "S\t" << b.branch_name << '\t' << b.seq
               << "\tLN:i:" << b.seq.size()
               << "\tCF:f:" << b.confidence
               << "\tHP:i:" << b.anchor_hap << '\n';
            ss << "L\thap" << (b.anchor_hap + 1)
               << "\t+\t" << b.branch_name << "\t+\t0M"
               << "\tAT:Z:hap" << (b.anchor_hap + 1) << ':' << b.anchor_pos << '\n';
        }
        for (std::size_t i = 0; i < out.haplotype_seqs.size(); ++i) {
            ss << "P\thap" << (i + 1) << "\thap" << (i + 1) << "+\t*\n";
        }

        std::string gfa_text = ss.str();
        // zstd compression deferred (system zstd dev pkg not present)
        const std::string& compressed = gfa_text;

        std::string path = prefix_ + ".gfa";
        std::ofstream f(path, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open " + path);
        f.write(compressed.data(), static_cast<std::streamsize>(compressed.size()));
    }
}

}  // namespace branch::wg
