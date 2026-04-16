#include "io/bam_reader.hpp"

#include <array>
#include <cstdint>
#include <cstring>
#include <utility>

#include <htslib/hts.h>
#include <htslib/sam.h>

namespace branch::io {

namespace {

// htslib's 4-bit packed-nibble encoding. seq_nt16_str is a global in
// htslib (char[16] -> {"=ACMGRSVTWYHKDBN"}), but we use only the four
// canonical bases and collapse everything else to 'N' to simplify
// reverse-complementation.
constexpr char kNibbleToBase[16] = {
    'N', 'A', 'C', 'N',  // 0=?, 1=A, 2=C, 3=M
    'G', 'N', 'N', 'N',  // 4=G, 5=R, 6=S, 7=V
    'T', 'N', 'N', 'N',  // 8=T, 9=W, 10=Y, 11=H
    'N', 'N', 'N', 'N',  // 12=K, 13=D, 14=B, 15=N
};

// ASCII -> complement table covering A/C/G/T/N in both cases; everything
// else maps to 'N' so downstream code never sees a non-ACGTN character.
constexpr std::array<char, 256> make_complement_table() {
    std::array<char, 256> t{};
    for (auto& c : t) c = 'N';
    t[static_cast<unsigned char>('A')] = 'T';
    t[static_cast<unsigned char>('T')] = 'A';
    t[static_cast<unsigned char>('C')] = 'G';
    t[static_cast<unsigned char>('G')] = 'C';
    t[static_cast<unsigned char>('N')] = 'N';
    t[static_cast<unsigned char>('a')] = 'T';
    t[static_cast<unsigned char>('t')] = 'A';
    t[static_cast<unsigned char>('c')] = 'G';
    t[static_cast<unsigned char>('g')] = 'C';
    t[static_cast<unsigned char>('n')] = 'N';
    return t;
}
constexpr auto kComplement = make_complement_table();

// In-place reverse-complement. Safe for both even and odd lengths.
void reverse_complement_inplace(std::string& s) {
    const std::size_t n = s.size();
    for (std::size_t i = 0, j = n == 0 ? 0 : n - 1; i < j; ++i, --j) {
        const char a = kComplement[static_cast<unsigned char>(s[i])];
        const char b = kComplement[static_cast<unsigned char>(s[j])];
        s[i] = b;
        s[j] = a;
    }
    if (n % 2 == 1) {
        const std::size_t mid = n / 2;
        s[mid] = kComplement[static_cast<unsigned char>(s[mid])];
    }
}

}  // namespace

BamReader::~BamReader() { close(); }

BamReader::BamReader(BamReader&& other) noexcept
    : file_(other.file_),
      header_(other.header_),
      record_(other.record_),
      n_read_(other.n_read_) {
    other.file_ = nullptr;
    other.header_ = nullptr;
    other.record_ = nullptr;
    other.n_read_ = 0;
}

BamReader& BamReader::operator=(BamReader&& other) noexcept {
    if (this != &other) {
        close();
        file_ = other.file_;
        header_ = other.header_;
        record_ = other.record_;
        n_read_ = other.n_read_;
        other.file_ = nullptr;
        other.header_ = nullptr;
        other.record_ = nullptr;
        other.n_read_ = 0;
    }
    return *this;
}

bool BamReader::open(const std::string& path) {
    close();
    htsFile* f = sam_open(path.c_str(), "r");
    if (f == nullptr) return false;
    sam_hdr_t* h = sam_hdr_read(f);
    if (h == nullptr) {
        sam_close(f);
        return false;
    }
    bam1_t* b = bam_init1();
    if (b == nullptr) {
        sam_hdr_destroy(h);
        sam_close(f);
        return false;
    }
    file_ = f;
    header_ = h;
    record_ = b;
    return true;
}

bool BamReader::next(std::string& name, std::string& seq) {
    if (file_ == nullptr || header_ == nullptr || record_ == nullptr) {
        return false;
    }
    auto* f = static_cast<htsFile*>(file_);
    auto* h = static_cast<sam_hdr_t*>(header_);
    auto* b = static_cast<bam1_t*>(record_);

    for (;;) {
        const int rc = sam_read1(f, h, b);
        if (rc < 0) return false;  // EOF (-1) or decode error (<-1).

        // Skip records that would cause a read to be yielded more than once.
        const std::uint16_t flag = b->core.flag;
        if ((flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) != 0) {
            continue;
        }

        // Name.
        const char* qn = bam_get_qname(b);
        name.assign(qn == nullptr ? "" : qn);

        // Sequence: 4-bit nibble-packed, two bases per byte.
        const std::int32_t l_qseq = b->core.l_qseq;
        seq.resize(l_qseq < 0 ? 0 : static_cast<std::size_t>(l_qseq));
        if (l_qseq > 0) {
            const std::uint8_t* packed = bam_get_seq(b);
            for (std::int32_t i = 0; i < l_qseq; ++i) {
                const std::uint8_t nib = bam_seqi(packed, i);
                seq[static_cast<std::size_t>(i)] = kNibbleToBase[nib & 0x0F];
            }
        }

        // BAM stores the sequence reference-forward. Restore original
        // sequencing orientation for reverse-strand reads so the minimizer
        // sketcher sees the genuine read.
        if ((flag & BAM_FREVERSE) != 0) {
            reverse_complement_inplace(seq);
        }

        ++n_read_;
        return true;
    }
}

void BamReader::close() {
    if (record_ != nullptr) {
        bam_destroy1(static_cast<bam1_t*>(record_));
        record_ = nullptr;
    }
    if (header_ != nullptr) {
        sam_hdr_destroy(static_cast<sam_hdr_t*>(header_));
        header_ = nullptr;
    }
    if (file_ != nullptr) {
        sam_close(static_cast<htsFile*>(file_));
        file_ = nullptr;
    }
}

}  // namespace branch::io
