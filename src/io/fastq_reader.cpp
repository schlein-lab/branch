#include "io/fastq_reader.hpp"

namespace branch::io {

FastqReader::FastqReader(const std::string& path)
    : owned_file_(std::make_unique<std::ifstream>(path)),
      stream_(owned_file_.get()) {}

FastqReader::FastqReader(std::istream& stream)
    : owned_file_(), stream_(&stream) {}

bool FastqReader::ok() const noexcept {
    return stream_ != nullptr && static_cast<bool>(*stream_);
}

bool FastqReader::next_record(FastqRecord& out) {
    if (!ok()) return false;

    std::string header;
    if (!std::getline(*stream_, header)) return false;
    while (!header.empty() && (header.back() == '\r')) header.pop_back();
    if (header.empty() || header[0] != '@') return false;

    std::string seq, plus, qual;
    if (!std::getline(*stream_, seq)) return false;
    if (!std::getline(*stream_, plus)) return false;
    if (!std::getline(*stream_, qual)) return false;

    while (!seq.empty() && seq.back() == '\r') seq.pop_back();
    while (!plus.empty() && plus.back() == '\r') plus.pop_back();
    while (!qual.empty() && qual.back() == '\r') qual.pop_back();

    if (plus.empty() || plus[0] != '+') return false;
    if (seq.size() != qual.size()) return false;

    // Strip '@' and anything after first whitespace for the name.
    std::string name = header.substr(1);
    const auto ws = name.find_first_of(" \t");
    if (ws != std::string::npos) name.resize(ws);

    out.name = std::move(name);
    out.sequence = std::move(seq);
    out.qualities = std::move(qual);
    ++n_read_;
    return true;
}

}  // namespace branch::io
