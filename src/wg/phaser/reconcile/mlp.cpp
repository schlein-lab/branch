#include "mlp.hpp"

#include <cmath>
#include <fstream>
#include <sstream>

namespace branch::wg::phaser {

namespace {

bool parse_activation(const std::string& s, MLP::Activation& out) {
    if (s == "identity") { out = MLP::Activation::Identity; return true; }
    if (s == "relu")     { out = MLP::Activation::ReLU;     return true; }
    if (s == "softmax")  { out = MLP::Activation::Softmax;  return true; }
    return false;
}

// Read the next whitespace-delimited token, skipping '#' comment lines and
// blank lines. Returns false at EOF.
bool next_token(std::istream& in, std::string& tok) {
    while (in >> tok) {
        if (!tok.empty() && tok[0] == '#') {
            std::string discard;
            std::getline(in, discard);  // skip rest of comment line
            continue;
        }
        return true;
    }
    return false;
}

bool next_float(std::istream& in, float& v) {
    std::string tok;
    if (!next_token(in, tok)) return false;
    try {
        v = std::stof(tok);
    } catch (...) {
        return false;
    }
    return true;
}

bool next_uint(std::istream& in, std::uint32_t& v) {
    std::string tok;
    if (!next_token(in, tok)) return false;
    try {
        v = static_cast<std::uint32_t>(std::stoul(tok));
    } catch (...) {
        return false;
    }
    return true;
}

}  // namespace

bool MLP::load(const std::string& path, std::string* err) {
    layers_.clear();
    loaded_ = false;

    std::ifstream f(path);
    if (!f) {
        if (err) *err = "cannot open weights file: " + path;
        return false;
    }

    // Header: "BRANCH_MLP v1".
    std::string magic, version;
    if (!next_token(f, magic) || magic != "BRANCH_MLP" ||
        !next_token(f, version) || version != "v1") {
        if (err) *err = "bad magic/version (expected 'BRANCH_MLP v1')";
        return false;
    }

    std::uint32_t n_layers = 0;
    if (!next_uint(f, n_layers) || n_layers == 0) {
        if (err) *err = "missing or zero layer count";
        return false;
    }

    std::vector<Layer> layers;
    layers.reserve(n_layers);
    for (std::uint32_t li = 0; li < n_layers; ++li) {
        Layer L;
        std::string act_str;
        if (!next_uint(f, L.in_dim) || !next_uint(f, L.out_dim) ||
            !next_token(f, act_str) || !parse_activation(act_str, L.act)) {
            if (err) *err = "bad layer header at layer " + std::to_string(li);
            return false;
        }
        if (L.in_dim == 0 || L.out_dim == 0) {
            if (err) *err = "zero dimension at layer " + std::to_string(li);
            return false;
        }
        // Successive layers must chain.
        if (li > 0 && L.in_dim != layers.back().out_dim) {
            if (err) *err = "layer " + std::to_string(li) +
                            " in_dim does not match previous out_dim";
            return false;
        }
        const std::size_t n_w = static_cast<std::size_t>(L.in_dim) * L.out_dim;
        L.weights.resize(n_w);
        for (std::size_t i = 0; i < n_w; ++i) {
            if (!next_float(f, L.weights[i])) {
                if (err) *err = "truncated weights at layer " + std::to_string(li);
                return false;
            }
        }
        L.bias.resize(L.out_dim);
        for (std::uint32_t i = 0; i < L.out_dim; ++i) {
            if (!next_float(f, L.bias[i])) {
                if (err) *err = "truncated bias at layer " + std::to_string(li);
                return false;
            }
        }
        layers.push_back(std::move(L));
    }

    layers_ = std::move(layers);
    loaded_ = true;
    return true;
}

std::vector<float> MLP::forward(const std::vector<float>& input) const {
    if (!loaded_ || layers_.empty()) return {};
    if (input.size() != layers_.front().in_dim) return {};

    std::vector<float> x = input;
    for (const auto& L : layers_) {
        std::vector<float> y(L.out_dim, 0.0f);
        for (std::uint32_t o = 0; o < L.out_dim; ++o) {
            float acc = L.bias[o];
            const float* w = &L.weights[static_cast<std::size_t>(o) * L.in_dim];
            for (std::uint32_t i = 0; i < L.in_dim; ++i) acc += w[i] * x[i];
            y[o] = acc;
        }
        switch (L.act) {
            case Activation::Identity:
                break;
            case Activation::ReLU:
                for (float& v : y) if (v < 0.0f) v = 0.0f;
                break;
            case Activation::Softmax: {
                float mx = y[0];
                for (float v : y) if (v > mx) mx = v;
                float sum = 0.0f;
                for (float& v : y) { v = std::exp(v - mx); sum += v; }
                if (sum > 0.0f) for (float& v : y) v /= sum;
                break;
            }
        }
        x = std::move(y);
    }
    return x;
}

}  // namespace branch::wg::phaser
