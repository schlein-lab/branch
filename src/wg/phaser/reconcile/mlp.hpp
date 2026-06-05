#pragma once
//
// MLP — a tiny, self-contained feed-forward neural network for inference.
//
// This is the runtime side of the dual-phaser neural voter (v0.9). It carries
// NO training framework and NO external runtime: weights are trained offline
// (PyTorch on GPU) and exported to a plain-text file that this class loads and
// evaluates with a hand-rolled forward pass. That keeps the BRANCH binary a
// dependency-free static artifact — matching the project's "self-contained C++"
// convention — while the heavy lifting stays in the offline trainer.
//
// Weights file format (text, see scripts/train_neural_voter.py for the writer):
//
//   BRANCH_MLP v1
//   <n_layers>
//   # then, for each layer:
//   <in_dim> <out_dim> <activation>          activation ∈ {relu, softmax, identity}
//   <out_dim*in_dim floats>                  weight matrix W, row-major (out × in)
//   <out_dim floats>                         bias vector b
//
// Lines beginning with '#' and blank lines are ignored. A layer computes
//   y = activation(W · x + b).

#include <cstdint>
#include <string>
#include <vector>

namespace branch::wg::phaser {

class MLP {
public:
    enum class Activation : std::uint8_t { Identity, ReLU, Softmax };

    struct Layer {
        std::uint32_t in_dim  = 0;
        std::uint32_t out_dim = 0;
        Activation    act     = Activation::Identity;
        std::vector<float> weights;  // out_dim × in_dim, row-major
        std::vector<float> bias;     // out_dim
    };

    MLP() = default;

    // Load weights from `path`. Returns false (and leaves the MLP unloaded)
    // on any I/O or format error, writing a reason to `err` if non-null.
    bool load(const std::string& path, std::string* err = nullptr);

    [[nodiscard]] bool loaded() const noexcept { return loaded_; }

    // Input dimension expected by the first layer (0 if unloaded).
    [[nodiscard]] std::uint32_t input_dim() const noexcept {
        return layers_.empty() ? 0u : layers_.front().in_dim;
    }
    // Output dimension of the last layer (0 if unloaded).
    [[nodiscard]] std::uint32_t output_dim() const noexcept {
        return layers_.empty() ? 0u : layers_.back().out_dim;
    }

    // Forward pass. `input.size()` must equal input_dim(); returns the final
    // activations. Returns an empty vector on dimension mismatch / unloaded.
    [[nodiscard]] std::vector<float> forward(const std::vector<float>& input) const;

private:
    std::vector<Layer> layers_;
    bool loaded_ = false;
};

}  // namespace branch::wg::phaser
