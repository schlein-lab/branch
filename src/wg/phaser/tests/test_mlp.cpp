// Unit tests for the self-contained MLP voter and its wiring into the
// Reconciler. Standalone (no GoogleTest); assertions kept live via -UNDEBUG.

#include "wg/phaser/reconcile/mlp.hpp"
#include "wg/phaser/reconcile/neural_features.hpp"
#include "wg/phaser/reconcile/reconciler.hpp"
#include "wg/phaser/types.hpp"
#include "wg/phaser/types_dual.hpp"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace branch::wg::phaser;

namespace {

std::string tmp_path(const std::string& name) {
    if (const char* p = std::getenv("TMPDIR")) {
        if (*p) return std::string(p) + "/" + name;
    }
    return "/tmp/" + name;
}

// Write a single-softmax-layer BRANCH_MLP model (in=kFeatureDim, out=3) whose
// biases force a chosen class; all weights are zero so the bias dominates.
std::string write_constant_model(const std::string& name,
                                 float bA, float bB, float bFlag) {
    const std::string path = tmp_path(name);
    std::ofstream o(path);
    assert(o.good());
    o << "BRANCH_MLP v1\n";
    o << "1\n";
    o << "# layer 0\n";
    o << kFeatureDim << " " << kVoterClasses << " softmax\n";
    for (std::uint32_t i = 0; i < kVoterClasses * kFeatureDim; ++i) o << "0 ";
    o << "\n" << bA << " " << bB << " " << bFlag << "\n";
    o.close();
    return path;
}

void test_load_and_forward() {
    const std::string path = write_constant_model("test_mlp_pickA.txt",
                                                  10.0f, 0.0f, 0.0f);
    MLP m;
    std::string err;
    assert(m.load(path, &err));
    assert(m.loaded());
    assert(m.input_dim() == kFeatureDim);
    assert(m.output_dim() == kVoterClasses);

    std::vector<float> in(kFeatureDim, 0.3f);
    auto out = m.forward(in);
    assert(out.size() == kVoterClasses);
    // Softmax: sums to 1, class 0 dominates (bias 10 vs 0).
    float sum = 0.0f;
    for (float v : out) sum += v;
    assert(std::fabs(sum - 1.0f) < 1e-4f);
    assert(out[0] > 0.9f);
    std::remove(path.c_str());
}

void test_bad_model_does_not_load() {
    MLP m;
    assert(!m.load("/nonexistent/path/model.txt"));
    assert(!m.loaded());
    assert(m.forward(std::vector<float>(kFeatureDim, 0.0f)).empty());
}

ReadAssignment mk(ReadId id, ReadTag t, float c) {
    ReadAssignment a;
    a.read_id = id; a.tag = t; a.confidence = c; a.home_node = 0;
    return a;
}

void test_reconcile_neural_path() {
    // Disagreement that the heuristic would resolve to A (higher conf);
    // load a model that forces class B and confirm the model overrides it.
    const std::string path = write_constant_model("test_mlp_pickB.txt",
                                                  0.0f, 10.0f, 0.0f);
    std::vector<ReadAssignment> a = { mk(0, ReadTag::H1_ONLY, 0.95f) };
    std::vector<ReadAssignment> b = { mk(0, ReadTag::H2_ONLY, 0.30f) };

    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts opts;
    opts.verbose = false;
    opts.neural_model_path = path;
    r.reconcile(a, b, out, opts);

    assert(out.size() == 1);
    assert(r.stats().model_used == true);
    assert(out[0].source == PhaseSource::NEURAL_B);   // model, not heuristic
    assert(out[0].tag_final == ReadTag::H2_ONLY);
    assert(r.stats().n_neural_b == 1);
    assert(r.stats().n_heuristic_a == 0);
    std::remove(path.c_str());
}

void test_reconcile_neural_flag() {
    const std::string path = write_constant_model("test_mlp_flag.txt",
                                                  0.0f, 0.0f, 10.0f);
    std::vector<ReadAssignment> a = { mk(0, ReadTag::H1_ONLY, 0.95f) };
    std::vector<ReadAssignment> b = { mk(0, ReadTag::H2_ONLY, 0.95f) };

    Reconciler r;
    std::vector<ReconciledAssignment> out;
    ReconcilerOpts opts;
    opts.verbose = false;
    opts.neural_model_path = path;
    r.reconcile(a, b, out, opts);

    assert(out[0].source == PhaseSource::FLAGGED);
    assert(out[0].tag_final == ReadTag::UNCERTAIN);
    assert(r.stats().n_flagged == 1);
    std::remove(path.c_str());
}

void test_feature_vector_layout() {
    auto f = make_feature_vector(mk(0, ReadTag::H1_ONLY, 0.8f),
                                 mk(0, ReadTag::H2_ONLY, 0.2f));
    assert(f.size() == kFeatureDim);
    assert(f[read_tag_index(ReadTag::H1_ONLY)] == 1.0f);
    assert(f[kNumReadTags + read_tag_index(ReadTag::H2_ONLY)] == 1.0f);
    assert(std::fabs(f[2 * kNumReadTags + 0] - 0.8f) < 1e-6f);
    assert(std::fabs(f[2 * kNumReadTags + 1] - 0.2f) < 1e-6f);
    assert(std::fabs(f[2 * kNumReadTags + 2] - 0.6f) < 1e-6f);
    assert(f[2 * kNumReadTags + 3] == 1.0f);  // same home_node
}

}  // namespace

int main() {
    test_load_and_forward();
    test_bad_model_does_not_load();
    test_feature_vector_layout();
    test_reconcile_neural_path();
    test_reconcile_neural_flag();
    std::printf("test_mlp: all tests passed\n");
    return 0;
}
