#!/usr/bin/env python3
"""Train the dual-phaser neural voter and export BRANCH_MLP v1 weights.

The runtime side is a self-contained C++ MLP (src/wg/phaser/reconcile/mlp.cpp)
that loads the text format this script writes. The FEATURE LAYOUT and CLASS
ORDER here MUST match src/wg/phaser/reconcile/neural_features.hpp exactly:

  feature vector (dim 22):
    [0..8]   one-hot(tag_a)              # 9 ReadTag classes
    [9..17]  one-hot(tag_b)
    [18]     conf_a
    [19]     conf_b
    [20]     conf_a - conf_b
    [21]     home_node_a == home_node_b ? 1 : 0
  classes (3): 0 = pick A (de-novo), 1 = pick B (anchored), 2 = flag

Input: a TSV with `dim` feature columns followed by one integer label column
(0/1/2). Generate it by running the dual phaser on data with known truth
(pbsim3 simulated reads or HG002-trio-phased reads): for each DISAGREEMENT,
emit the feature row and the label = which track matched truth (or 2 = flag
when neither did). See scripts/make_voter_trainset.py (data generation).

Usage:
  train_neural_voter.py --data examples.tsv --out voter.mlp.txt [--hidden 32]
  train_neural_voter.py --selftest --out /tmp/voter.mlp.txt   # synthetic round-trip
"""
import argparse
import sys

NUM_TAGS = 9
FEATURE_DIM = 2 * NUM_TAGS + 4   # 22 — keep in sync with neural_features.hpp
NUM_CLASSES = 3


def load_tsv(path):
    import numpy as np
    rows = np.loadtxt(path, dtype=float)
    if rows.ndim == 1:
        rows = rows.reshape(1, -1)
    if rows.shape[1] != FEATURE_DIM + 1:
        raise SystemExit(
            f"expected {FEATURE_DIM + 1} columns ({FEATURE_DIM} features + label), "
            f"got {rows.shape[1]}")
    X = rows[:, :FEATURE_DIM].astype("float32")
    y = rows[:, FEATURE_DIM].astype("int64")
    return X, y


def make_synthetic(n=4000, seed=0):
    """Toy separable trainset: label = argmax(conf_a, conf_b, both-weak->flag)."""
    import numpy as np
    rng = np.random.default_rng(seed)
    X = np.zeros((n, FEATURE_DIM), dtype="float32")
    y = np.zeros(n, dtype="int64")
    for i in range(n):
        ta, tb = rng.integers(0, NUM_TAGS, size=2)
        ca, cb = rng.random(), rng.random()
        X[i, ta] = 1.0
        X[i, NUM_TAGS + tb] = 1.0
        X[i, 2 * NUM_TAGS + 0] = ca
        X[i, 2 * NUM_TAGS + 1] = cb
        X[i, 2 * NUM_TAGS + 2] = ca - cb
        X[i, 2 * NUM_TAGS + 3] = float(rng.integers(0, 2))
        if ca < 0.5 and cb < 0.5:
            y[i] = 2                      # flag
        else:
            y[i] = 0 if ca >= cb else 1   # pick stronger
    return X, y


def train(X, y, hidden, epochs, lr, seed):
    import torch
    import torch.nn as nn
    torch.manual_seed(seed)
    model = nn.Sequential(
        nn.Linear(FEATURE_DIM, hidden),
        nn.ReLU(),
        nn.Linear(hidden, NUM_CLASSES),
    )
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(dev)
    Xt = torch.from_numpy(X).to(dev)
    yt = torch.from_numpy(y).to(dev)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    loss_fn = nn.CrossEntropyLoss()
    for ep in range(epochs):
        opt.zero_grad()
        out = model(Xt)
        loss = loss_fn(out, yt)
        loss.backward()
        opt.step()
        if (ep + 1) % max(1, epochs // 10) == 0:
            acc = (out.argmax(1) == yt).float().mean().item()
            print(f"epoch {ep+1}/{epochs} loss={loss.item():.4f} acc={acc:.4f}",
                  file=sys.stderr)
    return model.cpu()


def export_branch_mlp(model, path):
    """Write nn.Sequential(Linear,ReLU,Linear) as BRANCH_MLP v1.

    Layer activations: hidden Linear -> relu, output Linear -> softmax (the C++
    side applies softmax; CrossEntropyLoss trained on logits is consistent).
    """
    import torch.nn as nn
    layers = [m for m in model if isinstance(m, nn.Linear)]
    acts = ["relu"] * (len(layers) - 1) + ["softmax"]
    with open(path, "w") as o:
        o.write("BRANCH_MLP v1\n")
        o.write(f"{len(layers)}\n")
        for li, (lin, act) in enumerate(zip(layers, acts)):
            W = lin.weight.detach().numpy()   # (out, in)
            b = lin.bias.detach().numpy()      # (out,)
            in_dim, out_dim = W.shape[1], W.shape[0]
            o.write(f"# layer {li}\n{in_dim} {out_dim} {act}\n")
            o.write(" ".join(f"{v:.8g}" for v in W.reshape(-1)) + "\n")
            o.write(" ".join(f"{v:.8g}" for v in b) + "\n")
    print(f"wrote {path} ({len(layers)} layers)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", help="labeled TSV (features + label)")
    ap.add_argument("--out", required=True, help="output BRANCH_MLP weights file")
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--epochs", type=int, default=300)
    ap.add_argument("--lr", type=float, default=1e-2)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--selftest", action="store_true",
                    help="train on synthetic separable data (no --data needed)")
    args = ap.parse_args()

    if args.selftest:
        X, y = make_synthetic(seed=args.seed)
    elif args.data:
        X, y = load_tsv(args.data)
    else:
        raise SystemExit("provide --data or --selftest")

    model = train(X, y, args.hidden, args.epochs, args.lr, args.seed)
    export_branch_mlp(model, args.out)


if __name__ == "__main__":
    main()
