#!/usr/bin/env python3
"""Train the dual-phaser neural voter and export BRANCH_MLP v1 weights.

Pure-NumPy implementation (no torch / no GPU runtime needed) of a small
2-layer MLP: Linear(FEATURE_DIM, hidden) -> ReLU -> Linear(hidden, 3) -> softmax.
Trained with class-weighted cross-entropy and Adam, exported to the text format
that src/wg/phaser/reconcile/mlp.cpp loads. The C++ binary stays dependency-free.

The FEATURE LAYOUT and CLASS ORDER MUST match
src/wg/phaser/reconcile/neural_features.hpp exactly:
  features (22): one-hot(tag_a)[9] | one-hot(tag_b)[9] | conf_a | conf_b |
                 conf_a-conf_b | (home_node_a==home_node_b)
  classes (3):   0 = pick A (de-novo), 1 = pick B (anchored), 2 = flag

Input TSV: `dim` feature columns + 1 integer label column (see
make_voter_trainset.py). Usage:
  train_neural_voter.py --data examples.tsv --out voter.mlp.txt [--hidden 32]
  train_neural_voter.py --selftest --out /tmp/voter.mlp.txt
"""
import argparse
import sys

import numpy as np

NUM_TAGS = 9
FEATURE_DIM = 2 * NUM_TAGS + 4   # 22
NUM_CLASSES = 3


def load_tsv(path):
    rows = np.loadtxt(path, dtype=np.float64)
    if rows.ndim == 1:
        rows = rows.reshape(1, -1)
    if rows.shape[1] != FEATURE_DIM + 1:
        raise SystemExit(
            f"expected {FEATURE_DIM + 1} columns, got {rows.shape[1]}")
    return rows[:, :FEATURE_DIM], rows[:, FEATURE_DIM].astype(np.int64)


def make_synthetic(n=4000, seed=0):
    rng = np.random.default_rng(seed)
    X = np.zeros((n, FEATURE_DIM))
    y = np.zeros(n, dtype=np.int64)
    for i in range(n):
        ta, tb = rng.integers(0, NUM_TAGS, size=2)
        ca, cb = rng.random(), rng.random()
        X[i, ta] = 1.0
        X[i, NUM_TAGS + tb] = 1.0
        X[i, 2 * NUM_TAGS:2 * NUM_TAGS + 4] = [ca, cb, ca - cb,
                                               float(rng.integers(0, 2))]
        y[i] = 2 if (ca < 0.5 and cb < 0.5) else (0 if ca >= cb else 1)
    return X, y


def relu(z):
    return np.maximum(z, 0.0)


def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def train(X, y, hidden, epochs, lr, seed, val_frac=0.15):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(y))
    X, y = X[perm], y[perm]
    n_val = max(1, int(len(y) * val_frac))
    Xva, yva = X[:n_val], y[:n_val]
    Xtr, ytr = X[n_val:], y[n_val:]

    counts = np.bincount(ytr, minlength=NUM_CLASSES).astype(np.float64)
    counts[counts == 0] = 1.0
    cw = counts.sum() / (NUM_CLASSES * counts)        # inverse-frequency weights
    print(f"train n={len(ytr)} val n={len(yva)} counts={counts.astype(int)} "
          f"weights={np.round(cw, 3)}", file=sys.stderr)

    # He-init.
    W1 = rng.standard_normal((hidden, FEATURE_DIM)) * np.sqrt(2.0 / FEATURE_DIM)
    b1 = np.zeros(hidden)
    W2 = rng.standard_normal((NUM_CLASSES, hidden)) * np.sqrt(2.0 / hidden)
    b2 = np.zeros(NUM_CLASSES)
    params = [W1, b1, W2, b2]
    m = [np.zeros_like(p) for p in params]
    v = [np.zeros_like(p) for p in params]
    b1m, b2m, eps = 0.9, 0.999, 1e-8

    Yoh = np.eye(NUM_CLASSES)[ytr]
    sample_w = cw[ytr]
    n = len(ytr)

    t = 0
    for ep in range(epochs):
        # forward
        z1 = Xtr @ W1.T + b1
        a1 = relu(z1)
        p = softmax(a1 @ W2.T + b2)
        # weighted CE gradient on logits
        dlogits = (p - Yoh) * sample_w[:, None] / n
        gW2 = dlogits.T @ a1
        gb2 = dlogits.sum(axis=0)
        da1 = dlogits @ W2
        dz1 = da1 * (z1 > 0)
        gW1 = dz1.T @ Xtr
        gb1 = dz1.sum(axis=0)
        grads = [gW1, gb1, gW2, gb2]
        # Adam
        t += 1
        for i in range(4):
            m[i] = b1m * m[i] + (1 - b1m) * grads[i]
            v[i] = b2m * v[i] + (1 - b2m) * grads[i] ** 2
            mhat = m[i] / (1 - b1m ** t)
            vhat = v[i] / (1 - b2m ** t)
            params[i] -= lr * mhat / (np.sqrt(vhat) + eps)
        W1, b1, W2, b2 = params
        if (ep + 1) % max(1, epochs // 10) == 0:
            pv = softmax(relu(Xva @ W1.T + b1) @ W2.T + b2).argmax(1)
            print(f"epoch {ep+1}/{epochs} val_acc={(pv == yva).mean():.4f}",
                  file=sys.stderr)

    pv = softmax(relu(Xva @ W1.T + b1) @ W2.T + b2).argmax(1)
    names = {0: "pickA(denovo)", 1: "pickB(anchored)", 2: "flag"}
    print("=== validation per-class ===", file=sys.stderr)
    for c in range(NUM_CLASSES):
        mask = (yva == c)
        nc = int(mask.sum())
        rec = float((pv[mask] == c).mean()) if nc else 0.0
        print(f"  class {c} {names[c]:16s} n={nc:7d} recall={rec:.3f}",
              file=sys.stderr)
    print(f"  overall val_acc={(pv == yva).mean():.4f}", file=sys.stderr)
    return [("relu", W1, b1), ("softmax", W2, b2)]


def export_branch_mlp(layers, path):
    with open(path, "w") as o:
        o.write("BRANCH_MLP v1\n")
        o.write(f"{len(layers)}\n")
        for li, (act, W, b) in enumerate(layers):
            out_dim, in_dim = W.shape
            o.write(f"# layer {li}\n{in_dim} {out_dim} {act}\n")
            o.write(" ".join(f"{x:.8g}" for x in W.reshape(-1)) + "\n")
            o.write(" ".join(f"{x:.8g}" for x in b) + "\n")
    print(f"wrote {path} ({len(layers)} layers)", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hidden", type=int, default=32)
    ap.add_argument("--epochs", type=int, default=400)
    ap.add_argument("--lr", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        X, y = make_synthetic(seed=args.seed)
    elif args.data:
        X, y = load_tsv(args.data)
    else:
        raise SystemExit("provide --data or --selftest")

    layers = train(X, y, args.hidden, args.epochs, args.lr, args.seed)
    export_branch_mlp(layers, args.out)


if __name__ == "__main__":
    main()
