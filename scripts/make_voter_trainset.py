#!/usr/bin/env python3
"""Build a labelled trainset for the neural voter from a dual-phaser run.

Inputs (one synthetic diploid run):
  --reads      reads.fastq.gz  (headers carry the TRUE hap: "... hap=1|2 ...")
  --reconciled <prefix>.reconciled.tsv  (read_id, tag_a, tag_b, ..., conf_a,
               conf_b, home_node_a, home_node_b)

Output: a TSV with kFeatureDim feature columns + 1 integer label column, in the
EXACT layout of src/wg/phaser/reconcile/neural_features.hpp. Append many runs
to assemble a large trainset, then feed to train_neural_voter.py.

Truth labelling
---------------
The anchored track maps reads against hap1.fa / hap2.fa, so its h1_only/h2_only
calls are directly truth-aligned. The de-novo track's h1/h2 labels are only
defined up to a global swap, so we resolve that swap by the orientation that
best agrees with the known true haps, then compare. For each DISAGREEMENT:
  label 0 (pick de-novo) if de-novo got the true hap and anchored did not
  label 1 (pick anchored) if anchored got the true hap and de-novo did not
  label 2 (flag)         if neither track recovered the true hap
A track "got the true hap" only when it emitted a hap-bearing tag (h1/h2_only)
matching the read's true source hap. Reads where both tracks agree are skipped
(those are AGREED, not the voter's job).
"""
import argparse
import gzip

# ReadTag string -> index, matching the enum order in types.hpp / the C++
# read_tag_index(). Keep in lockstep.
TAG_INDEX = {
    "untagged": 0, "h1_only": 1, "h2_only": 2, "shared": 3, "branch": 4,
    "uncertain": 5, "branch_bridge": 6, "branch_novel": 7,
    "branch_extinct_parent": 8,
}
NUM_TAGS = 9
FEATURE_DIM = 2 * NUM_TAGS + 4


def read_truth(fastq_gz):
    """read-name -> (true hap, het_overlap) from '@<id> hap=N het=M ...'."""
    truth = {}
    op = gzip.open if fastq_gz.endswith(".gz") else open
    with op(fastq_gz, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            f.readline(); f.readline(); f.readline()  # seq, +, qual
            if not h.startswith("@"):
                continue
            parts = h[1:].split()
            name = parts[0]
            hap = het = None
            for p in parts[1:]:
                if p.startswith("hap="):
                    hap = int(p[4:])
                elif p.startswith("het="):
                    het = int(p[4:])
            if hap is not None:
                truth[name] = (hap, het if het is not None else 0)
    return truth


def tag_hap(tag):
    """The hap a tag asserts, or None if the tag is not hap-bearing."""
    if tag == "h1_only":
        return 1
    if tag == "h2_only":
        return 2
    return None


def load_reconciled(path):
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows.append({
                "read_id": p[col["read_id"]],
                "tag_a": p[col["tag_a"]],
                "tag_b": p[col["tag_b"]],
                "source": p[col["source"]],
                "conf_a": float(p[col["conf_a"]]),
                "conf_b": float(p[col["conf_b"]]),
                "home_a": int(p[col["home_node_a"]]),
                "home_b": int(p[col["home_node_b"]]),
            })
    return rows


def resolve_denovo_flip(rows, truth):
    """Pick the de-novo->true mapping (identity vs swap) with higher concordance.
    Returns a function tag_a -> asserted true hap (or None)."""
    ident = swap = 0
    for r in rows:
        h = tag_hap(r["tag_a"])
        tt = truth.get(r["read_id"])
        if h is None or tt is None or tt[1] == 0:  # need an informative read
            continue
        t = tt[0]
        if h == t:
            ident += 1
        else:
            swap += 1
    flip = swap > ident
    def denovo_hap(tag):
        h = tag_hap(tag)
        if h is None:
            return None
        return (3 - h) if flip else h   # 1<->2 under flip
    return denovo_hap


def feature_row(r):
    f = [0.0] * FEATURE_DIM
    ia = TAG_INDEX.get(r["tag_a"], 0)
    ib = TAG_INDEX.get(r["tag_b"], 0)
    f[ia] = 1.0
    f[NUM_TAGS + ib] = 1.0
    f[2 * NUM_TAGS + 0] = r["conf_a"]
    f[2 * NUM_TAGS + 1] = r["conf_b"]
    f[2 * NUM_TAGS + 2] = r["conf_a"] - r["conf_b"]
    f[2 * NUM_TAGS + 3] = 1.0 if r["home_a"] == r["home_b"] else 0.0
    return f


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True)
    ap.add_argument("--reconciled", required=True)
    ap.add_argument("--out", required=True, help="append-mode TSV output")
    args = ap.parse_args()

    truth = read_truth(args.reads)
    rows = load_reconciled(args.reconciled)
    denovo_hap = resolve_denovo_flip(rows, truth)

    n_emit = [0, 0, 0]
    with open(args.out, "a") as o:
        for r in rows:
            if r["source"] == "agreed":
                continue  # not a disagreement
            tt = truth.get(r["read_id"])
            if tt is None or tt[1] == 0:
                continue  # only label reads that overlap a het SNV (hap is real)
            t = tt[0]
            a_hap = denovo_hap(r["tag_a"])
            b_hap = tag_hap(r["tag_b"])          # anchored is truth-aligned
            a_ok = (a_hap == t)
            b_ok = (b_hap == t)
            if a_ok and not b_ok:
                label = 0
            elif b_ok and not a_ok:
                label = 1
            elif not a_ok and not b_ok:
                label = 2
            else:
                continue  # both recovered truth → not an informative disagreement
            f = feature_row(r)
            o.write("\t".join(f"{v:.6g}" for v in f) + f"\t{label}\n")
            n_emit[label] += 1

    print(f"{args.reconciled}: emitted A={n_emit[0]} B={n_emit[1]} flag={n_emit[2]}")


if __name__ == "__main__":
    main()
