# Dual-phaser neural voter

When the de-novo and pangenome-anchored phasers disagree on a read, the
reconciler can consult a trained neural voter to decide which track to trust
(or to FLAG the read as genuinely uncertain). The voter is a small MLP run by
a dependency-free C++ forward pass (`src/wg/phaser/reconcile/mlp.cpp`); weights
are trained offline and shipped as a text file.

## Using it

```sh
branch wg --reads sample.fastq.gz --read-tech hifi --use-phaser \
          --ref-hap1 hap1.fa --ref-hap2 hap2.fa \
          --neural-voter models/voter.mlp.txt \
          --out-prefix sample.wg
```

Without `--neural-voter`, disagreements fall back to a confidence heuristic,
honestly labelled `heuristic_*` in the `.reconciled.tsv` `source` column. With
the model, decisions are labelled `neural_denovo` / `neural_anchored` /
`flagged`.

## Feature layout (the contract)

`src/wg/phaser/reconcile/neural_features.hpp` defines the 22-dim feature vector
(one-hot tag_a | one-hot tag_b | conf_a | conf_b | conf_a-conf_b |
home_node_a==home_node_b) and the 3 classes (pick A / pick B / flag). The
trainer must mirror this exactly.

## Training

Pipeline (all under `scripts/`, runs on SLURM):

1. `make_synth_diploid.py` — synthetic diploids with known read origin, het-site
   annotation, and `--ref-divergence` (anchored refs differ from the read
   source, simulating mapping to a panel rather than to self).
2. `slurm_trainset_array.sh` — SLURM array spanning genome size / coverage /
   het rate / ref-divergence; runs `branch wg --use-phaser` per diploid.
3. `make_voter_trainset.py` — labels each disagreement by truth (anchored is
   truth-aligned; de-novo flip resolved by concordance; only het-overlapping
   reads are used).
4. `train_neural_voter.py` — pure-NumPy MLP, class-weighted cross-entropy,
   train/val split, per-class recall report; exports `BRANCH_MLP v1`.

`models/voter.mlp.txt` is a reference model trained on 831k labelled
disagreements. Held-out validation (124,659 examples), class-weighted to favour
minority-class recall over raw accuracy:

| class | n | recall |
|---|---|---|
| pick de-novo | 1,354 | 1.000 |
| pick anchored | 109,846 | 0.775 |
| flag | 13,459 | 0.690 |

The de-novo class is ~1% of the data yet fully recalled — the voter learns
*feature-conditional* trust, not a constant "always anchored" rule.

> Note: the reference model is trained on synthetic diploids. The anchored
> track's reliability there is bounded by `--ref-divergence`; on real cohorts
> the model should be retrained on cohort-derived disagreements (e.g. trio- or
> assembly-truth). The pipeline above is the production path for that.
