"""Simulated pathogenicity scores for missense variants (practice B7, seed 20261027 + 8).

16,000 missense variants in real MANE Select genes (gene symbols only; the variants themselves are simulated),
3.5 % pathogenic (560), split 60/20/20 into train / calibration / test, stratified by label (exactly 3.5 % each).

Generative model (one latent "damage" value d per variant):
* benign:     d ~ N(0, 1); a few "hard negatives" (benign variants that look damaging) d ~ N(MU_HARD_NEG, 1)
* pathogenic: most d ~ N(MU_HI, 1); some "hard positives" (look harmless) d ~ N(MU_LO, 1)
The true probability of being pathogenic given d follows from Bayes' rule (prevalence 3.5 %): p_true(d).

* Model A (single conservation-like score): a noisy copy of d, reported as a percentile (0-1). It is a RANK,
  not a probability.
* Model B (ensemble): ranks variants by p_true(d), but its raw score is mis-calibrated on purpose:
  logit(p_true) = K * logit(score_B) + C, i.e. score_B = 0.6 means p_true = 0.25 (over-confident).
  On the TRAIN split model B has seen the labels (memorisation): d is shifted by +/- DELTA_TRAIN, so train
  metrics are too optimistic. Calibration and test splits are untouched by the model.

Writes data/classifier_scores.tsv, data/truth/classifier_truth.tsv (variant_id, p_true) and
data/truth/classifier_truth.json (design and reference values).
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logit, logsumexp
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 20261027 + 8
rng = np.random.default_rng(SEED)

N, PREVALENCE = 16_000, 0.035
N_POS = int(round(N * PREVALENCE))                       # 560 pathogenic variants
SPLITS = [("train", 0.6), ("calibration", 0.2), ("test", 0.2)]
Q_HI, MU_HI, MU_LO = 0.95, 4.2, 0.8                      # pathogenic: 95 % clearly damaging, 5 % hard
H_NEG, MU_HARD_NEG = 0.045, 2.3                          # benign: 4.5 % hard negatives
K = 2.0                                                  # mis-calibration slope (logit scale)
C = float(logit(0.25) - K * logit(0.6))                  # anchored so that raw 0.6 -> true 0.25
DELTA_TRAIN = 1.0                                        # memorisation of training labels by model B
W_A, NOISE_A = 1.0, 1.3                                  # model A = W_A * d + NOISE_A * noise
TARGET_PRECISION = 0.90


def p_true_logit(d):
    """log-odds of pathogenic given d (Bayes rule with the mixtures above and prevalence 3.5 %)."""
    d = np.asarray(d, dtype=float)
    lpos = logsumexp(np.stack([np.log(Q_HI) + MU_HI * d - MU_HI ** 2 / 2,
                               np.log(1 - Q_HI) + MU_LO * d - MU_LO ** 2 / 2]), axis=0)
    lneg = logsumexp(np.stack([np.full_like(d, np.log(1 - H_NEG)),
                               np.log(H_NEG) + MU_HARD_NEG * d - MU_HARD_NEG ** 2 / 2]), axis=0)
    return logit(PREVALENCE) + lpos - lneg


def raw_b_from_logit(lt):
    return expit((lt - C) / K)


# ---- labels and stratified splits ---------------------------------------------------------------
y = np.zeros(N, dtype=int)
y[rng.permutation(N)[:N_POS]] = 1
split = np.empty(N, dtype=object)
for label in (0, 1):
    idx = rng.permutation(np.flatnonzero(y == label))
    n_tr, n_ca = int(round(0.6 * idx.size)), int(round(0.2 * idx.size))
    split[idx[:n_tr]], split[idx[n_tr:n_tr + n_ca]], split[idx[n_tr + n_ca:]] = "train", "calibration", "test"

# ---- latent damage, true probability and the two model scores -----------------------------------
hard = rng.random(N)
d = rng.standard_normal(N) + np.where(y == 1, np.where(hard < Q_HI, MU_HI, MU_LO),
                                      np.where(hard < H_NEG, MU_HARD_NEG, 0.0))
p_true = expit(p_true_logit(d))
x_b = d + np.where(split == "train", DELTA_TRAIN * (2 * y - 1), 0.0)
score_b = raw_b_from_logit(p_true_logit(x_b))
x_a = W_A * d + NOISE_A * rng.standard_normal(N)
score_a = (pd.Series(x_a).rank(method="first").to_numpy() - 0.5) / N   # conservation percentile, 0-1

# ---- genes (real MANE Select symbols, sampled uniformly) and ids ---------------------------------
genes = pd.read_csv(ROOT / "build" / "ref" / "genes_hg38_mane_select.tsv.gz", sep="\t", usecols=["gene", "chrom"])
genes = genes[~genes["chrom"].isin(["chrY", "chrM"]) & ~genes["gene"].str.startswith("LOC")]
gene_list = np.sort(genes["gene"].unique())
gene = rng.choice(gene_list, size=N)
order = rng.permutation(N)                                  # ids in random order (no information in the id)
variant_id = np.empty(N, dtype=object)
variant_id[order] = [f"mv{i + 1:05d}" for i in range(N)]

cs = pd.DataFrame({"variant_id": variant_id, "gene": gene, "split": split, "label_pathogenic": y,
                   "score_modelA": np.round(score_a, 6), "score_modelB": np.round(score_b, 6)})
cs = cs.sort_values("variant_id").reset_index(drop=True)
truth = pd.DataFrame({"variant_id": variant_id, "p_true": np.round(p_true, 6), "latent_d": np.round(d, 4)})
truth = truth.sort_values("variant_id").reset_index(drop=True)
(ROOT / "data" / "truth").mkdir(parents=True, exist_ok=True)
cs.to_csv(ROOT / "data" / "classifier_scores.tsv", sep="\t", index=False, na_rep="")
truth.to_csv(ROOT / "data" / "truth" / "classifier_truth.tsv", sep="\t", index=False, na_rep="")


# ---- reference values (computed exactly as in the notebook) --------------------------------------
def r6(x):
    return round(float(x), 6)


def threshold_for_precision(yy, ss, target):
    order_ = np.argsort(-ss, kind="mergesort")
    tp = np.cumsum(yy[order_])
    prec = tp / np.arange(1, yy.size + 1)
    ok = np.flatnonzero((prec >= target) & (yy[order_] == 1))
    return float(ss[order_][ok.max()])


cs = cs.merge(truth[["variant_id", "p_true"]], on="variant_id")
per_split = {}
for s, _ in SPLITS:
    t = cs[cs["split"] == s]
    per_split[s] = {"n": int(len(t)), "n_pathogenic": int(t.label_pathogenic.sum()),
                    "prevalence": r6(t.label_pathogenic.mean()),
                    "AUROC_A": r6(roc_auc_score(t.label_pathogenic, t.score_modelA)),
                    "AUROC_B": r6(roc_auc_score(t.label_pathogenic, t.score_modelB)),
                    "AP_A": r6(average_precision_score(t.label_pathogenic, t.score_modelA)),
                    "AP_B": r6(average_precision_score(t.label_pathogenic, t.score_modelB))}
cal, test = cs[cs["split"] == "calibration"], cs[cs["split"] == "test"]
iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(cal.score_modelB, cal.label_pathogenic)
b_cal = iso.predict(test.score_modelB)
cut = threshold_for_precision(cal.label_pathogenic.to_numpy(), cal.score_modelB.to_numpy(), TARGET_PRECISION)


def pr_at(t, c):
    called = t.score_modelB.to_numpy() >= c
    yy = t.label_pathogenic.to_numpy().astype(bool)
    return {"n_called": int(called.sum()), "TP": int((called & yy).sum()),
            "precision": r6((called & yy).sum() / called.sum()), "recall": r6((called & yy).sum() / yy.sum())}


grid = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
truth_json = {
    "description": "Planted design and reference values for data/classifier_scores.tsv (practice B7). SIMULATED data.",
    "seed": SEED, "n": N, "n_pathogenic": N_POS, "prevalence": PREVALENCE,
    "design": {"Q_HI": Q_HI, "MU_HI": MU_HI, "MU_LO": MU_LO, "H_NEG": H_NEG, "MU_HARD_NEG": MU_HARD_NEG,
               "K": K, "C": round(C, 6), "DELTA_TRAIN": DELTA_TRAIN, "W_A": W_A, "NOISE_A": NOISE_A,
               "model_A": "percentile of W_A*d + NOISE_A*noise (a rank, not a probability)",
               "model_B": "expit((logit p_true(x) - C) / K); x = d on calibration/test, d +/- DELTA_TRAIN on train"},
    "planted_calibration_map": {"formula": "true probability = expit(K * logit(score_B) + C)",
                                "points": {str(g): r6(expit(K * logit(g) + C)) for g in grid}},
    "per_split": per_split,
    "test": {"brier_raw_B": r6(brier_score_loss(test.label_pathogenic, test.score_modelB)),
             "brier_calibrated_B": r6(brier_score_loss(test.label_pathogenic, b_cal)),
             "brier_true_probability": r6(brier_score_loss(test.label_pathogenic, test.p_true)),
             "AUROC_calibrated_B": r6(roc_auc_score(test.label_pathogenic, b_cal)),
             "mae_vs_p_true_raw_B": r6(np.mean(np.abs(test.score_modelB - test.p_true))),
             "mae_vs_p_true_calibrated_B": r6(np.mean(np.abs(b_cal - test.p_true))),
             "isotonic_steps": int(np.unique(iso.predict(cal.score_modelB)).size)},
    "threshold": {"target_precision": TARGET_PRECISION, "chosen_on": "calibration", "cut_score_B": r6(cut),
                  "calibration": pr_at(cal, cut), "test": pr_at(test, cut),
                  "calibrated_probability_at_cut": r6(iso.predict([cut])[0])},
}
(ROOT / "data" / "truth" / "classifier_truth.json").write_text(json.dumps(truth_json, indent=1) + "\n")

if __name__ == "__main__":
    import dviz  # noqa: F401  (only to show the calibration table the notebook will see)
    print(f"rows {len(cs):,}  pathogenic {int(cs.label_pathogenic.sum())}  genes {cs.gene.nunique():,}")
    for s, v in per_split.items():
        print(f"{s:12s} n {v['n']:5d} pos {v['n_pathogenic']:3d} prev {v['prevalence']:.4f} | "
              f"AUROC A {v['AUROC_A']:.3f} B {v['AUROC_B']:.3f} | AP A {v['AP_A']:.3f} B {v['AP_B']:.3f}")
    print("test:", truth_json["test"])
    print("threshold:", truth_json["threshold"])
    print(dviz.calibration_table(test.label_pathogenic, test.score_modelB).round(3).to_string())
    print(dviz.calibration_table(test.label_pathogenic, b_cal).round(3).to_string())
