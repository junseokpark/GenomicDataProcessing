"""Build electives/B7_classifier.ipynb (practice B7 · Classifiers: ROC, precision-recall and calibration).

Numbers quoted in the markdown are computed here from data/classifier_scores.tsv and data/truth/classifier_truth.json
(written by build/make_classifier.py) with the same functions the notebook uses, so text and CHECK cells agree.
usage:  python build/build_B7_notebook.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score, roc_curve

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
from nbbuild import NB  # noqa: E402
import dviz  # noqa: E402

T = json.loads((ROOT / "data" / "truth" / "classifier_truth.json").read_text())
CS = pd.read_csv(ROOT / "data" / "classifier_scores.tsv", sep="\t")
PS = T["per_split"]
TE, TH = T["test"], T["threshold"]
test = CS[CS.split == "test"]
calib = CS[CS.split == "calibration"]
y = test.label_pathogenic.to_numpy()
n_benign_test = int((y == 0).sum())
fp_at_80 = {}
for m in ["A", "B"]:
    fpr, tpr, _ = roc_curve(y, test[f"score_model{m}"])
    i = int(np.argmax(tpr >= 0.80))
    fp_at_80[m] = (fpr[i], fpr[i] * n_benign_test)
rng = np.random.default_rng(20261027 + 8)                 # bootstrap of the test split (for one answer below)
boot = []
for _ in range(1000):
    idx = rng.integers(0, len(test), len(test))
    boot.append(roc_auc_score(y[idx], test.score_modelB.to_numpy()[idx]))
ci_lo, ci_hi = np.percentile(boot, [2.5, 97.5])
rec_at_085 = {}
for m in ["A", "B"]:
    pp, rr, _ = precision_recall_curve(y, test[f"score_model{m}"])
    rec_at_085[m] = float(rr[pp >= 0.85].max())
ct_raw = dviz.calibration_table(test.label_pathogenic, test.score_modelB)
row6 = ct_raw[ct_raw["bin"] == 5].iloc[0]                 # raw scores 0.5-0.6
row7 = ct_raw[ct_raw["bin"] == 6].iloc[0]                 # raw scores 0.6-0.7
iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(calib.score_modelB, calib.label_pathogenic)
b_cal = iso.predict(test.score_modelB)
top_step = float(np.max(iso.predict(calib.score_modelB)))
n_top = int(np.sum(b_cal == top_step))
n_top_benign = int(np.sum((b_cal == top_step) & (y == 0)))
below = TH['test']['precision'] < TH['target_precision']
ct_cal = dviz.calibration_table(test.label_pathogenic, b_cal)
n_bins_cal = len(ct_cal)
gap = (ct_cal['mean_predicted'] - ct_cal['observed']).abs()
worst = ct_cal.loc[gap.idxmax()]
pm = T["planted_calibration_map"]["points"]
prev = PS["test"]["prevalence"]
d_auc = TE["AUROC_calibrated_B"] - PS["test"]["AUROC_B"]


def p3(x):
    return f"{x:.3f}"


nb = NB()
# ------------------------------------------------------------------------------------------------ title
nb.md(f"""
# Practice B7 · Classifiers: ROC, precision–recall and calibration

**Prevalence per split → ROC and AUROC → precision–recall and AP against the prevalence baseline → calibration with isotonic regression → a cut-off for 90% precision**

*Genomic Data Processing · Day 5 (10/27) · Data visualization in genetics and genomics · about 18 minutes*

> **All data in this practice are simulated (computer-generated). No real patient data are used.**
> The {T['n']:,} missense variants carry real gene symbols from the MANE Select list (hg38), so the table looks like a
> real benchmark; the variants, labels and both model scores are simulated with a known answer (the true probability of
> every variant is stored in `data/truth/`).

**How to use this notebook.** Run the cells from top to bottom (Shift + Enter). Every step follows the same four beats:
**Question → Test → Check → Visualize**. A cell that starts with `# CHECK ·` compares your number with the reference
value and prints ✓ or ✗; if you see ✗, stop and find out why. Each step ends with **What changed?** and a short
**Your turn** task with a folded answer (click *Answer*). Prompt card 8 (`prompts/prompt_cards.md`) is one way to ask an
AI assistant for the plotting code; the provided cells are the other. Either way, the CHECK cells decide.

| File | What is inside |
|---|---|
| `data/classifier_scores.tsv` | one row per variant: gene, split, label (1 = pathogenic), score of model A, raw score of model B |
| `data/truth/classifier_truth.tsv`, `.json` | the planted true probability per variant and the reference values – CHECK cells only |

| Step | Topic | min |
|---|---|---|
| B7.1 | prevalence per split; why we judge on the test split | 3 |
| B7.2 | ROC curve and AUROC | 3 |
| B7.3 | precision–recall, AP and the prevalence baseline | 3 |
| B7.4 | calibration: does a score of 0.6 mean 60%? | 5 |
| B7.5 | a cut-off for 90% precision | 4 |
""")
nb.code("""
# Step 0 · Set-up: imports, folders and figure style
import os, sys, time
from pathlib import Path
ROOT = Path.cwd() if (Path.cwd() / "dviz.py").exists() else Path.cwd().parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))   # needed by umap-learn on some systems
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dviz
dviz.use_style()
DATA, FIG = ROOT / "data", ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
T_START = time.time()
import json
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score, brier_score_loss
from sklearn.isotonic import IsotonicRegression
""")

# ------------------------------------------------------------------------------------------------ B7.1
nb.md(f"""
## Step B7.1 · Prevalence per split, and why we judge on the test split  *(3 min)*

### Question

Two computer models give every missense variant a score for "how likely is this variant to cause disease?".
**Which model should a diagnostic lab use when only about 3.5% of the variants are pathogenic?**

**Why this step.** A *missense* variant changes one amino acid in a protein. Most are harmless; a few break the
protein and cause disease. Model A is a single conservation-like score (amino acids that never change during evolution
are probably important); model B is an *ensemble* that combines many features and was **trained** on labelled variants.
A model that has seen a variant's label during training scores it too well, so we keep variants apart: **train**
(model B learned from them), **calibration** (we use them in Steps B7.4–B7.5 to adjust the scores and to choose a
cut-off) and **test** (never used for anything, so it gives an honest answer).

**New terms**

* **Label** – the known answer per variant: `label_pathogenic` = 1 (pathogenic) or 0 (benign), e.g. from expert curation.
* **Prevalence** – fraction of variants that are pathogenic (here {prev:.1%}).
* **Stratified split** – the split keeps the same prevalence in every part.
* **Leakage / optimism** – judging a model on data it (or a step tuned on it) has already seen; the result looks too good.

### Test
""")
nb.code("""
# Step B7.1 · load the scores and count variants, pathogenic variants and prevalence per split
cs = pd.read_csv(DATA / "classifier_scores.tsv", sep="\\t")
print(f"{len(cs):,} missense variants in {cs['gene'].nunique():,} genes")
SPLITS = ["train", "calibration", "test"]
per_split = cs.groupby("split")["label_pathogenic"].agg(variants="size", pathogenic="sum").loc[SPLITS]
per_split["prevalence"] = (per_split["pathogenic"] / per_split["variants"]).round(4)
for s in SPLITS:                                      # model B judged on each split
    part = cs[cs["split"] == s]
    per_split.loc[s, "AUROC_B"] = round(roc_auc_score(part["label_pathogenic"], part["score_modelB"]), 3)
    per_split.loc[s, "AP_B"] = round(average_precision_score(part["label_pathogenic"], part["score_modelB"]), 3)
test = cs[cs["split"] == "test"].copy()               # the only split used for the final numbers
calib = cs[cs["split"] == "calibration"].copy()       # used to calibrate (B7.4) and to choose a cut-off (B7.5)
y = test["label_pathogenic"].to_numpy()
prevalence = float(y.mean())
per_split
""")
nb.md("### Check")
nb.code("""
# CHECK · same prevalence in every split (stratified), and the train split looks too good for model B
ref = json.loads((DATA / "truth" / "classifier_truth.json").read_text())
dviz.check("B7.prevalence_test", round(prevalence, 4), ref["per_split"]["test"]["prevalence"], tol=1e-4)
dviz.check("B7.n_pathogenic_test", int(y.sum()), ref["per_split"]["test"]["n_pathogenic"])
dviz.check("B7.same_prevalence_all_splits", bool(per_split["prevalence"].max() - per_split["prevalence"].min() < 1e-3), True)
dviz.check("B7.train_AUROC_B_inflated", bool(per_split.loc["train", "AUROC_B"] > per_split.loc["test", "AUROC_B"] + 0.01), True)
""")
nb.md("### Visualize")
nb.code("""
# Step B7.1 · figure: model B raw scores on the test split, benign vs pathogenic (y axis logarithmic)
bins = np.linspace(0, 1, 41)
fig, ax = plt.subplots(figsize=(5.4, 3.1))
for lab, name, col in [(0, "benign", dviz.SLATE), (1, "pathogenic", dviz.BLUE)]:
    vals = test.loc[test["label_pathogenic"] == lab, "score_modelB"]
    ax.hist(vals, bins=bins, histtype="step", lw=1.6, color=col, label=f"{name} (n = {len(vals):,})")
ax.set_yscale("log")
ax.set(xlabel="Model B raw score", ylabel="Variants per bin (log scale)", xlim=(0, 1),
       title=f"Test split: {prevalence:.1%} of variants are pathogenic")
ax.legend(loc="upper center")
dviz.savefig(fig, "B7_scores")
plt.show()
""")
nb.md(f"""
**What changed?** We counted the labels per split: every split has {prev:.1%} pathogenic variants
({PS['train']['n_pathogenic']} / {PS['calibration']['n_pathogenic']} / {PS['test']['n_pathogenic']} of
{PS['train']['n']:,} / {PS['calibration']['n']:,} / {PS['test']['n']:,}). Model B looks almost perfect on the train
split (AUROC {p3(PS['train']['AUROC_B'])}, AP {p3(PS['train']['AP_B'])}) but not on the test split (AUROC
{p3(PS['test']['AUROC_B'])}, AP {p3(PS['test']['AP_B'])}): the train numbers are inflated because the model saw those
labels. The histogram shows the imbalance – {n_benign_test:,} benign against {PS['test']['n_pathogenic']} pathogenic
variants – which is why the y axis is logarithmic. From now on, every number comes from the test split.

**Your turn.** Model B never saw the calibration split either, yet its AUROC there is {p3(PS['calibration']['AUROC_B'])}
against {p3(PS['test']['AUROC_B'])} on the test split. Why can two untouched splits disagree?

<details><summary>Answer</summary>

Chance. Each split holds only {PS['test']['n_pathogenic']} pathogenic variants, so the metrics move by a few hundredths
between random splits. A bootstrap of the test split (1,000 resamples) gives a 95% interval for model B's AUROC of
{ci_lo:.3f}–{ci_hi:.3f}. That is why the test split is chosen before looking, used once, and reported with its n
(and ideally an interval) – never "the split that looks best".
</details>
""")

# ------------------------------------------------------------------------------------------------ B7.2
nb.md("""
## Step B7.2 · ROC curve and AUROC  *(3 min)*

### Question

Over all possible cut-offs, how well does each score separate pathogenic from benign variants?

**Why this step.** A model gives scores, not decisions. If we call every variant above a cut-off "pathogenic", each
cut-off finds some pathogenic variants (true positives) and wrongly flags some benign ones (false positives). The **ROC
curve** shows the trade-off for every cut-off at once, and the area under it summarises the curve in one number.

**New terms**

* **Cut-off (threshold)** – the score above which we call a variant pathogenic.
* **TPR** (true-positive rate, sensitivity, recall) – pathogenic variants above the cut-off / all pathogenic variants.
* **FPR** (false-positive rate, 1 − specificity) – benign variants above the cut-off / all benign variants.
* **AUROC** – area under the ROC curve = probability that a random pathogenic variant scores higher than a random
  benign one (0.5 = coin flip, 1 = perfect).

### Test
""")
nb.code("""
# Step B7.2 · ROC curve and AUROC on the test split for both models
roc, auroc = {}, {}
for m in ["A", "B"]:
    fpr, tpr, thr = roc_curve(y, test[f"score_model{m}"])
    roc[m] = (fpr, tpr)
    auroc[m] = roc_auc_score(y, test[f"score_model{m}"])
    i = np.argmax(tpr >= 0.80)                                 # first cut-off that finds 80% of pathogenic variants
    print(f"model {m}: AUROC {auroc[m]:.3f}; to find 80% of pathogenic variants FPR = {fpr[i]:.3f}"
          f" -> {fpr[i] * (y == 0).sum():.0f} of {(y == 0).sum():,} benign variants flagged")
""")
nb.md("### Check")
nb.code("""
# CHECK · AUROC of both models (reference values and the expected ranges of the design)
dviz.check("B7.AUROC_A", round(auroc["A"], 4), round(ref["per_split"]["test"]["AUROC_A"], 4), tol=1e-4)
dviz.check("B7.AUROC_B", round(auroc["B"], 4), round(ref["per_split"]["test"]["AUROC_B"], 4), tol=1e-4)
dviz.check("B7.AUROC_A_in_design_range", auroc["A"], between=(0.92, 0.94))
dviz.check("B7.AUROC_B_in_design_range", auroc["B"], between=(0.96, 0.98))
""")
nb.md("### Visualize")
nb.code("""
# Step B7.2 · figure: ROC curves of model A (reference, slate) and model B (blue)
fig, ax = plt.subplots(figsize=(3.9, 3.7))
ax.plot([0, 1], [0, 1], ls="--", color=dviz.GREY, lw=1.2, label="random (AUROC 0.50)")
for m, col in [("A", dviz.SLATE), ("B", dviz.BLUE)]:
    ax.plot(*roc[m], color=col, lw=1.8, label=f"model {m}: AUROC {auroc[m]:.3f}")
ax.set(xlabel="False-positive rate (1 − specificity)", ylabel="True-positive rate (recall)", xlim=(0, 1), ylim=(0, 1.01),
       title=f"ROC, test split (n = {len(y):,}; {int(y.sum())} pathogenic)")
ax.set_aspect("equal")
ax.legend(loc="lower right")
dviz.savefig(fig, "B7_roc")
plt.show()
""")
nb.md(f"""
**What changed?** Each score became a curve of {len(y):,} possible cut-offs and one number: AUROC
{p3(PS['test']['AUROC_A'])} for model A and {p3(PS['test']['AUROC_B'])} for model B. Both look excellent. But read the
printed line: to find 80% of the pathogenic variants, model A also flags about {fp_at_80['A'][1]:.0f} benign variants
(FPR {p3(fp_at_80['A'][0])}) and model B about {fp_at_80['B'][1]:.0f} (FPR {p3(fp_at_80['B'][0])}) – a small FPR times
{n_benign_test:,} benign variants is a large number next to {PS['test']['n_pathogenic']} pathogenic ones.

**Your turn.** Predict: if we added 10,000 more benign variants that look like the ones we have, would the ROC curve
and AUROC change? Would the number of false alarms change?

<details><summary>Answer</summary>

The ROC curve and AUROC would stay (almost) the same: FPR is a *rate* among benign variants, so it does not depend on
how many there are. The number of false alarms at the same cut-off would grow about 4-fold. ROC ignores prevalence;
the precision–recall curve in the next step does not.
</details>
""")

# ------------------------------------------------------------------------------------------------ B7.3
nb.md(f"""
## Step B7.3 · Precision–recall, AP and the prevalence baseline  *(3 min)*

### Question

Of the variants a model flags, what fraction is really pathogenic – and how much better is that than guessing?

**Why this step.** With {prev:.1%} prevalence there are about 28 benign variants for every pathogenic one. A lab that
reviews flagged variants cares about **precision**: how many of the flagged variants are real. The precision–recall
(PR) curve shows precision against recall for every cut-off. A model that guesses at random has precision equal to
the prevalence at every recall, so the prevalence is the **baseline** of the PR plot (Saito & Rehmsmeier, 2015).

**New terms**

* **Precision** (positive predictive value) – pathogenic variants above the cut-off / all variants above the cut-off.
* **AP (average precision)** – area under the PR curve (the mean precision over all recall steps); baseline = prevalence.

### Test
""")
nb.code("""
# Step B7.3 · precision-recall curve and average precision (AP) on the test split
pr, ap = {}, {}
for m in ["A", "B"]:
    precision, recall, thr = precision_recall_curve(y, test[f"score_model{m}"])
    pr[m] = (recall, precision)
    ap[m] = average_precision_score(y, test[f"score_model{m}"])
    print(f"model {m}: AP {ap[m]:.3f} = {ap[m] / prevalence:.0f} x the baseline (prevalence {prevalence:.3f})")
""")
nb.md("### Check")
nb.code("""
# CHECK · both models beat the prevalence baseline, and model B has the higher AP
dviz.check("B7.AP_A", round(ap["A"], 4), round(ref["per_split"]["test"]["AP_A"], 4), tol=1e-4)
dviz.check("B7.AP_B", round(ap["B"], 4), round(ref["per_split"]["test"]["AP_B"], 4), tol=1e-4)
dviz.check("B7.AP_A_above_prevalence", bool(ap["A"] > prevalence), True)
dviz.check("B7.AP_B_above_prevalence", bool(ap["B"] > prevalence), True)
dviz.check("B7.AP_B_above_AP_A", bool(ap["B"] > ap["A"]), True)
""")
nb.md("### Visualize")
nb.code("""
# Step B7.3 · figure: precision-recall curves with the prevalence baseline (dashed)
fig, ax = plt.subplots(figsize=(3.9, 3.7))
for m, col in [("A", dviz.SLATE), ("B", dviz.BLUE)]:
    ax.plot(*pr[m], color=col, lw=1.8, drawstyle="steps-post", label=f"model {m}: AP {ap[m]:.3f}")
ax.axhline(prevalence, ls="--", color=dviz.GREY, lw=1.2)
ax.text(0.02, prevalence + 0.02, f"baseline = prevalence {prevalence:.1%}", color=dviz.SLATE, fontsize=8)
ax.set(xlabel="Recall (true-positive rate)", ylabel="Precision", xlim=(0, 1), ylim=(0, 1.02),
       title=f"Precision–recall, test split (n = {len(y):,})")
ax.set_aspect("equal")
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.3))
dviz.savefig(fig, "B7_pr")
plt.show()
""")
nb.md(f"""
**What changed?** The same scores, plotted as precision against recall, separate the models clearly: AP
{p3(PS['test']['AP_A'])} for model A and {p3(PS['test']['AP_B'])} for model B, against a baseline of {prev:.3f}
({PS['test']['AP_A'] / prev:.0f}× and {PS['test']['AP_B'] / prev:.0f}× better than guessing). With precision still ≥ 0.85,
model B reaches a recall of {rec_at_085['B']:.2f}, model A only {rec_at_085['A']:.2f}.

**Your turn.** At 3.5% prevalence a model has AP 0.20. Is that a bad model? What would AP 0.20 mean at 50% prevalence?

<details><summary>Answer</summary>

Compare AP with its baseline, not with 1: at 3.5% prevalence, AP 0.20 is almost 6× better than guessing and can be
useful. At 50% prevalence the baseline is 0.50, so AP 0.20 would be worse than random. Always print the prevalence next
to an AP and draw it as the dashed baseline.
</details>
""")

# ------------------------------------------------------------------------------------------------ B7.4
nb.md(f"""
## Step B7.4 · Calibration: does a score of 0.6 mean 60%?  *(5 min)*

### Question

Can we read model B's score as a probability, and if not, can we fix it without changing the ranking?

**Why this step.** A clinician who reads "0.6" thinks "6 in 10 such variants are pathogenic". That is only true if the
score is **calibrated**. A *calibration curve* groups variants by predicted value (10 bins) and compares the mean
prediction with the observed fraction of pathogenic variants; a calibrated model lies on the diagonal. Model B's raw
score is **over-confident**. We correct it with **isotonic regression**, fitted on the calibration split, and judge the
result on the test split (Niculescu-Mizil & Caruana, 2005).

**New terms**

* **Calibration (reliability) curve** – observed fraction pathogenic per bin of predicted probability; bins with fewer
  than 20 variants are skipped (`dviz.calibration_table`).
* **Over-confident** – the score claims a higher probability than is observed.
* **Isotonic regression** – fits the best non-decreasing step function from score to observed fraction; it never
  reverses the order of two variants, but it can merge neighbouring scores into one step (a tie).
* **Brier score** – mean of (predicted probability − label)²; 0 = perfect, lower is better.

### Test
""")
nb.code("""
# Step B7.4 · calibration table of model B's RAW score on the test split (10 bins; bins with < 20 variants dropped)
cal_raw = dviz.calibration_table(test["label_pathogenic"], test["score_modelB"], n_bins=10, min_count=20)
cal_raw.round(3)
""")
nb.code("""
# Step B7.4 · fit isotonic regression on the CALIBRATION split, then apply it to the TEST split
iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip")
iso.fit(calib["score_modelB"], calib["label_pathogenic"])
test["score_modelB_calibrated"] = iso.predict(test["score_modelB"])
cal_iso = dviz.calibration_table(test["label_pathogenic"], test["score_modelB_calibrated"])
brier = {"raw": brier_score_loss(y, test["score_modelB"]),
         "calibrated": brier_score_loss(y, test["score_modelB_calibrated"])}
auroc_cal = roc_auc_score(y, test["score_modelB_calibrated"])
n_steps = len(np.unique(iso.predict(calib["score_modelB"])))
print(f"Brier score: raw {brier['raw']:.4f} -> calibrated {brier['calibrated']:.4f}")
print(f"AUROC: raw {auroc['B']:.4f} -> calibrated {auroc_cal:.4f}  (isotonic step function with {n_steps} steps)")
cal_iso.round(3)
""")
nb.md("### Check")
nb.code("""
# CHECK · lower Brier score, the order of the variants is kept, AUROC moves only through ties
dviz.check("B7.brier_raw_B", round(brier["raw"], 4), round(ref["test"]["brier_raw_B"], 4), tol=1e-4)
dviz.check("B7.brier_calibrated_B", round(brier["calibrated"], 4), round(ref["test"]["brier_calibrated_B"], 4), tol=1e-4)
dviz.check("B7.brier_improved", bool(brier["calibrated"] < brier["raw"]), True)
order = np.argsort(test["score_modelB"].to_numpy(), kind="mergesort")
kept = bool(np.all(np.diff(test["score_modelB_calibrated"].to_numpy()[order]) >= 0))
dviz.check("B7.isotonic_never_reverses_order", kept, True)
dviz.check("B7.AUROC_change_after_calibration", round(auroc_cal - auroc["B"], 4), 0, tol=0.01, note="ties only")
p_true = pd.read_csv(DATA / "truth" / "classifier_truth.tsv", sep="\\t").set_index("variant_id")["p_true"]
p_true = p_true.loc[test["variant_id"]].to_numpy()             # planted probability of each test variant
err = {k: float(np.mean(np.abs(test[c].to_numpy() - p_true))) for k, c in
       [("raw", "score_modelB"), ("calibrated", "score_modelB_calibrated")]}
dviz.check("B7.mean_error_vs_planted_raw", round(err["raw"], 4), round(ref["test"]["mae_vs_p_true_raw_B"], 4), tol=1e-4)
dviz.check("B7.mean_error_vs_planted_calibrated", round(err["calibrated"], 4),
           round(ref["test"]["mae_vs_p_true_calibrated_B"], 4), tol=1e-4)
""")
nb.md("### Visualize")
nb.code("""
# Step B7.4 · figure: calibration curve of model B before (blue) and after (green) isotonic calibration
fig, ax = plt.subplots(figsize=(3.9, 3.7))
ax.plot([0, 1], [0, 1], ls="--", color=dviz.GREY, lw=1.2, label="perfectly calibrated")
for tab, name, col in [(cal_raw, "raw", dviz.BLUE), (cal_iso, "calibrated", dviz.GREEN)]:
    ax.plot(tab["mean_predicted"], tab["observed"], "-o", color=col, ms=4.5,
            label=f"model B {name} (Brier {brier[name]:.4f})")
ax.set(xlabel="Predicted probability (mean per bin)", ylabel="Observed fraction pathogenic", xlim=(-0.01, 1.03), ylim=(-0.01, 1.03),
       title="Calibration, test split (10 bins, ≥ 20 variants)")
ax.set_aspect("equal")
ax.legend(loc="upper left", fontsize=7.5)
dviz.savefig(fig, "B7_calibration")
plt.show()
""")
nb.md(f"""
**What changed?** We added one column, `score_modelB_calibrated`. Before calibration, test variants with raw scores
0.5–0.6 (mean {row6['mean_predicted']:.2f}) were pathogenic only {row6['observed']:.0%} of the time and those with
0.6–0.7 (mean {row7['mean_predicted']:.2f}) only {row7['observed']:.0%}: model B is over-confident (the planted rule is
raw 0.6 → true {pm['0.6']:.2f}). After isotonic regression the {n_bins_cal} remaining bins move toward the diagonal
(the largest gap left – predicted {worst['mean_predicted']:.2f}, observed {worst['observed']:.2f}, {int(worst['n'])} variants –
comes from steps fitted on only {PS['calibration']['n_pathogenic']} pathogenic calibration variants), the Brier score fell from
{TE['brier_raw_B']:.4f} to {TE['brier_calibrated_B']:.4f} (the planted true probabilities score {TE['brier_true_probability']:.4f},
the best possible), and the mean distance to the planted probability fell from {TE['mae_vs_p_true_raw_B']:.3f} to
{TE['mae_vs_p_true_calibrated_B']:.3f}. The order of the variants is unchanged; AUROC moved by {d_auc:+.4f} only because
isotonic regression puts many variants on the same step (ties).

**Your turn.** The highest isotonic step is {top_step:.2f}, and {n_top} test variants received it. Would you tell a
family "this variant is pathogenic with probability {top_step:.2f}"? Why might isotonic regression be too sure at the
extremes?

<details><summary>Answer</summary>

No. The top step is {top_step:.2f} because every calibration variant above some raw score happened to be pathogenic – a
handful of variants out of {PS['calibration']['n_pathogenic']} pathogenic ones. With so few variants at the extremes,
isotonic regression copies their observed fraction exactly; on the test split {n_top_benign} of the {n_top} variants on
that step {'is' if n_top_benign == 1 else 'are'} benign.
Report calibrated probabilities with their uncertainty (or cap them), and remember that calibration is only valid for
data like the calibration split (same population, same prevalence).
</details>
""")

# ------------------------------------------------------------------------------------------------ B7.5
nb.md(f"""
## Step B7.5 · A cut-off for 90% precision  *(4 min)*

### Question

Which variants should go to the expert panel if at least 9 out of 10 flagged variants should be pathogenic, and how
many pathogenic variants will we find that way?

**Why this step.** In practice we need a yes/no decision. We choose the cut-off on the **calibration** split (going
down the ranked list until precision would drop below 0.90) and then **report** its precision and recall on the test
split. Choosing the cut-off on the test split would make the result look better than it will be on new data.

**New terms**

* **Operating point** – one chosen cut-off, with its precision and recall.
* **Optimism** – a cut-off tuned on one data set usually performs a little worse on new data.

### Test
""")
nb.code("""
# Step B7.5 · walk down the ranked calibration split; keep the lowest cut-off whose precision is still >= 0.90
def threshold_for_precision(labels, scores, target=0.90):
    order = np.argsort(-scores, kind="mergesort")          # highest score first
    hits = np.cumsum(labels[order])                          # pathogenic variants found so far
    precision = hits / np.arange(1, len(labels) + 1)
    ok = (precision >= target) & (labels[order] == 1)        # stop only at a pathogenic variant
    return float(scores[order][np.flatnonzero(ok).max()])


def precision_recall_at(labels, scores, cut):
    called = scores >= cut
    tp = int(np.sum(called & (labels == 1)))
    return tp / called.sum(), tp / labels.sum(), int(called.sum())


yc, sc = calib["label_pathogenic"].to_numpy(), calib["score_modelB"].to_numpy()
cut = threshold_for_precision(yc, sc)
p_cal, r_cal, n_cal = precision_recall_at(yc, sc, cut)
p_test, r_test, n_test = precision_recall_at(y, test["score_modelB"].to_numpy(), cut)
print(f"cut-off on the raw score of model B: {cut:.4f} (calibrated probability there: {iso.predict([cut])[0]:.2f})")
print(f"calibration split: {n_cal} flagged, precision {p_cal:.3f}, recall {r_cal:.3f}")
print(f"test split:        {n_test} flagged, precision {p_test:.3f}, recall {r_test:.3f}")
""")
nb.md("### Check")
nb.code("""
# CHECK · the cut-off matches the reference; on new data precision is close to, but below, the target
dviz.check("B7.cut_score_B", round(cut, 6), ref["threshold"]["cut_score_B"], tol=1e-6)
dviz.check("B7.precision_calibration_at_cut", round(p_cal, 4), round(ref["threshold"]["calibration"]["precision"], 4), tol=1e-4)
dviz.check("B7.precision_test_at_cut", round(p_test, 4), round(ref["threshold"]["test"]["precision"], 4), tol=1e-4)
dviz.check("B7.recall_test_at_cut", round(r_test, 4), round(ref["threshold"]["test"]["recall"], 4), tol=1e-4)
dviz.check("B7.test_precision_near_target", p_test, between=(0.80, 1.0))
""")
nb.md("### Visualize")
nb.code("""
# Step B7.5 · figure: PR curves of model B on the calibration and test splits, with the chosen operating point
fig, ax = plt.subplots(figsize=(4.4, 3.7))
prec_c, rec_c, _ = precision_recall_curve(yc, sc)
ax.plot(rec_c, prec_c, color=dviz.GREY, lw=1.6, drawstyle="steps-post", label="calibration split (cut-off chosen here)")
ax.plot(*pr["B"], color=dviz.BLUE, lw=1.8, drawstyle="steps-post", label="test split (reported)")
ax.axhline(0.90, ls="--", color=dviz.SLATE, lw=1)
ax.text(0.02, 0.915, "target precision 0.90", fontsize=8, color=dviz.SLATE)
ax.plot(r_cal, p_cal, "o", color=dviz.SLATE, ms=7)
ax.plot(r_test, p_test, "o", color=dviz.BLUE, ms=7)
ax.annotate(f"cut-off {cut:.2f}\\ntest: precision {p_test:.2f}, recall {r_test:.2f}", (r_test, p_test),
            xytext=(-120, -45), textcoords="offset points", fontsize=8, color=dviz.INK,
            arrowprops=dict(arrowstyle="-", color=dviz.SLATE, lw=0.8))
ax.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0.4, 1.02),
       title="Model B: one cut-off, two splits")
ax.legend(loc="lower left", fontsize=7.5)
dviz.savefig(fig, "B7_threshold")
plt.show()
""")
nb.md(f"""
**What changed?** One number now turns scores into decisions: raw score ≥ {TH['cut_score_B']:.4f}. On the calibration
split it flags {TH['calibration']['n_called']} variants with precision {p3(TH['calibration']['precision'])}; on the
untouched test split it flags {TH['test']['n_called']} variants, {TH['test']['TP']} of them pathogenic (precision
{p3(TH['test']['precision'])}, recall {p3(TH['test']['recall'])}). {'Precision fell below the 0.90 target on new data – the optimism of a cut-off chosen on ' + str(PS['calibration']['n_pathogenic']) + ' pathogenic variants.' if below else 'Precision stayed at or above the 0.90 target on new data, but with ' + str(PS['calibration']['n_pathogenic']) + ' pathogenic variants this is partly luck.'}

**Your turn.** The calibrated probability at the cut-off is about {TH['calibrated_probability_at_cut']:.2f}, yet the
precision of the flagged list is about 0.9. How can both be true?

<details><summary>Answer</summary>

The cut-off is the *lowest* score that is still flagged. Variants right at the cut-off are pathogenic about
{TH['calibrated_probability_at_cut']:.0%} of the time, but most flagged variants score much higher (calibrated
probabilities near 1), so the average over the whole list – the precision – is about 0.9. A cut-off probability and the
precision of a list are different numbers; say which one you report.
</details>
""")

# ------------------------------------------------------------------------------------------------ summary
nb.md(f"""
## Summary

| Step | What you did | Number |
|---|---|---|
| B7.1 | prevalence per split; train vs test | {prev:.1%} in every split; model B AUROC {p3(PS['train']['AUROC_B'])} (train) vs {p3(PS['test']['AUROC_B'])} (test) |
| B7.2 | ROC and AUROC on the test split | A {p3(PS['test']['AUROC_A'])} · B {p3(PS['test']['AUROC_B'])} |
| B7.3 | PR and AP with the prevalence baseline | A {p3(PS['test']['AP_A'])} · B {p3(PS['test']['AP_B'])} · baseline {prev:.3f} |
| B7.4 | isotonic calibration (fit: calibration split) | Brier {TE['brier_raw_B']:.4f} → {TE['brier_calibrated_B']:.4f}; AUROC change {d_auc:+.4f} (ties) |
| B7.5 | cut-off for precision ≥ 0.90 (chosen on calibration) | test: precision {p3(TH['test']['precision'])}, recall {p3(TH['test']['recall'])} |

**Your presentation sentence** (fill in from your outputs):

> **Question.** Which pathogenicity model should we use when {prev:.1%} of missense variants are pathogenic?
> **Takeaway.** On the held-out test split both models have a high AUROC (A {p3(PS['test']['AUROC_A'])}, B
> {p3(PS['test']['AUROC_B'])}), but precision–recall separates them (AP {p3(PS['test']['AP_A'])} vs
> {p3(PS['test']['AP_B'])}, baseline {prev:.3f}); after isotonic calibration model B's Brier score fell from
> {TE['brier_raw_B']:.3f} to {TE['brier_calibrated_B']:.3f}, and a cut-off for 90% precision gave precision
> {p3(TH['test']['precision'])} and recall {p3(TH['test']['recall'])} on new data.
> **Limitation.** Simulated scores with {PS['test']['n_pathogenic']} pathogenic test variants (wide uncertainty); real
> benchmarks must also split by gene and check for circularity with the training labels.
""")
nb.code("""
# Step B7 · save the numbers used on the slides to results/B7.json (every slide number comes from here)
slide = {"prevalence_test": (round(prevalence, 4), "fraction", "fraction pathogenic in the test split (3,200 variants)"),
         "n_pathogenic_test": (int(y.sum()), "variants", "pathogenic variants in the test split"),
         "AUROC_A": (round(auroc["A"], 4), "AUROC", "model A, test split"),
         "AUROC_B": (round(auroc["B"], 4), "AUROC", "model B (raw), test split"),
         "AP_A": (round(ap["A"], 4), "AP", "model A average precision, test split"),
         "AP_B": (round(ap["B"], 4), "AP", "model B average precision, test split"),
         "AUROC_B_train": (float(per_split.loc["train", "AUROC_B"]), "AUROC", "model B on its own training split"),
         "brier_raw_B": (round(brier["raw"], 4), "Brier score", "model B raw, test split"),
         "brier_calibrated_B": (round(brier["calibrated"], 4), "Brier score", "model B after isotonic calibration"),
         "AUROC_B_calibrated": (round(auroc_cal, 4), "AUROC", "model B after isotonic calibration (ties)"),
         "cut_score_B": (round(cut, 4), "raw score B", "cut-off for precision >= 0.90 chosen on the calibration split"),
         "precision_calibration_at_cut": (round(p_cal, 4), "fraction", "precision at the cut-off, calibration split"),
         "precision_test_at_cut": (round(p_test, 4), "fraction", "precision at the cut-off, test split"),
         "recall_test_at_cut": (round(r_test, 4), "fraction", "recall at the cut-off, test split"),
         "n_flagged_test_at_cut": (n_test, "variants", "test variants at or above the cut-off")}
slide_numbers = {k: {"value": v, "unit": u, "meaning": m} for k, (v, u, m) in slide.items()}
tables = {"per_split": per_split.reset_index(names="split").to_dict(orient="records"),
          "calibration_raw": cal_raw.round(6).to_dict(orient="records"),
          "calibration_isotonic": cal_iso.round(6).to_dict(orient="records")}
dviz.save_results("B7", extra={"slide_numbers": slide_numbers, **tables})
print(f"Notebook time: {time.time() - T_START:.0f} s")
""")

out = ROOT / "electives" / "B7_classifier.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
nb.write(str(out))
print("wrote", out.relative_to(ROOT), len(nb.cells), "cells")
