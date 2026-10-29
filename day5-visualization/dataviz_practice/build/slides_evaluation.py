"""Slide figures for the Evaluation track (B6 variant callers, B7 classifiers) at the exact pixel size of the deck.

Targets (build/picture_frames.tsv; 288 dpi, no tight bbox):
  image-20-1  four-set Venn with ellipses and counts (deliberately hard to read; every count in its true region)
  image-20-2  UpSet (slide 20, shown 6.25 in wide)       image-59-1  UpSet (demo slide, shown 5.28 in wide)
  image-73-1  UpSet with large text (prompt-card thumbnail, shown 3.61 in wide)
  image-46-1, image-60-1  ROC | PR with prevalence baseline | calibration raw vs calibrated
  image-74-1  the same three panels with large text (prompt-card thumbnail, shown 4.44 in wide)
  practice_B6, practice_B7  (1800 x 1210 px) one figure per exercise
All numbers are computed here from data/ with the same methods as the notebooks (B6, B7).
usage:  python build/slides_evaluation.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Ellipse  # noqa: E402
from scipy import ndimage  # noqa: E402
from sklearn.isotonic import IsotonicRegression  # noqa: E402
from sklearn.metrics import (average_precision_score, brier_score_loss, precision_recall_curve,  # noqa: E402
                             roc_auc_score, roc_curve)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
import dviz  # noqa: E402
from slidefig import new_fig, save_exact  # noqa: E402

CALLERS = ["GATK", "DeepVariant", "FreeBayes", "Strelka2"]

# ------------------------------------------------------------------------------------------------ B6 data
vc = pd.read_csv(ROOT / "data" / "variant_callsets.tsv", sep="\t")
vc["n_callers"] = vc[CALLERS].sum(axis=1)
vc["key"] = vc[CALLERS].astype(int).astype(str).agg("".join, axis=1)
counts = dviz.upset_counts(vc, CALLERS)
counts["key"] = counts[CALLERS].astype(int).astype(str).agg("".join, axis=1)
counts["real"] = counts["key"].map(vc.groupby("key")["in_truth_set"].sum()).astype(int)
SIZES = vc[CALLERS].sum()
truth = vc["in_truth_set"].astype(bool)


def prf(mask):
    m = np.asarray(mask, bool)
    tp, fp, fn = int((m & truth).sum()), int((m & ~truth).sum()), int((~m & truth).sum())
    p, r = tp / (tp + fp), tp / (tp + fn)
    return {"TP": tp, "FP": fp, "precision": p, "recall": r, "F1": 2 * p * r / (p + r)}


RULES = {"≥ 1": vc.n_callers >= 1, "≥ 2": vc.n_callers >= 2, "≥ 3": vc.n_callers >= 3, "= 4": vc.n_callers == 4}
CONS = {k: prf(m) for k, m in RULES.items()}
DV_ALONE = prf(vc["DeepVariant"] == 1)

# ------------------------------------------------------------------------------------------------ B7 data
cs = pd.read_csv(ROOT / "data" / "classifier_scores.tsv", sep="\t")
test, calib = cs[cs.split == "test"].copy(), cs[cs.split == "calibration"]
y = test["label_pathogenic"].to_numpy()
PREV = float(y.mean())
iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip").fit(calib.score_modelB, calib.label_pathogenic)
test["score_modelB_calibrated"] = iso.predict(test.score_modelB)
AUROC = {m: roc_auc_score(y, test[f"score_model{m}"]) for m in "AB"}
AP = {m: average_precision_score(y, test[f"score_model{m}"]) for m in "AB"}
ROC = {m: roc_curve(y, test[f"score_model{m}"])[:2] for m in "AB"}
PR = {m: precision_recall_curve(y, test[f"score_model{m}"])[:2] for m in "AB"}
BRIER = {"raw": brier_score_loss(y, test.score_modelB), "calibrated": brier_score_loss(y, test.score_modelB_calibrated)}
CAL = {"raw": dviz.calibration_table(y, test.score_modelB), "calibrated": dviz.calibration_table(y, test.score_modelB_calibrated)}


def threshold_for_precision(labels, scores, target=0.90):
    order = np.argsort(-scores, kind="mergesort")
    hits = np.cumsum(labels[order])
    precision = hits / np.arange(1, len(labels) + 1)
    ok = (precision >= target) & (labels[order] == 1)
    return float(scores[order][np.flatnonzero(ok).max()])


CUT = threshold_for_precision(calib.label_pathogenic.to_numpy(), calib.score_modelB.to_numpy())
called = test.score_modelB.to_numpy() >= CUT
P_TEST, R_TEST = (called & (y == 1)).sum() / called.sum(), (called & (y == 1)).sum() / y.sum()


# ================================================================================================ UpSet
def upset_axes(fig, left, right, bar_bottom, bar_top, dot_bottom, dot_top):
    ax_bar = fig.add_axes([left, bar_bottom, right - left, bar_top - bar_bottom])
    ax_dot = fig.add_axes([left, dot_bottom, right - left, dot_top - dot_bottom], sharex=ax_bar)
    return ax_bar, ax_dot


def draw_upset(ax_bar, ax_dot, top, label_sizes=True, bar_fs=None, dot_s=60, lw=2.0, yticks=None, head=1.12):
    """UpSet like the original slide: bars sorted by size (all four in blue), counts on top, dot matrix below."""
    x = np.arange(len(top))
    all4 = top[CALLERS].all(axis=1).to_numpy()
    ax_bar.bar(x, top["count"], width=0.6, color=np.where(all4, dviz.BLUE, dviz.SLATE))
    fs = bar_fs or plt.rcParams["font.size"] * 0.9
    for xi, c in zip(x, top["count"]):
        ax_bar.text(xi, c, f"{c:,}", ha="center", va="bottom", fontsize=fs, color=dviz.INK)
    ax_bar.set_ylabel("Variants")
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.spines["bottom"].set_visible(False)
    if yticks is not None:
        ax_bar.set_yticks(yticks)
    ax_bar.set_ylim(0, top["count"].max() * head)
    rows = list(reversed(CALLERS))                                   # Strelka2 on top, GATK at the bottom
    for xi, (_, r) in zip(x, top.iterrows()):
        on = [j for j, c in enumerate(rows) if bool(r[c])]
        ax_dot.scatter([xi] * len(rows), range(len(rows)), s=dot_s, color=dviz.GREY, zorder=1, linewidths=0)
        ax_dot.scatter([xi] * len(on), on, s=dot_s, color=dviz.INK, zorder=3, linewidths=0)
        if len(on) > 1:
            ax_dot.plot([xi, xi], [min(on), max(on)], color=dviz.INK, lw=lw, zorder=2)
    labels = [f"{c} ({SIZES[c]:,})" if label_sizes else c for c in rows]
    ax_dot.set_yticks(range(len(rows)), labels)
    ax_dot.tick_params(axis="y", left=False)
    ax_dot.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["left", "bottom"]:
        ax_dot.spines[sp].set_visible(False)
    ax_dot.set_ylim(-0.6, len(rows) - 0.4)
    ax_bar.set_xlim(-0.6, len(top) - 0.4)


def slide_upset(target, base, n_bars=10, label_sizes=True, title="UpSet: each intersection is one bar",
                left=0.30, bar_fs=None, dot_s=60, lw=2.0, yticks=None, head=1.12, title_x=None):
    fig, ax0 = new_fig(target, base_size=base)
    ax0.remove()
    ax_bar, ax_dot = upset_axes(fig, left, 0.985, 0.43, 0.86, 0.03, 0.39)
    draw_upset(ax_bar, ax_dot, counts.head(n_bars), label_sizes=label_sizes, bar_fs=bar_fs, dot_s=dot_s, lw=lw,
               yticks=yticks, head=head)
    fig.text(left - 0.005 if title_x is None else title_x, 0.975, title, ha="left", va="top", fontsize=base * 1.15,
             fontweight="bold", color=dviz.INK)
    return save_exact(fig, target)


# ================================================================================================ Venn
VENN_ELLIPSES = [(0.350, 0.400, 0.72, 0.45, 140.0), (0.450, 0.500, 0.72, 0.45, 140.0),
                 (0.544, 0.500, 0.72, 0.45, 40.0), (0.644, 0.400, 0.72, 0.45, 40.0)]
VENN_COLORS = [dviz.ORANGE, dviz.OI_SKY, dviz.GREEN, dviz.PURPLE]


def inside(e, X, Y):
    cx, cy, w, h, ang = e
    t = np.deg2rad(ang)
    dx, dy = X - cx, Y - cy
    u, v = dx * np.cos(t) + dy * np.sin(t), -dx * np.sin(t) + dy * np.cos(t)
    return (u / (w / 2)) ** 2 + (v / (h / 2)) ** 2 <= 1


def slide_venn(target="image-20-1.png", base=12.5):
    fig, ax = new_fig(target, base_size=base)
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.0, top=0.86)
    step = 0.001
    xs, ys = np.arange(-0.02, 1.02, step), np.arange(0.0, 0.92, step)
    X, Y = np.meshgrid(xs, ys)
    member = np.stack([inside(e, X, Y) for e in VENN_ELLIPSES])         # 4 x ny x nx
    code = (member * (2 ** np.arange(3, -1, -1))[:, None, None]).sum(axis=0)   # bit 3 = first caller
    for e, col in zip(VENN_ELLIPSES, VENN_COLORS):
        cx, cy, w, h, ang = e
        ax.add_patch(Ellipse((cx, cy), w, h, angle=ang, facecolor=col, alpha=0.28, edgecolor="none"))
        ax.add_patch(Ellipse((cx, cy), w, h, angle=ang, facecolor="none", edgecolor=col, lw=1.6, alpha=0.9))
    n_by_key = dict(zip(counts["key"], counts["count"]))
    placed = {}
    for k in range(1, 16):
        key = format(k, "04b")                                          # e.g. '0100' = DeepVariant only
        mask = code == k
        assert mask.any(), f"region {key} missing from the ellipse layout"
        dist = ndimage.distance_transform_edt(mask)
        iy, ix = np.unravel_index(np.argmax(dist), dist.shape)          # deepest point inside the region
        placed[key] = (xs[ix], ys[iy], dist[iy, ix] * step)
        ax.text(xs[ix], ys[iy], f"{n_by_key.get(key, 0):,}", ha="center", va="center", fontsize=base * 0.82,
                color=dviz.INK)
    centre = np.array([0.497, 0.45])
    for (cx, cy, w, h, ang), name, col in zip(VENN_ELLIPSES, CALLERS, VENN_COLORS):
        d = np.array([np.cos(np.deg2rad(ang)), np.sin(np.deg2rad(ang))])
        tips = [np.array([cx, cy]) + s * (w / 2) * d for s in (1, -1)]
        tip = max(tips, key=lambda p: np.linalg.norm(p - centre))       # the tip that points outwards
        out = (tip - centre) / np.linalg.norm(tip - centre)
        pos = tip + 0.035 * out
        ax.text(pos[0], pos[1], name, ha="right" if out[0] < 0 else "left", va="bottom", color=col,
                fontsize=base * 1.0, fontweight="bold")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(0.0, 0.92)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.text(0.04, 0.975, "Four-set Venn", ha="left", va="top", fontsize=base * 1.3, fontweight="bold", color=dviz.INK)
    return save_exact(fig, target), placed


# ================================================================================================ ROC | PR | calibration
def three_panels(target, base, short=False, lw=1.8, ms=5):
    fig, axs = new_fig(target, 1, 3, base_size=base)
    if short:
        fig.subplots_adjust(left=0.085, right=0.985, bottom=0.30, top=0.80, wspace=0.45)
    else:
        fig.subplots_adjust(left=0.055, right=0.985, bottom=0.16, top=0.88, wspace=0.30)
    ax_roc, ax_pr, ax_cal = axs
    lab = (lambda m, v: f"{m} {v:.2f}") if short else (lambda m, v: f"Model {m}  {v:.2f}")
    ax_roc.plot([0, 1], [0, 1], ls="--", color=dviz.GREY, lw=lw * 0.8)
    for m, col in [("A", dviz.SLATE), ("B", dviz.BLUE)]:
        ax_roc.plot(*ROC[m], color=col, lw=lw, label=lab(m, AUROC[m]) if short else f"Model {m}  AUROC {AUROC[m]:.2f}")
        rec, prec = PR[m][1], PR[m][0]
        ax_pr.plot(rec, prec, color=col, lw=lw, drawstyle="steps-post",
                   label=lab(m, AP[m]) if short else f"Model {m}  AP {AP[m]:.2f}")
    ax_pr.axhline(PREV, ls="--", color=dviz.GREY, lw=lw * 0.8)
    ax_pr.text(0.36 if short else 0.02, PREV + 0.03, f"baseline {PREV:.1%}" if short else f"baseline = prevalence {PREV:.1%}",
               color=dviz.SLATE, fontsize=base * (0.6 if short else 0.85))
    ax_cal.plot([0, 1], [0, 1], ls="--", color=dviz.GREY, lw=lw * 0.8)
    for name, col in [("raw", dviz.BLUE), ("calibrated", dviz.GREEN)]:
        t = CAL[name]
        ax_cal.plot(t["mean_predicted"], t["observed"], "-o", color=col, lw=lw, ms=ms,
                    label=name if short else f"{name} (Brier {BRIER[name]:.3f})")
    titles = ["ROC (AUROC)", "PR (AP)", "Calibration"] if short else ["ROC", "Precision–recall", "Calibration"]
    xl = ["FPR", "Recall", "Predicted"] if short else ["False positive rate", "Recall", "Predicted probability"]
    yl = ["TPR", "Precision", "Observed"] if short else ["True positive rate", "Precision", "Observed fraction pathogenic"]
    for ax, ti, a, b in zip(axs, titles, xl, yl):
        ax.set(title=ti, xlabel=a, ylabel=b, xlim=(0, 1.0), ylim=(0, 1.03))
        if short:
            ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
            ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
        else:
            ax.set_xticks(np.linspace(0, 1, 6))
            ax.set_yticks(np.linspace(0, 1, 6))
    ax_cal.set_xlim(-0.01, 1.03)
    leg_fs = base * (0.75 if short else 0.9)
    kw = dict(fontsize=leg_fs, handlelength=1.0 if short else 2, labelspacing=0.3 if short else 0.5)
    ax_roc.legend(loc="lower right", **kw)
    ax_pr.legend(loc="lower left", bbox_to_anchor=(0.0, 0.08 if short else 0.12), **kw)
    cal_kw = dict(kw, fontsize=leg_fs if short else base * 0.8, handlelength=1.0 if short else 1.4)
    ax_cal.legend(loc="upper left", title=None if short else "Model B", title_fontsize=cal_kw["fontsize"], **cal_kw)
    return save_exact(fig, target)


# ================================================================================================ practice figures
def practice_b6(base=10.5):
    target = "practice_B6.png"
    fig, ax0 = new_fig(target, base_size=base)
    ax0.remove()
    top = counts.head(8)
    ax_bar, ax_dot = upset_axes(fig, 0.175, 0.585, 0.52, 0.85, 0.20, 0.46)
    draw_upset(ax_bar, ax_dot, top, label_sizes=False, bar_fs=base * 0.8, dot_s=34, lw=1.4, yticks=[0, 2000, 4000])
    ax_bar.set_ylabel("Variants called")
    n_false4 = int(top.iloc[0]["count"] - top.iloc[0]["real"])
    ax_bar.annotate(f"{n_false4} false", (0, top.iloc[0]["count"]), xytext=(14, -12), textcoords="offset points",
                    fontsize=base * 0.8, color=dviz.VERMILION, ha="left")
    for xi, (_, r) in enumerate(top.iterrows()):                       # % of calls in the truth set, per bar
        share = 100 * r["real"] / r["count"]
        ax_dot.text(xi, -1.25, f"{share:.0f}%", ha="center", va="center", fontsize=base * 0.78,
                    color=dviz.VERMILION if share < 50 else dviz.SLATE, fontweight="bold")
    ax_dot.text(-0.75, -1.25, "% real", ha="right", va="center", fontsize=base * 0.8, color=dviz.INK)
    ax_dot.set_ylim(-1.7, len(CALLERS) - 0.4)
    fig.text(0.015, 0.975, "A  Agreement is not truth: the 8 largest intersections", ha="left", va="top",
             fontsize=base * 1.05, fontweight="bold", color=dviz.INK)
    ax = fig.add_axes([0.715, 0.17, 0.27, 0.66])
    names = list(CONS)
    xx = np.arange(len(names))
    ax.plot(xx, [CONS[n]["precision"] for n in names], "-o", color=dviz.BLUE, ms=4, label="precision")
    ax.plot(xx, [CONS[n]["recall"] for n in names], "-o", color=dviz.SLATE, ms=4, label="recall")
    ax.axhline(DV_ALONE["precision"], color=dviz.BLUE, ls=":", lw=1)
    ax.axhline(DV_ALONE["recall"], color=dviz.SLATE, ls=":", lw=1)
    for xi, n in zip(xx, names):
        p, r = CONS[n]["precision"], CONS[n]["recall"]
        for v, col, dy in [(p, dviz.BLUE, 5 if p >= r else -12), (r, dviz.SLATE, -12 if p >= r else 5)]:
            ax.annotate(f"{v:.3f}", (xi, v), xytext=(0, dy), textcoords="offset points", ha="center",
                        fontsize=base * 0.7, color=col, zorder=5,
                        bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none", alpha=0.85))
    ax.set_xticks(xx, names)
    ax.set(ylim=(0.88, 1.012), xlim=(-0.4, 3.4), xlabel="callers that agree", ylabel="fraction")
    ax.plot([], [], ls=":", color=dviz.SLATE, lw=1, label="DeepVariant alone")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 0.0), fontsize=base * 0.72, handlelength=1.4)
    fig.text(0.635, 0.975, "B  Consensus voting", ha="left", va="top", fontsize=base * 1.05, fontweight="bold",
             color=dviz.INK)
    fig.text(0.015, 0.03, f"Simulated GIAB-style benchmark, chr20 10–16 Mb (hg38), {int(truth.sum()):,} truth variants",
             fontsize=base * 0.72, color=dviz.SLATE)
    return save_exact(fig, target)


def practice_b7(base=10.5):
    target = "practice_B7.png"
    fig, (ax_pr, ax_cal) = new_fig(target, 1, 2, base_size=base)
    fig.subplots_adjust(left=0.10, right=0.985, bottom=0.19, top=0.90, wspace=0.32)
    for m, col in [("A", dviz.SLATE), ("B", dviz.BLUE)]:
        ax_pr.plot(PR[m][1], PR[m][0], color=col, lw=1.6, drawstyle="steps-post",
                   label=f"{m}: AUROC {AUROC[m]:.2f} · AP {AP[m]:.2f}")
    ax_pr.axhline(PREV, ls="--", color=dviz.GREY, lw=1.1)
    ax_pr.text(0.02, PREV + 0.025, f"baseline = prevalence {PREV:.1%}", fontsize=base * 0.75, color=dviz.SLATE)
    ax_pr.plot(R_TEST, P_TEST, "o", color=dviz.INK, ms=5.5, zorder=4, ls="none",
               label=f"cut-off {CUT:.2f}: precision {P_TEST:.2f}, recall {R_TEST:.2f}")
    ax_pr.set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1.03), title="A  Precision–recall, test split")
    ax_pr.legend(loc="lower left", bbox_to_anchor=(0.0, 0.07), fontsize=base * 0.7, handlelength=1.2)
    ax_cal.plot([0, 1], [0, 1], ls="--", color=dviz.GREY, lw=1.1)
    for name, col in [("raw", dviz.BLUE), ("calibrated", dviz.GREEN)]:
        t = CAL[name]
        ax_cal.plot(t["mean_predicted"], t["observed"], "-o", color=col, lw=1.6, ms=4,
                    label=f"model B {name} (Brier {BRIER[name]:.3f})")
    ax_cal.set(xlabel="Predicted probability", ylabel="Observed fraction pathogenic", xlim=(-0.01, 1.03), ylim=(0, 1.03),
               title="B  Calibration of model B")
    ax_cal.legend(loc="upper left", fontsize=base * 0.75, handlelength=1.4)
    fig.text(0.015, 0.012, f"Simulated scores; test split n = {len(y):,} ({int(y.sum())} pathogenic); isotonic "
             f"regression fitted on the calibration split", fontsize=base * 0.72, color=dviz.SLATE)
    return save_exact(fig, target)


if __name__ == "__main__":
    out = []
    venn_path, placed = slide_venn()
    out.append(venn_path)
    out.append(slide_upset("image-20-2.png", base=12.5))
    out.append(slide_upset("image-59-1.png", base=14.0, left=0.33, dot_s=72))
    out.append(slide_upset("image-73-1.png", base=21.0, label_sizes=False, title="UpSet: each intersection is one bar",
                           left=0.25, bar_fs=16, dot_s=170, lw=3.2, yticks=[0, 2000, 4000], head=1.3, title_x=0.03))
    out.append(three_panels("image-46-1.png", base=12.5))
    out.append(three_panels("image-60-1.png", base=13.5))
    out.append(three_panels("image-74-1.png", base=23.0, short=True, lw=3.2, ms=8))
    out.append(practice_b6())
    out.append(practice_b7())
    for p in out:
        print("wrote", p.relative_to(ROOT))
    print("Venn label clearance (data units):", {k: round(v[2], 3) for k, v in sorted(placed.items(), key=lambda t: t[1][2])[:4]})
    print(f"numbers: all four {int(counts.iloc[0]['count']):,} ({int(counts.iloc[0]['count'] - counts.iloc[0]['real'])} false);"
          f" DV only {int(counts.loc[counts.key == '0100', 'count'].iloc[0])}; AUROC A {AUROC['A']:.4f} B {AUROC['B']:.4f};"
          f" AP A {AP['A']:.4f} B {AP['B']:.4f}; Brier {BRIER['raw']:.4f} -> {BRIER['calibrated']:.4f};"
          f" cut {CUT:.4f}: P {P_TEST:.4f} R {R_TEST:.4f}")
