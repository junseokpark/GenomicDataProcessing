#!/usr/bin/env python3
"""Draw the six slide figures of Practice 1 from the notebook's results.

Run from the package folder AFTER the notebook (it reads results/p1_results.json,
results/coverage/depth.tsv and data/provenance/ucsc_hg38_mane_chr11_5200001-5300000.json):

    python scripts/make_slide_figures.py [--font-dir DIR]

The lecture deck uses the SUIT font. If DIR contains SUIT-Regular.ttf and SUIT-Bold.ttf they
are registered with matplotlib; otherwise the default sans-serif font is used.
Output: results/figures/p1_fig_*.png (300 dpi, white background, exact physical size).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
import numpy as np                                   # noqa: E402
import pandas as pd                                  # noqa: E402
from matplotlib import font_manager                  # noqa: E402
from matplotlib.patheffects import AbstractPathEffect  # noqa: E402
from matplotlib.transforms import Affine2D           # noqa: E402

BLUE, DARK_BLUE, INK, GRAY = "#1F5FEA", "#173FAF", "#17202A", "#5B6573"
LIGHT_GRAY, LIGHT_BLUE, PALE, ORANGE = "#8A93A0", "#B7CEF5", "#EEF3FE", "#B25E09"
SAMPLES = ["S01", "S02", "S03"]
ROLE = {"S01": "child (affected)", "S02": "father", "S03": "mother"}
SAMPLE_COLOR = {"S01": BLUE, "S02": DARK_BLUE, "S03": GRAY}
BASE, SMALL = 10, 9          # role-mapped sizes: axis labels / series identity; ticks, legend, notes
OUT = "results/figures"


class Oblique(AbstractPathEffect):
    """Slanted text for gene symbols (the deck font has no italic face)."""

    def __init__(self, angle=11):
        super().__init__()
        self.angle = angle

    def draw_path(self, renderer, gc, tpath, affine, rgbFace=None):
        renderer.draw_path(gc, tpath, Affine2D().skew_deg(self.angle, 0) + affine, rgbFace)


def setup(font_dir):
    family = "sans-serif"
    if font_dir:
        files = [os.path.join(font_dir, f) for f in ("SUIT-Regular.ttf", "SUIT-Bold.ttf")]
        if all(os.path.exists(f) for f in files):
            for f in files:
                font_manager.fontManager.addfont(f)
            family = "SUIT"
    matplotlib.rcParams.update({
        "font.family": [family, "DejaVu Sans"], "axes.unicode_minus": False, "font.size": BASE, "axes.labelsize": BASE, "axes.titlesize": BASE,
        "legend.fontsize": SMALL, "xtick.labelsize": SMALL, "ytick.labelsize": SMALL,
        "axes.linewidth": 0.6, "axes.edgecolor": GRAY, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": GRAY, "ytick.color": GRAY, "xtick.labelcolor": INK, "ytick.labelcolor": INK,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.dpi": 300, "lines.linewidth": 1.2,
        "patch.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    return family


def check_layout(fig, name):
    """Render-then-verify: report overlapping text boxes and text outside the canvas."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip() and t.get_visible()]
    boxes = [(t, t.get_window_extent(r)) for t in texts]
    problems = []
    for i, (a, ba) in enumerate(boxes):
        if not fig.bbox.contains(ba.x0 + 0.5, ba.y0 + 0.5) or not fig.bbox.contains(ba.x1 - 0.5, ba.y1 - 0.5):
            problems.append(f"outside canvas: {a.get_text()!r}")
        for b, bb in boxes[i + 1:]:
            if ba.overlaps(bb):
                problems.append(f"overlap: {a.get_text()!r} / {b.get_text()!r}")
    print(f"{name}: {'OK' if not problems else problems}")
    return problems


def save(fig, name, size):
    path = os.path.join(OUT, name)
    check_layout(fig, name)
    fig.savefig(path, dpi=300)
    plt.close(fig)
    from PIL import Image
    w, h = Image.open(path).size
    assert (w, h) == (round(size[0] * 300), round(size[1] * 300)), (name, w, h)
    return path


def fig_quality(R):
    q1 = np.array(R["quality_curve_S01"]["R1_mean_by_cycle"])
    q2 = np.array(R["quality_curve_S01"]["R2_mean_by_cycle"])
    size = (6.0, 3.2)
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    x = np.arange(1, len(q1) + 1)
    ax.axhline(30, color=LIGHT_GRAY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(4, 29.75, "Q30 = 1 error in 1,000 bases", color=GRAY, fontsize=SMALL, va="top")
    ax.plot(x, q1, color=BLUE, lw=1.4, zorder=3)
    ax.plot(x, q2, color=ORANGE, lw=1.4, zorder=3)
    ax.text(153, q1[-8:].mean(), "R1", color=BLUE, fontsize=BASE, va="center", fontweight="bold")
    ax.text(153, q2[-8:].mean(), "R2", color=ORANGE, fontsize=BASE, va="center", fontweight="bold")
    n_pairs = R["read_pairs"]["S01"]
    ax.text(0.99, 0.97, f"Child S01 · {n_pairs:,} read pairs · 2 × 150 bp", transform=ax.transAxes,
            ha="right", va="top", fontsize=SMALL, color=GRAY)
    ax.set_xlim(0, 165)
    ax.set_xticks([1, 25, 50, 75, 100, 125, 150])
    ax.set_ylim(27, 38)
    ax.set_xlabel("Position in read (sequencing cycle)")
    ax.set_ylabel("Mean base quality (Phred)")
    return save(fig, "p1_fig_quality.png", size)


def mane_genes():
    m = json.load(open("data/provenance/ucsc_hg38_mane_chr11_5200001-5300000.json"))["mane"]
    genes = []
    for r in m:
        starts = [int(v) for v in r["chromStarts"].strip(",").split(",")]
        sizes = [int(v) for v in r["blockSizes"].strip(",").split(",")]
        exons = [(r["chromStart"] + s + 1, r["chromStart"] + s + z) for s, z in zip(starts, sizes)]
        genes.append(dict(name=r["geneName2"], start=r["chromStart"] + 1, end=r["chromEnd"],
                          strand=r["strand"], exons=exons))
    return genes


def fig_depth(R):
    depth = pd.read_csv("results/coverage/depth.tsv", sep="\t")
    depth.columns = ["chrom", "pos"] + SAMPLES
    lc = [int(v) for v in R["coverage_interval"]["low_coverage_interval"].split(":")[1].split("-")]
    bins = depth.groupby((depth.pos - 5_200_001) // 500)[SAMPLES].mean()
    xc = (5_200_001 + bins.index * 500 + 249.5) / 1e6
    size = (12.0, 3.0)
    fig, (ax, gx) = plt.subplots(2, 1, figsize=size, sharex=True, layout="constrained",
                                 gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.02})
    ax.axvspan(lc[0] / 1e6, lc[1] / 1e6, color=ORANGE, alpha=0.18, lw=0, zorder=0)
    ax.axhline(20, color=LIGHT_GRAY, lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.text(1.004, 20, "20×", color=GRAY, fontsize=SMALL, va="center", ha="left",
            transform=ax.get_yaxis_transform())
    for s in SAMPLES:
        ax.plot(xc, bins[s], color=SAMPLE_COLOR[s], lw=1.1, label=f"{s} {ROLE[s]}", zorder=3)
    lowmean = np.mean([R["coverage"][s]["mean_depth_lowcov_interval"] for s in SAMPLES])
    ax.annotate(f"Low-coverage interval\n1 kb, ~{lowmean:.0f}× (planted)",
                xy=(lc[1] / 1e6, 10), xytext=(lc[1] / 1e6 + 0.0012, 4.5), color=ORANGE,
                fontsize=SMALL, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8, shrinkA=0, shrinkB=2))
    ax.set_ylim(0, 54)
    ax.set_yticks([0, 10, 20, 30, 40])
    ax.set_ylabel("Depth (reads)")
    ax.legend(loc="upper left", ncol=3, bbox_to_anchor=(0.0, 1.04), handlelength=1.6,
              columnspacing=1.6, fontsize=SMALL)
    means = ", ".join(f"{s} {R['coverage'][s]['mean_depth']:.0f}×" for s in SAMPLES)
    ax.text(5.2875, 53, f"Mean depth after duplicate marking: {means} · 500-bp bins",
            ha="right", va="top", fontsize=SMALL, color=GRAY)
    ax.tick_params(axis="x", which="both", length=0)
    ax.spines["bottom"].set_visible(False)
    for g in mane_genes():
        a, b = max(g["start"], 5_200_001), min(g["end"], 5_300_000)
        gx.plot([a / 1e6, b / 1e6], [0, 0], color=GRAY, lw=0.8, solid_capstyle="butt")
        for ea, eb in g["exons"]:
            ea, eb = max(ea, 5_200_001), min(eb, 5_300_000)
            if ea < eb:
                gx.add_patch(plt.Rectangle((ea / 1e6, -0.32), (eb - ea) / 1e6, 0.64, color=INK, lw=0))
        xm = (a + b) / 2e6 if g["name"] != "OR51V1" else a / 1e6 + 0.0002
        t = gx.text(xm, 0.55, g["name"], ha="center" if g["name"] != "OR51V1" else "left", va="bottom",
                    fontsize=SMALL, color=INK)
        t.set_path_effects([Oblique()])
    assert all(g["strand"] == "-" for g in mane_genes())
    gx.set_ylim(-0.6, 1.6)
    gx.set_yticks([])
    gx.spines["left"].set_visible(False)
    gx.set_xlim(5.2, 5.3)
    gx.set_xticks(np.round(np.arange(5.20, 5.301, 0.01), 2))
    gx.set_xlabel("chr11 position (GRCh38, Mb)")
    gx.text(-0.005, 0.15, "MANE genes\n(minus strand)", transform=gx.get_yaxis_transform(), ha="right",
            va="center", fontsize=SMALL, color=GRAY, linespacing=1.0)
    return save(fig, "p1_fig_depth.png", size)


def fig_duplicates(R):
    d = R["duplicates"]
    jp = d["jackpot"]
    before, after = jp["counts_before_markdup"], jp["counts_after_markdup"]
    size = (6.0, 3.4)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=size, layout="constrained",
                                 gridspec_kw={"width_ratios": [1.1, 1.0], "wspace": 0.08})
    xs = [0, 1]
    refc = [before.get("G", 0), after.get("G", 0)]
    altc = [before.get("C", 0), after.get("C", 0)]
    a1.bar(xs, refc, width=0.55, color=GRAY, label="G (reference)")
    a1.bar(xs, altc, width=0.55, bottom=refc, color=ORANGE, label="C (copied PCR error)")
    for x, r_, a_ in zip(xs, refc, altc):
        a1.text(x, r_ / 2, str(r_), ha="center", va="center", color="white", fontsize=SMALL)
        if a_ >= 4:
            a1.text(x, r_ + a_ / 2, str(a_), ha="center", va="center", color="white", fontsize=SMALL)
        else:
            a1.text(x + 0.31, r_ + a_ / 2, str(a_), ha="left", va="center", color=ORANGE, fontsize=SMALL)
        a1.text(x, r_ + a_ + 1.2, f"C = {100 * a_ / (r_ + a_):.0f}%", ha="center", va="bottom",
                color=ORANGE, fontsize=SMALL)
    a1.set_xticks(xs, ["Before\nmarkdup", "After\nmarkdup"])
    a1.set_xlim(-0.55, 1.6)
    a1.set_ylim(0, 66)
    a1.set_yticks([0, 10, 20, 30, 40, 50, 60])
    a1.set_ylabel(f"S03 reads at chr11:{jp['pos']:,}")
    a1.legend(loc="upper right", bbox_to_anchor=(1.03, 1.03), fontsize=SMALL, handlelength=1.0)
    pct = [d["duplicate_pct"][s] for s in SAMPLES]
    a2.bar(range(3), pct, width=0.55, color=BLUE)
    for i, p in enumerate(pct):
        a2.text(i, p + 0.08, f"{p:.2f}%", ha="center", va="bottom", fontsize=SMALL)
    extra = pct[2] - np.mean(pct[:2])
    a2.annotate(f"includes the {jp['copies_with_identical_ends'] - 1}\njackpot copies (+{extra:.2f} points)",
                xy=(2, pct[2] + 0.45), xytext=(1.25, 4.95), fontsize=SMALL, color=ORANGE, ha="center",
                va="bottom", arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8, shrinkA=1, shrinkB=1))
    a2.set_xticks(range(3), [f"{s}\n{ROLE[s].split(' ')[0]}" for s in SAMPLES])
    a2.set_ylim(0, 6.2)
    a2.set_yticks([0, 1, 2, 3, 4, 5, 6])
    a2.set_ylabel("Reads flagged as duplicates (%)")
    for ax, letter in ((a1, "a"), (a2, "b")):
        ax.text(-0.02, 1.0, letter, transform=ax.transAxes, fontsize=BASE + 1, fontweight="bold",
                ha="right", va="bottom")
    return save(fig, "p1_fig_duplicates.png", size)


def fig_pileup(R):
    pile = R["hbb_pileup"]
    site = R["hbb_site"]["samples"]
    size = (6.0, 3.4)
    fig, ax = plt.subplots(figsize=size, layout="constrained")
    ys = np.arange(3)[::-1]
    for y, s in zip(ys, SAMPLES):
        t, a = pile[s]["T (ref)"], pile[s]["A (alt)"]
        ax.barh(y, t, height=0.55, color=GRAY)
        ax.barh(y, a, left=t, height=0.55, color=BLUE)
        if t:
            ax.text(t / 2, y, str(t), ha="center", va="center", color="white", fontsize=SMALL)
        ax.text(t + a / 2, y, str(a), ha="center", va="center", color="white", fontsize=SMALL)
        ax.text(t + a + 1.0, y, f"VAF {a / (t + a):.2f}  ·  GT {site[s]['GT']}", ha="left", va="center",
                fontsize=BASE, color=INK)
    ax.set_yticks(ys, [f"{s}\n{ROLE[s]}" for s in SAMPLES])
    ax.set_xlim(0, 56)
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_ylim(-0.6, 3.05)
    ax.set_xlabel("Reads at chr11:5,227,002 (HBB c.20A>T), after duplicate marking")
    h1 = plt.Rectangle((0, 0), 1, 1, color=GRAY)
    h2 = plt.Rectangle((0, 0), 1, 1, color=BLUE)
    ax.legend([h1, h2], ["T  reference allele", "A  sickle (HbS) allele"], loc="upper left",
              bbox_to_anchor=(0.0, 1.0), ncol=2, fontsize=SMALL, handlelength=1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    return save(fig, "p1_fig_pileup_hbb.png", size)


def fig_gl(R):
    g = R["genotype_likelihood_example"]
    gts = ["0/0", "0/1", "1/1"]
    ll = [g["log10_likelihood_by_hand"][k] for k in gts]
    plh = [g["PL_by_hand"][k] for k in gts]
    plb = [g["PL_bcftools"][k] for k in gts]
    size = (6.0, 3.0)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=size, sharey=True, layout="constrained",
                                 gridspec_kw={"width_ratios": [1.0, 1.25], "wspace": 0.04})
    ys = np.arange(3)[::-1]
    a1.barh(ys, ll, height=0.5, color=BLUE)
    for y, v in zip(ys, ll):
        a1.text(v - 1.5, y, f"{v:.1f}", ha="right", va="center", fontsize=SMALL)
    a1.set_xlim(-66, 0)
    a1.set_ylim(-0.5, 2.75)
    a1.set_xticks([-60, -40, -20, 0])
    a1.set_yticks(ys, ["0/0  (T/T)", "0/1  (T/A)", "1/1  (A/A)"])
    a1.set_xlabel("log10 likelihood, by hand")
    a1.text(0.02, 1.04, f"Father S02: {g['n_ref_T']} T + {g['n_alt_A']} A reads, ε = {g['epsilon']}",
            transform=a1.transAxes, ha="left", va="bottom", fontsize=SMALL, color=GRAY)
    h = 0.34
    a2.barh(ys + h / 2, plh, height=h, color=BLUE, label="by hand (binomial)")
    a2.barh(ys - h / 2, plb, height=h, color=GRAY, label="bcftools mpileup")
    for y, v1, v2 in zip(ys, plh, plb):
        for yy, v, c, inside in ((y + h / 2, v1, BLUE, False), (y - h / 2, v2, GRAY, True)):
            if v == 0:
                a2.plot([0], [yy], marker="o", ms=4, color=c, clip_on=False, zorder=4)
                a2.text(12, yy, "0  best", ha="left", va="center", fontsize=SMALL, color=c)
            elif inside:
                a2.text(v - 7, yy, str(v), ha="right", va="center", fontsize=SMALL, color="white")
            else:
                a2.text(v + 8, yy, str(v), ha="left", va="center", fontsize=SMALL, color=c)
    a2.axvline(255, color=LIGHT_GRAY, lw=0.8, ls=(0, (4, 3)))
    a2.text(262, 1.0, "bcftools stores\nPL > 255 as 255", color=GRAY, fontsize=SMALL, ha="left", va="center")
    a2.set_xlim(0, 500)
    a2.set_xticks([0, 100, 200, 300, 400])
    a2.set_xlabel("PL (Phred-scaled; 0 = most likely)")
    a2.legend(loc="lower right", bbox_to_anchor=(1.02, 0.99), ncol=2, fontsize=SMALL, handlelength=1.0,
              columnspacing=1.0)
    a2.tick_params(axis="y", length=0)
    a1.tick_params(axis="y", length=0)
    a1.spines["left"].set_visible(False)
    a2.spines["left"].set_visible(False)
    return save(fig, "p1_fig_genotype_likelihood.png", size)


def fig_benchmark(R):
    b = pd.DataFrame(R["benchmark"])
    b = b[b["sample"] == "all"].set_index("scenario")
    order = ["main pipeline", "no duplicate marking", "no base-quality minimum (-Q 0)", "main pipeline + filters"]
    labels = ["Main pipeline", "No duplicate marking", "No base-quality\nminimum (-Q 0)", "Main pipeline\n+ filters"]
    err = pd.DataFrame(R["benchmark_errors"])
    why = {"ARTIFACT_JACKPOT": "PCR jackpot", "ARTIFACT_LOWBQ": "low base quality", "LOWCOV_HET": "low coverage"}
    size = (6.0, 3.4)
    fig, axes = plt.subplots(1, 2, figsize=size, sharey=True, layout="constrained",
                             gridspec_kw={"wspace": 0.05})
    fig.get_layout_engine().set(rect=(0, 0, 1, 0.92))
    ys = np.arange(len(order))[::-1]
    for ax, metric, kind in ((axes[0], "precision", "FP"), (axes[1], "recall", "FN")):
        for y, sc in zip(ys, order):
            v = b.loc[sc, metric]
            c = BLUE if v >= 0.9995 else ORANGE
            if v < 0.9995:
                ax.plot([v, 1.0], [y, y], color=ORANGE, lw=0.8, zorder=2)
            ax.plot([v], [y], marker="o", ms=6, color=c, zorder=3, clip_on=False)
            if v < 0.9995:
                e = err[(err.scenario == sc) & (err.type == kind)].iloc[0]
                ax.text(v - 0.004, y, f"{v:.3f}\n1 {kind}: {why.get(e.planted_event, e.planted_event)} ({e['sample']})",
                        ha="right", va="center", fontsize=SMALL, color=ORANGE, linespacing=1.1)
            else:
                ax.text(v - 0.004, y + 0.0, f"{v:.3f}", ha="right", va="center", fontsize=SMALL, color=BLUE)
        ax.set_xlim(0.90, 1.004)
        ax.set_xticks([0.90, 0.95, 1.00], ["0.90", "0.95", "1.00"])
        ax.set_xlabel(metric.capitalize())
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.axvline(1.0, color=PALE, lw=6, zorder=0)
    axes[0].set_yticks(ys, labels)
    axes[0].set_ylim(-0.6, len(order) - 0.4)
    n_true = int(b.loc["main pipeline", "TP"] + b.loc["main pipeline", "FN"])
    fig.text(0.995, 0.995, f"Pooled over S01, S02, S03 · {n_true} true variant genotypes · higher = better",
             ha="right", va="top", fontsize=SMALL, color=GRAY)
    return save(fig, "p1_fig_benchmark.png", size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--font-dir", default=os.environ.get("P1_FONT_DIR"))
    args = ap.parse_args()
    family = setup(args.font_dir)
    print("font family:", family)
    os.makedirs(OUT, exist_ok=True)
    R = json.load(open("results/p1_results.json"))
    for f in (fig_quality, fig_depth, fig_duplicates, fig_pileup, fig_gl, fig_benchmark):
        print("wrote", f(R))


if __name__ == "__main__":
    main()
