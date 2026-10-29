"""Slide figures for the Tracks electives: B5 (copy number on one coordinate axis) and B8 (TP53 lollipop).

Reads data/ and repeats the notebooks' analysis steps; results/B5.json and results/B8.json are read only to
assert that every number drawn on a slide equals the executed notebook. Writes build/slide_figs/<target>.png at
the exact pixel size of the slide frame (288 dpi, build/picture_frames.tsv):
    image-43-1.png, image-58-1.png  tracks: depth, log2 ratio + segments, genes (real MANE names), LCR22 repeats
    image-72-1.png                  the same tracks, simplified, with larger text (thumbnail on prompt card 6)
    image-44-1.png                  genome-wide 1-Mb view: patient row and tumour row
    image-45-1.png                  TP53 lollipop
    practice_B5.png, practice_B8.png  practice-slide figures (1800 x 1210 px)
usage: python build/slides_tracks.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
import dviz  # noqa: E402
import slidefig as sf  # noqa: E402

DATA = ROOT / "data"
R5 = json.loads((ROOT / "results" / "B5.json").read_text())["slide_numbers"]
R8 = json.loads((ROOT / "results" / "B8.json").read_text())["slide_numbers"]
MINUS = "\u2212"


def num(x, fmt):
    s = fmt.format(x)
    return MINUS + s[1:] if s.startswith("-") else s


def sci(p):
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10{str(int(e)).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"


# ================================================================================================
# B5: the same steps as electives/B5_cnv_tracks.ipynb
# ================================================================================================
cn = pd.read_csv(DATA / "cnv_chr22_10kb.tsv", sep="\t")
cn["mappable"] = cn["mappable"].astype(str).eq("True")
ok = cn["mappable"]
ratio = (cn["depth_patient"] / cn.loc[ok, "depth_patient"].median()) / (
    cn["depth_reference"] / cn.loc[ok, "depth_reference"].median())
cn["log2"] = np.log2(ratio.where(ok))
seg = dviz.segment(cn["log2"], min_size=5, z=6)
seg["start"] = cn["start"].to_numpy()[seg["first"]]
seg["end"] = cn["end"].to_numpy()[seg["last"]]
seg["copy_number"] = 2 * 2 ** seg["mean"]
dele = seg.loc[seg["copy_number"] < 1.5].iloc[0]
D0, D1 = int(dele["start"]), int(dele["end"])
genes = pd.read_csv(DATA / "cnv_genes_chr22.tsv", sep="\t")
sd = pd.read_csv(DATA / "cnv_segdups_chr22.tsv", sep="\t")
lcr = sd.loc[sd["lcr22"].notna()].groupby("lcr22").agg(start=("start", "min"), end=("end", "max"))
truth = pd.read_csv(DATA / "truth" / "cnv_truth.tsv", sep="\t")
T = truth.loc[truth["sample"] == "patient"].iloc[0]
gw = pd.read_csv(DATA / "cnv_genome_1mb.tsv", sep="\t")
gx, ticks = dviz.manhattan_x(gw["chrom"], (gw["start"] + gw["end"]) / 2, gap=1e7)
gw["x"] = gx / 1e6
for col in ["log2_patient", "log2_tumour"]:
    gw[col + "_seg"] = np.nan
    for chrom, rows in gw.groupby("chrom", sort=False).groups.items():
        for _, s in dviz.segment(gw.loc[rows, col], min_size=2, z=6).iterrows():
            gw.loc[rows[int(s["first"]):int(s["last"]) + 1], col + "_seg"] = s["mean"]
p22 = gw.loc[gw["chrom"] == "chr22"]
N_LOW = int((p22["log2_patient"] < -0.5).sum())
away = ok & ((cn["end"] <= 18_000_000) | (cn["start"] >= 22_000_000))
inside = genes.loc[(genes["start"] < D1) & (genes["end"] > D0)]

# every number drawn below must equal the executed notebook
assert D0 + 1 == R5["del_start_1based"]["value"] and D1 == R5["del_end"]["value"]
assert abs(float(dele["mean"]) - json.loads((ROOT / "results" / "B5.json").read_text())["B5.mean_log2_del"]["value"]) < 1e-9
assert round(float(dele["copy_number"]), 2) == R5["copy_number"]["value"]
assert N_LOW == R5["bins_1mb_low"]["value"] and len(inside) == R5["n_genes_in_call"]["value"]
assert int(T["start"]) + 1 == R5["truth_start_1based"]["value"]
MEAN_DEL = float(dele["mean"])
LEN_MB = (D1 - D0) / 1e6


def state(m):
    return np.where(m > 0.15, "gain", np.where(m < -0.15, "loss", "none"))


STATE_COL = {"gain": dviz.VERMILION, "loss": dviz.OI_BLUE, "none": dviz.GREY}


def draw_depth(ax, w, fs):
    x = (w["start"] + 5_000) / 1e6
    ax.fill_between(x, w["depth_patient"].where(w["mappable"]), step="mid", color=dviz.SKY, lw=0)
    ax.fill_between(x, w["depth_patient"].where(~w["mappable"]), step="mid", color=dviz.GREY, lw=0)
    ax.set_ylim(0, 52)
    ax.set_yticks([0, 20, 40])


def draw_log2(ax, w, lw_seg=2.2, s=3, label=True, fs=9):
    x = (w["start"] + 5_000) / 1e6
    ax.scatter(x, w["log2"], s=s, color=dviz.SLATE, lw=0, rasterized=True)
    for _, r in seg.iterrows():
        is_del = r["copy_number"] < 1.5
        ax.hlines(r["mean"], r["start"] / 1e6, r["end"] / 1e6, color=dviz.BLUE if is_del else dviz.INK,
                  lw=lw_seg * (1.5 if is_del else 1.0), zorder=3)
    ax.axhline(0, color=dviz.GREY, lw=0.8, zorder=1)
    ax.axhline(-1, color=dviz.SLATE, lw=0.9, ls=(0, (4, 3)), zorder=1)
    ax.set_ylim(-2.2, 1.1)
    ax.set_yticks([-2, -1, 0, 1])
    if label:
        ax.text(22.97, -1.07, f"1 copy expected: log2(1/2) = {MINUS}1", ha="right", va="top", fontsize=fs,
                color=dviz.SLATE)
        ax.text((D0 + D1) / 2e6, MEAN_DEL - 0.62,
                f"DEL · 1 copy · mean log2 {num(MEAN_DEL, '{:.2f}')} · {LEN_MB:.2f} Mb",
                ha="center", va="top", fontsize=fs, color=dviz.BLUE, fontweight="bold")


def draw_genes(ax, key, fs, y_lab=(0.30, 0.70)):
    g = genes.loc[(genes["end"] > 17_000_000) & (genes["start"] < 23_000_000)]
    ax.broken_barh(list(zip(g["start"] / 1e6, (g["end"] - g["start"]) / 1e6)), (0.42, 0.16), color=dviz.SLATE, lw=0)
    lab = g.loc[g["gene"].isin(key)].sort_values("start")
    for j, (_, r) in enumerate(lab.iterrows()):
        up = j % 2 == 1
        ax.text((r["start"] + r["end"]) / 2e6, y_lab[1] if up else y_lab[0], r["gene"], ha="center",
                va="bottom" if up else "top", fontsize=fs, style="italic", color=dviz.INK)
    ax.set_ylim(-0.25, 1.25)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def draw_lcr(ax, fs, truth_marks=False):
    hi = sd.loc[sd["frac_match"] >= 0.98]
    ax.broken_barh(list(zip(hi["start"] / 1e6, (hi["end"] - hi["start"]) / 1e6)), (0.1, 0.4), color=dviz.GREY, lw=0)
    for j, (name, r) in enumerate(lcr.iterrows()):
        ax.text((r["start"] + r["end"]) / 2e6, 0.58 if j % 2 == 0 else 1.0, f"LCR22{name}", ha="center",
                va="bottom", fontsize=fs, color=dviz.SLATE)
    if truth_marks:
        for pos in (T["start"], T["end"]):
            ax.plot(pos / 1e6, 0.3, marker="v", ms=6, color=dviz.INK, zorder=4)
    ax.set_ylim(0, 1.55)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def shade(axs):
    for ax in axs:
        ax.axvspan(D0 / 1e6, D1 / 1e6, color=dviz.HIGHLIGHT, zorder=0)


def side_label(ax, text, fs):
    ax.set_ylabel(text, rotation=0, ha="right", va="center", fontsize=fs, labelpad=8)


def tracks_figure(target):
    """Lecture slides 43 and 58: four aligned tracks (same layout idea as the original image)."""
    fig, axs = sf.new_fig(target, 4, 1, sharex=True, gridspec_kw={"height_ratios": [1.0, 1.55, 0.62, 0.42]})
    fs = 10.5
    w = cn.loc[(cn["start"] >= 17_000_000) & (cn["end"] <= 23_000_000)]
    draw_depth(axs[0], w, fs)
    draw_log2(axs[1], w, fs=fs)
    draw_genes(axs[2], ["USP18", "PRODH", "HIRA", "TBX1", "COMT", "DGCR8", "CRKL", "MAPK1"], fs)
    draw_lcr(axs[3], fs)
    shade(axs)
    axs[0].set_ylabel("Depth (x)")
    axs[1].set_ylabel("log2 ratio")
    side_label(axs[2], "Genes", 11)
    side_label(axs[3], "Repeats", 11)
    axs[0].set_title(f"chr22 · 17–23 Mb · 10-kb bins · patient vs panel of 10 normals (simulated)   "
                     f"shaded: called deletion {D0 / 1e6:.2f}–{D1 / 1e6:.2f} Mb", fontsize=13)
    axs[3].set_xlim(17, 23)
    axs[3].set_xlabel("chr22 position (Mb, hg38) · grey depth = masked bins (repeats, gaps) · "
                      "segments: binary segmentation")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.9, bottom=0.15, hspace=0.18)
    return sf.save_exact(fig, target)


def tracks_thumbnail(target):
    """Prompt card 6 (small frame): depth + log2 + three gene names, text about twice as large."""
    fig, axs = sf.new_fig(target, 3, 1, base_size=24, sharex=True, gridspec_kw={"height_ratios": [0.9, 1.6, 0.6]})
    w = cn.loc[(cn["start"] >= 17_000_000) & (cn["end"] <= 23_000_000)]
    draw_depth(axs[0], w, 20)
    draw_log2(axs[1], w, lw_seg=4, s=9, label=False)
    axs[1].text(22.95, -1.12, f"{MINUS}1 = one copy", ha="right", va="top", fontsize=20, color=dviz.SLATE)
    draw_genes(axs[2], [], 20)
    for gname in ["TBX1", "CRKL"]:                                    # two labels, both above the gene boxes
        g = genes.loc[genes["gene"] == gname].iloc[0]
        axs[2].text((g["start"] + g["end"]) / 2e6, 0.7, gname, ha="center", va="bottom", fontsize=20,
                    style="italic", color=dviz.INK)
    axs[2].set_ylim(0.1, 1.9)
    shade(axs)
    axs[0].set_ylabel("Depth")
    axs[0].set_yticks([20, 40])
    axs[1].set_ylabel("log2")
    axs[1].set_yticks([-2, -1, 0])
    axs[2].set_xlim(17, 23)
    axs[2].set_xticks([17, 18, 19, 20, 21, 22, 23])
    axs[2].set_xlabel("chr22 position (Mb, hg38)")
    fig.subplots_adjust(left=0.085, right=0.985, top=0.97, bottom=0.25, hspace=0.18)
    return sf.save_exact(fig, target)


def genome_figure(target):
    """Lecture slide 44: one genome axis, two rows (patient, tumour), colour = segment state."""
    fig, axs = sf.new_fig(target, 2, 1, sharex=True, sharey=True)
    titles = ["Patient · constitutional DNA (22q11.2 deletion)", "Colorectal tumour · purity 0.6 (another person)"]
    for ax, col, title in zip(axs, ["log2_patient", "log2_tumour"], titles):
        st = state(gw[col + "_seg"])
        for name, c in STATE_COL.items():
            m = st == name
            ax.scatter(gw.loc[m, "x"], gw.loc[m, col], s=4 if name == "none" else 6, color=c, lw=0, zorder=2)
        for b in gw.groupby("chrom", sort=False)["x"].min().iloc[1:] - 5:
            ax.axvline(b, color=dviz.LIGHT, lw=0.8, zorder=0)
        ax.axhline(0, color=dviz.SLATE, lw=0.6, zorder=1)
        ax.text(0.003, 0.97, title, transform=ax.transAxes, ha="left", va="top", fontsize=11, color=dviz.INK)
        ax.set_ylabel("log2 ratio")
    axs[0].set_ylim(-1.35, 1.05)
    axs[0].set_yticks([-1, 0, 1])
    low = gw.loc[(gw["chrom"] == "chr22") & (gw["log2_patient"] < -0.5)]
    axs[0].annotate(f"22q11.2 deletion = {len(low)} full 1-Mb bins", xy=(low["x"].mean() - 8, -0.95),
                    xytext=(low["x"].mean() - 70, -0.95), ha="right", va="center", fontsize=10.5,
                    arrowprops={"arrowstyle": "->", "color": dviz.INK, "lw": 0.8})
    q8 = gw.loc[(gw["chrom"] == "chr8") & (gw["start"] >= 50_000_000)]
    axs[1].text(q8["x"].mean(), 0.78, "8q gain", ha="center", fontsize=10.5, color=dviz.VERMILION)
    p17 = gw.loc[(gw["chrom"] == "chr17") & (gw["log2_tumour_seg"] < -0.15)]
    axs[1].text(p17["x"].mean(), -1.0, "17p loss (TP53)", ha="center", fontsize=10.5, color=dviz.OI_BLUE)
    axs[1].set_xticks(ticks["centre"] / 1e6, [c.replace("chr", "") if c not in ("chr19", "chr21") else ""
                                              for c in ticks["chrom"]], fontsize=10)
    axs[1].set_xlabel("Chromosome (hg38, 1-Mb bins) · colour = segment state: gain (vermilion), loss (blue)")
    axs[1].set_xlim(gw["x"].min() - 10, gw["x"].max() + 10)
    fig.subplots_adjust(left=0.065, right=0.99, top=0.97, bottom=0.2, hspace=0.12)
    return sf.save_exact(fig, target)


# ================================================================================================
# B8: the same steps as electives/B8_tp53_lollipop.ipynb
# ================================================================================================
var = pd.read_csv(DATA / "tp53_variants.tsv", sep="\t")
dom = pd.read_csv(DATA / "tp53_domains.tsv", sep="\t")
LENGTH = 393
dbd = dom.loc[dom["short"] == "DBD"].iloc[0]
P_DBD = (dbd["aa_end"] - dbd["aa_start"] + 1) / LENGTH
var["in_dbd"] = var["position"].between(dbd["aa_start"], dbd["aa_end"])
per_pos = var.groupby(["position", "variant_class"]).size().rename("n").reset_index()
mine = var.groupby("position").size()
top = var.groupby(["position", "protein_change"]).size().rename("n").reset_index()
top = top.sort_values(["position", "n", "protein_change"], ascending=[True, False, True])
HOT = {}
for pos in mine[mine >= 30].index:
    a = top.loc[top["position"] == pos]
    keep = a.loc[a["n"] >= 10, "protein_change"].tolist() or a["protein_change"].tolist()[:1]
    HOT[pos] = keep[0][:-1] + "/".join(c[-1] for c in keep)
TESTS = {}
for cls in ["missense", "truncating"]:
    sub = var.loc[var["variant_class"] == cls]
    k, n = int(sub["in_dbd"].sum()), len(sub)
    TESTS[cls] = {"k": k, "n": n, "f": k / n, "p": float(stats.binomtest(k, n, P_DBD).pvalue)}
ODDS, P_FISHER = stats.fisher_exact([[TESTS["missense"]["k"], TESTS["missense"]["n"] - TESTS["missense"]["k"]],
                                     [TESTS["truncating"]["k"], TESTS["truncating"]["n"] - TESTS["truncating"]["k"]]])
assert TESTS["missense"]["p"] == R8["missense_binom_p"]["value"] and TESTS["truncating"]["p"] == R8["truncating_binom_p"]["value"]
assert float(P_FISHER) == R8["fisher_p"]["value"] and len(var) == R8["n_tumours"]["value"]
assert [HOT[p] for p in mine[mine >= 30].sort_values(ascending=False, kind="stable").index] == \
    [r["label"] for r in R8["top_residues"]["value"]]
CLASS_COL = {"missense": dviz.OI_BLUE, "truncating": dviz.VERMILION, "splice/other": dviz.GREEN}


def draw_lollipop(ax, fs, size_k=2.4, label_fs=None, legend=True, nudge=None):
    label_fs = label_fs or fs
    nudge = nudge or {}
    ax.axvspan(dbd["aa_start"] - 0.5, dbd["aa_end"] + 0.5, color=dviz.HIGHLIGHT, zorder=0)
    for _, r in per_pos.sort_values("n", ascending=False).iterrows():
        ax.vlines(r["position"], 0, r["n"], color=dviz.GREY, lw=0.7, zorder=1)
        ax.scatter(r["position"], r["n"], s=10 + size_k * r["n"], color=CLASS_COL[r["variant_class"]],
                   edgecolor="white", lw=0.5, zorder=2)
    ax.hlines(-7.5, 1, LENGTH, color=dviz.GREY, lw=4, zorder=1)
    for _, d in dom.iterrows():
        ax.add_patch(plt.Rectangle((d["aa_start"] - 0.5, -12.5), d["length_aa"], 10, color=dviz.SLATE, lw=0, zorder=2))
        ax.text((d["aa_start"] + d["aa_end"]) / 2, -7.5, d["domain"] if d["length_aa"] > 60 else d["short"],
                ha="center", va="center", fontsize=fs * 0.9, color="white", zorder=3)
    for pos, lab in HOT.items():
        dx, dy, ha = nudge.get(int(pos), (0, 8, "center"))          # crowded labels can be moved sideways
        ax.annotate(lab, (pos, mine[pos]), xytext=(dx, dy), textcoords="offset points", ha=ha,
                    va="bottom" if dy > 0 else "center", fontsize=label_fs, color=dviz.INK)
    if legend:
        for cls, c in CLASS_COL.items():
            ax.scatter([], [], color=c, s=40, label=f"{cls.capitalize()} ({int((var['variant_class'] == cls).sum())})")
        ax.legend(loc="upper right", fontsize=fs, title="Variant class (n tumours)", title_fontsize=fs,
                  borderaxespad=0.2)
    ax.set_xlim(0, LENGTH + 5)
    ax.set_ylim(-15, mine.max() * 1.22)
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.spines["left"].set_bounds(0, 80)


def lollipop_figure(target):
    """Lecture slide 45."""
    fig, ax = sf.new_fig(target)
    draw_lollipop(ax, 10.5, size_k=2.6)
    ax.text(4, mine.max() * 1.12,
            f"Inside the DNA-binding domain ({P_DBD:.0%} of the length):\n"
            f"missense {TESTS['missense']['f']:.0%} (binomial p = {sci(TESTS['missense']['p'])})\n"
            f"truncating {TESTS['truncating']['f']:.0%} (p = {TESTS['truncating']['p']:.2f})",
            ha="left", va="top", fontsize=10.5, color=dviz.INK, linespacing=1.35)
    ax.set_title(f"TP53 lollipop · {len(var):,} tumours · counts simulated; domains (IARC TP53 Database) and hotspot "
                 f"residues real", fontsize=13)
    ax.set_xlabel("Amino-acid position (TP53, 393 aa) · stem height and dot size = tumours")
    ax.set_ylabel("Tumours")
    fig.subplots_adjust(left=0.06, right=0.99, top=0.9, bottom=0.17)
    return sf.save_exact(fig, target)


def practice_b5(target):
    """Practice slide B5: test (log2 + segments) → check (truth inside the LCR22 repeats) → overview (1 Mb)."""
    fig = plt.figure(figsize=(sf.PRACTICE_PX[0] / sf.DPI, sf.PRACTICE_PX[1] / sf.DPI), dpi=sf.DPI)
    dviz.use_style(base_size=12.5)
    outer = fig.add_gridspec(2, 1, height_ratios=[2.75, 1.0], hspace=0.5)
    inner = outer[0].subgridspec(3, 1, height_ratios=[1.75, 0.55, 0.45], hspace=0.16)
    ax1 = fig.add_subplot(inner[0])
    ax2 = fig.add_subplot(inner[1], sharex=ax1)
    ax3 = fig.add_subplot(inner[2], sharex=ax1)
    ax4 = fig.add_subplot(outer[1])
    w = cn.loc[(cn["start"] >= 17_000_000) & (cn["end"] <= 23_000_000)]
    draw_log2(ax1, w, lw_seg=2.4, s=3, label=False)
    ax1.text(22.95, -1.08, f"{MINUS}1 = 1 copy", ha="right", va="top", fontsize=9.5, color=dviz.SLATE)
    ax1.text((D0 + D1) / 2e6, 0.72, f"called: {D0 / 1e6:.2f}–{D1 / 1e6:.2f} Mb · log2 {num(MEAN_DEL, '{:.2f}')} "
             f"= {dele['copy_number']:.2f} copies", ha="center", va="center", fontsize=9.5, color=dviz.BLUE,
             fontweight="bold")
    ax1.set_ylabel("log2 ratio")
    ax1.set_title("B5 · 22q11.2 deletion: test, check, visualize", fontsize=13)
    draw_genes(ax2, ["HIRA", "TBX1", "COMT", "CRKL", "MAPK1"], 9, y_lab=(0.2, 0.8))
    draw_lcr(ax3, 8.5, truth_marks=True)
    ax3.text(22.98, 0.95, "▼ planted breakpoints", fontsize=8.5, ha="right", va="bottom", color=dviz.INK)
    shade([ax1, ax2, ax3])
    for ax in (ax1, ax2):
        ax.tick_params(labelbottom=False)
    ax2.spines["bottom"].set_visible(False)                         # gene labels below the boxes need the room
    ax2.tick_params(bottom=False)
    ax3.set_xlim(17, 23)
    ax3.set_xlabel("chr22 position (Mb, hg38) · 10-kb bins", fontsize=10.5, labelpad=2)
    st = state(gw["log2_patient_seg"])
    for name, c in STATE_COL.items():
        m = st == name
        ax4.scatter(gw.loc[m, "x"], gw.loc[m, "log2_patient"], s=3 if name == "none" else 7, color=c, lw=0)
    ax4.axhline(0, color=dviz.SLATE, lw=0.6)
    low = gw.loc[(gw["chrom"] == "chr22") & (gw["log2_patient"] < -0.5)]
    ax4.annotate(f"the same deletion at 1 Mb: {len(low)} bins", xy=(low["x"].mean() - 10, -0.95),
                 xytext=(low["x"].mean() - 90, -0.95), ha="right", va="center", fontsize=9.5,
                 arrowprops={"arrowstyle": "->", "color": dviz.INK, "lw": 0.8})
    ax4.set_ylim(-1.3, 0.6)
    ax4.set_yticks([-1, 0])
    ax4.set_ylabel("log2 ratio")
    ax4.set_xticks(ticks["centre"] / 1e6, [c[3:] if (int(c[3:]) % 2 == 1 and c != "chr21") or c == "chr22" else ""
                                           for c in ticks["chrom"]], fontsize=8.5)
    ax4.set_xlabel("Patient genome-wide, 1-Mb bins (chromosome)", fontsize=10.5, labelpad=2)
    fig.subplots_adjust(left=0.115, right=0.98, top=0.93, bottom=0.125)
    return sf.save_exact(fig, target)


def practice_b8(target):
    """Practice slide B8: lollipop (visualize) + share inside the DBD vs expectation (test, check)."""
    fig = plt.figure(figsize=(sf.PRACTICE_PX[0] / sf.DPI, sf.PRACTICE_PX[1] / sf.DPI), dpi=sf.DPI)
    dviz.use_style(base_size=12.5)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.75, 1.0], hspace=0.55)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    draw_lollipop(ax1, 9, size_k=1.4, label_fs=8.5, legend=False, nudge={273: (9, 5, "left")})
    for i, (cls, c) in enumerate(CLASS_COL.items()):                   # class key in the empty top-right corner
        y = 50 - 10 * i
        ax1.scatter(332, y, s=40, color=c)
        ax1.text(340, y, cls, va="center", fontsize=9, color=dviz.INK)
    ax1.set_ylabel("Tumours")
    ax1.set_xlabel("Amino-acid position (TP53, 393 aa)", fontsize=10.5, labelpad=2)
    ax1.set_title(f"B8 · TP53 in {len(var):,} tumours (simulated): where on the protein?", fontsize=13)
    rows = [("missense", TESTS["missense"]), ("truncating", TESTS["truncating"])]
    for i, (cls, t) in enumerate(rows):
        ax2.barh(i, t["f"] * 100, color=CLASS_COL[cls], height=0.55)
        ax2.text(t["f"] * 100 - 1.5, i, f"{t['k']}/{t['n']} = {t['f']:.0%}", ha="right", va="center",
                 fontsize=9, color="white", fontweight="bold")
        p_txt = f"p = {sci(t['p'])}" if t["p"] < 0.001 else f"p = {t['p']:.2f}"
        ax2.text(t["f"] * 100 + 1.5, i, p_txt, ha="left", va="center", fontsize=9, color=dviz.INK)
    ax2.axvline(P_DBD * 100, color=dviz.INK, lw=1, ls=(0, (4, 3)))
    ax2.text(P_DBD * 100 + 1.2, -0.52, f"expected from length: 199/393 = {P_DBD:.0%}", ha="left", va="center",
             fontsize=8.5, color=dviz.INK)
    ax2.set_yticks([0, 1], ["Missense", "Truncating"])
    ax2.set_ylim(1.4, -0.7)
    ax2.set_xlim(0, 125)
    ax2.set_xticks([0, 25, 50, 75, 100])
    ax2.spines["bottom"].set_bounds(0, 100)
    ax2.set_xlabel(f"% inside the DNA-binding domain · binomial test vs 51% · Fisher OR {ODDS:.1f}",
                   fontsize=10, labelpad=2)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.93, bottom=0.135)
    return sf.save_exact(fig, target)


def contact_sheet(paths, out):
    imgs = [plt.imread(str(p)) for p in paths]
    fig, axs = plt.subplots(len(imgs), 1, figsize=(9, 2.9 * len(imgs)), dpi=110)
    for ax, im, p in zip(axs, imgs, paths):
        ax.imshow(im)
        ax.set_title(f"{p.name}  ({im.shape[1]} × {im.shape[0]} px)", fontsize=10, loc="left")
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(out, dpi=110, metadata={"Software": None})
    plt.close(fig)
    return out


if __name__ == "__main__":
    made = [tracks_figure("image-43-1.png"), tracks_figure("image-58-1.png"), tracks_thumbnail("image-72-1.png"),
            genome_figure("image-44-1.png"), lollipop_figure("image-45-1.png"),
            practice_b5("practice_B5.png"), practice_b8("practice_B8.png")]
    for p in made:
        print(p.relative_to(ROOT), sf.size_of(p.name))
    print(contact_sheet(made, HERE / "contact_sheet_tracks.png").relative_to(ROOT))
