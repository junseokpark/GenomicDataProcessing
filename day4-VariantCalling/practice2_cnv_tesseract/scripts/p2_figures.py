"""Slide figures for Practice 2 (deck style: SUIT font, course palette, no titles, 300 dpi, white).

Every function draws one figure at the exact physical size used on the slide, saves it with
fig.savefig(out, dpi=300) and returns the figure. Data are computed in the notebook and passed in.
The notebook displays the saved PNGs; open this file if you want to see how a plot is built.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Rectangle

BLUE, DARK, INK, GRAY, LGRAY, LIGHT, PALE, ORANGE = ("#1F5FEA", "#173FAF", "#17202A", "#5B6573",
                                                     "#8A93A0", "#B7CEF5", "#EEF3FE", "#B25E09")
SEGDUP_FILL = (0.541, 0.576, 0.627, 0.22)    # LGRAY at 22% opacity
BASE, SMALL, TICK = 11, 10, 9                 # three-size ladder: labels / annotations / ticks


def setup(font_dir="fonts"):
    """Register the SUIT font (if present) and set the deck rcParams."""
    family = "DejaVu Sans"
    for f in ("SUIT-Regular.ttf", "SUIT-Bold.ttf"):
        p = os.path.join(font_dir, f)
        if os.path.exists(p):
            font_manager.fontManager.addfont(p)
            family = "SUIT"
    mpl.rcParams.update({
        # a font LIST gives per-glyph fallback: SUIT has no minus sign, en dash, >= or infinity glyphs
        "font.family": [family, "DejaVu Sans"] if family != "DejaVu Sans" else ["DejaVu Sans"],
        "font.size": BASE, "axes.labelsize": BASE, "axes.titlesize": BASE, "legend.fontsize": SMALL,
        "xtick.labelsize": TICK, "ytick.labelsize": TICK, "axes.linewidth": 0.6,
        "axes.edgecolor": GRAY, "axes.labelcolor": INK, "xtick.color": GRAY, "ytick.color": GRAY,
        "xtick.labelcolor": INK, "ytick.labelcolor": INK, "text.color": INK,
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "figure.dpi": 100, "savefig.dpi": 300,
        "savefig.facecolor": "white", "figure.facecolor": "white", "axes.facecolor": "white",
        "lines.linewidth": 1.2, "patch.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.unicode_minus": True})
    return family


def _save(fig, out):
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=300, facecolor="white")
    return out


def layout_problems(fig):
    """Render-then-verify: overlapping visible texts and texts outside the figure (should be [])."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    fb = fig.bbox
    skip = set()
    for ax in fig.axes:
        for labs, lim, k in ((ax.get_xticklabels(), sorted(ax.get_xlim()), 0), (ax.get_yticklabels(), sorted(ax.get_ylim()), 1)):
            for lab in labs:
                v = lab.get_position()[k]
                if not (lim[0] - 1e-9 <= v <= lim[1] + 1e-9):
                    skip.add(id(lab))
    texts = [t for t in fig.findobj(mpl.text.Text) if t.get_visible() and t.get_text().strip() and id(t) not in skip]
    boxes = [(t, t.get_window_extent(r)) for t in texts]
    probs = []
    for i, (a, ba) in enumerate(boxes):
        if ba.x0 < fb.x0 - 1 or ba.y0 < fb.y0 - 1 or ba.x1 > fb.x1 + 1 or ba.y1 > fb.y1 + 1:
            probs.append(("outside", a.get_text()))
        for b, bb in boxes[i + 1:]:
            if ba.overlaps(bb) and ba.width > 0 and bb.width > 0:
                probs.append(("overlap", a.get_text(), b.get_text()))
    return probs


def _mb(ax, axis="x", decimals=1):
    fmt = mpl.ticker.FuncFormatter(lambda v, _: f"{v / 1e6:.{decimals}f}")
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_formatter(fmt)


# ------------------------------------------------------------------ 1. CNV profiles (Step 3)
def fig_cnv_profiles(bins, log2c, masked, segs, out, samples=("P1", "P2", "P3"),
                     windows=(("chr22", "22q11.2"), ("chr17", "17p12"), ("chr4", "4q13.2"))):
    fig = plt.figure(figsize=(12.0, 3.6))
    widths = [3.3, 1.35, 1.15]
    gs = fig.add_gridspec(len(samples), len(windows), width_ratios=widths, left=0.088, right=0.885,
                          top=0.90, bottom=0.165, wspace=0.07, hspace=0.16)
    ylim = (-3.35, 1.45)
    for i, s in enumerate(samples):
        for j, (chrom, band) in enumerate(windows):
            ax = fig.add_subplot(gs[i, j])
            sel = ((bins.chrom == chrom) & (bins.n_frac <= 0.1)).values     # assembly gaps have no data: not drawn
            x = ((bins.start + bins.end) / 2).values[sel]
            y = log2c[s].values[sel]
            m = masked.values[sel]
            col = np.array([GRAY] * len(x), dtype=object)
            sg = segs[(segs["sample"] == s) & (segs.chrom == chrom)]
            for _, g in sg.iterrows():
                inside = (x >= g.start) & (x <= g.end)
                col[inside] = BLUE if g.cn < 2 else ORANGE
            ax.axhline(0, color=LGRAY, lw=0.6, zorder=1)
            for yv in (-1, np.log2(1.5)):
                ax.axhline(yv, color=LGRAY, lw=0.6, ls=(0, (3, 3)), zorder=1)
            ax.scatter(x[m], np.clip(np.nan_to_num(y[m], nan=0), *ylim), s=2.5, color=LGRAY, alpha=0.45,
                       lw=0, zorder=2, rasterized=True)
            ax.scatter(x[~m], y[~m], s=3.0, c=list(col[~m]), lw=0, zorder=3, rasterized=True)
            for _, g in sg.iterrows():
                c = BLUE if g.cn < 2 else ORANGE
                ax.plot([g.start, g.end], [g.log2, g.log2], color=DARK if g.cn < 2 else ORANGE, lw=2.4,
                        solid_capstyle="butt", zorder=4)
                xm = (g.start + g.end) / 2
                lab = f"CN {g.cn}"
                yt = g.log2 + 0.62 if g.cn >= 2 else g.log2 - 0.95
                if g.cn == 0:
                    yt = g.log2 + 0.62
                ax.text(xm, yt, lab, ha="center", va="center", fontsize=SMALL, color=c, zorder=5)
            allw = (bins.chrom == chrom).values
            w0, w1 = bins.start[allw].min(), bins.end[allw].max()
            ax.set_xlim(w0, w1)
            ax.set_ylim(*ylim)
            ax.set_yticks([-3, -2, -1, 0, 1])
            if j > 0:
                ax.tick_params(labelleft=False)
            if i < len(samples) - 1:
                ax.tick_params(labelbottom=False)
            _mb(ax, decimals=0 if chrom == "chr22" else 1)
            if chrom == "chr22":
                ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([19e6, 20e6, 21e6, 22e6]))
            elif chrom == "chr17":
                ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([14.8e6, 15.2e6, 15.6e6]))
            else:
                ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([69.1e6, 69.4e6, 69.7e6]))
            if i == 0:
                ax.text(0.0, 1.06, f"{chrom} · {band}", transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=BASE, color=INK)
            if j == 0:
                ax.text(-0.068, 0.5, s, transform=ax.transAxes, ha="right", va="center", fontsize=BASE,
                        color=INK, fontweight="bold")
            if j == len(windows) - 1:
                ax.text(1.03, -1, "−1  (CN 1)", transform=ax.get_yaxis_transform(), ha="left", va="center",
                        fontsize=SMALL, color=BLUE)
                ax.text(1.03, np.log2(1.5) + 0.12, "+0.58  (CN 3)", transform=ax.get_yaxis_transform(),
                        ha="left", va="center", fontsize=SMALL, color=ORANGE)
                if i == 0:
                    ax.text(1.03, -3.0, "● masked bin", transform=ax.get_yaxis_transform(), ha="left",
                            va="center", fontsize=SMALL, color=LGRAY)
    fig.text(0.012, 0.53, "log2 ratio to panel of normals", rotation=90, ha="left", va="center", fontsize=BASE)
    fig.text(0.47, 0.02, "hg19 position (Mb)", ha="center", va="bottom", fontsize=BASE)
    _save(fig, out)
    return fig


# ------------------------------------------------------------------ 2. caller steps (Step 3)
def fig_caller_steps(x0, x1, xc, raw, pm, panel_pm, ratio, log2c, masked, seg, truth_iv, out,
                     ctrl_pm_min=None, ctrl_pm_max=None):
    fig, axs = plt.subplots(4, 1, figsize=(6.4, 3.6), sharex=True,
                            gridspec_kw=dict(left=0.23, right=0.80, top=0.975, bottom=0.135, hspace=0.24))
    m = masked
    labels = ["1  raw\nfragments", "2  per million\nfragments", "3  ratio to\npanel", "4  log2 ratio\n+ HMM state"]
    for ax in axs:
        ax.axvspan(truth_iv[0], truth_iv[1], color=PALE, zorder=0, lw=0)
    # a) raw counts
    ax = axs[0]
    ax.scatter(xc[~m], raw[~m], s=7, color=BLUE, lw=0, zorder=3)
    ax.scatter(xc[m], raw[m], s=7, color=LGRAY, alpha=0.6, lw=0, zorder=2)
    ax.set_ylim(0, max(raw) * 1.18)
    ax.text(1.02, 0.72, "P2", transform=ax.transAxes, color=BLUE, fontsize=SMALL, va="center")
    # b) per million + panel median (+ control range)
    ax = axs[1]
    if ctrl_pm_min is not None:
        ax.fill_between(xc, ctrl_pm_min, ctrl_pm_max, step="mid", color=LGRAY, alpha=0.25, lw=0, zorder=1)
    ax.step(xc, panel_pm, where="mid", color=GRAY, lw=1.0, zorder=2)
    ax.scatter(xc[~m], pm[~m], s=7, color=BLUE, lw=0, zorder=3)
    ax.set_ylim(0, max(np.nanmax(pm), np.nanmax(panel_pm)) * 1.2)
    ax.text(1.02, 0.86, "P2", transform=ax.transAxes, color=BLUE, fontsize=SMALL, va="center")
    ax.text(1.02, 0.38, "controls:\nmedian, range", transform=ax.transAxes, color=GRAY,
            fontsize=SMALL, va="center", linespacing=1.0)
    # c) ratio
    ax = axs[2]
    ax.axhline(1, color=LGRAY, lw=0.6)
    ax.axhline(1.5, color=ORANGE, lw=0.7, ls=(0, (3, 3)))
    ax.scatter(xc[~m], ratio[~m], s=7, color=BLUE, lw=0, zorder=3)
    ax.set_ylim(0.3, 2.1)
    ax.set_yticks([0.5, 1.0, 1.5, 2.0])
    ax.text(1.02, (1.5 - 0.3) / 1.8, "1.5 = 3 copies / 2", transform=ax.transAxes, color=ORANGE,
            fontsize=SMALL, va="center")
    # d) log2 + HMM segment
    ax = axs[3]
    ax.axhline(0, color=LGRAY, lw=0.6)
    ax.axhline(np.log2(1.5), color=ORANGE, lw=0.7, ls=(0, (3, 3)))
    inside = (xc >= seg[0]) & (xc <= seg[1])
    ax.scatter(xc[~m & ~inside], log2c[~m & ~inside], s=7, color=GRAY, lw=0, zorder=3)
    ax.scatter(xc[~m & inside], log2c[~m & inside], s=7, color=ORANGE, lw=0, zorder=3)
    ax.plot([seg[0], seg[1]], [seg[2], seg[2]], color=ORANGE, lw=2.4, solid_capstyle="butt", zorder=4)
    ax.plot([x0, seg[0]], [0, 0], color=GRAY, lw=1.6, zorder=2)
    ax.plot([seg[1], x1], [0, 0], color=GRAY, lw=1.6, zorder=2)
    ax.set_ylim(-1.2, 1.3)
    ax.set_yticks([-1, 0, 1])
    ax.text(1.02, (np.log2(1.5) + 1.2) / 2.5, "HMM: CN 3", transform=ax.transAxes, color=ORANGE,
            fontsize=SMALL, va="center")
    ax.text(1.02, 1.2 / 2.5 - 0.17, "HMM: CN 2", transform=ax.transAxes, color=GRAY, fontsize=SMALL,
            va="center")
    for ax, lab in zip(axs, labels):
        ax.set_ylabel(lab, fontsize=SMALL, rotation=0, ha="left", va="center", labelpad=8, linespacing=1.0)
        ax.yaxis.set_label_coords(-0.39, 0.5)
        ax.tick_params(axis="y", labelsize=TICK)
    axs[-1].set_xlim(x0, x1)
    _mb(axs[-1], decimals=2)
    axs[-1].xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.1e6))
    axs[-1].set_xlabel("hg19 position on chr17 (Mb)   ·   shaded: planted PMP22 duplication", fontsize=SMALL)
    _save(fig, out)
    return fig


# ------------------------------------------------------------------ 3. liftover (Step 4)
def fig_liftover(hg19_iv, pieces, bp_iv, outer_iv, clinvar38_iv, sd19, sd38, out, xlim19, xlim38,
                 off_chrom_note="", n_pieces_total=None):
    fig = plt.figure(figsize=(6.4, 3.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 2.6], left=0.255, right=0.975, top=0.955,
                          bottom=0.15, hspace=0.62)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    for ax, sd, xl in ((ax1, sd19, xlim19), (ax2, sd38, xlim38)):
        for a, b in sd:
            if b > xl[0] and a < xl[1]:
                ax.axvspan(max(a, xl[0]), min(b, xl[1]), color=SEGDUP_FILL, lw=0, zorder=0)
        ax.set_xlim(*xl)
        ax.spines["left"].set_visible(False)
        ax.tick_params(axis="y", length=0)
        _mb(ax, decimals=1)
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5e6))

    def bar(ax, y, iv, color, size_txt, hollow=False):
        if hollow:
            ax.add_patch(Rectangle((iv[0], y - 0.18), iv[1] - iv[0], 0.36, facecolor="white", edgecolor=color,
                                   lw=1.4, zorder=3))
        else:
            ax.add_patch(Rectangle((iv[0], y - 0.18), iv[1] - iv[0], 0.36, facecolor=color, edgecolor="none",
                                   zorder=3))
        ax.text(iv[1] + (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.012, y, size_txt, ha="left", va="center",
                fontsize=SMALL, color=color, zorder=4)

    ax1.set_ylim(-0.6, 0.6)
    bar(ax1, 0, hg19_iv, BLUE, f"{(hg19_iv[1] - hg19_iv[0]) / 1e6:.2f} Mb")
    ax1.set_yticks([0])
    ax1.set_yticklabels(["hg19 deletion\n(ClinVar GRCh37)"], fontsize=SMALL, linespacing=1.0)
    ax1.set_xlabel("hg19 chr22 (Mb)", fontsize=SMALL, labelpad=2)
    ax2.set_ylim(-0.6, 3.75)
    ys = {"pieces": 3, "bp": 2, "outer": 1, "clinvar": 0}
    for a, b in pieces:                                   # each lifted piece: a thin bar ...
        ax2.add_patch(Rectangle((a, ys["pieces"] - 0.13), max(b - a, 1500), 0.26, facecolor=DARK,
                                edgecolor="none", zorder=3))
    starts = sorted(a for a, b in pieces)                 # ... and a tick where each piece starts
    ax2.vlines(starts, ys["pieces"] + 0.17, ys["pieces"] + 0.42, color=DARK, lw=0.5, zorder=3)
    if off_chrom_note:
        ax2.text(xlim38[1], ys["pieces"] + 0.55, off_chrom_note, ha="right", va="bottom", fontsize=SMALL, color=GRAY)
    bar(ax2, ys["bp"], bp_iv, BLUE, f"{(bp_iv[1] - bp_iv[0]) / 1e6:.2f} Mb")
    bar(ax2, ys["outer"], outer_iv, GRAY, f"{(outer_iv[1] - outer_iv[0]) / 1e6:.2f} Mb")
    bar(ax2, ys["clinvar"], clinvar38_iv, DARK, f"{(clinvar38_iv[1] - clinvar38_iv[0]) / 1e6:.2f} Mb", hollow=True)
    ax2.set_yticks([ys["pieces"], ys["bp"], ys["outer"], ys["clinvar"]])
    n = n_pieces_total if n_pieces_total is not None else len(pieces)
    ax2.set_yticklabels([f"CrossMap pieces ({n})", "breakpoint lift", "outer span of pieces",
                         "ClinVar GRCh38 record"], fontsize=SMALL)
    ax2.set_xlabel("hg38 chr22 (Mb)   ·   gray shading: segmental duplications", fontsize=SMALL, labelpad=2)
    _save(fig, out)
    return fig


# ------------------------------------------------------------------ 4. RNA dosage (Step 5)
def fig_rna_dosage(panels, out):
    fig = plt.figure(figsize=(12.0, 3.4))
    widths = [p.get("width", 1.0) for p in panels]
    gs = fig.add_gridspec(1, len(panels), width_ratios=widths, left=0.068, right=0.995, top=0.80,
                          bottom=0.175, wspace=0.10)
    ylim = (-3.4, 1.5)
    for k, p in enumerate(panels):
        ax = fig.add_subplot(gs[k])
        g = p["genes"]
        loss = p["cn"] < 2
        col = BLUE if loss else ORANGE
        ax.axvspan(p["cnv_start"], p["cnv_end"], color=PALE, lw=0, zorder=0)
        ax.axhline(0, color=LGRAY, lw=0.6, zorder=1)
        exp_y = -1.0 if p["cn"] == 1 else (np.log2(1.5) if p["cn"] == 3 else None)
        if exp_y is not None:
            ax.plot([p["cnv_start"], p["cnv_end"]], [exp_y, exp_y], color=col, lw=0.9, ls=(0, (3, 2)), zorder=2)
        xm = (g.start0 + g.end) / 2
        y = g.log2fc.clip(ylim[0] + 0.25, ylim[1] - 0.1)
        ins = g.inside.values
        ex = g.expressed.values
        clipped = (g.log2fc < ylim[0] + 0.25).values
        ax.scatter(xm[~ins & ex], y[~ins & ex], s=12, color=LGRAY, lw=0, zorder=3)
        ax.scatter(xm[~ins & ~ex], y[~ins & ~ex], s=12, facecolor="white", edgecolor=LGRAY, lw=0.7, zorder=3)
        ax.scatter(xm[ins & ex & ~clipped], y[ins & ex & ~clipped], s=18, color=col, lw=0, zorder=4)
        ax.scatter(xm[ins & ex & clipped], y[ins & ex & clipped], s=30, color=col, marker="v", lw=0, zorder=4)
        ax.scatter(xm[ins & ~ex], y[ins & ~ex], s=18, facecolor="white", edgecolor=col, lw=0.8, zorder=4)
        for gene in p.get("label_genes", []):
            r = g[g.gene == gene]
            if len(r):
                r = r.iloc[0]
                yy = float(np.clip(r.log2fc, ylim[0] + 0.25, ylim[1] - 0.1))
                dy = p.get("label_dy", {}).get(gene, 0.45 if yy < 0 else -0.45)
                dx = p.get("label_dx", {}).get(gene, 0.0) * (p["xlim"][1] - p["xlim"][0])
                ax.annotate(gene, ((r.start0 + r.end) / 2, yy), xytext=((r.start0 + r.end) / 2 + dx, yy + dy),
                            fontsize=SMALL, color=INK, ha="center", va="center", style="italic",
                            arrowprops=dict(arrowstyle="-", color=LGRAY, lw=0.6, shrinkA=0, shrinkB=2.5))
        ax.set_ylim(*ylim)
        ax.set_xlim(p["xlim"])
        if k > 0:
            ax.tick_params(labelleft=False)
        else:
            ax.set_ylabel("RNA log2 fold change\n(patient vs control mean)", fontsize=BASE, linespacing=1.0)
        _mb(ax, decimals=p.get("decimals", 1))
        ax.xaxis.set_major_locator(mpl.ticker.MaxNLocator(p.get("nticks", 3)))
        ax.text(0.0, 1.17, p["title"], transform=ax.transAxes, ha="left", va="bottom", fontsize=BASE,
                color=INK, fontweight="bold")
        ax.text(0.0, 1.03, p["subtitle"], transform=ax.transAxes, ha="left", va="bottom", fontsize=SMALL,
                color=col)
        if p.get("note"):
            ax.text(0.5, 0.96, p["note"], transform=ax.transAxes, ha="center", va="top", fontsize=SMALL,
                    color=GRAY)
    fig.text(0.5, 0.015, "hg38 position (Mb)   ·   shaded: CNV   ·   filled: expressed gene   ·   "
             "hollow: not expressed (CPM < 1)   ·   gray: flanking genes   ·   dashed: expected for copy number",
             ha="center", va="bottom", fontsize=SMALL, color=GRAY)
    _save(fig, out)
    return fig


# ------------------------------------------------------------------ 5. track panels (Step 7)
HI_COL = {3: DARK, 2: BLUE, 1: LIGHT}          # every other score (0, 30, 40, not yet evaluated): LGRAY


def _hicol(v):
    return HI_COL.get(int(v), LGRAY) if str(v).isdigit() else LGRAY


def fig_tracks(region, cnv_iv, cnv_type, segments, tr, out, figsize, label_genes=(), bin_bp=20000, xtick=0.25e6,
               decimals=1, seg_label=None, n_similar=None):
    """Genome-browser style panel. tr: dict with DataFrames (hg38, already restricted to the region)."""
    rows = ["CNV · 60-kb segments", "genes (MANE)", "ClinGen HI score", "ClinGen TS score", "ClinGen gene–disease",
            "gnomAD pLI", "ClinVar P/LP CNVs", "ClinVar P/LP SNVs", "conserved elements", "ENCODE cCREs",
            "segmental dups", "common SVs (gnomAD)"]
    n = len(rows)
    W, H = figsize
    fig = plt.figure(figsize=figsize)
    left = 1.62 / W
    bottom = (0.78 if W > 8 else 0.98) / H
    ax = fig.add_axes([left, bottom, 0.985 - left - (0.0 if W > 8 else 0.0), 0.985 - bottom - 0.01])
    x0, x1 = region
    ax.set_xlim(x0, x1)
    ax.set_ylim(-0.6, n - 0.4)
    ax.invert_yaxis()
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_yticks(range(n))
    ax.set_yticklabels(rows, fontsize=SMALL)
    for i in range(1, n):
        ax.axhline(i - 0.5, color=PALE, lw=0.8, zorder=0)
    col = BLUE if cnv_type == "DEL" else ORANGE
    grid_col = LIGHT if cnv_type == "DEL" else "#E8C9A6"
    span = x1 - x0
    ax.axvspan(cnv_iv[0], cnv_iv[1], color=PALE, lw=0, zorder=0)
    # row 0: CNV bar with the 60-kb segment grid
    ax.add_patch(Rectangle((cnv_iv[0], -0.22), cnv_iv[1] - cnv_iv[0], 0.44, facecolor=col, edgecolor="none", zorder=3))
    for a, b in segments:
        ax.plot([a, a], [-0.5, n - 0.5], color=grid_col, lw=0.5, zorder=1)
        ax.plot([a, a], [-0.22, 0.22], color="white", lw=0.8, zorder=4)
    if segments:
        ax.plot([segments[-1][1]] * 2, [-0.5, n - 0.5], color=grid_col, lw=0.5, zorder=1)
    if seg_label:
        right_room = (x1 - cnv_iv[1]) / span
        if right_room > 0.16:
            ax.text(cnv_iv[1] + span * 0.006, 0, seg_label, ha="left", va="center", fontsize=SMALL, color=col)
        else:
            ax.text(cnv_iv[0] - span * 0.006, 0, seg_label, ha="right", va="center", fontsize=SMALL, color=col)
    # row 1: genes (gene bodies); a few labels
    g = tr["genes"]
    minw = span * 0.0015
    for _, r in g.iterrows():
        ax.add_patch(Rectangle((max(r.start0, x0), 1 - 0.16), max(min(r.end, x1) - max(r.start0, x0), minw), 0.32,
                               facecolor=INK, edgecolor="none", zorder=3))
    for gene in label_genes:
        r = g[g.gene == gene]
        if len(r):
            r = r.iloc[0]
            xm = (max(r.start0, x0) + min(r.end, x1)) / 2
            ax.text(xm, 1 - 0.23, gene, ha="center", va="bottom", fontsize=SMALL - 1, style="italic", color=INK, zorder=5)
    # rows 2-3: ClinGen HI / TS (regions: thin bar at the bottom of the row; genes: boxes)
    for row, key, scol in ((2, "hi", "haploScore"), (3, "ts", "triploScore")):
        d = tr[key]
        for _, r in d.sort_values("chromEnd").iterrows():
            is_region = str(r["name"]).startswith("ISCA")
            h = 0.16 if is_region else 0.42
            yy = row + (0.24 if is_region else -0.06)
            a, b = max(r.chromStart, x0), min(r.chromEnd, x1)
            ax.add_patch(Rectangle((a, yy - h / 2), max(b - a, minw), h, facecolor=_hicol(r[scol]), edgecolor="none",
                                   zorder=3))
    # row 4: gene-disease validity
    for _, r in tr["gd"].iterrows():
        c = DARK if r.Classification in ("Definitive", "Strong") else LIGHT
        ax.add_patch(Rectangle((max(r.chromStart, x0), 4 - 0.2), max(min(r.chromEnd, x1) - max(r.chromStart, x0), minw),
                               0.4, facecolor=c, edgecolor="none", zorder=3))
    # row 5: pLI (bar height = pLI; dark if >= 0.9)
    for _, r in tr["pli"].iterrows():
        h = 0.8 * float(r["_pli"])
        ax.add_patch(Rectangle((max(r.chromStart, x0), 5 + 0.4 - h), max(min(r.chromEnd, x1) - max(r.chromStart, x0), minw),
                               h, facecolor=DARK if r["_pli"] >= 0.9 else LIGHT, edgecolor="none", zorder=3))
    # binned rows
    edges = np.arange(x0, x1 + bin_bp, bin_bp)
    mids = (edges[:-1] + edges[1:]) / 2

    def binned(starts, ends, mode):
        v = np.zeros(len(mids))
        if mode == "cover":
            for a, b in zip(starts, ends):
                i0 = max(0, int((a - x0) // bin_bp)); i1 = min(len(mids) - 1, int((b - x0) // bin_bp))
                for i in range(i0, i1 + 1):
                    v[i] += max(0, min(b, edges[i + 1]) - max(a, edges[i])) / bin_bp
        else:
            idx = ((np.asarray(starts) - x0) // bin_bp).astype(int)
            np.add.at(v, idx[(idx >= 0) & (idx < len(mids))], 1)
        return v

    def bars(row, v, color, sqrt=False):
        if v.max() > 0:
            vv = np.sqrt(v) if sqrt else v
            hh = 0.8 * vv / vv.max()
            ax.bar(mids, hh, width=bin_bp * 0.92, bottom=row + 0.4 - hh, color=color, edgecolor="none", zorder=3)

    cv = tr["clinvar_cnv"]                     # similar records only (same type, >= 50% reciprocal overlap)
    bars(6, binned(cv.chromStart.values, cv.chromEnd.values, "cover"), col)
    if n_similar is not None:
        ax.text(x1 - span * 0.004, 6 - 0.30, f"n = {n_similar} similar records", ha="right", va="center",
                fontsize=SMALL - 1, color=GRAY, zorder=5)
    sn = tr["clinvar_snv"]
    bars(7, binned(sn.chromStart.values, sn.chromEnd.values, "count"), DARK, sqrt=True)
    pc = tr["phastcons"]
    bars(8, binned(pc.chromStart.values, pc.chromEnd.values, "cover"), GRAY)
    cc = tr["ccre"]
    bars(9, binned(cc.chromStart.values, cc.chromEnd.values, "count"), GRAY)
    sd = tr["segdup"]
    _intervals(ax, 10, list(zip(np.maximum(sd.chromStart.values, x0), np.minimum(sd.chromEnd.values, x1))), LGRAY, h=0.5)
    for row, d in ((10, sd), (11, tr["sv"])):
        if len(d) == 0:
            ax.text(x0 + span * 0.01, row, "none in this window", ha="left", va="center", fontsize=SMALL - 1,
                    color=LGRAY, zorder=5)
    for _, r in tr["sv"].iterrows():
        c = BLUE if r.svtype == "DEL" else ORANGE
        a, b = max(r.chromStart, x0), min(r.chromEnd, x1)
        ax.add_patch(Rectangle((a, 11 - 0.2), max(b - a, span * 0.002), 0.4, facecolor=c, edgecolor="none", zorder=3))
    _mb(ax, decimals=decimals)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(xtick))
    ax.set_xlabel("hg38 position (Mb)", fontsize=SMALL, labelpad=1)
    # legend (bottom; two lines on narrow figures)
    line1 = [("ClinGen score:", GRAY), ("■ 3", DARK), ("■ 2", BLUE), ("■ 1", LIGHT), ("■ 0 / 30 / 40", LGRAY)]
    line2 = [("gnomAD SVs (AF ≥ 1%, ≥ 5 kb):", GRAY), ("■ deletion", BLUE), ("■ duplication", ORANGE),
             ("   bars: per " + (f"{bin_bp // 1000}-kb" if bin_bp >= 1000 else f"{bin_bp}-bp") + " bin", GRAY)]
    lines = [line1 + [("   ", GRAY)] + line2] if W > 8 else [line1, line2]
    r_ = fig.canvas.get_renderer()
    for k, items in enumerate(lines):
        xpos, y_leg = left, (0.10 + 0.21 * (len(lines) - 1 - k)) / H
        for txt, c in items:
            t = fig.text(xpos, y_leg, txt, ha="left", va="bottom", fontsize=SMALL - 1, color=c)
            xpos += t.get_window_extent(r_).width / fig.bbox.width + 0.006
    _save(fig, out)
    return fig


def _intervals(ax, y, ivs, color, h=0.55, minw=0):
    for a, b in ivs:
        ax.add_patch(Rectangle((a, y - h / 2), max(b - a, minw), h, facecolor=color, edgecolor="none", zorder=3))


# ------------------------------------------------------------------ 6. TESSERACT-style image (Step 7)
def fig_tesseract_image(M, channel_labels, group_labels, x_edges, out):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("pale_blue", ["#FFFFFF", PALE, LIGHT, BLUE, DARK])
    fig = plt.figure(figsize=(6.4, 3.0))
    ax = fig.add_axes([0.30, 0.165, 0.585, 0.72])
    cax = fig.add_axes([0.905, 0.165, 0.018, 0.72])
    im = ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="nearest",
                   extent=(0, M.shape[1], M.shape[0], 0))
    ax.set_yticks(np.arange(M.shape[0]) + 0.5)
    ax.set_yticklabels(channel_labels, fontsize=TICK)
    for y in range(1, M.shape[0]):
        ax.axhline(y, color="white", lw=0.6)
    prev = None
    for i, gname in enumerate(group_labels):
        if gname != prev and i > 0:
            ax.axhline(i, color=GRAY, lw=0.8)
        prev = gname
    ax.tick_params(axis="y", length=0)
    nt = M.shape[1]
    ticks = [0, nt // 4, nt // 2, 3 * nt // 4, nt]
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel(f"500-bp window inside the duplication (n = {nt})", fontsize=SMALL, labelpad=2)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_visible(False)
    cb.set_ticks([0, 0.5, 1])
    cb.ax.tick_params(labelsize=TICK, length=2)
    cb.set_label("scaled 0–1 per channel", fontsize=SMALL, labelpad=3)
    fig.text(0.30, 0.975, f"Illustration from {M.shape[0]} public tracks — TESSERACT uses 78 channels",
             ha="left", va="top", fontsize=SMALL, color=ORANGE)
    _save(fig, out)
    return fig


# ------------------------------------------------------------------ 7. wrong build (Step 5)
def fig_wrong_build(rows, out):
    fig = plt.figure(figsize=(6.4, 3.0))
    ax = fig.add_axes([0.255, 0.17, 0.72, 0.70])
    n = len(rows)
    for i, r in enumerate(rows):
        y = n - 1 - i
        a, mi, sp = r["shared"], len(r["missed"]), len(r["spurious"])
        ax.barh(y, a, color=LGRAY, height=0.56, edgecolor="none")
        ax.barh(y, mi, left=a, color=BLUE, height=0.56, edgecolor="none")
        ax.barh(y, sp, left=a + mi, color=ORANGE, height=0.56, edgecolor="none")
        txt = []
        if mi:
            txt.append("missed: " + ", ".join(r["missed"][:3]) + (f" +{mi - 3}" if mi > 3 else ""))
        if sp:
            txt.append("added: " + ", ".join(r["spurious"][:3]) + (f" +{sp - 3}" if sp > 3 else ""))
        ax.text(a + mi + sp + 0.8, y, "   ".join(txt) if txt else "same genes", va="center", ha="left",
                fontsize=SMALL - 1, color=INK)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r["label"] for r in rows][::-1], fontsize=SMALL, linespacing=1.0)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, max(r["shared"] + len(r["missed"]) + len(r["spurious"]) for r in rows) * 1.9)
    ax.set_ylim(-0.55, n - 0.45)
    ax.set_xlabel("protein-coding genes (hg38 MANE)", fontsize=SMALL, labelpad=2)
    r_ = fig.canvas.get_renderer()
    xpos = 0.255
    for txt, c in [("■ found by both", GRAY), ("■ missed (hg19 numbers used on hg38)", BLUE),
                   ("■ wrongly added", ORANGE)]:
        t = fig.text(xpos, 0.955, txt, ha="left", va="top", fontsize=SMALL - 1, color=c)
        xpos += t.get_window_extent(r_).width / fig.bbox.width + 0.025
    _save(fig, out)
    return fig
