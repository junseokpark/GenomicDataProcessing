"""Slide figures for elective B3 (GWAS), drawn from data/ at the exact pixel size of their slide frames.

Targets (build/picture_frames.tsv; 288 dpi, no tight bbox):
  image-36-1.png   slide 36  Manhattan plot, like the original (3328 x 960 px, shown 11.6 in wide)
  image-37-1.png   slide 37  QQ plot of both p-value columns, identity line, lambda for both in the legend (1280 x 1080)
  image-57-1.png   slide 57  Manhattan for the demo slide (3328 x 960, shown 9.1 in wide -> text x 1.25)
  image-71-1.png   slide 71  the same Manhattan as a thumbnail (3328 x 960, shown 4.4 in wide -> text x 2.3)
  practice_B3.png  practice slide: PCA by case | QQ with and without PCs | small Manhattan (1800 x 1210)
plus build/contact_sheet_B3.png: all five at their on-slide sizes (not a slide target).

usage:  python build/slides_gwas.py      (about 15 s)
The helpers below repeat the notebook's code (electives/B3_gwas.ipynb) so the slides show the same numbers.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
import dviz  # noqa: E402
import slidefig  # noqa: E402

DATA = ROOT / "data"
warnings.filterwarnings("error", message=".*Glyph.*missing.*")   # a missing glyph must stop the build

# ------------------------------------------------------------------------------------------ data (as in the notebook)
gwas = pd.read_csv(DATA / "gwas_sumstats.tsv.gz", sep="\t")
samples = pd.read_csv(DATA / "gwas_samples.tsv", sep="\t")
snps = pd.read_csv(DATA / "gwas_snps_subset.tsv", sep="\t")
G = np.load(DATA / "gwas_genotypes_subset.npz")["genotypes"]
genes = pd.read_csv(DATA / "gwas_genes_hg38.tsv", sep="\t")
genes["CHR"] = genes["chrom"].str.replace("chr", "").astype(int)
y = samples["case"].to_numpy(dtype=float)


def nearest_gene(chrom, bp):
    g = genes[genes["CHR"] == chrom]
    first, last = g["start0"].to_numpy() + 1, g["end"].to_numpy()
    dist = np.where(bp < first, first - bp, np.where(bp > last, bp - last, 0))
    return g["gene"].iloc[np.lexsort((np.abs((first + last) / 2 - bp), dist))[0]]


def define_loci(hits, pcol="P", merge_bp=1_000_000):
    new = (hits["CHR"].diff() != 0) | (hits["BP"].diff() > merge_bp)
    hits = hits.assign(locus=np.cumsum(new.to_numpy()))
    lead = hits.loc[hits.groupby("locus")[pcol].idxmin()].reset_index(drop=True)
    lead["n_snps"] = hits.groupby("locus").size().to_numpy()
    lead["nearest_gene"] = [nearest_gene(c, b) for c, b in zip(lead["CHR"], lead["BP"])]
    return lead.copy()


def standardise(g):
    g = g.astype(np.float64)
    return (g - g.mean(axis=0)) / g.std(axis=0)


loci = define_loci(gwas[gwas["P"] < 5e-8].sort_values(["CHR", "BP"]))
lam_pc = dviz.genomic_lambda(gwas["P"])
lam_no = dviz.genomic_lambda(gwas["P_no_ancestry_PCs"])
G_pca = G[:, snps["pca_snp"].to_numpy() == 1]
pca = PCA(n_components=10, svd_solver="randomized", random_state=0).fit(standardise(G_pca))
pcs = pca.transform(standardise(G_pca))
var_pct = 100 * pca.explained_variance_ratio_
X_CUM, TICKS = dviz.manhattan_x(gwas["CHR"], gwas["BP"])
MLP = -np.log10(gwas["P"].to_numpy())


# ------------------------------------------------------------------------------------------ drawing helpers
def chrom_ticks(ax, min_gap_in):
    """Label chromosome centres at least min_gap_in inches apart; chromosomes with a locus are always labelled."""
    x0, x1 = ax.get_xlim()
    width_in = ax.get_position().width * ax.figure.get_figwidth()
    pos_in = {int(c): (xc - x0) / (x1 - x0) * width_in for c, xc in zip(TICKS["chrom"], TICKS["centre"])}
    keep = {int(c) for c in loci["CHR"]}
    for c in sorted(pos_in):
        if c not in keep and all(abs(pos_in[c] - pos_in[k]) >= min_gap_in for k in keep):
            keep.add(c)
    centre = dict(zip(TICKS["chrom"].astype(int), TICKS["centre"]))
    ax.set_xticks([centre[c] for c in sorted(keep)], [str(c) for c in sorted(keep)])


def label_loci(ax, fontsize, gap_pts=3):
    """Gene name above each locus; a label that would overlap another one is moved up until it is free."""
    fig = ax.figure
    fig.canvas.draw()
    ren = fig.canvas.get_renderer()
    placed = []
    xmax = X_CUM.max()
    order = loci.assign(mlp=-np.log10(loci["P"])).sort_values("mlp", ascending=False)
    for r in order.itertuples():
        i = int(np.flatnonzero(gwas["SNP"].to_numpy() == r.SNP)[0])
        ha = "right" if X_CUM[i] > 0.96 * xmax else "center"
        dy = gap_pts
        t = ax.annotate(r.nearest_gene, (X_CUM[i], MLP[i]), xytext=(0, dy), textcoords="offset points",
                        ha=ha, va="bottom", fontsize=fontsize, color=dviz.INK)
        for _ in range(60):
            bb = t.get_window_extent(ren).expanded(1.05, 1.1)
            if not any(bb.overlaps(p) for p in placed):
                break
            dy += 2
            t.set_position((0, dy))
        placed.append(t.get_window_extent(ren))


def manhattan(ax, size=2.0, label_fs=10, show_threshold_text=True, line_w=1.0, min_gap_in=0.32):
    is_sig = gwas["P"].to_numpy() < 5e-8
    colour = np.where(gwas["CHR"].to_numpy() % 2 == 1, dviz.DARKBLUE, dviz.SKY)
    ax.scatter(X_CUM[~is_sig], MLP[~is_sig], s=size, c=colour[~is_sig], lw=0, rasterized=True)
    ax.scatter(X_CUM[is_sig], MLP[is_sig], s=size * 1.8, color=dviz.VERMILION, lw=0, rasterized=True)
    ax.axhline(-np.log10(5e-8), color=dviz.VERMILION, ls="--", lw=line_w)
    ax.axhline(5, color=dviz.SLATE, ls=":", lw=line_w * 0.7)
    if show_threshold_text:
        ax.text(X_CUM[0], -np.log10(5e-8) + 1.0, "p = 5 × 10⁻⁸", color=dviz.VERMILION,
                fontsize=plt.rcParams["font.size"] * 0.85, va="bottom")
    ax.set_xlim(X_CUM.min() - 1.2e7, X_CUM.max() + 1.2e7)
    ax.set_ylim(0, MLP.max() * 1.16)
    ax.set_ylabel("−log10 p")
    chrom_ticks(ax, min_gap_in)
    ax.tick_params(axis="x", length=3)
    label_loci(ax, label_fs)


def qq(ax, cap, size=6, legend_fs=None, short=False, ms=7):
    """QQ plot of both p-value columns; values above `cap` are drawn as triangles on the top edge.

    The legend (upper left) carries lambda for both columns and explains the triangles; the curves only
    rise far from the diagonal on the right, so the upper-left corner stays free of data."""
    from matplotlib.lines import Line2D
    beyond = {}
    for col, colour in [("P_no_ancestry_PCs", dviz.ORANGE), ("P", dviz.BLUE)]:
        e, o = dviz.qq_points(gwas[col])
        m = o <= cap
        ax.scatter(e[m], o[m], s=size, color=colour, lw=0, rasterized=True)
        ax.scatter(e[~m], np.full((~m).sum(), cap), s=size * 2.2, marker="^", color=colour, lw=0, clip_on=False)
        beyond[col] = int((~m).sum())
    ax.plot([0, e.max()], [0, e.max()], color=dviz.INK, lw=1.0)
    ax.set_ylim(0, cap)
    ax.set_yticks(np.arange(0, cap + 0.1, 5))
    ax.set_xlim(0, e.max() * 1.03)
    ax.set_xlabel("Expected −log10 p")
    ax.set_ylabel("Observed −log10 p")
    if short:
        labels = [f"no PCs · λ = {lam_no:.2f}", f"4 PCs · λ = {lam_pc:.2f}", f"> {cap:g}, up to {MLP.max():.0f}"]
    else:
        labels = [f"no ancestry PCs\nλ = {lam_no:.2f}", f"with 4 ancestry PCs\nλ = {lam_pc:.2f}",
                  f"−log10 p > {cap:g} (max {MLP.max():.0f})"]
    handles = [Line2D([], [], ls="", marker="o", ms=ms, color=dviz.ORANGE),
               Line2D([], [], ls="", marker="o", ms=ms, color=dviz.BLUE),
               Line2D([], [], ls="", marker="^", ms=ms, color=dviz.SLATE)]
    return ax.legend(handles, labels, loc="upper left", fontsize=legend_fs, handletextpad=0.4, borderaxespad=0.2,
                     labelspacing=0.5)


# ------------------------------------------------------------------------------------------ 1. slide 36 Manhattan
fig, ax = slidefig.new_fig("image-36-1.png", base_size=12.5)
fig.subplots_adjust(left=0.052, right=0.995, bottom=0.2, top=0.87)
manhattan(ax, size=2.2, label_fs=10.5)
ax.set_xlabel("Chromosome (hg38)")
ax.set_title(f"{len(gwas):,} SNPs · {len(loci)} loci pass genome-wide significance (simulated)")
paths = [slidefig.save_exact(fig, "image-36-1.png")]

# ------------------------------------------------------------------------------------------ 2. slide 37 QQ
fig, ax = slidefig.new_fig("image-37-1.png", base_size=12.5)
fig.subplots_adjust(left=0.13, right=0.96, bottom=0.15, top=0.88)
qq(ax, cap=20, size=7, legend_fs=10.5)
ax.set_title("QQ plot and genomic inflation λ")
paths.append(slidefig.save_exact(fig, "image-37-1.png"))

# ------------------------------------------------------------------------------------------ 3. slide 57 demo Manhattan
fig, ax = slidefig.new_fig("image-57-1.png", base_size=15.5)
fig.subplots_adjust(left=0.06, right=0.995, bottom=0.235, top=0.85)
manhattan(ax, size=2.4, label_fs=13, min_gap_in=0.4)
ax.set_xlabel("Chromosome (hg38)")
ax.set_title(f"{len(gwas):,} SNPs · λ = {lam_pc:.2f} with 4 PCs · {len(loci)} loci pass 5 × 10⁻⁸ (simulated)")
paths.append(slidefig.save_exact(fig, "image-57-1.png"))

# ------------------------------------------------------------------------------------------ 4. slide 71 thumbnail
fig, ax = slidefig.new_fig("image-71-1.png", base_size=27)
fig.subplots_adjust(left=0.105, right=0.99, bottom=0.33, top=0.84)
manhattan(ax, size=4.0, label_fs=22, show_threshold_text=False, line_w=2.2, min_gap_in=0.9)
ax.set_ylim(0, MLP.max() * 1.42)                                  # room for the gene labels under the title
ax.set_yticks([0, 20, 40])
ax.set_xlabel("Chromosome (hg38)", fontsize=22, labelpad=2)
ax.set_title(f"Manhattan · {len(loci)} loci · λ = {lam_pc:.2f}", pad=6)
paths.append(slidefig.save_exact(fig, "image-71-1.png"))

# ------------------------------------------------------------------------------------------ 5. practice slide
fig, axs = slidefig.new_fig("practice_B3.png", 2, 2, base_size=12.5,
                            gridspec_kw={"height_ratios": [1.3, 1], "hspace": 0.66, "wspace": 0.36})
fig.subplots_adjust(left=0.11, right=0.985, bottom=0.135, top=0.93)
gs = axs[1, 0].get_gridspec()
for a in axs[1, :]:
    a.remove()
ax_pca, ax_qq, ax_man = axs[0, 0], axs[0, 1], fig.add_subplot(gs[1, :])
left = pcs[:, 0] < 0
for case, colour, lab in [(0, dviz.SLATE, "control"), (1, dviz.BLUE, "case")]:
    m = y == case
    ax_pca.scatter(pcs[m, 0], pcs[m, 1], s=3, color=colour, alpha=0.45, lw=0, label=lab, rasterized=True)
for m in (left, ~left):
    ax_pca.text(np.median(pcs[m, 0]), pcs[:, 1].max() * 1.08, f"{100 * y[m].mean():.0f}% cases",
                ha="center", va="bottom", fontsize=10.5, color=dviz.INK)
ax_pca.set_ylim(pcs[:, 1].min() * 1.75, pcs[:, 1].max() * 1.6)        # free bands for the labels and the legend
ax_pca.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(5))
ax_pca.set_xlabel(f"PC1 ({var_pct[0]:.2f}% of variance)")
ax_pca.set_ylabel(f"PC2 ({var_pct[1]:.2f}%)")
ax_pca.set_title("Genotype PCA by case status")
leg = ax_pca.legend(loc="lower center", ncol=2, fontsize=10, markerscale=3, handletextpad=0.1, columnspacing=0.8,
                    borderaxespad=0.1)
for h in leg.legend_handles:
    h.set_alpha(1)
qq(ax_qq, cap=15, size=4, legend_fs=9.5, short=True, ms=5)
ax_qq.set_title("QQ: with vs without PCs")
manhattan(ax_man, size=1.2, label_fs=9.5, show_threshold_text=False, line_w=0.8, min_gap_in=0.26)
ax_man.set_title(f"Manhattan · {len(gwas):,} SNPs · {len(loci)} loci (4 PCs)")
ax_man.set_xlabel("Chromosome (hg38)")
paths.append(slidefig.save_exact(fig, "practice_B3.png"))

# ------------------------------------------------------------------------------------------ contact sheet
frames = slidefig.FRAMES.set_index("media")
show = []
for p in paths:
    w_in = float(frames.loc[p.name, "w"]) if p.name in frames.index else 6.25     # on-slide width in inches
    show.append((p, w_in))
from PIL import ImageDraw, ImageFont  # noqa: E402

scale = 110                                                                     # px per slide inch
rows = [[show[0]], [show[2]], [show[3], show[1], show[4]]]
pad, cap_h = 24, 52
font = ImageFont.truetype(str(ROOT / "fonts" / "SUIT-Regular.ttf"), 18)
def img_size(path):
    with Image.open(path) as im:
        return im.size


size_on_sheet = lambda p, w: (int(w * scale), int(img_size(p)[1] * int(w * scale) / img_size(p)[0]))
canvas_w = max(sum(size_on_sheet(p, w)[0] for p, w in r) + pad * (len(r) + 1) for r in rows)
heights = [max(size_on_sheet(p, w)[1] for p, w in r) + cap_h + pad for r in rows]
sheet = Image.new("RGB", (canvas_w, sum(heights) + pad), "white")
draw = ImageDraw.Draw(sheet)
yy = pad
for r, h in zip(rows, heights):
    xx = pad
    for p, w in r:
        with Image.open(p) as im0:
            im = im0.convert("RGB")
        tw, th = size_on_sheet(p, w)
        sheet.paste(im.resize((tw, th), Image.LANCZOS), (xx, yy + cap_h))
        draw.rectangle([xx - 1, yy + cap_h - 1, xx + tw, yy + cap_h + th], outline=dviz.GREY)
        draw.text((xx, yy + 4), p.name, fill=dviz.INK, font=font)
        draw.text((xx, yy + 26), f"{im.size[0]} × {im.size[1]} px · {w:.2f} in wide on the slide",
                  fill=dviz.SLATE, font=font)
        xx += tw + pad
    yy += h
sheet.save(HERE / "contact_sheet_B3.png")
for p in paths:
    print(f"{p.relative_to(ROOT)}  {img_size(p)}")
print("build/contact_sheet_B3.png", img_size(HERE / "contact_sheet_B3.png"))
print(json.dumps({"lambda_noPC": round(lam_no, 3), "lambda_PC": round(lam_pc, 3), "n_loci": len(loci),
                  "loci": loci["nearest_gene"].tolist()}))
