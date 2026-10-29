"""Lecture-slide figures for the single-cell track (Practice B4 data) at exact pixel sizes.

  image-17-1.png   NKG7 feature plot: zeros drawn on top vs sorted, with a colour bar
  image-18-1.png   the same cells four ways: PCA, t-SNE, UMAP 5 and 200 neighbours, labelled directly
  image-19-1.png   exploded pie + legend vs sorted bars with direct % labels (composition, truth file)
  image-21-1.png   violin of NKG7 in CD8 T, NK, B, CD14 Mono: raw values vs tiny noise added
  image-41-1.png   UMAP with direct cluster labels and the settings in a footnote
  image-42-1.png   dot plot of marker genes x clusters
  practice_B4.png  the story of Practice B4 (1800 x 1210 px)

The analysis is recomputed from data/ with the notebook's helpers (build/scrna_helpers.py). Every number
written on a figure is read from results/B4.json (written by the executed notebook), and the script
asserts that its own recomputed values agree with that file.
usage:  python build/slides_single_cell.py      (run the notebook first)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.cluster import KMeans  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(HERE), str(ROOT)]
import dviz  # noqa: E402
import slidefig as SF  # noqa: E402
import scrna_helpers as H  # noqa: E402

DATA = ROOT / "data"
RES = json.loads((ROOT / "results" / "B4.json").read_text())
SN = {k: v["value"] for k, v in RES["slide_numbers"].items()}
CT = H.CELL_TYPES
COL = H.CT_COLORS
NAMES = {"pca2": "PCA (PC1, PC2)", "tsne30": "t-SNE, perplexity 30", "umap5": "UMAP, 5 neighbours",
         "umap15": "UMAP, 15 neighbours", "umap200": "UMAP, 200 neighbours"}


def analysis():
    """The notebook's pipeline, step by step, with the same settings and seeds."""
    counts = sp.load_npz(DATA / "scrna_counts.npz").tocsr()
    genes = pd.read_csv(DATA / "scrna_genes.tsv", sep="\t")
    cells = pd.read_csv(DATA / "scrna_cells.tsv", sep="\t")
    mt = pd.read_csv(DATA / "scrna_markers.tsv", sep="\t")
    markers = {r.cell_type: r.markers.split(",") for r in mt.itertuples()}
    truth = cells[["cell_id"]].merge(pd.read_csv(DATA / "truth" / "scrna_truth.tsv", sep="\t"), on="cell_id",
                                     how="left", validate="one_to_one")["cell_type"].to_numpy()
    keep = np.asarray((counts > 0).sum(axis=0)).ravel() >= 3
    C, g = counts[:, keep], genes["gene"].to_numpy()[keep]
    lognorm, scaled = H.normalise_log1p(C)
    hvg, _ = H.highly_variable_genes(scaled, n_top=1000)
    pcs, _ = H.run_pca(lognorm[:, hvg], n_pcs=30, seed=0)
    cluster = KMeans(n_clusters=8, n_init=10, random_state=0).fit_predict(pcs)
    dot_genes = list(dict.fromkeys(x for v in markers.values() for x in v))
    gi = [int(np.flatnonzero(g == x)[0]) for x in dot_genes]
    pct, mean = H.dotplot_table(C[:, gi].toarray(), lognorm[:, gi], cluster, dot_genes)
    names, _ = H.annotate_clusters(mean, markers)
    label = pd.Series(cluster).map(names).to_numpy()
    emb = {"pca2": pcs[:, :2], "tsne30": H.run_tsne(pcs, perplexity=30, seed=0)[0]}
    for k in (15, 5, 200):
        emb[f"umap{k}"] = H.run_umap(pcs, n_neighbors=k, min_dist=0.3, seed=0)[0]
    A = SimpleNamespace(C=C, genes=g, lognorm=lognorm, pcs=pcs, cluster=cluster, names=names, label=label,
                        truth=truth, emb=emb, cells=cells)
    # the figures must show exactly what the executed notebook reported
    assert round(adjusted_rand_score(truth, cluster), 3) == SN["ari"]
    assert round(100 * float((label == truth).mean()), 1) == SN["pct_correctly_labelled"]
    assert np.allclose(emb["umap15"], cells[["umap1", "umap2"]].to_numpy(), atol=1e-3)   # = stored UMAP
    ref = H.centroid_distances(pcs, label)
    for k, e in emb.items():
        assert round(100 * H.knn_overlap(pcs, e, k=15), 2) == SN[f"knn15_kept_{k}_pct"], k
        assert round(float(spearmanr(ref, H.centroid_distances(e, label)).statistic), 2) == SN[f"centre_rho_{k}"], k
    return A


def gene(A, name):
    return A.lognorm[:, int(np.flatnonzero(A.genes == name)[0])]


def blank_axes(ax):
    ax.set(xticks=[], yticks=[])


def subtitle(ax, text, size=10.5):
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(0, 4), textcoords="offset points",
                ha="left", va="bottom", fontsize=size, color=dviz.SLATE)


# PCA (PC1, PC2): the two monocyte clusters overlap, so these labels are placed by hand
# (x, y as fractions of the data range, horizontal alignment)
PCA_LABELS = {"NK": (0.32, 0.95, "left"), "CD8 T": (0.24, 0.76, "left"), "Memory CD4 T": (0.27, 0.53, "left"),
              "Naive CD4 T": (0.24, 0.29, "left"), "B": (0.44, 0.07, "left"), "FCGR3A Mono": (1.0, 0.74, "right"),
              "CD14 Mono": (1.0, 0.25, "right"), "DC": (0.66, 0.2, "right")}


# UMAP (15 neighbours) of slides 41 and practice_B4: labels placed by hand next to each island
UMAP15_LABELS = {"DC": (0.07, 0.86, "left"), "NK": (0.70, 0.97, "left"), "CD8 T": (0.47, 0.70, "center"),
                 "Memory CD4 T": (0.70, 0.39, "center"), "Naive CD4 T": (1.0, 0.68, "right"),
                 "B": (0.67, 0.04, "left"), "CD14 Mono": (0.21, 0.23, "left"), "FCGR3A Mono": (0.10, 0.5, "left")}


def manual_labels(ax, e, labels, spec, fontsize):
    import matplotlib.patheffects as pe
    mn, rg = e.min(axis=0), np.ptp(e, axis=0)
    for name, (fx, fy, ha) in spec.items():
        anchor = np.median(e[labels == name], axis=0)
        ax.annotate(name, xy=anchor, xytext=(mn[0] + fx * rg[0], mn[1] + fy * rg[1]), ha=ha, va="center",
                    fontsize=fontsize, color=dviz.INK, path_effects=[pe.withStroke(linewidth=3, foreground="white")],
                    arrowprops=dict(arrowstyle="-", color=dviz.SLATE, lw=0.6, shrinkA=1, shrinkB=0))
    ax.set(xlim=(mn[0] - 0.05 * rg[0], mn[0] + 1.05 * rg[0]), ylim=(mn[1] - 0.06 * rg[1], mn[1] + 1.06 * rg[1]))


def types_scatter(ax, e, labels, s=2.0):
    ax.scatter(e[:, 0], e[:, 1], s=s, c=[COL[t] for t in labels], linewidths=0, rasterized=False)


# ------------------------------------------------------------------------------------------------ 17
def fig17(A):
    fig, ax0 = SF.new_fig("image-17-1.png")
    ax0.remove()
    e, v = A.emb["umap15"], gene(A, "NKG7")
    orders = [("High values drawn first (hidden)", np.argsort(-v, kind="stable"), SN["visible_nkg7_zeros_on_top_pct"]),
              ("High values drawn last", np.argsort(v, kind="stable"), SN["visible_nkg7_sorted_pct"])]
    for i, (title, order, seen) in enumerate(orders):
        ax = fig.add_axes([0.065 + i * 0.43, 0.11, 0.37, 0.68])
        dots = ax.scatter(e[order, 0], e[order, 1], c=v[order], s=3.2, cmap="viridis", vmin=0, vmax=float(v.max()),
                          linewidths=0)
        ax.set_title(title, pad=24)
        subtitle(ax, f"NKG7-positive cells visible in {seen:.0f}% of the cell area", size=10.5)
        ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
        blank_axes(ax)
        for name in ["NK", "CD8 T", "Memory CD4 T"]:
            pts = e[A.label == name]
            if name == "NK":                                   # top cluster: label to its right
                xy, ha, va = (pts[:, 0].max() + 0.02 * np.ptp(e[:, 0]), np.median(pts[:, 1])), "left", "center"
            else:
                xy, ha, va = (np.median(pts[:, 0]), pts[:, 1].max() + 0.03 * np.ptp(e[:, 1])), "center", "bottom"
            ax.text(*xy, name, ha=ha, va=va, fontsize=9.5, color=dviz.SLATE)
        ax.margins(0.06)
    cax = fig.add_axes([0.905, 0.2, 0.014, 0.52])
    cb = fig.colorbar(dots, cax=cax)
    cb.set_label("NKG7, log(1 + UMIs per 10,000)")
    return SF.save_exact(fig, "image-17-1.png")


# ------------------------------------------------------------------------------------------------ 18
def fig18(A):
    fig, ax0 = SF.new_fig("image-18-1.png")
    ax0.remove()
    for i, k in enumerate(["pca2", "tsne30", "umap5", "umap200"]):
        ax = fig.add_axes([0.012 + i * 0.249, 0.03, 0.232, 0.74])
        types_scatter(ax, A.emb[k], A.label, s=1.8)
        if k == "pca2":
            manual_labels(ax, A.emb[k], A.label, PCA_LABELS, fontsize=8.5)
        else:
            H.label_clusters(ax, A.emb[k], A.label, fontsize=8.5)
        ax.set_title(NAMES[k], pad=22)
        subtitle(ax, f"15-NN kept {SN[f'knn15_kept_{k}_pct']:.0f}% · centre-distance ρ {SN[f'centre_rho_{k}']:.2f}")
        blank_axes(ax)
    return SF.save_exact(fig, "image-18-1.png")


# ------------------------------------------------------------------------------------------------ 19
def fig19(A):
    fig, ax0 = SF.new_fig("image-19-1.png")
    ax0.remove()
    share = pd.Series(A.truth).value_counts(normalize=True).mul(100)
    assert round(share["CD14 Mono"], 1) == SN["pct_cd14_mono"] and round(share["Memory CD4 T"], 1) == SN["pct_memory_cd4_t"]
    ax = fig.add_axes([0.055, 0.05, 0.235, 0.78])
    wedges, _ = ax.pie(share.loc[CT], colors=[COL[t] for t in CT], explode=[0.07] * len(CT), shadow=True,
                       startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=1))
    ax.set_title("Exploded pie + legend", pad=8)
    ax.legend(wedges, CT, loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=10.5, handlelength=1.4,
              labelspacing=0.55)
    bx = fig.add_axes([0.605, 0.17, 0.355, 0.66])
    s = share.sort_values()
    colors = [dviz.BLUE if t in ("CD14 Mono", "Memory CD4 T") else dviz.GREY for t in s.index]
    bx.barh(np.arange(len(s)), s.to_numpy(), color=colors, height=0.66)
    for y, x in enumerate(s.to_numpy()):
        bx.text(x + 0.4, y, f"{x:.1f}%", va="center", fontsize=11, color=dviz.INK)
    bx.set_yticks(np.arange(len(s)), s.index)
    bx.set(xlim=(0, 30), xlabel=f"Cells (%) · n = {len(A.truth):,}")
    bx.set_title("Sorted bars, direct labels", pad=8)
    return SF.save_exact(fig, "image-19-1.png")


# ------------------------------------------------------------------------------------------------ 21
def fig21(A):
    fig, ax0 = SF.new_fig("image-21-1.png")
    ax0.remove()
    groups = ["CD8 T", "NK", "B", "CD14 Mono"]
    v = gene(A, "NKG7")
    data = [v[A.truth == t] for t in groups]
    rng = np.random.default_rng(20261027 + 5)
    noisy = [d + rng.normal(0, 0.03, d.size) for d in data]         # visible stand-in for Seurat's tiny noise
    axes = [fig.add_axes([0.085, 0.2, 0.41, 0.64]), fig.add_axes([0.575, 0.2, 0.41, 0.64])]
    for ax, values, title in [(axes[0], data, "Raw values (all-zero groups vanish)"),
                              (axes[1], noisy, "Tiny noise added (shape invented)")]:
        for i, (t, d) in enumerate(zip(groups, values)):
            if np.ptp(d) > 0:
                part = ax.violinplot([d], positions=[i], widths=0.8, showextrema=False)
                for body in part["bodies"]:
                    body.set(facecolor=COL[t], edgecolor=COL[t], alpha=0.55, linewidth=1)
            else:
                ax.hlines(0, i - 0.4, i + 0.4, color=dviz.INK, lw=2.2)
            ax.scatter(i + dviz.jitter(d.size, 0.22, seed=i), d, s=1.2, color=dviz.INK, alpha=0.35, linewidths=0)
        ax.set_xticks(range(len(groups)), [f"{t}\nn = {np.sum(A.truth == t)}" for t in groups])
        ax.set_title(title, pad=8)
        ax.set_ylim(-0.3, float(v.max()) + 0.3)
    axes[0].set_ylabel("NKG7, log(1 + UMIs\nper 10,000)")
    axes[1].tick_params(labelleft=False)
    return SF.save_exact(fig, "image-21-1.png")


# ------------------------------------------------------------------------------------------------ 41
def fig41(A):
    fig, ax = SF.new_fig("image-41-1.png")
    fig.subplots_adjust(left=0.075, right=0.985, top=0.9, bottom=0.235)
    types_scatter(ax, A.emb["umap15"], A.label, s=3)
    manual_labels(ax, A.emb["umap15"], A.label, UMAP15_LABELS, fontsize=11)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    blank_axes(ax)
    ax.set_title(f"{SN['n_cells']:,} cells · 8 clusters named from marker genes", pad=8)
    fig.text(0.012, 0.015, f"{SN['umap_settings']} · {SN['umap_runtime_s']:.1f} s\n"
             f"clusters: k-means (k = 8) on the 30 PCs · {SN['pct_correctly_labelled']:.1f}% named correctly "
             f"(simulated truth)\ngaps, sizes and positions of the islands are not data", fontsize=8.5,
             color=dviz.SLATE, va="bottom", linespacing=1.35)
    return SF.save_exact(fig, "image-41-1.png")


# ------------------------------------------------------------------------------------------------ 42
def fig42(A):
    fig, ax = SF.new_fig("image-42-1.png")
    fig.subplots_adjust(left=0.2, right=0.79, top=0.87, bottom=0.2)
    genes = ["CCR7", "IL7R", "S100A4", "CD3D", "CD8A", "NKG7", "GNLY", "MS4A1", "CD79A", "CD14", "LYZ",
             "FCGR3A", "MS4A7", "FCER1A"]
    gi = [int(np.flatnonzero(A.genes == x)[0]) for x in genes]
    pct, mean = H.dotplot_table(A.C[:, gi].toarray(), A.lognorm[:, gi], A.cluster, genes)
    rows = [A.names[A.names == t].index[0] for t in CT]
    size, level = pct.loc[rows].to_numpy(), mean.loc[rows]
    level01 = ((level - level.min()) / (level.max() - level.min())).to_numpy()
    xx, yy = np.meshgrid(np.arange(len(genes)), np.arange(len(rows)))
    k = 2.0                                                            # dot area (pt^2) per % expressing
    dots = ax.scatter(xx.ravel(), yy.ravel(), s=k * size.ravel(), c=level01.ravel(), cmap="viridis", vmin=0, vmax=1,
                      edgecolors=dviz.SLATE, linewidths=0.4)
    ax.set_xticks(range(len(genes)), genes, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(rows)), CT)
    ax.set(xlim=(-0.6, len(genes) - 0.4), ylim=(len(rows) - 0.4, -0.6))
    ax.set_title("Dot plot · size = % cells expressing, colour = mean", pad=8)
    for p in (25, 50, 100):
        ax.scatter([], [], s=k * p, color=dviz.GREY, edgecolors=dviz.SLATE, linewidths=0.4, label=f"{p}%")
    ax.legend(title="expressing", loc="upper left", bbox_to_anchor=(1.03, 1.03), labelspacing=1.0, fontsize=10.5,
              title_fontsize=10.5, borderaxespad=0)
    cax = fig.add_axes([0.83, 0.2, 0.016, 0.3])
    cb = fig.colorbar(dots, cax=cax, ticks=[0, 0.5, 1])
    cb.set_label("mean, scaled\nper gene", fontsize=10.5)
    cb.ax.tick_params(labelsize=10)
    return SF.save_exact(fig, "image-42-1.png")


# ------------------------------------------------------------------------------------------------ practice
def practice(A):
    fig, ax0 = SF.new_fig("practice_B4.png")
    ax0.remove()
    ax = fig.add_axes([0.065, 0.19, 0.47, 0.63])
    types_scatter(ax, A.emb["umap15"], A.label, s=1.8)
    manual_labels(ax, A.emb["umap15"], A.label, UMAP15_LABELS, fontsize=8.5)
    blank_axes(ax)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    ax.set_title("Clusters are real", pad=22)
    subtitle(ax, f"k-means on 30 PCs: ARI {SN['ari']:.2f} · {SN['pct_correctly_labelled']:.1f}% named correctly",
             size=10)
    keys = ["pca2", "tsne30", "umap5", "umap15", "umap200"]
    short = {"pca2": "PC1-PC2", "tsne30": "t-SNE 30", "umap5": "UMAP 5", "umap15": "UMAP 15", "umap200": "UMAP 200"}
    colors = [dviz.SLATE if k == "pca2" else dviz.BLUE for k in keys]
    y = np.arange(len(keys))[::-1]
    for box, metric, title, xlabel, xlim, fmt in [
            ([0.72, 0.625, 0.26, 0.215], "knn15_kept_{}_pct", "Neighbours are kept", "15-NN kept (%)", (0, 36), "{:.0f}%"),
            ([0.72, 0.19, 0.26, 0.215], "centre_rho_{}", "Distances are not",
             "centre-distance ρ\n(1 = same order as 30 PCs)", (0, 1.15), "{:.2f}")]:
        bx = fig.add_axes(box)
        vals = [SN[metric.format(k)] for k in keys]
        bx.barh(y, vals, color=colors, height=0.66)
        for yi, x in zip(y, vals):
            bx.text(x + 0.015 * xlim[1], yi, fmt.format(x), va="center", fontsize=9.5, color=dviz.INK)
        bx.set_yticks(y, [short[k] for k in keys], fontsize=10)
        bx.set_xlim(*xlim)
        bx.tick_params(axis="x", labelsize=9.5)
        bx.set_xlabel(xlabel, fontsize=10.5)
        bx.set_title(title, fontsize=12.5, pad=6)
    fig.axes[-2].text(0.98, 0.97, f"chance {SN['knn15_chance_pct']:.1f}%", transform=fig.axes[-2].transAxes,
                      ha="right", va="top", fontsize=9, color=dviz.SLATE)
    fig.text(0.02, 0.015, f"{SN['n_cells']:,} simulated PBMC-like cells · {SN['n_genes']:,} genes\n"
             "UMAP of 30 PCs · n_neighbors 15 · min_dist 0.3 · seed 0", fontsize=8.5, color=dviz.SLATE,
             linespacing=1.3)
    return SF.save_exact(fig, "practice_B4.png")


if __name__ == "__main__":
    A = analysis()
    for f in (fig17, fig18, fig19, fig21, fig41, fig42, practice):
        p = f(A)
        print("wrote", p.relative_to(ROOT), SF.size_of(p.name))
