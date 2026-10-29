"""Slide figures drawn from the core practice data (Practice A, B1, B2) at the exact size of the deck images.

usage: python build/slides_core.py            -> build/slide_figs/image-<slide>-<k>.png
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
import dviz  # noqa: E402
import slidefig as S  # noqa: E402

D = HERE.parent / "data"
DPI = 288
BWO = LinearSegmentedColormap.from_list("bwo", [dviz.DARKBLUE, dviz.SKY, "white", "#F2C185", dviz.VERMILION])


def scale_of(target):
    if target.startswith("practice_"):
        return 1.0
    r = S.FRAMES.loc[S.FRAMES["media"] == target].iloc[0]
    return float(r["w"]) * DPI / S.size_of(target)[0]


def fig_for(target, nrows=1, ncols=1, on_slide=12.0, **kw):
    """Exact-size figure; base font chosen so that text is ~on_slide pt when shown in the deck frame."""
    base = min(max(on_slide / scale_of(target), on_slide), 34)
    return S.new_fig(target, nrows, ncols, base_size=base, **kw)


# ---------------------------------------------------------------------------------------------- data
ebg = pd.read_csv(D / "expression_by_group.tsv", sep="\t")
mr = pd.read_csv(D / "mapping_rates.tsv", sep="\t")
tc = pd.read_csv(D / "timecourse.tsv", sep="\t")
comp = pd.read_csv(D / "celltype_composition.tsv", sep="\t")
me = pd.read_csv(D / "methylation_cpg.tsv.gz", sep="\t")
de = pd.read_csv(D / "de_results.tsv", sep="\t")
meta = pd.read_csv(D / "rnaseq_metadata.tsv", sep="\t")
counts = pd.read_csv(D / "rnaseq_counts.tsv.gz", sep="\t", index_col="gene").drop(columns=["gene_id", "chrom"])[meta["sample"]]
expressed = (counts > 0).all(axis=1)
lg = np.log(counts[expressed])
sf = np.exp(lg.sub(lg.mean(axis=1), axis=0).median(axis=0))
norm = counts / sf
log2n = np.log2(norm + 1)
tested = de["padj"].notna()
up = tested & (de["padj"] < 0.05) & (de["log2FC_shrunk"] > 1)
dn = tested & (de["padj"] < 0.05) & (de["log2FC_shrunk"] < -1)


def groups(ds):
    d = ebg[ebg["dataset"] == ds]
    return [d.loc[d["group"] == g, "expression"].to_numpy() for g in ["Control", "Case"]]


# ---------------------------------------------------------------------------------------------- helpers
def box_points(ax, vals, violin=False, width=0.35, s=18, offset=0.0):
    for i, (v, c) in enumerate(zip(vals, [dviz.SLATE, dviz.BLUE])):
        if violin:
            vp = ax.violinplot([v], positions=[i], showextrema=False, widths=0.7)
            for b in vp["bodies"]:
                b.set_facecolor(dviz.PALEBLUE); b.set_edgecolor("none"); b.set_alpha(1)
        ax.boxplot([v], positions=[i], widths=0.14 if violin else width, showfliers=False,
                   medianprops={"color": dviz.INK, "lw": 1.6}, whiskerprops={"color": dviz.INK},
                   capprops={"color": dviz.INK}, patch_artist=True, boxprops={"facecolor": "white", "color": dviz.INK})
        ax.scatter(i + offset + dviz.jitter(len(v), 0.07, seed=i), v, s=s, color=c, zorder=3, edgecolor="white", lw=0.3)


def volcano(ax, labels=8, legend_counts=True, s=(5, 12), fs=None):
    x, y = de["log2FC_shrunk"], -np.log10(de["padj"])
    rest = tested & ~up & ~dn
    ax.scatter(x[rest], y[rest], s=s[0], color=dviz.GREY, rasterized=True, lw=0)
    ax.scatter(x[up], y[up], s=s[1], color=dviz.VERMILION, lw=0)
    ax.scatter(x[dn], y[dn], s=s[1], color=dviz.OI_BLUE, lw=0)
    ax.axhline(-np.log10(0.05), color=dviz.SLATE, ls="--", lw=0.9)
    for v in (-1, 1):
        ax.axvline(v, color=dviz.SLATE, ls="--", lw=0.9)
    lim = np.ceil(np.nanmax(np.abs(x[tested])))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-2, y[tested].max() * 1.12)
    if labels:
        top = de[tested].nsmallest(labels, "padj")
        xs, ys = top["log2FC_shrunk"].to_numpy(), -np.log10(top["padj"].to_numpy())
        texts = [ax.text(xx - 0.35 if xx > 0 else xx + 0.35, yy, g, fontsize=fs or mpl.rcParams["font.size"] * 0.72,
                         style="italic", color=dviz.INK, ha="right" if xx > 0 else "left", va="center")
                 for xx, yy, g in zip(xs, ys, top["gene"])]
        from adjustText import adjust_text
        sx = de.loc[up | dn, "log2FC_shrunk"].to_numpy(); sy = -np.log10(de.loc[up | dn, "padj"].to_numpy())
        adjust_text(texts, x=sx, y=sy, ax=ax, expand=(1.15, 1.5), force_text=(0.4, 0.6), force_static=(0.2, 0.4),
                    only_move={"text": "xy", "static": "xy", "explode": "xy", "pull": "xy"},
                    arrowprops=dict(arrowstyle="-", color=dviz.SLATE, lw=0.6))
    if legend_counts:
        f = mpl.rcParams["font.size"] * 0.85
        ax.text(0.02, 0.97, f"{dn.sum()} down", transform=ax.transAxes, color=dviz.OI_BLUE, fontsize=f, weight="bold", va="top")
        ax.text(0.02, 0.89, f"{up.sum()} up", transform=ax.transAxes, color=dviz.VERMILION, fontsize=f, weight="bold", va="top")
    ax.set_xlabel("Shrunken log2 fold change (treated vs control)")
    ax.set_ylabel("−log10 adjusted p (BH)")


def bad_volcano(ax, s=6):
    x, y = de["log2FC"], -np.log10(de["pvalue"].fillna(1.0))
    lim = np.nanpercentile(np.abs(x[tested]), 99.5)
    keep = np.abs(x) < 8
    ax.scatter(x[keep & (x < 0)], y[keep & (x < 0)], s=s, color="#00B000", rasterized=True, lw=0)
    ax.scatter(x[keep & (x >= 0)], y[keep & (x >= 0)], s=s, color="#E00000", rasterized=True, lw=0)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("−log10 p")


def pca_scatter(ax, s=110, legend=True, title=None, lfs=None):
    keep = norm.mean(axis=1) >= 10
    top = log2n[keep].var(axis=1).nlargest(1000).index
    X = log2n.loc[top].T
    X = X - X.mean()
    u, sv, vt = np.linalg.svd(X.to_numpy(), full_matrices=False)
    P = u * sv
    vr = sv ** 2 / np.sum(sv ** 2) * 100
    for cond, col in [("control", dviz.SLATE), ("treated", dviz.BLUE)]:
        for b, mk in [("A", "o"), ("B", "^")]:
            m = ((meta["condition"] == cond) & (meta["batch"] == b)).to_numpy()
            ax.scatter(P[m, 0], P[m, 1], color=col, marker=mk, s=s, edgecolor="white", lw=0.8, label=f"{cond}, batch {b}")
    ax.axhline(0, color=dviz.LIGHT, lw=1, zorder=0); ax.axvline(0, color=dviz.LIGHT, lw=1, zorder=0)
    ax.set_xlabel(f"PC1 ({vr[0]:.0f}% variance)"); ax.set_ylabel(f"PC2 ({vr[1]:.0f}% variance)")
    if title:
        ax.set_title(title)
    if legend:
        ax.legend(fontsize=lfs or mpl.rcParams["font.size"] * 0.8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2), handletextpad=0.3)
    return vr


def heatmap(fig, rect, labels=True, dendro=True, cbar=True, title=True):
    top40 = de[tested].nsmallest(40, "padj")["gene"]
    Z = log2n.loc[top40]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0)
    rl = linkage(Z.to_numpy(), "average", metric="correlation")
    cl = linkage(Z.T.to_numpy(), "average", metric="correlation")
    ro, co = leaves_list(rl), leaves_list(cl)
    Zo = Z.iloc[ro, co]
    x0, y0, w, h = rect
    dw = 0.07 * w if dendro else 0
    cw = 0.018 if cbar else 0
    fs = mpl.rcParams["font.size"]
    fig_w = fig.get_size_inches()[0]
    # right margin for the colour-bar tick labels ("−2") and its rotated label, from the font sizes in use
    lab_w = ((fs * 0.65 * 1.6 + fs * 0.75 * 1.5) / 72 + 0.12) / fig_w if cbar else 0
    row_w = (fs * 0.6 * 1.6 / 72 + 0.06) / fig_w if labels else 0     # rotated row label right of the heatmap
    gap = cw + lab_w + row_w + 0.01 if cbar else row_w
    hm = fig.add_axes([x0 + dw, y0, w - dw - gap, h * 0.74])
    im = hm.imshow(Zo, aspect="auto", cmap=BWO, norm=TwoSlopeNorm(0, -2.5, 2.5), interpolation="nearest")
    if labels:
        hm.set_yticks([])                                   # 40 names are unreadable at slide size
        hm.yaxis.set_label_position("right")                # left of the heatmap holds the row dendrogram
        hm.set_ylabel("40 genes (names in the notebook)", fontsize=fs * 0.6, color=dviz.SLATE, labelpad=4)
        hm.set_xticks(range(Zo.shape[1]), Zo.columns, rotation=90, fontsize=fs * 0.6)
    else:
        hm.set_yticks([]); hm.set_xticks([])
    hm.tick_params(length=0)
    m2 = meta.set_index("sample").loc[Zo.columns]
    for k, (col, cmap) in enumerate([("condition", {"control": dviz.SLATE, "treated": dviz.BLUE}),
                                     ("batch", {"A": dviz.GREY, "B": dviz.ORANGE})]):
        a = fig.add_axes([x0 + dw, y0 + h * (0.765 + 0.05 * k), w - dw - gap, h * 0.04])
        a.imshow([[mpl.colors.to_rgb(cmap[v]) for v in m2[col]]], aspect="auto")
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
        if labels:
            a.text(-0.6, 0, col, ha="right", va="center", fontsize=fs * 0.6, color=dviz.SLATE)
    if dendro:
        ad = fig.add_axes([x0, y0, dw * 0.95, h * 0.74])
        dendrogram(rl, orientation="left", ax=ad, color_threshold=0, above_threshold_color=dviz.SLATE, no_labels=True)
        ad.invert_yaxis(); ad.axis("off")
        at = fig.add_axes([x0 + dw, y0 + h * 0.87, w - dw - gap, h * 0.11])
        dendrogram(cl, ax=at, color_threshold=0, above_threshold_color=dviz.SLATE, no_labels=True)
        at.axis("off")
    if cbar:
        ca = fig.add_axes([x0 + w - cw - lab_w, y0, cw, h * 0.74])
        cb = fig.colorbar(im, cax=ca)
        cb.set_label("row z-score", fontsize=fs * 0.75)
        cb.ax.tick_params(labelsize=fs * 0.65)
    if title:
        fig.text(x0, y0 + h + 0.005, "Top 40 DE genes · rows z-scored · 1 − Pearson r, average linkage",
                 fontsize=fs * 0.85, weight="bold", color=dviz.INK, va="bottom")
    return Zo


# ---------------------------------------------------------------------------------------------- figures
def f14():
    t = "image-14-1.png"
    fig, axs = fig_for(t, 1, 3)
    small, large = groups("small_n"), groups("large_n")
    m = [v.mean() for v in small]; se = [v.std(ddof=1) / np.sqrt(len(v)) for v in small]
    axs[0].bar([0, 1], m, yerr=se, color=[dviz.GREY, dviz.BLUE], width=0.6, capsize=6, error_kw={"lw": 1.2})
    axs[0].text(1, m[1] + se[1] + 6, "*", ha="center", fontsize=mpl.rcParams["font.size"] * 1.6)
    axs[0].set_title("Bar ± SE (n = 4 + 4)")
    box_points(axs[1], small, s=34)
    axs[1].set_title("Same data, every point")
    box_points(axs[2], large, s=16)
    axs[2].set_title("Another study, n = 40 + 40")
    for ax in axs:
        ax.set_xticks([0, 1], ["Control", "Case"]); ax.set_ylim(0, 160)
    axs[0].set_ylabel("Expression (TPM)")
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.12, top=0.88, wspace=0.25)
    return S.save_exact(fig, t)


def f15():
    t = "image-15-1.png"
    fig, axs = fig_for(t, 1, 2)
    means = mr.groupby("pipeline")["pct_mapped"].mean()
    for ax, lo, title in [(axs[0], 95.5, "Axis starts at 95.5%"), (axs[1], 0, "Axis starts at 0")]:
        ax.bar(means.index, means.values, color=[dviz.GREY, dviz.GREY, dviz.BLUE], width=0.6)
        ax.set_ylim(lo, 100)
        for i, v in enumerate(means.values):
            if lo > 0:
                ax.text(i, v + (100 - lo) * 0.012, f"{v:.1f}", ha="center", va="bottom", fontsize=mpl.rcParams["font.size"] * 0.9)
            else:
                ax.text(i, v - 4, f"{v:.1f}", ha="center", va="top", fontsize=mpl.rcParams["font.size"] * 0.9,
                        color="white" if i == 2 else dviz.INK)
        ax.set_title(title); ax.set_ylabel("Reads mapped (%)")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.12, top=0.88, wspace=0.32)
    return S.save_exact(fig, t)


def f22():
    t = "image-22-1.png"
    fig, axs = fig_for(t, 1, 2)
    bad_volcano(axs[0]); axs[0].set_title("Before")
    volcano(axs[1], labels=6); axs[1].set_title("After")
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.17, top=0.9, wspace=0.22)
    return S.save_exact(fig, t)


def f27():
    t = "image-27-1.png"
    fig, axs = fig_for(t, 1, 4)
    case = groups("large_n")[1]
    for ax, b in zip(axs[:3], [3, 12, 60]):
        ax.hist(case, bins=b, range=(40, 130), color=dviz.SKY, edgecolor="white", density=True)
        ax.set_title(f"{b} bins")
    grid = np.linspace(30, 145, 300)
    for bw, col, ls in [(0.15, dviz.ORANGE, "-"), (0.4, dviz.BLUE, "-"), (1.2, dviz.SLATE, "--")]:
        axs[3].plot(grid, stats.gaussian_kde(case, bw_method=bw)(grid), color=col, ls=ls, lw=1.8, label=f"bw {bw}")
    axs[3].legend(fontsize=mpl.rcParams["font.size"] * 0.75); axs[3].set_title("KDE bandwidth")
    for ax in axs:
        ax.set_xlabel("Expression (TPM)")
    axs[0].set_ylabel("Density")
    fig.subplots_adjust(left=0.05, right=0.99, bottom=0.18, top=0.88, wspace=0.28)
    return S.save_exact(fig, t)


def f28():
    t = "image-28-1.png"
    fig, axs = fig_for(t, 1, 4, sharey=True)
    vals = groups("large_n")
    axs[0].boxplot(vals, positions=[0, 1], widths=0.5, patch_artist=True, medianprops={"color": dviz.INK, "lw": 1.8},
                   boxprops={"facecolor": dviz.LIGHT, "color": dviz.SLATE})
    vp = axs[1].violinplot(vals, positions=[0, 1], showmedians=True)
    for b in vp["bodies"]:
        b.set_facecolor(dviz.PALEBLUE); b.set_edgecolor(dviz.BLUE); b.set_alpha(1)
    for k in ["cbars", "cmins", "cmaxes", "cmedians"]:
        vp[k].set_color(dviz.INK)
    for i, (v, c) in enumerate(zip(vals, [dviz.SLATE, dviz.BLUE])):
        axs[2].scatter(i + dviz.jitter(len(v), 0.16, seed=i), v, s=16, color=c)
    box_points(axs[3], vals, violin=True, s=10, offset=0.3)
    for ax, ttl in zip(axs, ["Box", "Violin", "Strip (jitter)", "Violin + box + points"]):
        ax.set_title(ttl); ax.set_xticks([0, 1], ["Control", "Case"])
    axs[0].set_ylabel("Expression (TPM)")
    fig.subplots_adjust(left=0.05, right=0.99, bottom=0.12, top=0.88, wspace=0.12)
    return S.save_exact(fig, t)


def f31():
    t = "image-31-1.png"
    fig, axs = fig_for(t, 1, 3, sharey=True)
    x, y = me["beta_wgbs"], me["beta_nanopore"]
    axs[0].scatter(x, y, s=1.5, color=dviz.BLUE, rasterized=True, lw=0); axs[0].set_title(f"{len(me):,} points, opaque")
    axs[1].scatter(x, y, s=1.5, color=dviz.BLUE, alpha=0.04, rasterized=True, lw=0); axs[1].set_title("Transparency (alpha 0.04)")
    hb = axs[2].hexbin(x, y, gridsize=45, bins="log", cmap="viridis", mincnt=1)
    cb = fig.colorbar(hb, ax=axs[2], fraction=0.05, pad=0.02); cb.set_label("CpGs (log)")
    axs[2].set_title(f"Hexbin density · r = {x.corr(y):.3f}")
    for ax in axs:
        ax.set_xlabel("Bisulfite (WGBS) β"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    axs[0].set_ylabel("Nanopore β")
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.9, wspace=0.12)
    return S.save_exact(fig, t)


def f32a():
    t = "image-32-1.png"
    fig, ax = fig_for(t)
    for g, col in zip(["TSC22D3", "FKBP5", "PER1"], [dviz.BLUE, dviz.ORANGE, dviz.SLATE]):
        d = tc[tc["gene"] == g]
        m = d.groupby("time_h")["log2fc_vs_0h"].agg(["mean", "std"])
        ax.scatter(d["time_h"], d["log2fc_vs_0h"], color=col, s=14, alpha=0.45, lw=0)
        ax.errorbar(m.index, m["mean"], yerr=m["std"], color=col, marker="o", ms=6, capsize=3, lw=2)
        ax.text(24.8, m["mean"].iloc[-1], g, color=col, va="center", style="italic")
    ax.axhline(0, color=dviz.LIGHT, lw=1, zorder=0)
    ax.set_xticks([0, 1, 2, 4, 8, 24]); ax.set_xlim(-0.6, 30)
    ax.set_xlabel("Hours after dexamethasone"); ax.set_ylabel("log2 FC vs 0 h")
    ax.set_title("Line graph · mean ± SD, n = 3")
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.17, top=0.88)
    return S.save_exact(fig, t)


def f32b():
    t = "image-32-2.png"
    fig, ax = fig_for(t)
    order = ["Naive CD4 T", "CD14 Mono", "Memory CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono", "DC"]
    wide = comp.pivot(index="donor", columns="cell_type", values="percent")[order]
    left = np.zeros(len(wide))
    cols = [dviz.BLUE, dviz.ORANGE, dviz.GREEN, dviz.VERMILION, dviz.OI_SKY, dviz.PURPLE, "#8C6D31", dviz.SLATE]
    for ct, c in zip(order, cols):
        ax.barh(wide.index, wide[ct], left=left, color=c, label=ct, height=0.7)
        left += wide[ct].to_numpy()
    ax.set_xlim(0, 100); ax.set_xlabel("Cells (%)"); ax.set_title("Composition per donor")
    ax.legend(ncol=2, fontsize=mpl.rcParams["font.size"] * 0.72, loc="upper center", bbox_to_anchor=(0.5, -0.28),
              handlelength=1.2, columnspacing=1.0, frameon=False)
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.42, top=0.9)
    return S.save_exact(fig, t)


def f33():
    t = "image-33-1.png"
    fig, axs = fig_for(t, 1, 2, sharey=True)
    A = np.log2(de["baseMean"] + 1)
    sig = up | dn
    for ax, col, ttl, note in [(axs[0], "log2FC", "Raw log2 fold change", "wide scatter at low counts"),
                               (axs[1], "log2FC_shrunk", "Shrunken log2 fold change", "low-information genes pulled toward 0")]:
        ax.scatter(A[tested & ~sig], de.loc[tested & ~sig, col], s=3, color=dviz.GREY, rasterized=True, lw=0)
        ax.scatter(A[sig], de.loc[sig, col], s=7, color=dviz.BLUE, lw=0)
        ax.axhline(0, color=dviz.SLATE, lw=0.9)
        ax.set_title(ttl); ax.set_xlabel("A = log2 mean normalised count")
        ax.text(0.03, 0.93, note, transform=ax.transAxes, color=dviz.VERMILION, fontsize=mpl.rcParams["font.size"] * 0.85)
    axs[0].set_ylabel("M = log2 fold change"); axs[0].set_ylim(-5.5, 5.5)
    fig.subplots_adjust(left=0.06, right=0.99, bottom=0.16, top=0.9, wspace=0.06)
    return S.save_exact(fig, t)


def f_volcano(t, labels=8, on_slide=12.0):
    fig, ax = fig_for(t, on_slide=on_slide)
    volcano(ax, labels=labels)
    ax.set_title(f"FDR < 0.05 and |log2FC| > 1 (n = 6 vs 6)")
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.17, top=0.9)
    return S.save_exact(fig, t)


def f35():
    t = "image-35-1.png"
    fig, ax = fig_for(t)
    p = de.loc[tested & (de["baseMean"] >= 10), "pvalue"].to_numpy()
    pi0 = min(1.0, float((p > 0.5).mean() / 0.5))
    ax.hist(p, bins=40, range=(0, 1), color=dviz.SKY, edgecolor="white")
    floor = len(p) / 40 * pi0
    ax.axhline(floor, color=dviz.BLUE, ls="--", lw=1.4)
    ax.text(0.98, floor * 1.35, f"null genes ≈ uniform (π0 ≈ {pi0:.2f})", ha="right", va="bottom", color=dviz.BLUE,
            fontsize=mpl.rcParams["font.size"] * 0.8, bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    h0 = np.histogram(p, bins=40, range=(0, 1))[0][0]
    ax.annotate("true signal", xy=(0.03, h0 * 0.85), xytext=(0.12, h0 * 0.85), va="center",
                arrowprops=dict(arrowstyle="->", color=dviz.INK, lw=1), fontsize=mpl.rcParams["font.size"] * 0.8)
    ax.set_xlabel("Raw p-value"); ax.set_ylabel("Genes"); ax.set_title("p-value histogram (check before volcano)")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.17, top=0.88)
    return S.save_exact(fig, t)


def f_heatmap(t, labels=True, dendro=True, title=True, on_slide=12.0):
    fig, _ = fig_for(t, on_slide=on_slide)
    fig.axes[0].remove()
    heatmap(fig, (0.03 if dendro else 0.02, 0.18 if labels else 0.05, 0.96, 0.74 if labels else 0.86),
            labels=labels, dendro=dendro, title=title)
    return S.save_exact(fig, t)


def f40():
    t = "image-40-1.png"
    fig, axs = fig_for(t, 1, 2, gridspec_kw={"width_ratios": [1.5, 1]})
    vr = pca_scatter(axs[0], s=90, legend=False, title="12 samples · 1,000 most variable genes")
    axs[0].legend(fontsize=mpl.rcParams["font.size"] * 0.75, loc="center left", bbox_to_anchor=(1.0, 0.5), handletextpad=0.3)
    axs[1].bar(range(1, 7), vr[:6], color=[dviz.BLUE, dviz.BLUE] + [dviz.GREY] * 4, width=0.65)
    axs[1].set_xticks(range(1, 7)); axs[1].set_xlabel("Principal component"); axs[1].set_ylabel("Variance explained (%)")
    axs[1].set_title("Scree plot")
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.21, top=0.88, wspace=0.75)
    return S.save_exact(fig, t)


def f_pca(t, on_slide=12.0):
    fig, ax = fig_for(t, on_slide=on_slide)
    pca_scatter(ax, s=95 if on_slide <= 12 else 80, title="1,000 most variable genes")
    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.3, top=0.9)
    return S.save_exact(fig, t)


def f_distribution(t, from_zero=True, on_slide=12.0):
    fig, ax = fig_for(t, on_slide=on_slide)
    vals = groups("large_n")
    box_points(ax, vals, violin=True, s=14, offset=0.3)
    p = stats.mannwhitneyu(vals[1], vals[0], alternative="two-sided").pvalue
    ax.set_title(f"Mann–Whitney p = {p:.1g}")
    ax.set_xticks([0, 1], [f"Control\nn = {len(vals[0])}", f"Case\nn = {len(vals[1])}"])
    ax.set_ylabel("GENE_X expression (TPM)")
    if from_zero:
        ax.set_ylim(0, 130)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.9)
    return S.save_exact(fig, t)


def f75():
    t = "image-75-1.png"
    fig, ax = fig_for(t, on_slide=11)
    bad_volcano(ax, s=8)
    ax.set_title("Volcano")
    fig.subplots_adjust(left=0.13, right=0.97, bottom=0.17, top=0.9)
    return S.save_exact(fig, t)


if __name__ == "__main__":
    out = [f14(), f15(), f22(), f27(), f28(), f31(), f32a(), f32b(), f33(),
           f_volcano("image-34-1.png"), f35(), f_heatmap("image-38-1.png"), f40(),
           f_distribution("image-53-1.png"), f_volcano("image-54-1.png", labels=8), f_heatmap("image-55-1.png", on_slide=11),
           f_pca("image-56-1.png"),
           f_distribution("image-67-1.png", on_slide=10.5), f_volcano("image-68-1.png", labels=5, on_slide=10),
           f_heatmap("image-69-1.png", labels=False, dendro=True, title=False, on_slide=9), f_pca("image-70-1.png", on_slide=10),
           f75()]
    for p in out:
        print(p.name)
