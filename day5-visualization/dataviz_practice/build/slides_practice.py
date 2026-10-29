"""Practice-slide figures for the two core exercises (B1 two groups, B2 bulk RNA-seq), 1800 x 1210 px at 288 dpi."""
import numpy as np
import pandas as pd
import matplotlib as mpl
from scipy import stats

import slides_core as C           # loads the practice data and the shared plotting helpers
import slidefig as S
import dviz


def practice_b1():
    t = "practice_B1.png"
    fig, axs = C.fig_for(t, 1, 2, on_slide=11.0, sharey=True)
    for ax, ds, title in [(axs[0], "small_n", "Small study"), (axs[1], "large_n", "Larger study")]:
        ctrl, case = C.groups(ds)
        C.box_points(ax, [ctrl, case], s=26 if ds == "small_n" else 14)
        mw = stats.mannwhitneyu(case, ctrl, alternative="two-sided").pvalue
        we = stats.ttest_ind(case, ctrl, equal_var=False).pvalue
        ax.set_title(f"{title}, n = {len(ctrl)} per group", fontsize=mpl.rcParams["font.size"] * 0.95)
        f = mpl.rcParams["font.size"] * 0.78
        txt = f"Mann–Whitney p = {mw:.2g}\nWelch t-test p = {we:.2g}" if ds == "small_n" else \
              f"Mann–Whitney p = {mw:.1e}".replace("e-05", " × 10⁻⁵") + f"\nWelch t-test p = {we:.1e}".replace("e-05", " × 10⁻⁵")
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=f, color=dviz.INK)
        ax.set_xticks([0, 1], [f"Control\nn = {len(ctrl)}", f"Case\nn = {len(case)}"])
        ax.set_ylim(0, 175)
    axs[0].set_ylabel("GENE_X expression (TPM)")
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.17, top=0.91, wspace=0.12)
    return S.save_exact(fig, t)


def practice_b2():
    t = "practice_B2.png"
    fig, axs = C.fig_for(t, 2, 2, on_slide=9.5)
    f = mpl.rcParams["font.size"]
    # (a) PCA: biology or batch?
    C.pca_scatter(axs[0, 0], s=46, legend=False)
    axs[0, 0].set_title("1 · PCA: treatment, then batch", fontsize=f * 0.95)
    axs[0, 0].text(0.5, 0.97, "colour = condition\nshape = batch", transform=axs[0, 0].transAxes, ha="center",
                   va="top", fontsize=f * 0.72, color=dviz.SLATE)
    # (b) p-value histogram
    p = C.de.loc[C.tested & (C.de["baseMean"] >= 10), "pvalue"].to_numpy()
    pi0 = min(1.0, float((p > 0.5).mean() / 0.5))
    axs[0, 1].hist(p, bins=20, range=(0, 1), color=dviz.SKY, edgecolor="white")
    axs[0, 1].axhline(len(p) / 20 * pi0, color=dviz.BLUE, ls="--", lw=1.1)
    axs[0, 1].set_title("2 · p-values: flat floor + spike", fontsize=f * 0.95)
    axs[0, 1].set_xlabel("raw p-value"); axs[0, 1].set_ylabel("genes")
    # (c) volcano with both thresholds
    C.volcano(axs[1, 0], labels=0, legend_counts=True, s=(3, 7))
    axs[1, 0].set_xlabel("shrunken log2 fold change"); axs[1, 0].set_ylabel("−log10 padj")
    axs[1, 0].set_title("3 · Volcano: both thresholds", fontsize=f * 0.95)
    # (d) check against the planted truth
    truth = pd.read_csv(C.D / "truth" / "rnaseq_truth.tsv", sep="\t").set_index("gene")
    de = C.de.set_index("gene")
    is_de = truth.loc[de.index, "is_de"].astype(bool)
    sig = de["padj"].lt(0.05)
    both = sig & de["log2FC_shrunk"].abs().gt(1)
    rows = [("padj < 0.05", sig), ("padj < 0.05\nand |log2FC| > 1", both)]
    ax = axs[1, 1]
    for i, (lab, m) in enumerate(rows):
        tp, fp = int((m & is_de).sum()), int((m & ~is_de).sum())
        ax.barh(i, tp, color=dviz.BLUE, height=0.55)
        ax.barh(i, fp, left=tp, color=dviz.VERMILION, height=0.55)
        ax.text(tp + fp + 8, i, f"{tp + fp} genes\n{100 * fp / (tp + fp):.1f}% false", va="center", fontsize=f * 0.72)
    ax.set_yticks([0, 1], [r[0] for r in rows], fontsize=f * 0.8)
    ax.invert_yaxis()
    ax.set_xlim(0, 560)
    ax.set_xlabel("genes called: planted (blue), false (red)")
    ax.set_title("4 · Check with the truth", fontsize=f * 0.95)
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.135, top=0.95, wspace=0.45, hspace=0.66)
    return S.save_exact(fig, t)


if __name__ == "__main__":
    print(practice_b1())
    print(practice_b2())
