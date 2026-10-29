"""Reference code for the eight demonstrations (lecture Part 4) and prompt cards 1-8 (Part 5).

Open this after you have tried the prompt card yourself. Run from the package folder:
    python solutions/figures.py          # all eight
    python solutions/figures.py 2 5      # only demos 2 and 5
Each demo reads data/, saves results/figures/demo<k>_<name>.png and .svg, and PRINTS the numbers behind the figure,
so you can compare them with the CHECK cells of the notebooks (lecture slide "Verify the numbers").
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import dviz  # noqa: E402

DATA = ROOT / "data"
MM = 1 / 25.4                                     # inches per millimetre (journal widths: 89 mm, 183 mm)


def fmt_p(p):
    """0.029 -> '0.029'; 8e-05 -> '8.0 × 10⁻⁵' (scientific notation that a reader can parse)."""
    if p >= 0.001:
        return f"{p:.2g}"
    m, e = f"{p:.1e}".split("e")
    sup = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return f"{m} × 10{str(int(e)).translate(sup)}"


def save(fig, name):
    fig.savefig(ROOT / "results" / "figures" / f"{name}.png", dpi=300)
    fig.savefig(ROOT / "results" / "figures" / f"{name}.svg")      # vector, editable text
    plt.close(fig)


# ---------------------------------------------------------------------------------------- 1 · distribution
def demo1_distribution():
    """Card 1: box + every point, n under each label, Mann-Whitney U named, axis from 0, no bar of means."""
    from scipy.stats import mannwhitneyu
    e = pd.read_csv(DATA / "expression_by_group.tsv", sep="\t")
    fig, axs = plt.subplots(1, 2, figsize=(120 * MM, 55 * MM), sharey=True)
    for ax, ds in zip(axs, ["small_n", "large_n"]):
        s = e[e.dataset == ds]
        groups = []
        for i, g in enumerate(["Control", "Case"]):
            v = s[s.group == g].expression
            groups.append(v)
            if len(v) >= 10:                                      # a violin needs enough points to show a shape
                vp = ax.violinplot(v, positions=[i], showextrema=False)
                vp["bodies"][0].set_facecolor(dviz.PALEBLUE)
            ax.boxplot(v, positions=[i], widths=0.12, showfliers=False)
            ax.scatter(i + 0.3 + dviz.jitter(len(v)), v, s=14, color=[dviz.SLATE, dviz.OI_BLUE][i])
        ctrl, case = groups
        n0, n1 = len(ctrl), len(case)
        ax.set_xticks([0, 1], [f"Control n={n0}", f"Case n={n1}"])
        p = mannwhitneyu(ctrl, case).pvalue                       # two-sided
        ax.set_title(f"n = {n0} per group\nMann–Whitney p = {fmt_p(p)}", fontsize=8)
        print(f"demo 1 · {ds}: median control {ctrl.median():.1f}, case {case.median():.1f}; Mann–Whitney p = {p:.2g}")
    axs[0].set_ylim(0, None)
    axs[0].set_ylabel("GENE_X expression (TPM)")
    fig.tight_layout()
    save(fig, "demo1_distribution")


# ---------------------------------------------------------------------------------------- 2 · volcano + MA
def demo2_volcano():
    """Card 2: BH-adjusted p, shrunken log2FC, both thresholds, untested genes removed and counted."""
    from adjustText import adjust_text
    de = pd.read_csv(DATA / "de_results.tsv", sep="\t")
    ok = de.padj.notna()                                          # drop untested genes
    up = ok & (de.padj < .05) & (de.log2FC_shrunk > 1)
    dn = ok & (de.padj < .05) & (de.log2FC_shrunk < -1)
    x, y = de.log2FC_shrunk, -np.log10(de.padj)
    fig, (ax, ma) = plt.subplots(1, 2, figsize=(183 * MM, 70 * MM))
    ax.scatter(x[ok & ~up & ~dn], y[ok & ~up & ~dn], s=3, c="lightgrey", rasterized=True)   # grey first: draw order
    ax.scatter(x[up], y[up], s=6, c=dviz.VERMILION, label=f"{up.sum()} up")
    ax.scatter(x[dn], y[dn], s=6, c=dviz.OI_BLUE, label=f"{dn.sum()} down")
    ax.axhline(-np.log10(.05), ls="--", lw=0.8, color=dviz.SLATE)
    for v in (-1, 1):
        ax.axvline(v, ls="--", lw=0.8, color=dviz.SLATE)
    lim = np.ceil(x[ok].abs().max())
    ax.set_xlim(-lim, lim)                                        # symmetric x-axis
    top = de[ok].nsmallest(8, "padj")
    adjust_text([ax.text(r.log2FC_shrunk, -np.log10(r.padj), r.gene, fontsize=6) for r in top.itertuples()],
                ax=ax, arrowprops=dict(arrowstyle="-", lw=0.5, color=dviz.SLATE))
    ax.set(xlabel="shrunken log2 fold change", ylabel="−log10 adjusted p (BH)")
    ax.legend(title=f"{(~ok).sum()} genes not tested", fontsize=6, title_fontsize=6, loc="upper left")
    A = np.log2(de.baseMean + 1)
    ma.scatter(A[ok], de.log2FC[ok], s=2, c="lightgrey", label="raw", rasterized=True)
    ma.scatter(A[ok], de.log2FC_shrunk[ok], s=0.6, c=dviz.BLUE, label="shrunken", rasterized=True)
    ma.axhline(0, color=dviz.INK, lw=0.6)
    ma.set(xlabel="log2(baseMean + 1)", ylabel="log2 fold change")
    ma.legend(fontsize=6, markerscale=3)
    fig.tight_layout()
    save(fig, "demo2_volcano_ma")
    assert not (up & dn).any()                                    # 177 up, 106 down
    print(f"demo 2 · tested {ok.sum():,}, untested {(~ok).sum():,}, up {up.sum()}, down {dn.sum()}")


# ---------------------------------------------------------------------------------------- 3 · heatmap
def demo3_heatmap():
    """Card 3: top 40 DE genes, z-scored rows, 1 - Pearson r, average linkage, diverging map at 0, annotation bars."""
    from scipy.cluster.hierarchy import linkage, leaves_list
    from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap, ListedColormap
    from matplotlib.patches import Patch
    ln = pd.read_csv(DATA / "rnaseq_log2norm.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(DATA / "rnaseq_metadata.tsv", sep="\t").set_index("sample")
    de = pd.read_csv(DATA / "de_results.tsv", sep="\t")
    top = de.dropna().nsmallest(40, "padj").gene
    Z = ln.loc[top]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0)
    rows = linkage(Z, "average", metric="correlation")
    cols = linkage(Z.T, "average", metric="correlation")
    Z = Z.iloc[leaves_list(rows), leaves_list(cols)]
    blue_white_orange = LinearSegmentedColormap.from_list("bwo", [dviz.DARKBLUE, "white", dviz.VERMILION])
    fig = plt.figure(figsize=(120 * MM, 120 * MM))
    ann = fig.add_axes([0.20, 0.86, 0.62, 0.05])
    ax = fig.add_axes([0.20, 0.14, 0.62, 0.70])
    cax = fig.add_axes([0.85, 0.14, 0.03, 0.70])
    m = meta.loc[Z.columns]
    codes = [[0 if c == "control" else 1 for c in m.condition], [2 if b == "A" else 3 for b in m.batch]]
    ann_colours = [dviz.SLATE, dviz.BLUE, dviz.ORANGE, dviz.GREY]   # control, treated, batch A, batch B
    ann.imshow(codes, cmap=ListedColormap(ann_colours), vmin=0, vmax=3, aspect="auto")
    fig.legend(handles=[Patch(color=c, label=l) for c, l in zip(ann_colours, ["control", "treated", "batch A", "batch B"])],
               loc="lower center", ncol=4, fontsize=6, frameon=False)
    ann.set_yticks([0, 1], ["condition", "batch"], fontsize=6)
    ann.set_xticks([])
    im = ax.imshow(Z, aspect="auto", cmap=blue_white_orange, norm=TwoSlopeNorm(0, -2.5, 2.5))
    ax.set_xticks(range(Z.shape[1]), Z.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(Z.shape[0]), Z.index, fontsize=5)
    fig.colorbar(im, cax=cax, label="row z-score")
    save(fig, "demo3_heatmap")
    assert np.allclose(Z.mean(axis=1), 0, atol=1e-9) and np.allclose(Z.std(axis=1), 1)
    print("demo 3 · column order:", ", ".join(Z.columns))


# ---------------------------------------------------------------------------------------- 4 · PCA (+ UMAP)
def demo4_pca():
    """Card 4A: PCA on the 1,000 most variable genes; colour = condition, shape = batch; % variance; scree.
    Then card 4B (demo4b_umap) if the B4 notebook has been run."""
    from sklearn.decomposition import PCA
    ln = pd.read_csv(DATA / "rnaseq_log2norm.tsv.gz", sep="\t", index_col=0)
    meta = pd.read_csv(DATA / "rnaseq_metadata.tsv", sep="\t").set_index("sample").loc[ln.columns]
    X = ln.loc[ln.var(axis=1).nlargest(1000).index].T
    pca = PCA(5).fit(X - X.mean())
    S = pca.transform(X - X.mean())
    vr = pca.explained_variance_ratio_ * 100
    fig, (ax, scree) = plt.subplots(1, 2, figsize=(150 * MM, 60 * MM), gridspec_kw={"width_ratios": [1.6, 1]})
    for cond, col in [("control", dviz.SLATE), ("treated", dviz.BLUE)]:
        for b, mk in [("A", "o"), ("B", "^")]:
            m = ((meta.condition == cond) & (meta.batch == b)).to_numpy()
            ax.scatter(S[m, 0], S[m, 1], color=col, marker=mk, s=30, label=f"{cond}, batch {b}")
    ax.set_xlabel(f"PC1 ({vr[0]:.0f}%)")
    ax.set_ylabel(f"PC2 ({vr[1]:.0f}%)")
    ax.legend(fontsize=6)
    scree.bar(range(1, 6), vr, color=dviz.GREY)
    scree.set(xlabel="principal component", ylabel="% variance")
    fig.tight_layout()
    save(fig, "demo4_pca")
    print("demo 4 · % variance PC1-PC5:", ", ".join(f"{v:.1f}" for v in vr))
    demo4b_umap()


def demo4b_umap():
    """Card 4B: UMAP of the B4 cells with the cluster names saved by the B4 notebook; labels on the map, no axis
    numbers, settings stated."""
    from adjustText import adjust_text
    lab_file = ROOT / "results" / "B4_clusters.tsv"
    if not lab_file.exists():
        print("demo 4B · skipped: run electives/B4_single_cell.ipynb first (it writes results/B4_clusters.tsv)")
        return
    cells = pd.read_csv(DATA / "scrna_cells.tsv", sep="\t").merge(
        pd.read_csv(lab_file, sep="\t"), on="cell_id", validate="one_to_one")     # join on the barcode, not row order
    colours = {"Naive CD4 T": "#0072B2", "Memory CD4 T": "#E69F00", "CD8 T": "#009E73", "NK": "#D55E00",
               "B": "#56B4E9", "CD14 Mono": "#CC79A7", "FCGR3A Mono": "#8C6D1F", "DC": "#000000"}  # as in B4
    import matplotlib.patheffects as pe
    fig, ax = plt.subplots(figsize=(89 * MM, 80 * MM))
    labels, halo = [], [pe.withStroke(linewidth=2.5, foreground="white")]
    for name, d in cells.groupby("cluster_name"):
        ax.scatter(d.umap1, d.umap2, s=1.5, color=colours.get(name, dviz.SLATE), rasterized=True)
        labels.append(ax.text(d.umap1.median(), d.umap2.median(), name, fontsize=7, color=dviz.INK, path_effects=halo))
    pts = cells.sample(600, random_state=0)                             # labels also move off the points
    adjust_text(labels, x=pts.umap1.to_numpy(), y=pts.umap2.to_numpy(), ax=ax, expand=(1.3, 1.6),
                arrowprops=dict(arrowstyle="-", color=dviz.SLATE, lw=0.5))
    ax.set(xticks=[], yticks=[], xlabel="UMAP 1", ylabel="UMAP 2")       # UMAP coordinates have no units
    ax.set_title(f"{len(cells):,} cells · {cells.cluster_name.nunique()} clusters named from marker genes", loc="left")
    fig.text(0.01, 0.01, "UMAP of 30 PCs (1,000 HVGs) · n_neighbors 15 · min_dist 0.3 · seed 0\n"
                         "distances between clusters are not interpretable", fontsize=6, color=dviz.SLATE)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, "demo4b_umap")
    print("demo 4B · cells per cluster:", ", ".join(f"{k} {n:,}" for k, n in cells.cluster_name.value_counts().items()))


# ---------------------------------------------------------------------------------------- 5 · Manhattan + QQ
def demo5_manhattan_qq():
    """Card 5: cumulative position, alternating blues, 5e-8 line, rasterized points; QQ with lambda for both p columns."""
    from scipy.stats import chi2
    gwas = pd.read_csv(DATA / "gwas_sumstats.tsv.gz", sep="\t")
    g = gwas.sort_values(["CHR", "BP"]).reset_index(drop=True)
    off, centres = 0, {}
    for c in range(1, 23):                                        # cumulative position
        m = g.CHR == c
        g.loc[m, "x"] = g.BP[m] + off
        centres[c] = off + g.BP[m].max() / 2
        off += g.BP[m].max() + 2e7
    fig, (ax, qq) = plt.subplots(1, 2, figsize=(183 * MM, 60 * MM), gridspec_kw={"width_ratios": [3, 1]})
    ax.scatter(g.x, -np.log10(g.P), c=np.where(g.CHR % 2, dviz.BLUE, dviz.SKY), s=1, rasterized=True)
    ax.axhline(-np.log10(5e-8), ls="--", lw=0.8, color=dviz.VERMILION)
    ax.set_xticks(list(centres.values()), [str(c) if c < 13 or c % 2 else "" for c in centres], fontsize=6)
    ax.set(xlabel="chromosome (hg38)", ylabel="−log10 p")
    for col, colour, lab in [("P_no_ancestry_PCs", dviz.ORANGE, "no PCs"), ("P", dviz.BLUE, "4 PCs")]:
        lam = np.median(chi2.isf(g[col], 1)) / chi2.ppf(0.5, 1)  # 1.33 without PCs, 1.03 with them
        e, o = dviz.qq_points(g[col])
        qq.scatter(e, o, s=2, color=colour, label=f"{lab} · λ = {lam:.2f}", rasterized=True)
        print(f"demo 5 · λ ({col}) = {lam:.3f}")
    qq.plot([0, 6], [0, 6], color=dviz.INK, lw=0.8)
    qq.set(xlabel="expected −log10 p", ylabel="observed −log10 p")
    qq.legend(fontsize=6)
    fig.tight_layout()
    save(fig, "demo5_manhattan_qq")
    lam_p = np.median(chi2.isf(g.P, 1)) / chi2.ppf(0.5, 1)
    assert 0.95 < lam_p < 1.1, "wrong p-value column?"
    print(f"demo 5 · {len(g):,} SNPs; {(g.P < 5e-8).sum()} with P < 5e-8")


# ---------------------------------------------------------------------------------------- 6 · copy-number tracks
def demo6_cnv_tracks():
    """Card 6: depth and log2 ratio on one x-axis (Mb), segments, expected -1, genes; BED-style 0-based starts."""
    lc = pd.read_csv(DATA / "cnv_chr22_10kb.tsv", sep="\t")
    genes = pd.read_csv(DATA / "cnv_genes_chr22.tsv", sep="\t")
    ok = lc.mappable
    lc["log2_ratio"] = np.log2((lc.depth_patient / lc.depth_patient[ok].median())
                               / (lc.depth_reference / lc.depth_reference[ok].median())).where(ok)
    seg = dviz.segment(lc.log2_ratio)                             # one row per segment (first, last bin index)
    seg["start"], seg["end"] = lc.start.to_numpy()[seg["first"]], lc.end.to_numpy()[seg["last"]]
    win = (lc.start >= 17e6) & (lc.end <= 23e6)
    x = lc.start[win] / 1e6
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=(183 * MM, 90 * MM), gridspec_kw={"height_ratios": [1, 1.4, 0.5]})
    axs[0].fill_between(x, lc.depth_patient[win], step="mid", color=dviz.SKY)    # read depth
    axs[1].scatter(x, lc.log2_ratio[win], s=2, color=dviz.SLATE)                # bins
    s = seg[(seg.end > 17e6) & (seg.start < 23e6)]
    axs[1].hlines(s["mean"], s.start.clip(lower=17e6) / 1e6, s.end.clip(upper=23e6) / 1e6, lw=2, color=dviz.BLUE)
    axs[1].axhline(0, color=dviz.LIGHT, lw=0.8)
    axs[1].axhline(-1, color=dviz.SLATE, ls="--", lw=0.8)
    axs[1].text(22.9, -1.1, "1 copy expected: log2(1/2) = −1", ha="right", va="top", fontsize=6)
    g = genes[(genes.end > 17e6) & (genes.start < 23e6)]
    axs[2].hlines(np.zeros(len(g)), g.start / 1e6, g.end / 1e6, lw=6, color=dviz.SLATE)
    for k, name in enumerate(["HIRA", "TBX1", "COMT", "CRKL"]):
        r = g[g.gene == name].iloc[0]
        axs[2].text((r.start + r.end) / 2e6, 0.6 if k % 2 else -0.9, name, ha="center", fontsize=6)
    axs[2].set(ylim=(-1.2, 1.2), yticks=[], xlim=(17, 23))
    dl = s[s["mean"] < -0.5].iloc[0]
    for ax in axs:
        ax.axvspan(dl.start / 1e6, dl.end / 1e6, alpha=.08, color=dviz.VERMILION)
    axs[0].set_ylabel("depth (x)")
    axs[1].set_ylabel("log2 ratio")
    axs[2].set_xlabel("chr22 position (Mb, hg38)")
    fig.tight_layout()
    save(fig, "demo6_cnv_tracks")
    inside = lc.log2_ratio[(lc.start >= dl.start) & (lc.end <= dl.end)]
    print(f"demo 6 · deletion chr22:{dl.start + 1:,.0f}-{dl.end:,.0f} (1-based), {(dl.end - dl.start) / 1e6:.2f} Mb; "
          f"median log2 inside {inside.median():.2f}, outside {lc.log2_ratio[win & ~lc.start.between(dl.start, dl.end)].median():.2f}")


# ---------------------------------------------------------------------------------------- 7 · UpSet + lollipop
def demo7_upset_lollipop():
    """Card 7: UpSet sorted by size (no four-set Venn), precision vs the truth; TP53 lollipop with domains."""
    vc = pd.read_csv(DATA / "variant_callsets.tsv", sep="\t")
    callers = ["GATK", "DeepVariant", "FreeBayes", "Strelka2"]
    vc[callers] = vc[callers].astype(bool)                         # stored as 0/1
    fig, counts = dviz.upset(vc, callers, max_bars=10, title="Four callers: the 10 largest intersections")
    save(fig, "demo7_upset")
    for c in callers:                                              # check: precision of each caller vs the truth set
        print(f"demo 7 · {c}: precision {vc.loc[vc[c], 'in_truth_set'].mean():.3f}")
    var = pd.read_csv(DATA / "tp53_variants.tsv", sep="\t")
    dom = pd.read_csv(DATA / "tp53_domains.tsv", sep="\t")
    per = var.groupby(["position", "variant_class"]).size().rename("n").reset_index()
    colours = {"missense": dviz.OI_BLUE, "truncating": dviz.VERMILION, "splice/other": dviz.GREEN}
    fig, ax = plt.subplots(figsize=(183 * MM, 60 * MM))
    for cls, d in per.groupby("variant_class"):
        ax.vlines(d.position, 0, d.n, color=dviz.GREY, lw=0.6)
        ax.scatter(d.position, d.n, s=8 + d.n * 1.5, color=colours.get(cls, dviz.SLATE), label=cls, zorder=3)
    for r in dom.itertuples():                                     # domains as boxes (1-based, inclusive)
        ax.add_patch(plt.Rectangle((r.aa_start, -8), r.aa_end - r.aa_start + 1, 5, color=dviz.SLATE))
        ax.text((r.aa_start + r.aa_end) / 2, -5.5, r.short, ha="center", va="center", color="white", fontsize=6)
    hot = var.groupby("position").agg(n=("tumour_id", "size"), label=("protein_change", lambda s: s.mode()[0]))
    for pos, r in hot[hot.n >= 30].iterrows():
        ax.text(pos, r.n + 3, r.label, ha="center", fontsize=6)
    ax.set(xlim=(0, 394), ylim=(-10, None), xlabel="amino-acid position (TP53, 393 aa)", ylabel="tumours")
    ax.legend(fontsize=6)
    fig.tight_layout()
    save(fig, "demo7_lollipop")


# ---------------------------------------------------------------------------------------- 8 · ROC, PR, calibration
def demo8_classifier():
    """Card 8: test split only; ROC with AUROC, PR with AP and the prevalence baseline, calibration (10 bins, >= 20)."""
    from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score
    from sklearn.isotonic import IsotonicRegression
    cs = pd.read_csv(DATA / "classifier_scores.tsv", sep="\t")
    cal = cs[cs.split == "calibration"]
    iso = IsotonicRegression(out_of_bounds="clip").fit(cal.score_modelB, cal.label_pathogenic)
    t = cs[cs.split == "test"]
    y = t.label_pathogenic
    fig, (ax_roc, ax_pr, ax_cal) = plt.subplots(1, 3, figsize=(183 * MM, 60 * MM))
    for col, colour, lab in [("score_modelA", dviz.SLATE, "A"), ("score_modelB", dviz.BLUE, "B")]:
        fpr, tpr, _ = roc_curve(y, t[col])
        prec, rec, _ = precision_recall_curve(y, t[col])
        auc, ap = roc_auc_score(y, t[col]), average_precision_score(y, t[col])
        ax_roc.plot(fpr, tpr, color=colour, label=f"{lab} AUROC {auc:.2f}")
        ax_pr.plot(rec, prec, color=colour, label=f"{lab} AP {ap:.2f}")
        print(f"demo 8 · model {lab}: AUROC {auc:.3f}, AP {ap:.3f}")
    ax_roc.plot([0, 1], [0, 1], ls="--", color=dviz.GREY)
    ax_pr.axhline(y.mean(), ls="--", color=dviz.GREY)             # baseline = prevalence 3.5%
    ax_pr.text(0.02, y.mean() + 0.02, f"baseline = prevalence {y.mean():.1%}", fontsize=6, color=dviz.SLATE)
    for prob, colour, lab in [(t.score_modelB, dviz.BLUE, "raw"), (iso.predict(t.score_modelB), dviz.GREEN, "calibrated")]:
        ct = dviz.calibration_table(y, prob, n_bins=10, min_count=20)
        ax_cal.plot(ct.mean_predicted, ct.observed, marker="o", color=colour, label=lab)
    ax_cal.plot([0, 1], [0, 1], ls="--", color=dviz.GREY)
    ax_roc.set(xlabel="false-positive rate", ylabel="true-positive rate")
    ax_pr.set(xlabel="recall", ylabel="precision")
    ax_cal.set(xlabel="predicted probability", ylabel="observed fraction")
    for ax in (ax_roc, ax_pr, ax_cal):
        ax.legend(fontsize=6)
    fig.tight_layout()
    save(fig, "demo8_classifier")
    assert average_precision_score(y, t.score_modelB) > y.mean()
    print(f"demo 8 · prevalence {y.mean():.1%} in the test split ({len(t):,} variants)")


DEMOS = {1: demo1_distribution, 2: demo2_volcano, 3: demo3_heatmap, 4: demo4_pca,
         5: demo5_manhattan_qq, 6: demo6_cnv_tracks, 7: demo7_upset_lollipop, 8: demo8_classifier}

if __name__ == "__main__":
    dviz.use_style()
    (ROOT / "results" / "figures").mkdir(parents=True, exist_ok=True)
    for k in [int(a) for a in sys.argv[1:]] or DEMOS:
        DEMOS[k]()
