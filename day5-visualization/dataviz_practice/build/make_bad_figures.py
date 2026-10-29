"""Exercise 1: eight misleading figures (M1-M8, M = misleader) made from the practice data, plus the slide panel image-84-1.png.

They deliberately use naive defaults (matplotlib's default style, red/green, rainbow, 3-D-like pie) because these are
the mistakes the lecture discusses. Answer key: solutions/ANSWER_KEY_bad_figures.md.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import sparse

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT))
D, OUT = ROOT / "data", ROOT / "bad_figures"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261027 + 10)


def naive():
    mpl.rcdefaults()                       # the defaults many people never change
    plt.rcParams.update({"figure.dpi": 100, "savefig.dpi": 200})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=200, metadata={"Software": None})
    plt.close(fig)
    return OUT / f"{name}.png"


def b1():  # bar of means +- SE with a star, n = 4 per group
    e = pd.read_csv(D / "expression_by_group.tsv", sep="\t").query("dataset == 'small_n'")
    g = [e.loc[e["group"] == k, "expression"] for k in ["Control", "Case"]]
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.bar(["Control", "Case"], [v.mean() for v in g], yerr=[v.std(ddof=1) / 2 for v in g], color=["grey", "red"], capsize=4)
    ax.text(1, g[1].mean() + g[1].std(ddof=1) / 2 + 3, "*", ha="center", fontsize=16)
    ax.set_ylim(0, 135)
    ax.set_ylabel("GENE_X")
    ax.set_title("GENE_X is higher in cases")
    return save(fig, "M1")


def b2():  # truncated bar axis
    m = pd.read_csv(D / "mapping_rates.tsv", sep="\t").groupby("pipeline")["pct_mapped"].mean()
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.bar(m.index, m.values, color=["grey", "grey", "green"], edgecolor="black")
    ax.set_ylim(95.8, 98.1)
    ax.set_title("Pipeline C maps far more reads!")
    return save(fig, "M2")


def b3():  # rainbow heatmap, raw values, no colour bar, no labels that mean anything
    meta = pd.read_csv(D / "rnaseq_metadata.tsv", sep="\t")
    c = pd.read_csv(D / "rnaseq_counts.tsv.gz", sep="\t", index_col="gene").drop(columns=["gene_id", "chrom"])[meta["sample"]]
    de = pd.read_csv(D / "de_results.tsv", sep="\t")
    top = de[de["padj"].notna()].nsmallest(40, "padj")["gene"]
    X = np.log2(c.loc[top] + 1).to_numpy()
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.imshow(X, aspect="auto", cmap="jet")
    ax.set_title("Top genes")
    return save(fig, "M3")


def scrna():
    cells = pd.read_csv(D / "scrna_cells.tsv", sep="\t")
    truth = pd.read_csv(D / "truth" / "scrna_truth.tsv", sep="\t").set_index("cell_id")
    cells["cell_type"] = truth.loc[cells["cell_id"], "cell_type"].to_numpy()
    return cells


def b4():  # draw order hides the expressing cells
    cells = scrna()
    genes = pd.read_csv(D / "scrna_genes.tsv", sep="\t")
    X = sparse.load_npz(D / "scrna_counts.npz").tocsc()
    j = int(np.flatnonzero(genes["gene"].to_numpy() == "FCGR3A")[0])
    v = np.log1p(X[:, j].toarray().ravel() / np.asarray(X.sum(axis=1)).ravel() * 1e4)
    order = np.argsort(-v)                 # high values drawn FIRST, zeros on top
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.scatter(cells["umap1"].to_numpy()[order], cells["umap2"].to_numpy()[order], c=v[order], cmap="Reds", s=12)
    ax.set_title("FCGR3A is rarely expressed")
    return save(fig, "M4")


def b5():  # pie with many slices, shadow, explode, legend
    comp = pd.read_csv(D / "celltype_composition.tsv", sep="\t").groupby("cell_type")["n_cells"].sum()
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.pie(comp.values, explode=[0.08] * len(comp), shadow=True, startangle=20)
    ax.legend(comp.index, fontsize=5, loc="center left", bbox_to_anchor=(0.92, 0.5))
    ax.set_title("Cell types")
    return save(fig, "M5")


def b6():  # red-green volcano, raw p, no thresholds, untested genes drawn at p = 1
    de = pd.read_csv(D / "de_results.tsv", sep="\t")
    x = de["log2FC"].clip(-6, 6)
    y = -np.log10(de["pvalue"].fillna(1.0))
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.scatter(x[x < 0], y[x < 0], s=3, color="green")
    ax.scatter(x[x >= 0], y[x >= 0], s=3, color="red")
    ax.set_xlabel("FC"); ax.set_ylabel("p")
    return save(fig, "M6")


def b7():  # a distance claim on a UMAP
    cells = scrna()
    fig, ax = plt.subplots(figsize=(4, 3.4))
    cen = {}
    for t, d in cells.groupby("cell_type"):
        ax.scatter(d["umap1"], d["umap2"], s=4, label=t)
        cen[t] = (d["umap1"].median(), d["umap2"].median())
    others = np.array([v for k, v in cen.items() if k != "B"])
    far = others[np.argmax(((others - np.array(cen["B"])) ** 2).sum(axis=1))]
    ax.plot([cen["B"][0], far[0]], [cen["B"][1], far[1]], color="black", lw=1.5)
    ax.text(0.03, 0.04, "B cells are the most different\n(farthest from everything)", transform=ax.transAxes, fontsize=7)
    ax.legend(fontsize=5, markerscale=2, loc="upper right")
    return save(fig, "M7")


def b8():  # dual y-axes make two series look related
    tc = pd.read_csv(D / "timecourse.tsv", sep="\t").groupby(["gene", "time_h"])["log2fc_vs_0h"].mean().unstack(0)
    fig, ax = plt.subplots(figsize=(4, 3.4))
    ax.plot(tc.index, tc["TSC22D3"], "o-", color="blue")
    ax2 = ax.twinx()
    ax2.plot(tc.index, tc["FKBP5"], "s-", color="red")
    ax2.set_ylim(-0.4, 3.9)
    ax.set_title("TSC22D3 drives FKBP5")
    return save(fig, "M8")


def panel(paths):
    """Slide image for Exercise 1 (2 x 4 grid), exact size of the deck frame image-84-1.png."""
    import slidefig as S
    fig = plt.figure(figsize=(3328 / 288, 1180 / 288), dpi=288)
    for k, p in enumerate(paths):
        r, c = divmod(k, 4)
        ax = fig.add_axes([0.005 + c * 0.249, 0.51 - r * 0.5, 0.243, 0.47])
        if p is not None and Path(p).exists():
            ax.imshow(plt.imread(p))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color("#E4ECFB")
        ax.text(0.012, 0.975, f"M{k + 1}", transform=ax.transAxes, ha="left", va="top", fontsize=9,
                color="#1F5FEA", fontweight="bold", fontfamily="DejaVu Sans")
    S.OUT.mkdir(parents=True, exist_ok=True)
    out = S.OUT / "image-84-1.png"
    with plt.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(out, dpi=288, metadata={"Software": None})
    plt.close(fig)
    return out


if __name__ == "__main__":
    naive()
    have_sc = (D / "scrna_cells.tsv").exists() and (D / "scrna_counts.npz").exists()
    paths = [b1(), b2(), b3(), b4() if have_sc else None, b5(), b6(), b7() if have_sc else None, b8()]
    print("made:", [p.name for p in paths if p], "· single-cell figures:", "yes" if have_sc else "waiting for single-cell data (M4, M7)")
    print(panel(paths))
