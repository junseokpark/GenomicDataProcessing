"""Single-cell helpers for Practice B4 (one source of truth).

The notebook builder copies the source of these functions into visible notebook cells
(inspect.getsource), the data generator uses them for the fallback UMAP, and the slide script
uses them for the lecture figures. Plain numpy / pandas / scipy / scikit-learn only.
"""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

# Colours for the 8 cell types: Okabe-Ito categories (colour-blind safe). Okabe-Ito yellow is
# replaced by a dark gold because yellow points are almost invisible on white.
CELL_TYPES = ["Naive CD4 T", "Memory CD4 T", "CD8 T", "NK", "B", "CD14 Mono", "FCGR3A Mono", "DC"]
CT_COLORS = {"Naive CD4 T": "#0072B2", "Memory CD4 T": "#E69F00", "CD8 T": "#009E73", "NK": "#D55E00",
             "B": "#56B4E9", "CD14 Mono": "#CC79A7", "FCGR3A Mono": "#8C6D1F", "DC": "#000000"}


def normalise_log1p(counts, target=10_000):
    """Scale every cell to `target` UMIs, then take log(1 + x).

    counts: cells x genes (scipy sparse or numpy). Returns (lognorm, scaled) as dense arrays.
    """
    X = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    X = X.astype(np.float64)
    umis = X.sum(axis=1, keepdims=True)          # UMIs per cell (library size)
    scaled = X / umis * target                   # every row now sums to `target`
    return np.log1p(scaled), scaled


def highly_variable_genes(scaled, n_top=1000, n_bins=20):
    """Indices of the n_top genes that vary most compared with genes of similar mean.

    Dispersion = variance / mean of the normalised counts. Genes are put into n_bins bins of
    equal size by mean, and dispersions are z-scored inside each bin (Seurat / Scanpy idea).
    """
    mean = scaled.mean(axis=0)
    var = scaled.var(axis=0, ddof=1)
    log_disp = np.log(var / mean)
    bins = pd.qcut(np.log1p(mean), n_bins, labels=False, duplicates="drop")
    z = pd.Series(log_disp).groupby(bins).transform(lambda d: (d - d.mean()) / d.std())
    order = np.argsort(-z.to_numpy(), kind="stable")
    return np.sort(order[:n_top]), z.to_numpy()


def run_pca(lognorm_hvg, n_pcs=30, scale=False, seed=0):
    """PCA of the log-normalised HVG values (PCA subtracts each gene's mean itself).

    scale=True first z-scores every gene (clipped at 10), as Seurat's ScaleData does.
    Returns (cell scores: cells x n_pcs, fraction of variance explained by each PC).
    """
    from sklearn.decomposition import PCA
    X = lognorm_hvg
    if scale:
        X = np.clip((X - X.mean(axis=0)) / X.std(axis=0, ddof=1), -10, 10)
    pca = PCA(n_components=n_pcs, svd_solver="full", random_state=seed)
    scores = pca.fit_transform(X)
    signs = np.sign(scores[np.abs(scores).argmax(axis=0), np.arange(n_pcs)])
    return scores * signs, pca.explained_variance_ratio_      # fixed signs: same picture every run


def run_umap(pcs, n_neighbors=15, min_dist=0.3, seed=0):
    """UMAP of the PCs with fixed settings. Returns (embedding, seconds)."""
    import umap
    t0 = time.time()
    emb = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed,
                    n_jobs=1).fit_transform(pcs)
    return emb, time.time() - t0


def run_tsne(pcs, perplexity=30, seed=0):
    """t-SNE of the PCs (scikit-learn, Barnes-Hut, PCA initialisation). Returns (embedding, seconds)."""
    from sklearn.manifold import TSNE
    t0 = time.time()
    emb = TSNE(perplexity=perplexity, init="pca", learning_rate="auto",
               random_state=seed).fit_transform(pcs)
    return emb, time.time() - t0


def knn_overlap(high, low, k=15):
    """Average fraction of each cell's k nearest neighbours in `high` (e.g. 30 PCs)
    that are also among its k nearest neighbours in `low` (a 2-D embedding). 1 = all kept."""
    from sklearn.neighbors import NearestNeighbors
    nn_high = NearestNeighbors(n_neighbors=k).fit(high).kneighbors(return_distance=False)
    nn_low = NearestNeighbors(n_neighbors=k).fit(low).kneighbors(return_distance=False)
    shared = [np.intersect1d(a, b).size for a, b in zip(nn_high, nn_low)]
    return float(np.mean(shared) / k)


def centroid_distances(emb, labels):
    """Distance between the centres (means) of every pair of clusters in one embedding."""
    from scipy.spatial.distance import pdist
    labels = np.asarray(labels)
    names = sorted(set(labels))
    centres = np.array([emb[labels == c].mean(axis=0) for c in names])
    pairs = [f"{a} | {b}" for i, a in enumerate(names) for b in names[i + 1:]]
    return pd.Series(pdist(centres), index=pairs)


def dotplot_table(counts_genes, lognorm_genes, clusters, gene_names):
    """Per cluster and gene: % of cells with any UMI, and mean log expression (all cells)."""
    clusters = pd.Series(np.asarray(clusters), name="cluster")
    expressing = pd.DataFrame(np.asarray(counts_genes) > 0, columns=gene_names).astype(float)
    level = pd.DataFrame(np.asarray(lognorm_genes), columns=gene_names)
    pct = 100 * expressing.groupby(clusters).mean()
    mean = level.groupby(clusters).mean()
    return pct, mean


def annotate_clusters(mean_table, markers):
    """Give each cluster the cell type whose marker genes are highest in it.

    mean_table: clusters x genes (mean log expression); markers: {cell type: [genes]}.
    Each gene is z-scored across clusters, so every marker counts equally.
    """
    z = (mean_table - mean_table.mean()) / mean_table.std(ddof=0)
    score = pd.DataFrame({t: z[g].mean(axis=1) for t, g in markers.items()})
    return score.idxmax(axis=1), score


def label_clusters(ax, emb, labels, fontsize=None, color="#17202A", halo=3):
    """Write each cluster's name next to its cells (direct labels, no legend to look up).

    Each name starts at the cluster's median position; adjustText (if installed) then pushes it off
    the points and draws a thin line back to the cluster.
    """
    import matplotlib.patheffects as pe
    labels = np.asarray(labels)
    texts = []
    for name in pd.unique(labels):
        x, y = np.median(emb[labels == name], axis=0)
        texts.append(ax.text(x, y, name, ha="center", va="center", fontsize=fontsize, color=color,
                             path_effects=[pe.withStroke(linewidth=halo, foreground="white")]))
    ax.margins(0.1)                                    # room for labels at the edges
    try:
        from adjustText import adjust_text
        step = max(1, len(emb) // 1500)                # avoid (a sample of) the points as well
        adjust_text(texts, x=emb[::step, 0], y=emb[::step, 1], ax=ax, expand=(1.1, 1.3), iter_lim=25,
                    arrowprops=dict(arrowstyle="-", color="#5B6573", lw=0.6))   # fixed iterations: same result every run
    except ImportError:
        pass
    return texts


def isolation(dist):
    """Mean distance from each cluster centre to all other centres (input: centroid_distances)."""
    rows = [(name, d) for pair, d in dist.items() for name in pair.split(" | ")]
    table = pd.DataFrame(rows, columns=["cluster", "distance"])
    return table.groupby("cluster")["distance"].mean().sort_values(ascending=False)


def same_label_neighbours(emb, labels, k=15):
    """Share of each cell's k nearest neighbours (in `emb`) that carry the same label as the cell."""
    from sklearn.neighbors import NearestNeighbors
    labels = np.asarray(labels)
    nn = NearestNeighbors(n_neighbors=k).fit(emb).kneighbors(return_distance=False)
    return float((labels[nn] == labels[:, None]).mean())


def visible_share(xy, values, order, bins=60):
    """Put a bins x bins grid over the plot. In each square a reader sees only the cell drawn last.
    Returns the share of occupied squares whose visible (last-drawn) cell has a value > 0."""
    xy, v = np.asarray(xy)[order], np.asarray(values)[order]
    ij = np.floor((xy - xy.min(axis=0)) / np.ptp(xy, axis=0) * (bins - 1e-9)).astype(int)
    square = ij[:, 0] * bins + ij[:, 1]
    top = pd.Series(v).groupby(square).last()        # rows are in drawing order: last = on top
    return float((top > 0).mean())
