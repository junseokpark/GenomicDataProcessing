"""Numbers behind the folded "Your turn" answers of Practice B4 -> results/B4_yourturn.json.

Runs the notebook's pipeline (same helpers, same settings) and the small variations the students are
asked to try: UMAP with seed 1, k-means on the 2-D UMAP coordinates, PCA with every gene scaled.
usage:  python build/answers_B4.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import scrna_helpers as H  # noqa: E402

SEED = 0
counts = sp.load_npz(ROOT / "data" / "scrna_counts.npz").tocsr()
cells = pd.read_csv(ROOT / "data" / "scrna_cells.tsv", sep="\t")
truth = cells[["cell_id"]].merge(pd.read_csv(ROOT / "data" / "truth" / "scrna_truth.tsv", sep="\t"),
                                 on="cell_id", how="left", validate="one_to_one")["cell_type"].to_numpy()
keep = np.asarray((counts > 0).sum(axis=0)).ravel() >= 3
lognorm, scaled = H.normalise_log1p(counts[:, keep])
hvg, _ = H.highly_variable_genes(scaled, n_top=1000)
pcs, _ = H.run_pca(lognorm[:, hvg], n_pcs=30, seed=SEED)
cluster = KMeans(n_clusters=8, n_init=10, random_state=SEED).fit_predict(pcs)

umap0, _ = H.run_umap(pcs, n_neighbors=15, min_dist=0.3, seed=0)
umap1, _ = H.run_umap(pcs, n_neighbors=15, min_dist=0.3, seed=1)
nn0 = NearestNeighbors(n_neighbors=15).fit(umap0).kneighbors(return_distance=False)
nn1 = NearestNeighbors(n_neighbors=15).fit(umap1).kneighbors(return_distance=False)
same = np.mean([np.intersect1d(a, b).size / 15 for a, b in zip(nn0, nn1)])

pcs_scaled, _ = H.run_pca(lognorm[:, hvg], n_pcs=30, scale=True, seed=SEED)
out = {
    "description": "Practice B4 - numbers for the folded 'Your turn' answers (same pipeline and settings as the notebook)",
    "kmeans_on_pcs_ari": round(adjusted_rand_score(truth, cluster), 3),
    "seed0_knn15_kept_pct": round(100 * H.knn_overlap(pcs, umap0, k=15), 1),
    "seed1_knn15_kept_pct": round(100 * H.knn_overlap(pcs, umap1, k=15), 1),
    "seed_neighbours_same_pct": round(100 * same, 1),
    "seed0_same_cluster_neighbours_pct": round(100 * H.same_label_neighbours(umap0, cluster), 1),
    "seed1_same_cluster_neighbours_pct": round(100 * H.same_label_neighbours(umap1, cluster), 1),
    "kmeans_on_umap_ari": round(adjusted_rand_score(
        truth, KMeans(n_clusters=8, n_init=10, random_state=SEED).fit_predict(umap0)), 3),
    "kmeans_scaled_pca_ari": round(adjusted_rand_score(
        truth, KMeans(n_clusters=8, n_init=10, random_state=SEED).fit_predict(pcs_scaled)), 3),
}
# two common mistakes, measured: truth aligned by row instead of joined by barcode; PCA without log1p
truth_by_row = pd.read_csv(ROOT / "data" / "truth" / "scrna_truth.tsv", sep="\t")["cell_type"].to_numpy()
pcs_nolog, _ = H.run_pca(scaled[:, hvg], n_pcs=30, seed=SEED)
out["mistake_truth_by_row_order_ari"] = round(adjusted_rand_score(truth_by_row, cluster), 3)
out["mistake_no_log1p_ari"] = round(adjusted_rand_score(
    truth, KMeans(n_clusters=8, n_init=10, random_state=SEED).fit_predict(pcs_nolog)), 3)
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "results" / "B4_yourturn.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
