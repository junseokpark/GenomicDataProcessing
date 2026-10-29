"""data/rnaseq_log2norm.tsv.gz: log2(median-of-ratios normalised count + 1) for expressed genes (mean normalised count >= 10).

Same formula as notebook Step B2.2/B2.3, provided so that AI-written heatmap/PCA code (prompt cards 3 and 4) starts from
the same matrix as the CHECK cells. Deterministic.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "data"
meta = pd.read_csv(D / "rnaseq_metadata.tsv", sep="\t")
counts = pd.read_csv(D / "rnaseq_counts.tsv.gz", sep="\t", index_col="gene").drop(columns=["gene_id", "chrom"])[meta["sample"]]
expressed = (counts > 0).all(axis=1)
log_geo_mean = np.log(counts[expressed]).mean(axis=1)
size_factors = np.exp((np.log(counts[expressed]).sub(log_geo_mean, axis=0)).median(axis=0))
norm = counts / size_factors
keep = norm.mean(axis=1) >= 10
out = np.log2(norm[keep] + 1).round(4).reset_index()
out.to_csv(D / "rnaseq_log2norm.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
print(f"rnaseq_log2norm.tsv.gz: {len(out):,} expressed genes x {len(meta)} samples")
