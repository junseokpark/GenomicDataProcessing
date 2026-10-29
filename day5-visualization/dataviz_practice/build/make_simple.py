"""Practice A (and B1) data: small tables you can read by eye.

Writes to data/:
  anscombe.tsv            Anscombe's quartet (real published values, Anscombe 1973)
  mapping_rates.tsv       % reads mapped for 3 pipeline versions x 4 libraries
  expression_by_group.tsv one gene (TPM) in Control vs Case: n = 4 + 4 and n = 40 + 40
  timecourse.tsv          3 glucocorticoid-response genes x 6 time points x 3 replicates (log2 FC vs 0 h)
  methylation_cpg.tsv.gz  60,000 CpGs measured by bisulfite sequencing (WGBS) and nanopore (beta values)
(celltype_composition.tsv is written by make_cell_labels.py from the single-cell labels.)
All values except Anscombe are simulated (seed 20261027 + 1).
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
rng = np.random.default_rng(20261027 + 1)


def write(df, name, **kw):
    p = DATA / name
    if name.endswith(".gz"):
        df.to_csv(p, sep="\t", index=False, compression={"method": "gzip", "mtime": 0}, **kw)
    else:
        df.to_csv(p, sep="\t", index=False, **kw)
    print(f"{name:28s} {len(df):>7,} rows")


# 1. Anscombe's quartet (Anscombe FJ, 1973, Am Stat 27:17-21) -------------------------------------
x123 = [10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5]
ans = {
    "I": (x123, [8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68]),
    "II": (x123, [9.14, 8.14, 8.74, 8.77, 9.26, 8.10, 6.13, 3.10, 9.13, 7.26, 4.74]),
    "III": (x123, [7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73]),
    "IV": ([8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8], [6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.50, 5.56, 7.91, 6.89]),
}
write(pd.DataFrame([(k, x, y) for k, (xs, ys) in ans.items() for x, y in zip(xs, ys)],
                   columns=["dataset", "x", "y"]), "anscombe.tsv")

# 2. Mapping rates: three versions of one alignment pipeline on the same four libraries --------------
rates = {"A": [95.9, 96.3, 96.0, 96.2], "B": [97.2, 97.6, 97.3, 97.5], "C": [97.9, 98.2, 97.8, 98.1]}
write(pd.DataFrame([(p, f"L{i + 1}", v) for p, vs in rates.items() for i, v in enumerate(vs)],
                   columns=["pipeline", "library", "pct_mapped"]), "mapping_rates.tsv")

# 3. One gene in two groups -----------------------------------------------------------------------------
# small_n: 4 vs 4, chosen so the two tests disagree (one Case sample at 150 TPM)
small = [("Control", v) for v in [55.3, 59.6, 62.0, 66.4]] + [("Case", v) for v in [71.2, 84.1, 95.5, 150.0]]
# large_n: 40 vs 40; Cases are a mixture of non-responders (like controls) and responders
ctrl = rng.normal(61.0, 6.0, 40)
resp = rng.random(40) < 0.42
case = np.where(resp, rng.normal(106.0, 8.0, 40), rng.normal(62.5, 6.0, 40))
rows = [("small_n", g, f"{'C' if g == 'Control' else 'K'}{i % 4 + 1:02d}", v) for i, (g, v) in enumerate(small)]
rows += [("large_n", "Control", f"C{i + 1:02d}", v) for i, v in enumerate(ctrl)]
rows += [("large_n", "Case", f"K{i + 1:02d}", v) for i, v in enumerate(case)]
ebg = pd.DataFrame(rows, columns=["dataset", "group", "sample", "expression"])
ebg["expression"] = ebg["expression"].round(1)
write(ebg, "expression_by_group.tsv")
truth = pd.DataFrame({"sample": [f"K{i + 1:02d}" for i in range(40)], "responder": resp.astype(int)})
truth.to_csv(DATA / "truth" / "expression_by_group_truth.tsv", sep="\t", index=False)

# 4. Time course after dexamethasone (log2 fold change vs 0 h), 3 replicate wells per time ----------
times = [0, 1, 2, 4, 8, 24]
shape = {  # mean log2 FC: fast-sustained, slow-rising, fast-transient (illustrative kinetics)
    "TSC22D3": [0.0, 1.6, 2.4, 2.9, 3.1, 3.0],
    "FKBP5": [0.0, 0.4, 1.0, 1.9, 2.8, 3.6],
    "PER1": [0.0, 2.2, 1.9, 1.3, 0.9, 0.6],
}
rows = []
for g, mu in shape.items():
    for t, m in zip(times, mu):
        for r in range(3):
            rows.append((g, t, r + 1, round(m + rng.normal(0, 0.12 + 0.05 * m), 2) if t else round(rng.normal(0, 0.08), 2)))
write(pd.DataFrame(rows, columns=["gene", "time_h", "replicate", "log2fc_vs_0h"]), "timecourse.tsv")

# 5. CpG methylation: WGBS vs nanopore on the same sample ---------------------------------------------
n = 60_060
kind = rng.choice(3, n, p=[0.46, 0.44, 0.10])               # unmethylated, methylated, intermediate
m = np.where(kind == 0, rng.beta(0.6, 14, n), np.where(kind == 1, rng.beta(14, 0.7, n), rng.uniform(0.15, 0.85, n)))
d_wgbs = rng.poisson(18, n) + 4                               # read depth per CpG
d_ont = rng.poisson(24, n) + 4
b_wgbs = rng.binomial(d_wgbs, m) / d_wgbs                     # beta = methylated reads / all reads
m_ont = np.clip(0.015 + 0.97 * m + rng.normal(0, 0.035, n), 0, 1)   # small calling bias + model noise
b_ont = rng.binomial(d_ont, m_ont) / d_ont
mane = pd.read_csv(ROOT / "build" / "ref" / "genes_hg38_mane_select.tsv.gz", sep="\t")
mane = mane[mane["chrom"].isin([f"chr{i}" for i in range(1, 23)])].reset_index(drop=True)
g = mane.iloc[rng.integers(0, len(mane), n)]
tss = np.where(g["strand"].to_numpy() == "+", g["start0"].to_numpy(), g["end"].to_numpy())
pos = tss + rng.integers(-20_000, 20_000, n)                  # CpGs near genes (1-based positions)
meth = pd.DataFrame({"chrom": g["chrom"].to_numpy(), "pos": np.maximum(pos, 1),
                     "beta_wgbs": b_wgbs.round(4), "depth_wgbs": d_wgbs,
                     "beta_nanopore": b_ont.round(4), "depth_nanopore": d_ont})
meth["_k"] = meth["chrom"].str.replace("chr", "").astype(int)
meth = meth.drop_duplicates(["chrom", "pos"]).sample(n=60_000, random_state=11)          # exactly 60,000 CpGs
meth = meth.sort_values(["_k", "pos"]).drop(columns="_k").reset_index(drop=True)
write(meth, "methylation_cpg.tsv.gz")
