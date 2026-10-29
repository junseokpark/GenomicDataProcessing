"""Fixed cell labels for the single-cell data (B4) and the composition table (Practice A4).

2,940 PBMC-like cells from 4 donors, 8 cell types. The overall composition is fixed exactly; donors
differ a little (Dirichlet split per cell type), as in real blood samples.
Writes data/truth/scrna_truth.tsv (cell_id, donor, cell_type) and data/celltype_composition.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20261027 + 50)
TYPES = {"Naive CD4 T": 706, "CD14 Mono": 559, "Memory CD4 T": 479, "B": 353,
         "CD8 T": 338, "NK": 235, "FCGR3A Mono": 191, "DC": 79}
DONORS = ["D1", "D2", "D3", "D4"]
assert sum(TYPES.values()) == 2940

rows = []
for t, n in TYPES.items():
    share = rng.dirichlet(np.full(4, 40.0))                    # donor shares for this type (mild variation)
    raw = share * n
    k = np.floor(raw).astype(int)
    k[np.argsort(-(raw - k))[: n - k.sum()]] += 1               # largest remainder -> exact total
    for d, kk in zip(DONORS, k):
        rows += [(d, t)] * int(kk)
cells = pd.DataFrame(rows, columns=["donor", "cell_type"]).sample(frac=1, random_state=7).reset_index(drop=True)
bases = np.array(list("ACGT"))
codes = set()
ids = []
while len(ids) < len(cells):
    bc = "".join(rng.choice(bases, 16)) + "-1"                  # 10x-style barcode
    if bc not in codes:
        codes.add(bc); ids.append(bc)
cells.insert(0, "cell_id", ids)
(ROOT / "data" / "truth").mkdir(parents=True, exist_ok=True)
cells.to_csv(ROOT / "data" / "truth" / "scrna_truth.tsv", sep="\t", index=False)

comp = cells.groupby(["donor", "cell_type"]).size().rename("n_cells").reset_index()
comp["percent"] = (100 * comp["n_cells"] / comp.groupby("donor")["n_cells"].transform("sum")).round(1)
order = list(TYPES)
comp["cell_type"] = pd.Categorical(comp["cell_type"], order, ordered=True)
comp = comp.sort_values(["donor", "cell_type"]).reset_index(drop=True)
comp["cell_type"] = comp["cell_type"].astype(str)
comp.to_csv(ROOT / "data" / "celltype_composition.tsv", sep="\t", index=False)
print(cells.cell_type.value_counts().to_string())
print(cells.donor.value_counts().sort_index().to_string())
print(comp.pivot(index="cell_type", columns="donor", values="percent").loc[order].to_string())
