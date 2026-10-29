"""B2 data: bulk RNA-seq, 12 libraries (6 control, 6 dexamethasone-treated) in two sequencing batches.

Simulated with real anchors: 8,000 real protein-coding genes (MANE Select, hg38); the strongest induced genes
include known glucocorticoid-response genes (e.g. FKBP5, TSC22D3, ZBTB16, KLF15, PER1; Himes et al. 2014);
all other DE genes (including every down-regulated gene) are random. Counts are negative binomial with
realistic size factors and a mean-dispersion trend; batch shifts some genes but is balanced with treatment.

Writes data/rnaseq_counts.tsv.gz, data/rnaseq_metadata.tsv, data/de_results.tsv (PyDESeq2, design
~ batch + condition, apeGLM-shrunken log2FC), data/de_results_no_batch.tsv (design ~ condition) and
data/truth/rnaseq_truth.tsv. Seed 20261027 + 3.
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
rng = np.random.default_rng(20261027 + 3)
N_GENES = 8000
GR_UP = ["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "PER1", "DUSP1", "CRISPLD2", "KLF9", "ERRFI1",
         "CEBPD", "SGK1", "MT2A", "ALOX15B", "CDKN1C", "PDK4"]
P = dict(low_frac=0.08, low_max=0.03, mu_med=180.0, mu_sd=1.75, disp0=0.08, disp1=0.9, disp_sd=0.3,
         n_de=600, frac_up=0.57, scale=3.0, latent_sd=0.10, gr_mu=1000.0, de_mu_max=3000.0, large_max=3.4, batch_frac=0.35, batch_sd=0.65, strong_batch=150, strong_batch_lfc=1.1)
import os
P.update(json.loads(os.environ.get("RNASEQ_PARAMS", "{}")))          # tuning hook (defaults = shipped data)
FAST = os.environ.get("RNASEQ_FAST") == "1"

# ---- genes --------------------------------------------------------------------------------------------
mane = pd.read_csv(ROOT / "build" / "ref" / "genes_hg38_mane_select.tsv.gz", sep="\t")
mane = mane[mane["chrom"].isin([f"chr{i}" for i in range(1, 23)] + ["chrX"])].drop_duplicates("gene")
missing = [g for g in GR_UP if g not in set(mane["gene"])]
gr = [g for g in GR_UP if g not in missing]
others = mane[~mane["gene"].isin(gr)].sample(n=N_GENES - len(gr), random_state=int(rng.integers(1e9)))
genes = pd.concat([mane[mane["gene"].isin(gr)], others])
genes["_k"] = genes["chrom"].str.replace("chr", "").replace({"X": "23"}).astype(int)
genes = genes.sort_values(["_k", "start0"]).reset_index(drop=True)
G = len(genes)

# ---- samples: balanced design, 3 control + 3 treated in each batch ------------------------------------
meta = pd.DataFrame({
    "sample": [f"{c}_{b}{i}" for b in "AB" for c in ["ctrl", "trt"] for i in (1, 2, 3)],
    "condition": [c for b in "AB" for c in ["control", "treated"] for _ in range(3)],
    "batch": [b for b in "AB" for _ in range(6)]})
meta["RIN"] = np.round(rng.uniform(8.4, 9.8, len(meta)), 1)
J = len(meta)

# ---- baseline expression, dispersion, effects -----------------------------------------------------------
# genes that are not expressed in airway smooth-muscle cells (tissue-restricted families) -> barely any reads
TISSUE_RESTRICTED = r"^(OR\d|TAS[12]R|KRTAP|LCE\d|PRAMEF|SPANX|CT45|CT47|GAGE|PAGE\d|XAGE|SSX\d|MAGE[ABC]|DEFB|SPRR|IFNA|IFNW|TSPY|RBMY|VCX|VCY|USP17L|NBPF|FAM90A|HSFX|SPATA|TEX\d|ZP\d|OVOL3|CRYG|PRM\d|TNP\d)"
restricted = genes["gene"].str.match(TISSUE_RESTRICTED).to_numpy()
low = (rng.random(G) < P["low_frac"]) | restricted
mu0 = np.where(low, rng.uniform(0.0005, P["low_max"], G), np.exp(rng.normal(np.log(P["mu_med"]), P["mu_sd"], G)))
is_gr = genes["gene"].isin(gr).to_numpy()
mu0[is_gr] = np.exp(rng.normal(np.log(P["gr_mu"]), 0.6, is_gr.sum()))         # GR genes are well expressed
alpha = (P["disp0"] + P["disp1"] / np.maximum(mu0, 0.05)) * np.exp(rng.normal(0, P["disp_sd"], G))
alpha = np.clip(alpha, 0.005, 8)

lfc = np.zeros(G)
cand = np.flatnonzero(~low & ~is_gr & (mu0 >= 20) & (mu0 <= P["de_mu_max"]))   # random DE genes: moderate expression
de_idx = rng.choice(cand, P["n_de"] - is_gr.sum(), replace=False)
sign = np.where(rng.random(de_idx.size) < P["frac_up"], 1.0, -1.0)
size_class = rng.choice(3, de_idx.size, p=[0.40, 0.35, 0.25])
mag = np.where(size_class == 0, rng.uniform(0.3, 0.8, de_idx.size),
               np.where(size_class == 1, rng.uniform(0.8, 1.6, de_idx.size), rng.uniform(1.6, P["large_max"], de_idx.size)))
lfc[de_idx] = sign * mag
lfc[is_gr] = rng.uniform(2.6, 4.6, is_gr.sum())

batch_lfc = np.where(rng.random(G) < P["batch_frac"], rng.normal(0, P["batch_sd"], G), 0.0)
strong = rng.choice(np.flatnonzero(~low), P["strong_batch"], replace=False)   # e.g. library-prep sensitive genes
batch_lfc[strong] = rng.choice([-1, 1], strong.size) * rng.uniform(0.8, 1.5, strong.size) * P["strong_batch_lfc"]

sf = np.exp(rng.normal(0, 0.25, J))
sf = sf / np.exp(np.mean(np.log(sf)))
cond = (meta["condition"] == "treated").to_numpy().astype(float)
bat = (meta["batch"] == "B").to_numpy().astype(float)
# two latent "culture" factors: cultures differ coherently across many genes (realistic replicate scatter)
F = rng.normal(0, 1, (J, 2))
D = np.column_stack([np.ones(J), bat, cond])
F = F - D @ np.linalg.lstsq(D, F, rcond=None)[0]           # orthogonal to the design: within-group variation only
F = F / F.std(axis=0)
latent = F @ rng.normal(0, P["latent_sd"], (2, G))
log2mu = np.log2(mu0)[:, None] + lfc[:, None] * cond[None, :] + batch_lfc[:, None] * bat[None, :] + latent.T
mean = (2.0 ** log2mu) * sf[None, :] * P["scale"]                             # ~15-25 M reads per library
shape = 1.0 / alpha[:, None]
counts = rng.poisson(rng.gamma(shape, mean / shape))
cdf = pd.DataFrame(counts.T, index=meta["sample"], columns=genes["gene"])

# ---- write inputs ---------------------------------------------------------------------------------------
out = cdf.T.reset_index().rename(columns={"index": "gene"})
out.insert(1, "gene_id", genes["gene_id"].str.split(".").str[0].to_numpy())
out.insert(2, "chrom", genes["chrom"].to_numpy())
out.to_csv(DATA / "rnaseq_counts.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
meta.to_csv(DATA / "rnaseq_metadata.tsv", sep="\t", index=False)
truth = pd.DataFrame({"gene": genes["gene"], "is_de": (lfc != 0).astype(int), "true_log2fc": lfc.round(4),
                      "batch_log2_shift": batch_lfc.round(4), "baseline_mean": mu0.round(3),
                      "known_gr_target": is_gr.astype(int)})
truth.to_csv(DATA / "truth" / "rnaseq_truth.tsv", sep="\t", index=False)

# ---- DESeq2 (PyDESeq2) ----------------------------------------------------------------------------------
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
inf = DefaultInference(n_cpus=1)                       # 1 cpu -> deterministic
md = meta.set_index("sample")[["condition", "batch"]]


def run_deseq(design, shrink=True):
    dds = DeseqDataSet(counts=cdf, metadata=md, design=design, refit_cooks=True, inference=inf, quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=["condition", "treated", "control"], inference=inf, quiet=True)
    ds.summary()
    res = ds.results_df.copy()
    if shrink:
        coef = [c for c in dds.obsm["design_matrix"].columns if c.startswith("condition")][0]
        ds.lfc_shrink(coeff=coef)
        res["log2FC_shrunk"] = ds.results_df["log2FoldChange"]
    return res, dds


res, dds = run_deseq("~batch + condition")
res0 = run_deseq("~condition", shrink=False)[0] if not FAST else res.copy()


def table(r):
    t = pd.DataFrame({"gene": r.index, "gene_id": genes.set_index("gene").loc[r.index, "gene_id"].str.split(".").str[0].to_numpy(),
                      "baseMean": r["baseMean"].round(3), "log2FC": r["log2FoldChange"].round(5)})
    if "log2FC_shrunk" in r:
        t["log2FC_shrunk"] = r["log2FC_shrunk"].round(5)
    t["lfcSE"] = r["lfcSE"].round(5)
    t["pvalue"] = r["pvalue"].map(lambda v: float(f"{v:.6g}") if pd.notna(v) else np.nan)
    t["padj"] = r["padj"].map(lambda v: float(f"{v:.6g}") if pd.notna(v) else np.nan)
    return t.reset_index(drop=True)


de, de0 = table(res), table(res0)
de.to_csv(DATA / "de_results.tsv", sep="\t", index=False, na_rep="")
de0.to_csv(DATA / "de_results_no_batch.tsv", sep="\t", index=False, na_rep="")

# ---- summary for tuning ---------------------------------------------------------------------------------
ok = de["padj"].notna()
up = ok & (de["padj"] < 0.05) & (de["log2FC_shrunk"] > 1)
dn = ok & (de["padj"] < 0.05) & (de["log2FC_shrunk"] < -1)
sig = ok & (de["padj"] < 0.05)
sig0 = de0["padj"].notna() & (de0["padj"] < 0.05)
tr = truth.set_index("gene").loc[de["gene"]]
fdp = float((sig & (tr["is_de"].to_numpy() == 0)).sum() / max(sig.sum(), 1))
norm = cdf.div(dds.obs["size_factors"] if "size_factors" in dds.obs else dds.obsm["size_factors"], axis=0)
ln = np.log2(norm.T + 1)
keep = norm.mean(axis=0) >= 10
X = ln.loc[keep.to_numpy()]
X = X.loc[X.var(axis=1).nlargest(1000).index].T
Xc = X - X.mean()
u, s, vt = np.linalg.svd(Xc.to_numpy(), full_matrices=False)
var = s ** 2 / np.sum(s ** 2) * 100
pcs = u * s
pc_cond = [abs(np.corrcoef(pcs[:, k], cond)[0, 1]) for k in range(3)]
pc_bat = [abs(np.corrcoef(pcs[:, k], bat)[0, 1]) for k in range(3)]
# by-hand per-gene regression on log2 normalised counts (what students do in B2.4)
filt = (norm.mean(axis=0) >= 10).to_numpy()
Y = ln.loc[filt].to_numpy().T                                   # samples x genes
def ols_p(Xd):
    beta, *_ = np.linalg.lstsq(Xd, Y, rcond=None)
    resid = Y - Xd @ beta
    dfree = Xd.shape[0] - Xd.shape[1]
    s2 = (resid ** 2).sum(axis=0) / dfree
    se = np.sqrt(s2 * np.linalg.inv(Xd.T @ Xd)[-1, -1])
    from scipy import stats as st
    return 2 * st.t.sf(np.abs(beta[-1] / se), dfree)
one = np.ones(J)
p_b = ols_p(np.column_stack([one, bat, cond])); p_n = ols_p(np.column_stack([one, cond]))
sys.path.insert(0, str(ROOT)); import dviz
isde = truth.set_index("gene").loc[ln.index[filt], "is_de"].to_numpy()
def fdp_power(p):
    q = dviz.bh_adjust(p); sg = q < 0.05
    return int(sg.sum()), round(float((sg & (isde == 0)).sum() / max(sg.sum(), 1)), 4), int((sg & (isde == 1)).sum())
nullp = de["pvalue"].to_numpy()[(tr["is_de"].to_numpy() == 0) & de["pvalue"].notna().to_numpy()]
summary = {"ols_with_batch(sig,fdp,tp)": fdp_power(p_b), "ols_no_batch(sig,fdp,tp)": fdp_power(p_n),
           "ols_null_p05_with_batch": round(float((p_b[isde == 0] < 0.05).mean()), 4),
           "deseq_null_p05": round(float((nullp < 0.05).mean()), 4), "deseq_null_p01": round(float((nullp < 0.01).mean()), 4),
           "volcano_fdp": round(float(((up | dn) & (tr["is_de"].to_numpy() == 0)).sum() / max((up | dn).sum(), 1)), 4),
           "max_mlog10_padj": round(float(-np.log10(de["padj"].min())), 1), "top8": list(de[ok].nsmallest(8, "padj")["gene"]),
           "genes": G, "restricted_genes": int(restricted.sum()), "libsize_M": [round(v, 1) for v in (cdf.sum(axis=1) / 1e6)], "missing_gr": missing, "tested": int(ok.sum()), "untested": int((~ok).sum()),
           "sig_fdr05": int(sig.sum()), "sig_fdr05_no_batch": int(sig0.sum()), "up": int(up.sum()), "down": int(dn.sum()),
           "fdp_among_sig": round(fdp, 4), "planted_de": int((lfc != 0).sum()),
           "pc_var": [round(v, 1) for v in var[:4]], "pc_corr_condition": [round(v, 2) for v in pc_cond],
           "pc_corr_batch": [round(v, 2) for v in pc_bat], "size_factors": [round(v, 3) for v in sf]}
print(json.dumps(summary, indent=1))
json.dump(summary, open(ROOT / "build" / "rnaseq_build_summary.json", "w"), indent=1)
