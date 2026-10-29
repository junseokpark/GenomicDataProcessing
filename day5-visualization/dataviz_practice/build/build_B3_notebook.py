"""Write electives/B3_gwas.ipynb (elective B3 · GWAS: association along the genome).

usage:  python build/build_B3_notebook.py
The numbers quoted in the markdown cells ("What changed?", answers, summary) are read from results/B3.json,
which the executed notebook writes. Build twice: build -> execute -> build -> execute (the data are
deterministic, so the second execution reproduces the same JSON). Missing numbers are shown as '…'.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from nbbuild import NB  # noqa: E402

RES = ROOT / "results" / "B3.json"
R = json.loads(RES.read_text()) if RES.exists() else {}
NUM = R.get("numbers", {})


MISSING = "[?]"                                  # shown until results/B3.json exists


def n(key, fmt="{}"):
    """A number from results/B3.json ('numbers' block) formatted for the text."""
    v = NUM.get(key)
    return MISSING if v is None else fmt.format(v)


def p_sci(v):
    """p-value as 4.2 × 10⁻⁴⁶ (Unicode superscripts, as on the slides)."""
    m, e = f"{v:.1e}".split("e")
    return f"{m} × 10{str(int(e)).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"


def sci(key):
    v = NUM.get(key)
    return MISSING if v is None else p_sci(v)


def listing(key, fmt):
    """Join the records of a list in the 'numbers' block with a format function."""
    L = NUM.get(key)
    return MISSING if L is None else ("none" if not L else "; ".join(fmt(d) for d in L))


loci_text = lambda: listing("loci", lambda d: f"**{d['gene']}** (p = {p_sci(d['p'])}, OR {d['OR']:.2f})")
missed_text = lambda: listing("missed", lambda d: f"{d['locus']}, planted OR {d['OR']:.2f}, best p = {d['best_p']:.1e}")
false_text = lambda: listing("false_noPC", lambda d: f"chr{d['CHR']}:{d['BP']:,} (nearest gene {d['gene']}, "
                                                     f"p = {p_sci(d['p_noPC'])} without PCs)")
aim_pc_text = lambda: " and ".join(f"{d['p_PC']:.2f}" for d in NUM.get("false_noPC", [])) or MISSING
def cfh_block():
    if "cfh_block_kb" not in NUM:
        return MISSING
    a, b = NUM["cfh_block_kb"]
    if a >= 0:
        return f"on the right of the lead only, up to {b:.0f} kb away"
    if b <= 0:
        return f"on the left of the lead only, up to {-a:.0f} kb away"
    return f"from {-a:.0f} kb left to {b:.0f} kb right of the lead"
max_tower = lambda: MISSING if not NUM.get("loci") else str(max(d["n_snps"] for d in NUM["loci"]))
lead_vs_causal = lambda: (MISSING if "cfh_lead_is_causal" not in NUM else
                          "is the planted causal SNP itself" if NUM["cfh_lead_is_causal"] else
                          f"is a proxy of the planted causal SNP (r² = {NUM['cfh_causal_r2']:.2f})")
lead_here = lambda: (MISSING if "cfh_lead_is_causal" not in NUM else
                     "the lead happens to be the planted causal SNP" if NUM["cfh_lead_is_causal"] else
                     "the lead is not the planted causal SNP")
strongest = lambda: (MISSING if not NUM.get("loci") else " and ".join(
    f"{d['gene']} (OR {d['OR']:.2f})" for d in sorted(NUM["loci"], key=lambda d: d["p"])[:2]))
missed_or = lambda: (MISSING if not NUM.get("missed") else
                     f"{min(d['OR'] for d in NUM['missed']):.2f}–{max(d['OR'] for d in NUM['missed']):.2f}")
TRUTH_JSON = ROOT / "data" / "truth" / "gwas_truth.json"
THIN_KB = (f"{json.loads(TRUTH_JSON.read_text())['design']['thin_window_bp'] / 1e3:.0f}" if TRUTH_JSON.exists() else MISSING)
block_genes = lambda: (MISSING if "cfh_block_genes" not in NUM else
                       " and ".join(NUM["cfh_block_genes"]) if NUM["cfh_block_genes"] else "no gene body")


nb = NB()
# ------------------------------------------------------------------------------------------------ title
nb.md(f"""
# Practice B3 · GWAS: association along the genome

**Genotype PCA → trend test with and without ancestry PCs → λ and QQ plot → Manhattan plot → regional plot at CFH**

*Genomic Data Processing · Day 5 (10/27) · Data visualization in genetics and genomics · about 20 minutes*

> **All data in this practice are simulated (computer-generated). No real patient data are used.**
> The strong loci were placed inside real genes (hg38 MANE Select coordinates) that the largest GWAS of age-related
> macular degeneration reported (Fritsche et al., 2016, *Nat Genet* 48:134–143), so the gene names you meet are the ones in
> the AMD literature; all other SNPs, alleles and people are random.

**How to use this notebook.** Run the cells from top to bottom (Shift + Enter). Each step has four beats:
**Question** (what we want to know), **Test** (run the method and print numbers), **Check** (a `# CHECK ·` cell compares
your number with the planted truth or a reference: ✓ = as expected, ✗ = stop and fix before you go on) and **Visualize**
(the figure that shows the result honestly). **Your turn** asks you to *predict first* and then run the short cell below
it; the answer is folded (click **Answer**). The files in `data/truth/` hold the planted answers – open them only in
CHECK cells. The whole notebook runs in well under a minute.
""")
nb.code("""
# Step 0 · Set-up: imports, folders and figure style
import os, sys, time
from pathlib import Path
ROOT = Path.cwd() if (Path.cwd() / "dviz.py").exists() else Path.cwd().parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))   # needed by umap-learn on some systems
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dviz
dviz.use_style()
DATA, FIG = ROOT / "data", ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
T_START = time.time()
import json
from scipy.stats import chi2                      # chi-square distribution: p-values of the trend test
from sklearn.decomposition import PCA             # principal component analysis
""")

# ------------------------------------------------------------------------------------------------ B3.1
nb.md("""
## Step B3.1 · Question: what a GWAS asks, and why ancestry matters  *(3 min)*

**Why this step.** Age-related macular degeneration (AMD) damages the centre of the retina and is the commonest cause of
severe sight loss in older adults. A **genome-wide association study (GWAS)** asks one simple question at hundreds of
thousands of places in the genome: *is one allele of this SNP more common in people with the disease (cases) than in
people without it (controls)?* But allele frequencies also differ between **genetic-ancestry groups**. If cases come
more often from one group than controls do, almost every SNP looks a little "associated" – for reasons that have
nothing to do with the eye.

**New terms**
* **SNP** – single-nucleotide polymorphism: one position where people carry different letters (alleles).
* **A1 / A2** – A1 is the allele we count (0, 1 or 2 copies per person; here the alt allele), A2 is the other one.
* **Case / control** – a person with / without the disease.
* **Genetic-ancestry group** – people whose genomes share recent ancestry and therefore similar allele frequencies
  (a genetic description, not a social label).
* **Population stratification** – case–control differences in allele frequency caused by a different ancestry mix in
  cases and controls rather than by the disease.

### Question
*Which regions of the genome carry alleles that change the risk of this AMD-like disease – and how do we avoid being
fooled by ancestry?* First, who is in the study?
""")
nb.code("""
# Step B3.1 · load the table of people and count cases and controls per recruitment site
samples = pd.read_csv(DATA / "gwas_samples.tsv", sep="\\t")    # IID, case (1/0), site, sex, age
print(samples.head(3).to_string(index=False))
y = samples["case"].to_numpy(dtype=float)                       # 1 = case, 0 = control
tab = pd.crosstab(samples["site"], samples["case"].map({1: "case", 0: "control"}), margins=True)
tab["% cases"] = (100 * tab["case"] / tab["All"]).round(1)
tab
""")
nb.code("""
# CHECK · 3,000 people: 1,500 cases and 1,500 controls
dviz.check("B3.n_people", len(samples), 3000)
dviz.check("B3.n_cases", int(samples["case"].sum()), 1500)
dviz.check("B3.n_controls", int((samples["case"] == 0).sum()), 1500)
""")
nb.md(f"""
**What changed?** Nothing in the data yet – we only looked. The two recruitment sites differ: {n('pct_cases_S1', '{:.1f}')}%
of the people from site S1 are cases, but {n('pct_cases_S2', '{:.1f}')}% of those from site S2. If the sites also differ in
ancestry, we are in exactly the situation that creates stratification.

**Your turn.** A SNP has **no** effect on the disease. Its allele frequency is 0.30 in ancestry group A and 0.40 in
group B. Cases are 50% group A, controls 70% group A. *Predict* the allele frequency in cases and in controls, then run
the cell below.

<details><summary>Answer</summary>

Cases: 0.5 × 0.30 + 0.5 × 0.40 = 0.350. Controls: 0.7 × 0.30 + 0.3 × 0.40 = 0.330. A difference of
{n('yt_toy_difference', '{:+.3f}')} with no biology at all. Over ~300,000 SNPs such small differences lift *every* test a little; a few
SNPs whose frequencies differ a lot between groups can even pass genome-wide significance.
</details>
""")
nb.code("""
# Your turn · B3.1 · a SNP with NO effect on disease, but different frequencies in two ancestry groups
f_A, f_B = 0.30, 0.40                              # allele frequency in group A and in group B
share_A_cases, share_A_controls = 0.50, 0.70       # fraction of group-A people among cases / controls
f_cases = share_A_cases * f_A + (1 - share_A_cases) * f_B
f_controls = share_A_controls * f_A + (1 - share_A_controls) * f_B
print(f"cases {f_cases:.3f}   controls {f_controls:.3f}   difference {f_cases - f_controls:+.3f}")
""")

# ------------------------------------------------------------------------------------------------ B3.2
nb.md(f"""
## Step B3.2 · Genotype PCA: ancestry is written in the genotypes  *(4 min)*

**Why this step.** We are not told anybody's ancestry, but their genotypes contain it: across thousands of SNPs, people
from the same group share many small frequency differences. **Principal component analysis (PCA)** finds the directions
along which people differ most; in genotype data the first PCs usually measure ancestry. We use SNPs spread out along the
genome (**LD-thinned**), so that one region with many correlated SNPs cannot dominate a PC.

**New terms**
* **Genotype matrix** – a people × SNPs table of 0/1/2 (copies of A1).
* **Standardise** – for each SNP, subtract its mean and divide by its standard deviation, so every SNP counts equally.
* **Principal component (PC)** – a weighted sum of all SNPs; PC1 explains the most variance, PC2 the next most, …
* **LD (linkage disequilibrium)** – nearby SNPs are inherited together, so their genotypes are correlated.
* **LD thinning** – keep one SNP per window (here one common SNP per {THIN_KB} kb) so the kept SNPs are nearly independent.

### Question
*Do the genotypes show ancestry groups, and are cases more common in one of them?*

### Test
""")
nb.code("""
# Step B3.2 · load the genotype subset: people x SNPs, 0/1/2 copies of A1
z = np.load(DATA / "gwas_genotypes_subset.npz")
G = z["genotypes"]                                             # int8 array, rows = people, columns = SNPs
snps = pd.read_csv(DATA / "gwas_snps_subset.tsv", sep="\\t")    # one row per column of G
assert (z["iid"] == samples["IID"].to_numpy()).all() and (z["snp"] == snps["SNP"].to_numpy()).all()
print("genotype matrix:", G.shape, "| PCA SNPs:", int(snps["pca_snp"].sum()),
      "| SNPs within ±500 kb of the CFH lead:", int(snps["cfh_window"].sum()))
print("first person, first 12 SNPs:", G[0, :12])
""")
nb.code("""
# Step B3.2 · standardise each SNP, then PCA with 10 components
def standardise(g):
    \"\"\"people x SNPs genotype matrix -> z-scores per SNP (mean 0, SD 1).\"\"\"
    g = g.astype(np.float64)
    return (g - g.mean(axis=0)) / g.std(axis=0)

G_pca = G[:, snps["pca_snp"].to_numpy() == 1]                 # the LD-thinned genome-wide SNPs
X = standardise(G_pca)
pca = PCA(n_components=10, svd_solver="randomized", random_state=0).fit(X)
pcs = pca.transform(X)                                         # people x 10 PC scores
var_pct = 100 * pca.explained_variance_ratio_                  # % of the total variance per PC
print(f"{G_pca.shape[1]:,} SNPs | variance explained (%):", np.round(var_pct, 3))
left = pcs[:, 0] < 0                                           # the two clusters along PC1
for name, m in [("PC1 < 0", left), ("PC1 > 0", ~left)]:
    print(f"{name}: {m.sum():,} people, {100 * y[m].mean():.0f}% cases, "
          f"{100 * (samples['site'].to_numpy()[m] == 'S2').mean():.0f}% recruited at site S2")
""")
nb.md("""
### Check
""")
nb.code("""
# CHECK · PC1 separates the two hidden ancestry groups (truth labels: allowed here) and stands far above PC2
anc = pd.read_csv(DATA / "truth" / "gwas_ancestry.tsv", sep="\\t")
ref = json.loads((DATA / "truth" / "gwas_truth.json").read_text())["realized"]   # reference values
assert (anc["IID"] == samples["IID"]).all()
isB = (anc["ancestry_group"] == "B").to_numpy()
agreement = max(np.mean(left == isB), np.mean(left != isB))     # the sign of a PC is arbitrary
dviz.check("B3.n_pca_snps", G_pca.shape[1], ref["n_pca_snps"])
dviz.check("B3.pc1_ancestry_agreement", float(agreement), between=(0.99, 1.0))
dviz.check("B3.pc1_var_pct", float(var_pct[0]), ref["pc_var_pct"][0], tol=0.001)
dviz.check("B3.pc1_over_pc2", float(var_pct[0] / var_pct[1]), between=(2, 100))
""")
nb.md("""
### Visualize
""")
nb.code("""
# Step B3.2 · PC1 vs PC2: colour = case status, marker = recruitment site
fig, ax = plt.subplots(figsize=(4.8, 3.9))
for case, col, lab in [(0, dviz.SLATE, "control"), (1, dviz.BLUE, "case")]:
    for site, mk in [("S1", "o"), ("S2", "^")]:
        m = (y == case) & (samples["site"] == site).to_numpy()
        ax.scatter(pcs[m, 0], pcs[m, 1], s=7, marker=mk, color=col, alpha=0.5, lw=0, label=f"{lab}, site {site}")
for m in (left, ~left):                                         # case fraction of each cluster, written above it
    ax.text(np.median(pcs[m, 0]), pcs[:, 1].max() * 1.12, f"{100 * y[m].mean():.0f}% cases\\n(n = {m.sum():,})",
            ha="center", va="bottom", fontsize=8, color=dviz.INK)
ax.set_ylim(pcs[:, 1].min() * 1.1, pcs[:, 1].max() * 1.55)
ax.set_xlabel(f"PC1 ({var_pct[0]:.2f}% of variance)")
ax.set_ylabel(f"PC2 ({var_pct[1]:.2f}% of variance)")
ax.set_title(f"Genotype PCA · {G_pca.shape[1]:,} LD-thinned SNPs · n = {len(y):,}")
leg = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.42), ncol=2, markerscale=2, handletextpad=0.3)
for h in leg.legend_handles:
    h.set_alpha(1)                                              # legend symbols without transparency
dviz.savefig(fig, "B3_pca")
plt.show()
""")
nb.md(f"""
**What changed?** 3,000 people × {n('n_pca_snps', '{:,}')} genotypes → 10 numbers per person. PC1 ({n('pc1_var_pct', '{:.2f}')}% of
the variance) splits the people into two clusters of {n('cluster_sizes')} that match the hidden ancestry groups
(agreement {n('pc1_agreement', '{:.3f}')}); PC2 to PC10 each explain only ≈ {n('pc2_var_pct', '{:.2f}')}% (noise level). Cases are
{n('cluster_case_pct')} of the two clusters: **case status and ancestry are mixed up**. Most people from site S2 sit in the
cluster with more cases, but not all of them – site is an imperfect proxy for ancestry. A small % of variance is normal for PC1 in genotype data: most
genetic variation is between individuals, not between groups (within Europe, PC1 explains well under 1% of the variance;
Novembre et al., 2008).

**Your turn.** Redo the PCA with only the first 1,000 PCA SNPs. *Predict:* will PC1 still separate the two groups as
cleanly? Then run the cell below.

<details><summary>Answer</summary>

Agreement drops from {n('pc1_agreement', '{:.3f}')} to {n('yt_agreement_1000', '{:.3f}')}. Each SNP carries only a tiny bit of ancestry
information (the groups differ by Fst ≈ {n('fst', '{:.4f}')}), so PCA needs thousands of SNPs to see the groups clearly. This is
why GWAS compute PCs from about 10,000 or more LD-thinned SNPs.
</details>
""")
nb.code("""
# Your turn · B3.2 · PCA with only 1,000 SNPs: do the two groups still separate?
pcs_1k = PCA(n_components=2, svd_solver="randomized", random_state=0).fit_transform(standardise(G_pca[:, :1000]))
left_1k = pcs_1k[:, 0] < 0
agreement_1k = float(max(np.mean(left_1k == isB), np.mean(left_1k != isB)))
dviz.check("B3.yt_agreement_1000_snps", agreement_1k, note="Your turn")
""")

# ------------------------------------------------------------------------------------------------ B3.3
nb.md("""
## Step B3.3 · Trend test with and without ancestry PCs  *(5 min)*

**Why this step.** For every SNP we ask whether the number of A1 alleles goes up (or down) with case status. The
**Armitage trend test** measures the correlation *r* between genotype (0/1/2) and case status (1/0); with no association
N·r² follows a chi-square distribution with 1 degree of freedom. The **EIGENSTRAT** method (Price et al., 2006) removes
ancestry first: it regresses the PCs out of the genotypes *and* out of case status, then computes the same correlation on
what is left. A SNP that only tracks ancestry has nothing left to correlate.

**New terms**
* **Residual** – what is left of a variable after subtracting the part that the PCs predict (a least-squares fit).
* **χ² statistic** – here (N − K − 1)·r², with K = number of PCs; large = strong association. The **p-value** is the
  chance of a χ² at least this large if the SNP had no effect.
* **λ (genomic inflation)** – median χ² divided by 0.455, its expected median under no association. About 1 = well
  calibrated; above 1.1 suggests stratification (or relatedness); true polygenic signal also raises it slightly.
* **QQ plot** – sorted observed −log10 p against the values expected under no association; on the diagonal = as expected.

### Question
*Are the test statistics inflated when ancestry is ignored, and do 4 PCs fix it?*

### Test
""")
nb.code("""
# Step B3.3 · the trend test in EIGENSTRAT form: regress out the covariates, then chi2 = (N - K - 1) r^2
def trend_test(G, y, covariates=None):
    \"\"\"chi2 and p for every SNP (column of G). Covariates (e.g. PCs) are removed from G and y first.\"\"\"
    N = len(y)
    C = np.ones((N, 1)) if covariates is None else np.column_stack([np.ones(N), covariates])
    G = G.astype(np.float64)
    y_res = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]       # case status minus what the PCs predict
    G_res = G - C @ np.linalg.lstsq(C, G, rcond=None)[0]       # every SNP minus what the PCs predict
    r = (G_res.T @ y_res) / (np.linalg.norm(G_res, axis=0) * np.linalg.norm(y_res))
    chi2_stat = (N - C.shape[1]) * r ** 2                      # C has K PCs + 1 intercept column
    return chi2_stat, chi2.sf(chi2_stat, df=1)

PC4 = pcs[:, :4]
chi2_no, p_no = trend_test(G_pca, y)                           # K = 0: ancestry ignored
chi2_pc, p_pc = trend_test(G_pca, y, PC4)                      # K = 4 ancestry PCs
lam_no, lam_pc = dviz.genomic_lambda(p_no), dviz.genomic_lambda(p_pc)
print(f"{len(p_pc):,} thinned SNPs   λ without PCs = {lam_no:.3f}   λ with 4 PCs = {lam_pc:.3f}")
""")
nb.code("""
# Step B3.3 · two SNPs before and after removing the PCs
pca_ids = snps.loc[snps["pca_snp"] == 1, "SNP"].to_numpy()
j_anc = int(np.argmax(chi2_no - chi2_pc))          # the SNP that loses most when ancestry is removed
j_true = int(np.argmax(chi2_pc))                   # the SNP that is strongest after ancestry is removed
for what, j in [("largest drop with PCs", j_anc), ("strongest with PCs", j_true)]:
    r_pc1 = np.corrcoef(G_pca[:, j], pcs[:, 0])[0, 1]
    print(f"{what:22s} {pca_ids[j]:18s} chi2 {chi2_no[j]:5.1f} without PCs -> {chi2_pc[j]:5.1f} with 4 PCs"
          f"   (correlation with PC1: {r_pc1:+.2f})")
""")
nb.md("""
### Check
""")
nb.code("""
# CHECK · lambda without and with PCs, and your p-values equal those in the summary-statistics file
dviz.check("B3.lambda_subset_noPC", lam_no, between=(1.15, 1.5))
dviz.check("B3.lambda_subset_PC", lam_pc, between=(0.95, 1.1))
gwas = pd.read_csv(DATA / "gwas_sumstats.tsv.gz", sep="\\t")    # genome-wide results, used again in B3.4
file_p = gwas.set_index("SNP").loc[pca_ids]
diff = max(np.abs(np.log10(p_pc) - np.log10(file_p["P"])).max(),
           np.abs(np.log10(p_no) - np.log10(file_p["P_no_ancestry_PCs"])).max())
dviz.check("B3.max_abs_log10p_difference", float(diff), between=(0, 0.001), note="you vs the file")
""")
nb.md("""
### Visualize
""")
nb.code("""
# Step B3.3 · QQ plot of the thinned SNPs for both analyses
fig, ax = plt.subplots(figsize=(4.2, 3.8))
for p, col, lab in [(p_no, dviz.ORANGE, f"no ancestry PCs (λ = {lam_no:.2f})"),
                    (p_pc, dviz.BLUE, f"4 ancestry PCs (λ = {lam_pc:.2f})")]:
    e, o = dviz.qq_points(p)                                    # expected and observed -log10 p, sorted
    ax.scatter(e, o, s=5, color=col, lw=0, label=lab)
ax.plot([0, e.max()], [0, e.max()], color=dviz.INK, lw=0.8)     # identity line = no association
ax.set_xlabel("Expected −log10 p")
ax.set_ylabel("Observed −log10 p")
ax.set_title(f"QQ plot · trend test · {len(p_pc):,} thinned SNPs")
ax.legend(loc="lower right", markerscale=2)
dviz.savefig(fig, "B3_qq_subset")
plt.show()
""")
nb.md(f"""
**What changed?** Each SNP now has two p-values. Without PCs the whole QQ curve lies above the diagonal and
λ = {n('lambda_subset_noPC', '{:.2f}')}: almost every SNP looks slightly associated – the signature of stratification, not of
biology. With 4 PCs the curve follows the diagonal and λ = {n('lambda_subset_PC', '{:.2f}')}. Two SNPs show what the PCs do:
{n('anc_snp')} falls from χ² = {n('anc_chi2_no', '{:.1f}')} to {n('anc_chi2_pc', '{:.1f}')} – its genotype simply tracks PC1 (correlation
{n('anc_r_pc1', '{:+.2f}')}) – while {n('true_snp')}, {n('true_where')}, keeps its signal ({n('true_chi2_no', '{:.1f}')} → {n('true_chi2_pc', '{:.1f}')}).
**PCs remove ancestry, not biology.** Your p-values match the file to within {n('max_abs_log10p_difference', '{:.1e}')} on the −log10
scale, so the same test made the genome-wide summary statistics.

**Your turn.** Adjust for the recruitment **site** (S1/S2) instead of the genetic PCs. *Predict:* will λ be close to
{n('lambda_subset_noPC', '{:.2f}')}, close to {n('lambda_subset_PC', '{:.2f}')}, or in between? Then run the cell below.

<details><summary>Answer</summary>

λ = {n('yt_lambda_site', '{:.2f}')} – in between. Site is correlated with ancestry (most group-B people were recruited at S2) but not
perfectly, so part of the confounding remains. Genetic PCs measure ancestry directly; site, country or self-reported
labels are only proxies.
</details>
""")
nb.code("""
# Your turn · B3.3 · adjust for recruitment site instead of the genetic PCs
site_S2 = (samples["site"] == "S2").to_numpy(dtype=float)       # 1 = site S2, 0 = site S1
chi2_site, p_site = trend_test(G_pca, y, site_S2)
lam_site = dviz.genomic_lambda(p_site)
dviz.check("B3.yt_lambda_site", lam_site, note="Your turn: site instead of PCs")
""")

# ------------------------------------------------------------------------------------------------ B3.4
nb.md("""
## Step B3.4 · Genome-wide Manhattan plot  *(5 min)*

**Why this step.** The summary statistics contain the same test (with 4 PCs) for all SNPs. A **Manhattan plot** places
every SNP at its position along the genome (x) and its strength of association, −log10 p, on y. **Genome-wide
significance**, p < 5 × 10⁻⁸, is a Bonferroni correction for about one million independent common-variant tests
(0.05 / 10⁶). Neighbouring SNPs in LD share the signal, so a true locus looks like a tower. We merge significant SNPs
less than 1 Mb apart into one **locus** and label it with the nearest gene.

**New terms**
* **Summary statistics** – one row per SNP: position (hg38, 1-based), alleles, effect size, its SE and p.
* **BETA / odds ratio (OR)** – BETA is the (approximate) log odds ratio per copy of A1; OR = exp(BETA). OR > 1 raises the
  odds of disease, OR < 1 lowers them.
* **Locus / lead SNP** – a region with significant SNPs / the SNP with the smallest p in it.
* **Nearest gene** – the gene whose body is closest to the lead SNP (0 bp if inside). A label, not proof of causality.
* **hg38 vs BED** – the GWAS table uses 1-based positions and `CHR` = 1…22; the gene table uses `chr1` and 0-based
  starts (BED). Convert before you compare.

### Question
*Which regions pass genome-wide significance, and are they the planted ones?*

### Test
""")
nb.code("""
# Step B3.4 · genome-wide lambda for both p-value columns, and the SNPs with P < 5e-8
print(gwas.head(3).to_string(index=False))
lam_gw_no = dviz.genomic_lambda(gwas["P_no_ancestry_PCs"])
lam_gw_pc = dviz.genomic_lambda(gwas["P"])
sig = gwas[gwas["P"] < 5e-8].sort_values(["CHR", "BP"])
print(f"{len(gwas):,} SNPs | λ without PCs {lam_gw_no:.3f} | λ with 4 PCs {lam_gw_pc:.3f} | {len(sig)} SNPs with P < 5e-8")
""")
nb.code("""
# Step B3.4 · merge significant SNPs closer than 1 Mb into loci; label each lead SNP with its nearest gene
genes = pd.read_csv(DATA / "gwas_genes_hg38.tsv", sep="\\t")          # MANE Select genes: chrom, start0, end, gene
genes["CHR"] = genes["chrom"].str.replace("chr", "").astype(int)     # 'chr1' -> 1, as in the GWAS table

def nearest_gene(chrom, bp):
    g = genes[genes["CHR"] == chrom]
    first, last = g["start0"].to_numpy() + 1, g["end"].to_numpy()    # 1-based first and last base of each gene
    dist = np.where(bp < first, first - bp, np.where(bp > last, bp - last, 0))   # 0 = inside the gene
    best = np.lexsort((np.abs((first + last) / 2 - bp), dist))[0]      # inside two genes: closest midpoint
    return g["gene"].iloc[best]

def define_loci(hits, pcol="P", merge_bp=1_000_000):
    new = (hits["CHR"].diff() != 0) | (hits["BP"].diff() > merge_bp)  # a new locus after a gap > 1 Mb
    hits = hits.assign(locus=np.cumsum(new.to_numpy()))
    lead = hits.loc[hits.groupby("locus")[pcol].idxmin()].reset_index(drop=True)
    lead["n_snps"] = hits.groupby("locus").size().to_numpy()
    lead["OR"] = np.exp(lead["BETA"]).round(2)
    lead["nearest_gene"] = [nearest_gene(c, b) for c, b in zip(lead["CHR"], lead["BP"])]
    return lead[["CHR", "BP", "SNP", pcol, "OR", "n_snps", "nearest_gene"]].copy()

loci = define_loci(sig)
loci
""")
nb.md("""
### Check
""")
nb.code("""
# CHECK · every locus lies at a planted locus (truth file: allowed here) and the strong planted loci are found
truth = pd.read_csv(DATA / "truth" / "gwas_truth.tsv", sep="\\t")      # locus, gene, CHR, BP, SNP, A1, A1_FREQ, OR

def planted_near(chrom, bp, max_bp=500_000):
    t = truth[(truth["CHR"] == chrom) & ((truth["BP"] - bp).abs() <= max_bp)]
    return t["locus"].iloc[0] if len(t) else "none (false locus)"

loci["planted_locus"] = [planted_near(c, b) for c, b in zip(loci["CHR"], loci["BP"])]
strong = truth.loc[(truth["OR"] >= 1.4) | (truth["OR"] <= 1 / 1.4), "locus"]
same_side = loci.merge(truth, left_on="planted_locus", right_on="locus", suffixes=("", "_true"))
dviz.check("B3.n_snps", len(gwas), ref["n_snps"])
dviz.check("B3.lambda_noPC", lam_gw_no, between=(1.2, 1.4))
dviz.check("B3.lambda_PC", lam_gw_pc, between=(1.0, 1.08))
dviz.check("B3.n_snps_p_lt_5e-8", len(sig), ref["n_snps_p_lt_5e-8"])
dviz.check("B3.n_loci", len(loci), between=(4, 6))
dviz.check("B3.false_loci", int((loci["planted_locus"] == "none (false locus)").sum()), 0)
dviz.check("B3.strong_loci_found", int(strong.isin(loci["planted_locus"]).sum()), len(strong))
dviz.check("B3.or_direction_agrees", bool(((same_side["OR"] > 1) == (same_side["OR_true"] > 1)).all()), True)
loci
""")
nb.md("""
### Visualize
""")
nb.code("""
# Step B3.4 · Manhattan plot: two alternating blues, significant SNPs in vermilion, points rasterised
x, ticks = dviz.manhattan_x(gwas["CHR"], gwas["BP"])            # cumulative position, gaps between chromosomes
mlp = -np.log10(gwas["P"].to_numpy())
is_sig = gwas["P"].to_numpy() < 5e-8
colour = np.where(gwas["CHR"].to_numpy() % 2 == 1, dviz.DARKBLUE, dviz.SKY)
fig, ax = plt.subplots(figsize=(7.5, 3.0))
ax.scatter(x[~is_sig], mlp[~is_sig], s=1.2, c=colour[~is_sig], lw=0, rasterized=True)
ax.scatter(x[is_sig], mlp[is_sig], s=3, color=dviz.VERMILION, lw=0, rasterized=True)
ax.axhline(-np.log10(5e-8), color=dviz.VERMILION, ls="--", lw=0.8)
ax.axhline(5, color=dviz.SLATE, ls=":", lw=0.6)                   # suggestive line, p = 1e-5
ax.text(x[0], -np.log10(5e-8) + 0.8, "p = 5 × 10⁻⁸", color=dviz.VERMILION, fontsize=7)
for r in loci.itertuples():                                        # label each locus with its nearest gene
    i = int(np.flatnonzero(gwas["SNP"].to_numpy() == r.SNP)[0])
    ha = "right" if x[i] > 0.97 * x.max() else "center"            # keep labels inside the axes
    ax.text(x[i], mlp[i] + 1.0, r.nearest_gene, ha=ha, va="bottom", fontsize=7, color=dviz.INK)
shown = ticks[(ticks["chrom"] <= 12) | (ticks["chrom"] % 2 == 0)]
ax.set_xticks(shown["centre"], shown["chrom"])
ax.set_xlim(x.min() - 1e7, x.max() + 1e7)
ax.set_ylim(0, mlp.max() * 1.12)
ax.set_xlabel("Chromosome (hg38)")
ax.set_ylabel("−log10 p")
ax.set_title(f"{len(gwas):,} SNPs · trend test with 4 ancestry PCs · {len(loci)} loci pass 5 × 10⁻⁸ (simulated)")
dviz.savefig(fig, "B3_manhattan")
plt.show()
""")
nb.md(f"""
**What changed?** {n('n_snps', '{:,}')} rows → {n('n_loci')} loci. λ is {n('lambda_noPC', '{:.2f}')} without PCs and {n('lambda_PC', '{:.2f}')} with 4 PCs
(the subset gave {n('lambda_subset_noPC', '{:.2f}')} and {n('lambda_subset_PC', '{:.2f}')}). {n('n_snps_p_lt_5e-8')} SNPs pass 5 × 10⁻⁸; merged within 1 Mb they form
{n('n_loci')} loci: {loci_text()}. Every locus lies at a planted locus and there is no false locus. The planted loci that were
missed ({missed_text()}) are too weak to reach 5 × 10⁻⁸ with 3,000 people (the real AMD GWAS used about 34,000). The tallest
towers contain up to {max_tower()} significant SNPs, because neighbouring SNPs share the signal through LD.

**Your turn.** *Predict:* how many loci would pass if you used the `P_no_ancestry_PCs` column instead, and would they all
be planted loci? Then run the cell below.

<details><summary>Answer</summary>

{n('yt_n_loci_noPC')} loci, {n('yt_false_loci_noPC')} of them false: {false_text()}. They are not in the truth file. They sit in
regions that mimic LCT/MCM6 (lactase persistence) and OCA2/HERC2 (eye colour): in real European data these carry alleles
whose frequency differs strongly between populations – the textbook stratification artefacts. With 4 PCs the same lead
SNPs have p = {aim_pc_text()}: nothing is left.
</details>
""")
nb.code("""
# Your turn · B3.4 · the same loci from the p-values that ignore ancestry
loci_no = define_loci(gwas[gwas["P_no_ancestry_PCs"] < 5e-8].sort_values(["CHR", "BP"]), pcol="P_no_ancestry_PCs")
loci_no["planted_locus"] = [planted_near(c, b) for c, b in zip(loci_no["CHR"], loci_no["BP"])]
dviz.check("B3.yt_n_loci_noPC", len(loci_no), note="Your turn")
dviz.check("B3.yt_false_loci_noPC", int((loci_no["planted_locus"] == "none (false locus)").sum()), note="Your turn")
loci_no
""")

# ------------------------------------------------------------------------------------------------ B3.5
nb.md("""
## Step B3.5 · Regional plot at CFH: who travels with the lead SNP?  *(3 min)*

**Why this step.** A Manhattan plot shows *where*; a **regional plot** zooms in on one locus (here ±500 kb around the CFH
lead SNP) and colours every SNP by its **r²** with the lead SNP. SNPs with high r² rise together with the lead because
they share its signal through LD – a tall tower does not mean many causal variants. A gene track below shows which genes
the tower overlaps.

**New terms**
* **r²** – squared correlation between the genotypes of two SNPs across people: 1 = they always travel together,
  0 = independent. We compute it from the 3,000 genotyped people.
* **LD block** – a stretch of SNPs with high r² among themselves, bounded by recombination hotspots.
* **Fine-mapping** – methods that narrow a locus down to the likely causal variant(s); not done here.

### Question
*Which SNPs make up the CFH tower, and why do they rise together?*

### Test
""")
nb.code("""
# Step B3.5 · r² of every SNP in the CFH window with the lead SNP (genotypes from the subset)
cols = np.flatnonzero(snps["cfh_window"].to_numpy() == 1)        # columns of G inside the window
win = snps.iloc[cols].merge(gwas[["SNP", "P", "BETA"]], on="SNP")  # same order as cols
lead = win.loc[win["P"].idxmin()]
g_lead = G[:, cols[win["P"].to_numpy().argmin()]].astype(float)
win["r2"] = [np.corrcoef(g_lead, G[:, c].astype(float))[0, 1] ** 2 for c in cols]
print(f"lead SNP {lead['SNP']}   p = {lead['P']:.2e}   OR = {np.exp(lead['BETA']):.2f}   {len(win)} SNPs within ±500 kb")
print(win.groupby(pd.cut(win["r2"], [0, 0.2, 0.4, 0.6, 0.8, 1.0001], right=False), observed=False).size().to_string())
""")
nb.md("""
### Check
""")
nb.code("""
# CHECK · the lead is the reference lead, r² with itself is 1, and the planted causal SNP is part of the tower
causal = truth[truth["locus"] == "CFH"].iloc[0]
r2_causal = float(win.loc[win["BP"] == causal["BP"], "r2"].iloc[0])
dviz.check("B3.cfh_lead_snp", lead["SNP"], ref["cfh_lead_snp"])
dviz.check("B3.cfh_window_snps", len(win), ref["n_cfh_window_snps"])
dviz.check("B3.cfh_lead_r2_with_itself", float(win.loc[win["SNP"] == lead["SNP"], "r2"].iloc[0]), 1.0, tol=1e-9)
dviz.check("B3.cfh_causal_r2_with_lead", r2_causal, between=(0.5, 1.0), note=f"causal SNP {causal['SNP']}")
dviz.check("B3.cfh_n_r2_ge_0.6", int((win["r2"] >= 0.6).sum()), note="SNPs (incl. the lead) that travel with the lead")
""")
nb.md("""
### Visualize
""")
nb.code("""
# Step B3.5 · regional plot: -log10 p vs position, colour = r² with the lead SNP (binned viridis), genes below
edges, labels = [0, 0.2, 0.4, 0.6, 0.8, 1.0001], ["< 0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "≥ 0.8"]
colours = plt.cm.viridis(np.linspace(0.0, 0.92, 5))
fig, (ax, axg) = plt.subplots(2, 1, figsize=(6.2, 4.4), sharex=True, gridspec_kw={"height_ratios": [3.2, 1], "hspace": 0.06})
others = win[win["SNP"] != lead["SNP"]]
for k in range(5):
    m = (others["r2"] >= edges[k]) & (others["r2"] < edges[k + 1])
    ax.scatter(others.loc[m, "BP"] / 1e6, -np.log10(others.loc[m, "P"]), s=22, color=colours[k],
               edgecolor=dviz.INK, lw=0.3, label=f"r² {labels[k]}")
ax.scatter(lead["BP"] / 1e6, -np.log10(lead["P"]), s=70, marker="D", color=dviz.PURPLE, edgecolor=dviz.INK, lw=0.6,
           label=f"lead {lead['SNP']}")
ax.axhline(-np.log10(5e-8), color=dviz.VERMILION, ls="--", lw=0.8)
ax.set_ylabel("−log10 p")
ax.set_title("CFH region · ±500 kb around the lead SNP · r² from 3,000 people")
ax.legend(loc="upper right", fontsize=7, handletextpad=0.2, borderaxespad=0.2)
g = genes[(genes["CHR"] == 1) & (genes["end"] > win["BP"].min()) & (genes["start0"] < win["BP"].max())].reset_index(drop=True)
for k, r in g.iterrows():                                            # gene bodies in three staggered rows
    a, b = max(r["start0"], win["BP"].min()) / 1e6, min(r["end"], win["BP"].max()) / 1e6
    axg.plot([a, b], [-(k % 3)] * 2, lw=5, color=dviz.SLATE, solid_capstyle="butt")
    axg.text((a + b) / 2, -(k % 3) + 0.3, r["gene"], ha="center", va="bottom", fontsize=7, color=dviz.INK)
axg.set_ylim(-2.6, 1.2); axg.set_yticks([]); axg.spines["left"].set_visible(False)
axg.set_xlabel("Position on chromosome 1 (Mb, hg38)")
dviz.savefig(fig, "B3_regional_cfh")
plt.show()
""")
nb.md(f"""
**What changed?** {n('cfh_window_snps')} SNPs around CFH → one picture of LD. The lead SNP ({n('cfh_lead_snp')}, p = {sci('cfh_lead_p')},
OR {n('cfh_lead_or', '{:.2f}')}) {lead_vs_causal()}. {n('cfh_n_r2_ge_0.6')} SNPs (including the lead) have r² ≥ 0.6 with it; they lie
{cfh_block()} and form the top of the tower. Beyond a recombination hotspot r² falls below 0.2 and so does the
association. The high-LD block overlaps {block_genes()}: the "nearest gene" label (CFH) is plausible, but LD alone cannot
say which of these genes is causal.

**Your turn.** *Predict:* is the SNP with the smallest p always the causal SNP? Run the cell below and compare the lead
with its best proxies.

<details><summary>Answer</summary>

Not necessarily. Here {lead_here()}, but the second-best SNP ({n('yt_second_snp')}, r² = {n('yt_second_r2', '{:.2f}')},
{n('yt_second_dist_kb', '{:.0f}')} kb away) has p = {sci('yt_second_p')} – almost as strong; with slightly different noise it could have been the lead.
In real data the lead SNP is often only a proxy of the causal variant, which is why fine-mapping methods (for example
SuSiE or FINEMAP) give each SNP a probability of being causal.
</details>
""")
nb.code("""
# Your turn · B3.5 · how close are the lead's best proxies?
top = win.sort_values("P").head(5)[["SNP", "BP", "MAF", "r2", "P"]].reset_index(drop=True)
top["kb_from_lead"] = ((top["BP"] - lead["BP"]) / 1e3).round(1)
dviz.check("B3.yt_second_best_r2", float(top["r2"].iloc[1]), note="Your turn")
top
""")

# ------------------------------------------------------------------------------------------------ summary
nb.md(f"""
## Summary

| Step | What you did | Number |
|---|---|---|
| B3.1 Question | counted the people | {n('n_people', '{:,}')} people: {n('n_cases', '{:,}')} cases, {n('n_controls', '{:,}')} controls |
| B3.2 Genotype PCA | standardised {n('n_pca_snps', '{:,}')} LD-thinned SNPs, PCA | PC1 {n('pc1_var_pct', '{:.2f}')}%, PC2 {n('pc2_var_pct', '{:.2f}')}%; PC1 = ancestry (agreement {n('pc1_agreement', '{:.3f}')}) |
| B3.3 Trend test | EIGENSTRAT trend test on the subset, K = 0 and 4 PCs | λ {n('lambda_subset_noPC', '{:.2f}')} → {n('lambda_subset_PC', '{:.2f}')} |
| B3.4 Manhattan | {n('n_snps', '{:,}')} SNPs, 5 × 10⁻⁸, loci merged within 1 Mb | λ {n('lambda_noPC', '{:.2f}')} → {n('lambda_PC', '{:.2f}')}; {n('n_snps_p_lt_5e-8')} SNPs in {n('n_loci')} loci, 0 false |
| B3.5 Regional plot | r² with the CFH lead SNP, ±500 kb | lead p = {sci('cfh_lead_p')}, OR {n('cfh_lead_or', '{:.2f}')}; {n('cfh_n_r2_ge_0.6')} SNPs with r² ≥ 0.6 |

**Presentation sentence (fill in, then say it in 30 seconds).**
*Question:* which genomic regions are associated with an AMD-like disease in 1,500 cases and 1,500 controls?
*Takeaway with a number:* after adjusting for 4 ancestry PCs (λ {n('lambda_noPC', '{:.2f}')} → {n('lambda_PC', '{:.2f}')}), {n('n_loci')} loci pass
p < 5 × 10⁻⁸ – the strongest at {strongest()} – all of them planted and none false; without PCs,
{n('yt_false_loci_noPC')} extra false loci appear.
*Limitation:* 3,000 people detect only strong effects (planted ORs of {missed_or()} were missed), and the nearest gene is a
label, not proof of the causal gene.
""")
nb.code("""
# Summary · collect the numbers behind the text and the slides (1/2): Steps B3.1-B3.3 and the Your-turn cells
numbers = {
    "n_people": len(y), "n_cases": int(y.sum()), "n_controls": int(len(y) - y.sum()), "n_snps": len(gwas),
    "pct_cases_S1": float(tab.loc["S1", "% cases"]), "pct_cases_S2": float(tab.loc["S2", "% cases"]),
    "n_pca_snps": int(G_pca.shape[1]), "pc1_var_pct": float(var_pct[0]), "pc2_var_pct": float(var_pct[1]),
    "pc1_agreement": float(agreement), "fst": float(ref["fst_hudson"]),
    "cluster_sizes": f"{left.sum():,} and {(~left).sum():,}",
    "cluster_case_pct": f"{100 * y[left].mean():.0f}% and {100 * y[~left].mean():.0f}%",
    "lambda_subset_noPC": lam_no, "lambda_subset_PC": lam_pc, "lambda_noPC": lam_gw_no, "lambda_PC": lam_gw_pc,
    "anc_snp": str(pca_ids[j_anc]), "anc_chi2_no": float(chi2_no[j_anc]), "anc_chi2_pc": float(chi2_pc[j_anc]),
    "anc_r_pc1": float(np.corrcoef(G_pca[:, j_anc], pcs[:, 0])[0, 1]), "true_snp": str(pca_ids[j_true]),
    "true_chi2_no": float(chi2_no[j_true]), "true_chi2_pc": float(chi2_pc[j_true]),
    "max_abs_log10p_difference": float(diff), "n_snps_p_lt_5e-8": len(sig), "n_loci": len(loci),
    "yt_toy_difference": float(f_cases - f_controls), "yt_agreement_1000": agreement_1k, "yt_lambda_site": lam_site,
    "yt_n_loci_noPC": len(loci_no), "yt_second_snp": str(top["SNP"].iloc[1]), "yt_second_r2": float(top["r2"].iloc[1]),
    "yt_second_dist_kb": float(abs(top["kb_from_lead"].iloc[1])), "yt_second_p": float(top["P"].iloc[1]),
}
""")
nb.code("""
# Summary · collect the numbers behind the text and the slides (2/2): loci (B3.4) and the CFH region (B3.5)
p_of = gwas.set_index("SNP")["P"]
near = lambda t: (gwas["CHR"] == t.CHR) & ((gwas["BP"] - t.BP).abs() <= 500_000)
false_no = loci_no[loci_no["planted_locus"] == "none (false locus)"]
high = win.loc[win["r2"] >= 0.6, "BP"]                               # positions of the high-LD block
numbers["loci"] = [{"gene": r.nearest_gene, "planted_locus": r.planted_locus, "CHR": int(r.CHR), "BP": int(r.BP),
                    "SNP": r.SNP, "p": float(r.P), "OR": float(r.OR), "n_snps": int(r.n_snps)} for r in loci.itertuples()]
numbers["missed"] = [{"locus": t.locus, "OR": float(t.OR), "best_p": float(gwas.loc[near(t), "P"].min())}
                     for t in truth[~truth["locus"].isin(loci["planted_locus"])].itertuples()]
numbers["false_noPC"] = [{"CHR": int(r.CHR), "BP": int(r.BP), "gene": r.nearest_gene, "p_noPC": float(r.P_no_ancestry_PCs),
                          "p_PC": float(p_of[r.SNP])} for r in false_no.itertuples()]
t_snp = snps.loc[snps["SNP"] == pca_ids[j_true]].iloc[0]          # is the strongest thinned SNP in a planted locus?
t_loc = planted_near(t_snp["CHR"], t_snp["BP"])
numbers.update({
    "true_where": f"which lies in the planted {t_loc} locus" if t_loc in set(truth["locus"]) else "outside the planted loci",
    "yt_false_loci_noPC": len(false_no), "cfh_window_snps": len(win), "cfh_lead_snp": str(lead["SNP"]),
    "cfh_lead_p": float(lead["P"]), "cfh_lead_or": float(np.exp(lead["BETA"])), "cfh_causal_r2": r2_causal,
    "cfh_lead_is_causal": bool(lead["BP"] == causal["BP"]), "cfh_n_r2_ge_0.6": int((win["r2"] >= 0.6).sum()),
    "cfh_block_kb": [float((high.min() - lead["BP"]) / 1e3), float((high.max() - lead["BP"]) / 1e3)],
    "cfh_block_genes": g.loc[(g["end"] >= high.min()) & (g["start0"] < high.max()), "gene"].tolist(),
})
""")
nb.code("""
# Summary · numbers for the slides (value, unit, meaning), then save (copied into the lecture deck and the instructor notes)
S = lambda value, unit, meaning: {"value": value, "unit": unit, "meaning": meaning}
slide_numbers = {
    "n_snps": S(len(gwas), "SNPs", "autosomal SNPs tested (hg38, MAF >= 0.01)"),
    "n_cases": S(int(y.sum()), "people", "cases"), "n_controls": S(int(len(y) - y.sum()), "people", "controls"),
    "lambda_noPC": S(round(lam_gw_no, 3), "", "genomic inflation λ, all SNPs, trend test without ancestry PCs"),
    "lambda_PC": S(round(lam_gw_pc, 3), "", "genomic inflation λ, all SNPs, trend test with 4 ancestry PCs"),
    "lambda_subset_noPC": S(round(lam_no, 3), "", f"λ on the {G_pca.shape[1]:,} LD-thinned SNPs, no PCs (notebook B3.3)"),
    "lambda_subset_PC": S(round(lam_pc, 3), "", f"λ on the {G_pca.shape[1]:,} LD-thinned SNPs, 4 PCs (notebook B3.3)"),
    "n_snps_p_lt_5e-8": S(len(sig), "SNPs", "SNPs with P < 5e-8 (4 PCs)"),
    "n_loci": S(len(loci), "loci", "significant SNPs merged within 1 Mb; all planted, 0 false"),
    "loci": S([{k: d[k] for k in ("gene", "planted_locus", "SNP", "p", "OR")} for d in numbers["loci"]], "p; OR per copy of A1",
              "lead SNP per locus, labelled with the nearest MANE Select gene"),
    "pc1_var_pct": S(round(float(var_pct[0]), 3), "% of variance", "genotype PCA, PC1 (= ancestry)"),
    "pc2_var_pct": S(round(float(var_pct[1]), 3), "% of variance", "genotype PCA, PC2 (noise level)"),
    "n_loci_noPC": S(len(loci_no), "loci", f"loci without ancestry PCs, {len(false_no)} of them false"),
}
dviz.save_results("B3", extra={"numbers": numbers, "slide_numbers": slide_numbers})
print(f"Notebook time: {time.time() - T_START:.0f} s")
""")

out = ROOT / "electives" / "B3_gwas.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
nb.write(str(out))
missing = sum(c.source.count(MISSING) for c in nb.cells if c.cell_type == "markdown")
print(f"wrote {out.relative_to(ROOT)}: {len(nb.cells)} cells; results/B3.json {'found' if R else 'not found'}; "
      f"placeholders in markdown: {missing}")
