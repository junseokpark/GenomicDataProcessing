"""Simulated case-control GWAS of an age-related eye disease (AMD-like): the data for elective B3.

All data are computer-generated; no real people. Real anchors: hg38 chromosome lengths
(build/ref/hg38_chrom_sizes.tsv) and MANE Select gene coordinates (build/ref/genes_hg38_mane_select.tsv.gz).
The planted loci sit inside genes reported by the largest AMD GWAS (Fritsche et al., 2016, Nat Genet 48:134-143),
so the labels look familiar. Positions of all other SNPs, alleles, frequencies, LD and people are random.

Design (full data dictionary in docs/DATA_DICTIONARY.md)
--------------------------------------------------
People   3,000 = 1,800 of genetic-ancestry group A + 1,200 of group B. Cases: 756 of A (42 %) and 744 of B
         (62 %) -> 1,500 cases + 1,500 controls. Cases are over-represented in group B (recruitment), so a SNP
         whose frequency differs between A and B also differs between cases and controls: population
         stratification. Recruitment site follows ancestry imperfectly (A: 80 % S1, B: 75 % S2).
SNPs     300,000 autosomal positions (per chromosome proportional to its hg38 length), uniform between 1 Mb
         before the first and 1 Mb after the last MANE gene of the chromosome (so no SNPs on the unassembled short
         arms of chromosomes 13, 14, 15, 21, 22). Random ref/alt bases. QC: SNPs with sample MAF < 0.01 dropped.
Freqs    ancestral alt-allele frequency: MAF from 80 % Uniform(0.05, 0.5) + 20 % log-uniform(0.01, 0.05), alt
         allele minor or major at random; with probability Q_COPY a SNP copies (+- COPY_JITTER) the frequency of
         its left neighbour (SNPs that tag the same haplotype). Group frequencies: Balding-Nichols with F = FST.
LD       latent Gaussian AR(1) along each chromosome, one chain per haplotype:
         z_j = rho_j z_(j-1) + sqrt(1 - rho_j^2) e_j,  rho_j = exp(-d_j / L_BG - h_j),
         d_j = bp distance to the previous SNP, h_j = summed strength of the recombination hotspots between the
         two SNPs (Poisson, one per HS_SPACING bp; strength ~ Exponential(HS_MEAN)). Allele = alt if
         z_j < Phi^-1(p_group). Genotype = alt copies on the two haplotypes (0/1/2, A1 = alt).
Disease  liability = planted effects + polygenic background (N_BG random SNPs, tiny effects, total variance
         V_BG) + logistic noise; within each ancestry group the people with the highest liability are the cases
         (exact counts). Planted (conditional) effects are calibrated in a 1,000,000-person pilot so that a
         single-SNP logistic regression within ancestry groups returns the target OR.
Test     Armitage trend test in the EIGENSTRAT form: regress genotype and phenotype on [1, PC1..PCK], keep the
         residuals, chi2 = (N - K - 1) r^2, K = 0 or 4. PCs: randomized PCA (seed 0) of the standardised
         LD-thinned subset. BETA = slope of residual phenotype on residual genotype / (phi (1 - phi)),
         phi = case fraction (linear-to-logistic approximation of the log odds ratio); SE = score-test SE on the
         same scale, so BETA / SE = z and z^2 = chi2.
Run      python build/make_gwas.py      (about 1 minute, about 3 GB RAM; output is byte-identical between runs)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from scipy.special import ndtri
from sklearn.decomposition import PCA

T0 = time.time()
ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "build" / "ref"
DATA = ROOT / "data"
TRUTH = DATA / "truth"
SEED = 20261027 + 4
rng = np.random.default_rng(SEED)

# ------------------------------------------------------------------------------------------ design parameters
N_A, N_B = 1800, 1200                 # people per genetic-ancestry group
CASES_A, CASES_B = 756, 744           # 42 % and 62 % cases -> 1,500 cases in total
S2_A, S2_B = 360, 900                 # people recruited at site S2 (A: 20 %, B: 75 %); the rest at S1
N = N_A + N_B
M_TOTAL = 312_000                     # SNP positions before QC (about 300,000 remain after the MAF filter)
MAF_MIN = 0.01                        # QC filter
FST = 0.0028                          # Balding-Nichols F (see notes: Fst 0.01 would give lambda ~2.1 here)
L_BG, HS_SPACING, HS_MEAN = 3.0e6, 50e3, 0.3
Q_COPY, COPY_JITTER = 0.5, 0.02
N_BLOCK = 12                          # SNPs on each side of a planted / ancestry SNP that tag the same haplotype
N_BG, V_BG = 5000, 2.0                # polygenic background: SNPs and total liability variance (logit scale)
N_PCS = 4
AGE_CASE, AGE_CONTROL, AGE_SD, AGE_RANGE = 76, 71, 7, (55, 95)   # years; age depends on case status only
N_THIN = 10_000                       # target size of the LD-thinned PCA subset
THIN_MAF = 0.05
CFH_WINDOW = 500_000                  # subset keeps every SNP within +-500 kb of the CFH lead SNP
P_GWS, MERGE_BP = 5e-8, 1_000_000
PLANTED = [  # locus label, MANE gene that receives the causal SNP (gene midpoint), A1 (= alt) frequency, target OR
    ("CFH", "CFH", 0.38, 2.40),
    ("ARMS2/HTRA1", "HTRA1", 0.22, 2.60),
    ("C2/CFB", "CFB", 0.10, 0.55),
    ("C3", "C3", 0.28, 1.45),
    ("APOE", "APOE", 0.14, 0.75),
    ("LIPC", "LIPC", 0.30, 0.85),
    ("CETP", "CETP", 0.32, 1.20),
    ("TIMP3", "TIMP3", 0.35, 1.30),
]
# Ancestry-informative regions WITHOUT any disease effect: alt frequency in group A and group B. They mimic
# loci such as LCT (lactase persistence) and HERC2 (eye colour) whose frequencies differ strongly between
# European populations. They pass 5e-8 only when ancestry is ignored (a stratification artefact).
AIMS = [("LCT", "LCT", 0.28, 0.74), ("HERC2", "HERC2", 0.30, 0.76)]


def log(msg):
    print(f"[{time.time() - T0:5.1f}s] {msg}", flush=True)


# ------------------------------------------------------------------------------------------ reference tables
CHROMS = [f"chr{i}" for i in range(1, 23)]
sizes = pd.read_csv(REF / "hg38_chrom_sizes.tsv", sep="\t").set_index("chrom").loc[CHROMS, "length"]
LEN = sizes.to_numpy(dtype=np.int64)
genes = pd.read_csv(REF / "genes_hg38_mane_select.tsv.gz", sep="\t")
genes = genes[genes["chrom"].isin(CHROMS)].copy()
genes["CHR"] = genes["chrom"].str[3:].astype(int)
genes = genes.sort_values(["CHR", "start0", "end", "gene"]).reset_index(drop=True)
SPAN = []                              # SNP-bearing span per chromosome (1-based, inclusive)
for c in range(1, 23):
    g = genes[genes["CHR"] == c]
    SPAN.append((max(1, int(g["start0"].min()) - 1_000_000), min(int(LEN[c - 1]), int(g["end"].max()) + 1_000_000)))

# ------------------------------------------------------------------------------------------ SNP positions
raw = M_TOTAL * LEN / LEN.sum()
n_per = np.floor(raw).astype(int)
n_per[np.argsort(-(raw - n_per), kind="stable")[: M_TOTAL - n_per.sum()]] += 1
pos_by_chr = []
for c in range(1, 23):
    lo, hi = SPAN[c - 1]
    pos_by_chr.append(np.sort(rng.choice(hi - lo + 1, size=int(n_per[c - 1]), replace=False) + lo))

planted = pd.DataFrame(PLANTED, columns=["locus", "gene", "A1_FREQ_DESIGN", "OR"])
gi = genes.drop_duplicates("gene").set_index("gene")
planted["CHR"] = [int(gi.loc[g, "CHR"]) for g in planted["gene"]]
planted["BP"] = [int(gi.loc[g, "start0"]) + 1 + (int(gi.loc[g, "end"]) - int(gi.loc[g, "start0"])) // 2
                 for g in planted["gene"]]
aims = pd.DataFrame(AIMS, columns=["locus", "gene", "A1_FREQ_A", "A1_FREQ_B"])
aims["CHR"] = [int(gi.loc[g, "CHR"]) for g in aims["gene"]]
aims["BP"] = [int(gi.loc[g, "start0"]) + 1 + (int(gi.loc[g, "end"]) - int(gi.loc[g, "start0"])) // 2
              for g in aims["gene"]]
for c, bp in zip(list(planted["CHR"]) + list(aims["CHR"]), list(planted["BP"]) + list(aims["BP"])):
    p = pos_by_chr[c - 1]                                  # the special SNP replaces the nearest random SNP
    if bp not in set(p.tolist()):
        p[np.argmin(np.abs(p - bp))] = bp
        pos_by_chr[c - 1] = np.sort(p)

CHR = np.concatenate([np.full(len(p), c, dtype=np.int64) for c, p in zip(range(1, 23), pos_by_chr)])
BP = np.concatenate(pos_by_chr).astype(np.int64)
M = BP.size
chrom_start = np.r_[0, np.cumsum([len(p) for p in pos_by_chr])]
key = pd.Series(np.arange(M), index=pd.MultiIndex.from_arrays([CHR, BP]))
planted_idx = np.array([key.loc[(c, b)] for c, b in zip(planted["CHR"], planted["BP"])])
aim_idx = np.array([key.loc[(c, b)] for c, b in zip(aims["CHR"], aims["BP"])])
BASES = np.array(list("ACGT"))
ref_i = rng.integers(0, 4, M)
alt_i = (ref_i + rng.integers(1, 4, M)) % 4
log(f"{M:,} SNP positions on 22 autosomes; planted SNP rows {planted_idx.tolist()}")

# ------------------------------------------------------------------------------------------ allele frequencies
maf0 = np.where(rng.random(M) < 0.8, rng.uniform(0.05, 0.5, M), np.exp(rng.uniform(np.log(0.01), np.log(0.05), M)))
p_fresh = np.where(rng.random(M) < 0.5, maf0, 1 - maf0)
copy = rng.random(M) < Q_COPY
jit = rng.normal(0.0, COPY_JITTER, M)
special = np.r_[planted_idx, aim_idx]
special_freq = np.r_[planted["A1_FREQ_DESIGN"].to_numpy(), (aims["A1_FREQ_A"] + aims["A1_FREQ_B"]).to_numpy() / 2]
chr_of = np.searchsorted(chrom_start, np.arange(M), side="right") - 1
for j0 in special:                                      # a block of SNPs that tag the same haplotype
    a, b = max(chrom_start[chr_of[j0]] + 1, j0 - N_BLOCK + 1), min(chrom_start[chr_of[j0] + 1] - 1, j0 + N_BLOCK)
    copy[a:b + 1] = True
copy[chrom_start[:-1]] = False                          # the first SNP of a chromosome never copies
run_start = np.flatnonzero(~copy)
run = np.cumsum(~copy) - 1
cj = np.cumsum(np.where(copy, jit, 0.0))
p_anc = p_fresh[run_start][run] + cj - cj[run_start][run]
for j0, f0 in zip(special, special_freq):               # shift the block so the special SNP has its design frequency
    p_anc[run == run[j0]] += f0 - p_anc[j0]
p_anc = np.clip(p_anc, 0.005, 0.995)
a_bn, b_bn = p_anc * (1 - FST) / FST, (1 - p_anc) * (1 - FST) / FST
pA = np.clip(rng.beta(a_bn, b_bn), 1e-4, 1 - 1e-4)
pB = np.clip(rng.beta(a_bn, b_bn), 1e-4, 1 - 1e-4)
for j0, (fa0, fb0) in zip(aim_idx, zip(aims["A1_FREQ_A"], aims["A1_FREQ_B"])):
    blk = run == run[j0]                                # the whole tagging block differs between the groups
    pA[blk] = np.clip(p_anc[blk] - (fb0 - fa0) / 2, 0.01, 0.99)
    pB[blk] = np.clip(p_anc[blk] + (fb0 - fa0) / 2, 0.01, 0.99)
tA, tB = ndtri(pA).astype(np.float32), ndtri(pB).astype(np.float32)


# ------------------------------------------------------------------------------------------ haplotypes and LD
def simulate_chromosome(pos, tA_c, tB_c, lo, hi):
    """Latent-Gaussian AR(1) per haplotype; haplotype columns: [A hap1, A hap2, B hap1, B hap2]."""
    n_hs = rng.poisson((hi - lo) / HS_SPACING)
    hs_pos = np.sort(rng.uniform(lo, hi, n_hs))
    hs_sum = np.r_[0.0, np.cumsum(rng.exponential(HS_MEAN, n_hs))]
    k = np.searchsorted(hs_pos, pos)
    rho = np.exp(-np.diff(pos) / L_BG - (hs_sum[k[1:]] - hs_sum[k[:-1]])).astype(np.float32)
    s = np.sqrt(1 - rho.astype(np.float64) ** 2).astype(np.float32)
    m = pos.size
    haps = np.empty((m, 2 * N), dtype=np.bool_)
    z = rng.standard_normal(2 * N, dtype=np.float32)
    tmp = np.empty_like(z)
    a2 = 2 * N_A
    for b0 in range(0, m, 2048):
        eps = rng.standard_normal((min(2048, m - b0), 2 * N), dtype=np.float32)
        for jj in range(eps.shape[0]):
            j = b0 + jj
            if j > 0:
                z *= rho[j - 1]
                np.multiply(eps[jj], s[j - 1], out=tmp)
                z += tmp
            np.less(z[:a2], tA_c[j], out=haps[j, :a2])
            np.less(z[a2:], tB_c[j], out=haps[j, a2:])
    h8 = haps.view(np.int8)
    g = np.empty((m, N), dtype=np.int8)
    g[:, :N_A] = h8[:, :N_A] + h8[:, N_A:a2]
    g[:, N_A:] = h8[:, a2:a2 + N_B] + h8[:, a2 + N_B:]
    return g, rho


G = np.empty((M, N), dtype=np.int8)                     # SNP-major; people in internal order (A first, then B)
rho_all = np.zeros(M)                                   # correlation with the previous SNP (0 at chromosome starts)
for c in range(1, 23):
    sl = slice(chrom_start[c - 1], chrom_start[c])
    G[sl], rho_c = simulate_chromosome(BP[sl], tA[sl], tB[sl], *SPAN[c - 1])
    rho_all[chrom_start[c - 1] + 1:chrom_start[c]] = rho_c
log("haplotypes simulated")

# ------------------------------------------------------------------------------------------ QC (MAF >= 0.01)
f_alt = G.mean(axis=1, dtype=np.float64) / 2
maf = np.minimum(f_alt, 1 - f_alt)
keep = maf >= MAF_MIN
keep[special] = True
n_qc_removed = int((~keep).sum())
old_to_new = np.cumsum(keep) - 1
G, CHR, BP, ref_i, alt_i = G[keep], CHR[keep], BP[keep], ref_i[keep], alt_i[keep]
pA, pB, p_anc, f_alt, maf, rho_all = pA[keep], pB[keep], p_anc[keep], f_alt[keep], maf[keep], rho_all[keep]
planted_idx, aim_idx = old_to_new[planted_idx], old_to_new[aim_idx]
M = BP.size
log(f"QC removed {n_qc_removed:,} SNPs with MAF < {MAF_MIN}; {M:,} SNPs remain")

# ------------------------------------------------------------------------------------------ disease model
far = np.ones(M, dtype=bool)
for c, bp in zip(list(planted["CHR"]) + list(aims["CHR"]), list(planted["BP"]) + list(aims["BP"])):
    far &= ~((CHR == c) & (np.abs(BP - bp) < 2_000_000))
bg_idx = np.sort(rng.choice(np.flatnonzero(far & (maf >= 0.05)), N_BG, replace=False))
beta_bg = rng.normal(0.0, 1.0, N_BG)
var_g_bg = 2 * f_alt[bg_idx] * (1 - f_alt[bg_idx])
beta_bg *= np.sqrt(V_BG / np.sum(var_g_bg * beta_bg ** 2))     # total background variance exactly V_BG
grp = np.r_[np.zeros(N_A, dtype=np.int64), np.ones(N_B, dtype=np.int64)]   # 0 = group A, 1 = group B
case_frac = (CASES_A / N_A, CASES_B / N_B)


def top_k_cases(liability, group):
    y = np.zeros(liability.size, dtype=np.int64)
    for g_, f_ in enumerate(case_frac):
        idx = np.flatnonzero(group == g_)
        k = int(round(f_ * idx.size))
        y[idx[np.argsort(-liability[idx], kind="stable")[:k]]] = 1
    return y


def logit_slope(x, group, y):
    """Logistic regression y ~ 1 + group + x by Newton-Raphson; returns the coefficient of x."""
    X = np.column_stack([np.ones_like(x), group.astype(float), x])
    b = np.zeros(3)
    for _ in range(30):
        mu = 1 / (1 + np.exp(-(X @ b)))
        w = mu * (1 - mu)
        step = np.linalg.solve((X.T * w) @ X, X.T @ (y - mu))
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b[2]


# pilot: unlinked planted SNPs + Gaussian background with the same variance, same case fractions per group
target = np.log(planted["OR"].to_numpy())
n_pilot = 1_000_000
grp_p = np.r_[np.zeros(int(0.6 * n_pilot), dtype=np.int64), np.ones(n_pilot - int(0.6 * n_pilot), dtype=np.int64)]
fr = np.where(grp_p[:, None] == 0, pA[planted_idx][None, :], pB[planted_idx][None, :])
Gp = rng.binomial(2, fr).astype(np.float64)
noise_p = rng.normal(0.0, np.sqrt(V_BG), n_pilot) + rng.logistic(0.0, 1.0, n_pilot)
beta_cond = target.copy()
for _ in range(5):
    y_p = top_k_cases(Gp @ beta_cond + noise_p, grp_p)
    est = np.array([logit_slope(Gp[:, j], grp_p, y_p) for j in range(len(target))])
    beta_cond = beta_cond * target / est
pilot_or = np.exp(est)
del Gp, fr, noise_p
log("calibrated conditional OR: " + ", ".join(f"{l} {np.exp(b):.2f}" for l, b in zip(planted["locus"], beta_cond)))

liability = (G[planted_idx].T.astype(np.float64) @ beta_cond + G[bg_idx].T.astype(np.float64) @ beta_bg
             + rng.logistic(0.0, 1.0, N))
case_int = top_k_cases(liability, grp)

# ------------------------------------------------------------------------------------------ sample metadata
site_int = np.zeros(N, dtype=np.int64)                  # 0 = S1, 1 = S2
site_int[rng.choice(np.flatnonzero(grp == 0), S2_A, replace=False)] = 1
site_int[rng.choice(np.flatnonzero(grp == 1), S2_B, replace=False)] = 1
sex = np.where(rng.random(N) < 0.56, "F", "M")
age = np.where(case_int == 1, rng.normal(AGE_CASE, AGE_SD, N), rng.normal(AGE_CONTROL, AGE_SD, N))
age = np.clip(np.round(age), *AGE_RANGE).astype(int)

perm = rng.permutation(N)                               # output row r = internal person perm[r]
G = G[:, perm]
y = case_int[perm].astype(np.float64)
anc = np.where(grp[perm] == 0, "A", "B")
IID = np.array([f"ID{r + 1:04d}" for r in range(N)])
samples = pd.DataFrame({"IID": IID, "case": case_int[perm], "site": np.where(site_int[perm] == 1, "S2", "S1"),
                        "sex": sex[perm], "age": age[perm]})

# ------------------------------------------------------------------------------------------ LD-thinned subset + PCA
thin = []
total_span = sum(hi - lo + 1 for lo, hi in SPAN)
step = total_span / N_THIN
for c in range(1, 23):
    lo, hi = SPAN[c - 1]
    sl = np.arange(np.searchsorted(CHR, c), np.searchsorted(CHR, c, side="right"))
    ok = sl[maf[sl] >= THIN_MAF]
    for w0 in np.arange(lo, hi, step):
        centre = w0 + step / 2
        inwin = ok[(BP[ok] >= w0) & (BP[ok] < w0 + step)]
        if inwin.size:
            thin.append(int(inwin[np.argmin(np.abs(BP[inwin] - centre))]))
thin = np.array(sorted(set(thin)))


def standardise(g):
    """people x SNPs genotype matrix -> z-scores per SNP (mean 0, SD 1)."""
    g = g.astype(np.float64)
    return (g - g.mean(axis=0)) / g.std(axis=0)


X = standardise(np.ascontiguousarray(G[thin].T))        # people x SNPs, same layout as in the notebook
pca = PCA(n_components=10, svd_solver="randomized", random_state=0).fit(X)
scores = pca.transform(X)
var_pct = 100 * pca.explained_variance_ratio_
PCS = scores[:, :N_PCS]
log(f"PCA on {thin.size:,} thinned SNPs: PC1 {var_pct[0]:.3f} %, PC2 {var_pct[1]:.3f} %")
del X


# ------------------------------------------------------------------------------------------ association
def trend_test(g_snp_major, y, covariates):
    """EIGENSTRAT trend test for SNP-major genotypes (m x N). Returns chi2, slope and score-test SE."""
    C = np.column_stack([np.ones(y.size)] + ([covariates] if covariates is not None else []))
    Q, _ = np.linalg.qr(C)
    yr = y - Q @ (Q.T @ y)
    Xg = g_snp_major.T.astype(np.float64)
    Xr = Xg - Q @ (Q.T @ Xg)
    sxx = np.einsum("ij,ij->j", Xr, Xr)
    sxy = Xr.T @ yr
    syy = float(yr @ yr)
    df = y.size - C.shape[1]                      # N - K - 1
    r = sxy / np.sqrt(sxx * syy)
    return df * r ** 2, sxy / sxx, np.sqrt(syy / df / sxx)


chi_pc, chi_no = np.empty(M), np.empty(M)
slope, se_lin = np.empty(M), np.empty(M)
for a in range(0, M, 20_000):
    b = min(M, a + 20_000)
    chi_pc[a:b], slope[a:b], se_lin[a:b] = trend_test(G[a:b], y, PCS)
    chi_no[a:b], _, _ = trend_test(G[a:b], y, None)
phi = y.mean()
BETA = slope / (phi * (1 - phi))
SE = se_lin / (phi * (1 - phi))
P = chi2.sf(chi_pc, 1)
P_no = chi2.sf(chi_no, 1)
lam_pc = float(np.median(chi2.isf(P, 1)) / chi2.ppf(0.5, 1))
lam_no = float(np.median(chi2.isf(P_no, 1)) / chi2.ppf(0.5, 1))
log(f"lambda: no PCs {lam_no:.3f}, {N_PCS} PCs {lam_pc:.3f}")

SNP = np.array([f"{c}:{b}:{BASES[r]}:{BASES[t]}" for c, b, r, t in zip(CHR, BP, ref_i, alt_i)])
f_out = G.mean(axis=1, dtype=np.float64) / 2
MAF = np.minimum(f_out, 1 - f_out)
sumstats = pd.DataFrame({"SNP": SNP, "CHR": CHR, "BP": BP, "A1": BASES[alt_i], "A2": BASES[ref_i], "MAF": MAF,
                         "BETA": BETA, "SE": SE, "P": P, "P_no_ancestry_PCs": P_no})


def nearest_gene(c, bp):
    g = genes[genes["CHR"] == c]
    s1, e = g["start0"].to_numpy() + 1, g["end"].to_numpy()
    dist = np.where(bp < s1, s1 - bp, np.where(bp > e, bp - e, 0))
    order = np.lexsort((np.abs((s1 + e) / 2 - bp), dist))
    return str(g["gene"].iloc[order[0]]), int(dist[order[0]])


def define_loci(df, pcol):
    hits = df[df[pcol] < P_GWS].sort_values(["CHR", "BP"])
    rows, cur = [], []
    for r in hits.itertuples(index=False):
        if cur and (r.CHR != cur[-1].CHR or r.BP - cur[-1].BP > MERGE_BP):
            rows.append(cur)
            cur = []
        cur.append(r)
    if cur:
        rows.append(cur)
    out = []
    for grp_ in rows:
        lead = min(grp_, key=lambda r: getattr(r, pcol))
        gname, gdist = nearest_gene(lead.CHR, lead.BP)
        out.append({"CHR": int(lead.CHR), "start_bp": int(grp_[0].BP), "end_bp": int(grp_[-1].BP), "n_snps": len(grp_),
                    "lead_snp": lead.SNP, "lead_bp": int(lead.BP), "lead_p": float(getattr(lead, pcol)),
                    "lead_or": float(np.exp(lead.BETA)), "nearest_gene": gname, "gene_distance_bp": gdist})
    return out


def match_truth(loci):
    for L in loci:
        d = [(abs(L["lead_bp"] - b), loc) for c, b, loc in zip(planted["CHR"], planted["BP"], planted["locus"]) if c == L["CHR"]]
        best = min(d) if d else (None, None)
        L["planted_locus"] = best[1] if best[0] is not None and best[0] <= 500_000 else ""
        L["distance_to_causal_bp"] = int(best[0]) if best[0] is not None else None
    return loci


loci_pc = match_truth(define_loci(sumstats, "P"))
loci_no = match_truth(define_loci(sumstats, "P_no_ancestry_PCs"))
log(f"{int((P < P_GWS).sum())} SNPs < 5e-8 in {len(loci_pc)} loci (4 PCs); "
    f"{int((P_no < P_GWS).sum())} SNPs in {len(loci_no)} loci without PCs")

# ------------------------------------------------------------------------------------------ CFH window and subset
cfh_bp = int(planted.loc[planted["locus"] == "CFH", "BP"].iloc[0])
near = np.flatnonzero((CHR == 1) & (np.abs(BP - cfh_bp) <= 1_000_000))
cfh_lead = int(near[np.argmin(P[near])])
window = np.flatnonzero((CHR == 1) & (np.abs(BP - BP[cfh_lead]) <= CFH_WINDOW))
subset = np.union1d(thin, window)
pca_flag = np.isin(subset, thin).astype(int)
win_flag = np.isin(subset, window).astype(int)

# ------------------------------------------------------------------------------------------ write files
TRUTH.mkdir(parents=True, exist_ok=True)
gz = {"method": "gzip", "mtime": 0}
sumstats.to_csv(DATA / "gwas_sumstats.tsv.gz", sep="\t", index=False, float_format="%.6g", na_rep="", compression=gz)
np.savez_compressed(DATA / "gwas_genotypes_subset.npz", genotypes=np.ascontiguousarray(G[subset].T),
                    iid=IID, snp=SNP[subset])
snps_sub = sumstats.loc[subset, ["SNP", "CHR", "BP", "A1", "A2", "MAF"]].copy()
snps_sub["pca_snp"], snps_sub["cfh_window"] = pca_flag, win_flag
snps_sub.to_csv(DATA / "gwas_snps_subset.tsv", sep="\t", index=False, float_format="%.6g")
samples.to_csv(DATA / "gwas_samples.tsv", sep="\t", index=False)
genes[["chrom", "start0", "end", "gene", "strand"]].to_csv(DATA / "gwas_genes_hg38.tsv", sep="\t", index=False)

truth = planted[["locus", "gene", "CHR", "BP"]].copy()
truth["SNP"] = SNP[planted_idx]
truth["A1"] = BASES[alt_i[planted_idx]]
truth["A1_FREQ"] = np.round(f_out[planted_idx], 4)
truth["OR"] = planted["OR"].to_numpy()
truth.to_csv(TRUTH / "gwas_truth.tsv", sep="\t", index=False)
pd.DataFrame({"IID": IID, "ancestry_group": anc}).to_csv(TRUTH / "gwas_ancestry.tsv", sep="\t", index=False)

fa, fb = G[:, anc == "A"].mean(axis=1) / 2, G[:, anc == "B"].mean(axis=1) / 2   # sample group frequencies
na, nb = 2 * N_A, 2 * N_B
fst_hudson = float(np.mean((fa - fb) ** 2 - fa * (1 - fa) / (na - 1) - fb * (1 - fb) / (nb - 1))
                   / np.mean(fa * (1 - fb) + fb * (1 - fa)))
rnd = lambda v, k=6: float(f"{v:.{k}g}")
info = {
    "dataset": "gwas (elective B3)", "seed": SEED, "simulated": True, "genome_build": "hg38",
    "design": {"n_group_A": N_A, "n_group_B": N_B, "cases_group_A": CASES_A, "cases_group_B": CASES_B,
               "site_S2_group_A": S2_A, "site_S2_group_B": S2_B, "snp_positions_before_qc": M_TOTAL,
               "maf_min": MAF_MIN, "balding_nichols_F": FST, "ld_L_bg_bp": L_BG, "hotspot_spacing_bp": HS_SPACING,
               "hotspot_mean_strength": HS_MEAN, "freq_copy_prob": Q_COPY, "freq_copy_jitter": COPY_JITTER,
               "n_background_snps": N_BG, "background_variance_logit": V_BG, "n_pcs": N_PCS,
               "haplotype_block_snps_each_side": N_BLOCK, "pilot_n": n_pilot,
               "maf_spectrum": "80 % Uniform(0.05, 0.5) + 20 % log-uniform(0.01, 0.05)",
               "age_mean_case": AGE_CASE, "age_mean_control": AGE_CONTROL, "age_sd": AGE_SD, "age_range": list(AGE_RANGE),
               "ancestry_informative_design": {l: [fa0, fb0] for l, fa0, fb0 in zip(aims["locus"], aims["A1_FREQ_A"], aims["A1_FREQ_B"])},
               "thin_window_bp": rnd(step), "thin_maf_min": THIN_MAF, "cfh_window_bp": CFH_WINDOW,
               "planted_conditional_or": {l: rnd(np.exp(b), 4) for l, b in zip(planted["locus"], beta_cond)},
               "pilot_marginal_or": {l: rnd(o, 4) for l, o in zip(planted["locus"], pilot_or)}},
    "realized": {"n_snps": int(M), "n_snps_qc_removed": n_qc_removed, "n_cases": int(y.sum()),
                 "n_controls": int(N - y.sum()), "fst_hudson": rnd(fst_hudson, 4),
                 "lambda_noPC": rnd(lam_no, 5), "lambda_PC": rnd(lam_pc, 5),
                 "n_snps_p_lt_5e-8": int((P < P_GWS).sum()), "n_loci": len(loci_pc),
                 "n_snps_p_lt_5e-8_noPC": int((P_no < P_GWS).sum()), "n_loci_noPC": len(loci_no),
                 "pc_var_pct": [rnd(v, 4) for v in var_pct], "n_pca_snps": int(thin.size),
                 "n_subset_snps": int(subset.size), "n_cfh_window_snps": int(window.size),
                 "cfh_lead_snp": str(SNP[cfh_lead]), "cfh_lead_bp": int(BP[cfh_lead]), "cfh_lead_p": float(P[cfh_lead]),
                 "mean_chi2_PC": rnd(float(chi_pc.mean()), 5), "mean_chi2_noPC": rnd(float(chi_no.mean()), 5),
                 "ancestry_informative": [{"locus": l, "gene": g_, "CHR": int(CHR[j]), "BP": int(BP[j]), "SNP": str(SNP[j]),
                                           "A1_freq_A": rnd(fa[j], 4), "A1_freq_B": rnd(fb[j], 4),
                                           "P": float(P[j]), "P_no_ancestry_PCs": float(P_no[j])}
                                          for l, g_, j in zip(aims["locus"], aims["gene"], aim_idx)],
                 "loci_PC": loci_pc, "loci_noPC": loci_no,
                 "case_fraction_by_site": samples.groupby("site")["case"].mean().round(4).to_dict()},
}
(TRUTH / "gwas_truth.json").write_text(json.dumps(info, indent=1))
for f in ["gwas_sumstats.tsv.gz", "gwas_genotypes_subset.npz", "gwas_snps_subset.tsv", "gwas_samples.tsv",
          "gwas_genes_hg38.tsv", "truth/gwas_truth.tsv", "truth/gwas_ancestry.tsv", "truth/gwas_truth.json"]:
    log(f"wrote data/{f}  {(DATA / f).stat().st_size / 1e6:.2f} MB")
print(json.dumps({k: v for k, v in info["realized"].items() if not k.startswith("loci")}, indent=1))
for L in loci_pc:
    print(f"  locus chr{L['CHR']}:{L['start_bp']:,}-{L['end_bp']:,}  {L['n_snps']:3d} SNPs  lead {L['lead_snp']}  "
          f"p={L['lead_p']:.2e} OR={L['lead_or']:.2f}  nearest {L['nearest_gene']}  planted {L['planted_locus'] or '-'}")
print("without PCs:", [(L["CHR"], L["lead_bp"], f"{L['lead_p']:.1e}", L["planted_locus"] or "FALSE") for L in loci_no])
