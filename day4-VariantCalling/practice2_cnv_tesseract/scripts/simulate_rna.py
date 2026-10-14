#!/usr/bin/env python
"""Simulate a SYNTHETIC matched RNA-seq gene-count matrix (hg38, MANE Select genes, 8 samples).

This is NOT GTEx and NOT real data. The baseline is a synthetic 'fibroblast-like' profile:
  * every hg38 MANE Select protein-coding gene gets a baseline expression level
    (log10 CPM ~ Normal(1.0, 0.75), rescaled so that all genes sum to 1e6 CPM)
  * ~27% of genes are 'not expressed' (baseline 0.005-0.1 CPM, i.e. 0-3 reads)
  * chrY genes are only expressed in samples with sex M
  * a few named teaching genes have fixed baselines (PMP22, COMT, CRKL, TBX1, UGT2B17 ...) and a low
    biological noise (phi = 0.01) so that their expected dosage pattern is visible
Counts are drawn from a negative binomial (variance = mu + phi*mu^2, phi ~ 0.02*exp(N(0,0.4)); cultured cells vary little).
Library sizes differ between samples (18-32 M reads), so CPM normalisation matters.

Copy-number effect (patients and C2): every gene that overlaps a planted CNV (hg38 coordinates
from the breakpoint lift-over of the planted hg19 CNV) gets fold = (CN/2)^beta,
beta ~ U(0.85, 1.0) for ~85% of genes and U(0.2, 0.5) for ~15% ('buffered' genes);
CN 0 gives no expression. Genes only partly inside a DUPLICATION are left unchanged (an
incomplete extra copy does not make a full transcript); genes partly inside a DELETION are
treated as losing that allele.

Seed: numpy default_rng([20261015, 99]).  Output:
  data/rna/counts_hg38_mane.tsv.gz, data/rna/samples.tsv, data/truth/rna_truth_cnv_genes.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402
from fetch_public_data import lift_points  # noqa: E402

FIXED_CPM = {"PMP22": 400.0, "COMT": 150.0, "CRKL": 60.0, "TBX1": 20.0, "UGT2B17": 8.0,
             "HIRA": 25.0, "DGCR8": 15.0, "MAPK1": 120.0, "SNAP29": 30.0, "TXNRD2": 18.0,
             "GAPDH": 3000.0, "ACTB": 4000.0}
FIXED_BETA = {"COMT": 0.95, "CRKL": 0.92, "TBX1": 0.9, "PMP22": 0.95, "UGT2B17": 1.0}
NOT_EXPRESSED_FRACTION = 0.27
BUFFERED_FRACTION = 0.15
SEX_SPECIFIC = {"chrY": "M"}


def lifted_planted():
    rows = []
    for ev in C.PLANTED:
        a, b = lift_points([(ev["chrom"], ev["start_1based"] - 1), (ev["chrom"], ev["end_1based"] - 1)])
        rows.append(dict(cnv_id=ev["cnv_id"], sample=ev["sample"], chrom=a[0], start0=min(a[1], b[1]),
                         end=max(a[1], b[1]) + 1, cn=ev["cn"], type=ev["type"]))
    return pd.DataFrame(rows)


def main():
    rng = np.random.default_rng([C.SEED, 99])
    genes = pd.read_csv("data/annotation/genes_hg38_mane_select.tsv.gz", sep="\t")
    n = len(genes)
    log10cpm = rng.normal(1.0, 0.75, n)
    base = 10 ** log10cpm
    not_expr = rng.random(n) < NOT_EXPRESSED_FRACTION
    base[not_expr] = rng.uniform(0.005, 0.1, int(not_expr.sum()))
    sym = genes.gene.values
    for g, v in FIXED_CPM.items():
        base[sym == g] = v
        not_expr[sym == g] = False
    # rescale expressed genes so the profile sums to 1e6 CPM (fixed genes keep their values)
    fixed = np.isin(sym, list(FIXED_CPM))
    free = (~not_expr) & (~fixed)
    target = 1e6 - base[fixed].sum() - base[not_expr].sum()
    base[free] *= target / base[free].sum()
    phi = 0.02 * np.exp(rng.normal(0, 0.4, n))
    phi[fixed] = 0.01    # named teaching genes: lower biological noise so the expected dosage is visible
    beta = np.where(rng.random(n) < BUFFERED_FRACTION, rng.uniform(0.2, 0.5, n), rng.uniform(0.85, 1.0, n))
    for g, b in FIXED_BETA.items():
        beta[sym == g] = b

    cnvs = lifted_planted()
    samples = pd.DataFrame(C.SAMPLES)
    lib = {s: int(round(v)) for s, v in zip(samples["sample"], rng.uniform(18e6, 32e6, len(samples)))}
    fold = {s: np.ones(n) for s in samples["sample"]}
    truth = []
    for _, cv in cnvs.iterrows():
        on = (genes.chrom.values == cv.chrom) & (genes.start0.values < cv.end) & (genes.end.values > cv.start0)
        inside = on & (genes.start0.values >= cv.start0) & (genes.end.values <= cv.end)
        for i in np.where(on)[0]:
            partial = not inside[i]
            if cv.type == "DUP" and partial:
                f = 1.0
            else:
                f = 0.0 if cv.cn == 0 else (cv.cn / 2.0) ** beta[i]
            fold[cv["sample"]][i] *= f
            truth.append(dict(cnv_id=cv.cnv_id, sample=cv["sample"], gene=sym[i], chrom=cv.chrom,
                              gene_start0=int(genes.start0.values[i]), gene_end=int(genes.end.values[i]),
                              overlap="partial" if partial else "complete", cn=int(cv.cn),
                              baseline_cpm=round(float(base[i]), 3), expressed_baseline=bool(base[i] >= 1.0),
                              beta=round(float(beta[i]), 3), buffered=bool(beta[i] < 0.6),
                              expected_fold=round(float(f), 3),
                              expected_log2fc=(round(float(np.log2(f)), 3) if f > 0 else None)))
    counts = {}
    for _, srow in samples.iterrows():
        s = srow["sample"]
        mu = base * fold[s]
        for chrom, sex in SEX_SPECIFIC.items():
            if srow["sex"] != sex:
                mu[genes.chrom.values == chrom] = 0.0
        mu = mu / mu.sum() * lib[s] if mu.sum() > 0 else mu
        # negative binomial via gamma-Poisson mixture
        lam = rng.gamma(1.0 / phi, mu * phi)
        counts[s] = rng.poisson(lam)
    out = genes[["gene", "gene_id", "chrom", "start0", "end", "strand"]].copy()
    for s in samples["sample"]:
        out[s] = counts[s]
    os.makedirs("data/rna", exist_ok=True)
    out.to_csv("data/rna/counts_hg38_mane.tsv.gz", sep="\t", index=False)
    samples = samples.assign(tissue="fibroblast-like (SYNTHETIC baseline, not GTEx)",
                             library_reads=[int(out[s].sum()) for s in samples["sample"]],
                             data_type="synthetic")
    samples.to_csv("data/rna/samples.tsv", sep="\t", index=False)
    os.makedirs("data/truth", exist_ok=True)
    pd.DataFrame(truth).to_csv("data/truth/rna_truth_cnv_genes.tsv", sep="\t", index=False)
    print(f"genes={n} not_expressed={int(not_expr.sum())} ({not_expr.mean():.1%}) buffered={float((beta < 0.6).mean()):.1%}")
    print(samples[["sample", "library_reads"]].to_string(index=False))


if __name__ == "__main__":
    main()
