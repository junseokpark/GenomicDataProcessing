"""B5 data: read depth on hg38 chr22 for one patient with a heterozygous 22q11.2 (LCR22A-D) deletion.

All depths are SIMULATED. Real anchors (hg38): chromosome lengths, UCSC genomicSuperDups (segmental
duplications), the ClinGen region ISCA-37446 (22q11.2 recurrent DGS region, A-D) and MANE Select genes.

Writes (paths relative to dataviz_practice/)
  data/cnv_chr22_10kb.tsv        10-kb bins over all of chr22: chrom start end gc mappable depth_patient depth_reference
  data/cnv_genome_1mb.tsv        1-Mb bins chr1-22: log2_patient (same patient) and log2_tumour (colorectal tumour)
  data/cnv_genes_chr22.tsv       MANE Select protein-coding genes overlapping chr22 17-23 Mb (real)
  data/cnv_segdups_chr22.tsv     UCSC segmental duplications on chr22 17.98-22.41 Mb (real) + LCR22 block name
  data/truth/cnv_truth.tsv       planted events (copy number, purity, expected log2)

Coordinates: intervals are BED style (0-based start, end exclusive), genome build hg38.
Seed: np.random.default_rng(20261027 + 6). Running the script twice gives byte-identical files.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "build" / "ref"
TRK = REF / "d4_tracks_hg38"
DATA = ROOT / "data"
rng = np.random.default_rng(20261027 + 6)

CHROM, BIN = "chr22", 10_000
P_ARM_END = 15_410_000         # bins before this: short arm, centromere, pericentromeric satellites (no unique sequence)
SHORT_ARM_GAP_END = 10_510_000  # hg38 chr22 short-arm N gap (depth 0, GC undefined)
HI_ID = 0.98                    # "high identity" segmental duplication (reads cannot be placed uniquely)
READ_LEN = 150
N_NORMALS = 10
PURITY = 0.6

# ------------------------------------------------------------------------------------------------
# 1. Real anchors
# ------------------------------------------------------------------------------------------------
sizes = pd.read_csv(REF / "hg38_chrom_sizes.tsv", sep="\t").set_index("chrom")["length"]
L22 = int(sizes[CHROM])                                                    # 50,818,468

sd_raw = pd.read_csv(TRK / "genomicSuperDups.tsv.gz", sep="\t")
sd_raw = sd_raw[sd_raw["chrom"] == CHROM].drop(columns=["region_id"]).drop_duplicates()  # 2 windows overlap
sd = pd.DataFrame({"chrom": sd_raw["chrom"], "start": sd_raw["chromStart"], "end": sd_raw["chromEnd"],
                   "other_chrom": sd_raw["otherChrom"], "other_start": sd_raw["otherStart"],
                   "other_end": sd_raw["otherEnd"], "strand": sd_raw["strand"],
                   "frac_match": sd_raw["fracMatch"].round(6)})
sd = sd.sort_values(["start", "end", "other_chrom", "other_start"]).reset_index(drop=True)

haplo = pd.read_csv(TRK / "clinGenHaplo.tsv.gz", sep="\t")
isca = haplo.loc[haplo["name"] == "ISCA-37446"].iloc[0]
ISCA = (int(isca["chromStart"]), int(isca["chromEnd"]))                   # 18,924,717-21,111,383

mane = pd.read_csv(REF / "genes_hg38_mane_select.tsv.gz", sep="\t")


def merge(iv, max_gap=0):
    """Merge sorted (start, end) intervals closer than max_gap."""
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1] + max_gap:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


# Assembly gaps inside LCR22A: holes between segdup clusters of exactly 50 kb or 100 kb (the hg38 gap sizes)
# that contain no feature of any annotation track (phastCons elements checked here).
phast = pd.read_csv(TRK / "phastConsElements100way.tsv.gz", sep="\t")
phast = phast[phast["chrom"] == CHROM]
clusters = merge(zip(sd["start"], sd["end"]))
ASM_GAPS = []
for (a0, a1), (b0, b1) in zip(clusters[:-1], clusters[1:]):
    hole = b0 - a1
    if min(abs(hole - 50_000), abs(hole - 100_000)) <= 10:
        n_feat = int(((phast["chromStart"] < b0) & (phast["chromEnd"] > a1)).sum())
        if n_feat == 0:
            ASM_GAPS.append((int(a1), int(b0)))

# LCR22 blocks: high-identity segdups merged across holes < 20 kb or across an assembly gap
hi_iv = merge(zip(sd.loc[sd["frac_match"] >= HI_ID, "start"], sd.loc[sd["frac_match"] >= HI_ID, "end"]))
blocks = []
for a, b in hi_iv:
    if blocks and (a - blocks[-1][1] < 20_000 or any(g0 <= blocks[-1][1] + 10 and g1 >= a - 10 for g0, g1 in ASM_GAPS)):
        blocks[-1][1] = b
    else:
        blocks.append([a, b])
NAMED = {"LCR22A": 18_700_000, "LCR22B": 20_350_000, "LCR22C": 20_700_000, "LCR22D": 21_300_000}
LCR = {}
for name, probe in NAMED.items():
    hit = [blk for blk in blocks if blk[0] <= probe < blk[1]]
    assert len(hit) == 1, name
    LCR[name] = (int(hit[0][0]), int(hit[0][1]))
sd["lcr22"] = ""
for name, (a, b) in LCR.items():
    sd.loc[(sd["start"] >= a) & (sd["end"] <= b), "lcr22"] = name[-1]          # A, B, C or D

# Planted deletion: non-allelic homologous recombination between two directly oriented copies of one
# 99.7 %-identical segment, one copy in LCR22A and one in LCR22D; crossover at the same relative position.
pair = sd[(sd["start"] == 18_709_564) & (sd["other_chrom"] == CHROM) & (sd["other_start"] == 21_335_350)].iloc[0]
assert pair["strand"] == "+" and pair["frac_match"] > 0.99
F_CROSS = 0.485
DEL = (int(pair["start"] + round(F_CROSS * (pair["end"] - pair["start"]))),
       int(pair["other_start"] + round(F_CROSS * (pair["other_end"] - pair["other_start"]))))
assert DEL[0] <= ISCA[0] and DEL[1] >= ISCA[1], "deletion must contain ISCA-37446"
assert LCR["LCR22A"][0] <= DEL[0] < LCR["LCR22A"][1] and LCR["LCR22D"][0] <= DEL[1] < LCR["LCR22D"][1]


def ar1(n, rho, sd_, r):
    """Stationary AR(1) series: neighbouring values are correlated (rho), overall SD = sd_."""
    e = r.normal(0, sd_ * np.sqrt(1 - rho ** 2), n)
    x = np.empty(n)
    x[0] = r.normal(0, sd_)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + e[i]
    return x


def overlap_fraction(starts, ends, intervals):
    """Fraction of each bin covered by the union of intervals."""
    cov = np.zeros(len(starts))
    for a, b in merge(intervals):
        i0, i1 = a // BIN, (b - 1) // BIN
        for i in range(i0, min(i1, len(starts) - 1) + 1):
            cov[i] += max(0, min(b, ends[i]) - max(a, starts[i]))
    return cov / (ends - starts)


# ------------------------------------------------------------------------------------------------
# 2. Bins, masks and GC
# ------------------------------------------------------------------------------------------------
start = np.arange(0, L22, BIN, dtype=np.int64)
end = np.minimum(start + BIN, L22)
n = len(start)
blen = end - start
segdup_cov = overlap_fraction(start, end, list(zip(sd.loc[sd["frac_match"] >= HI_ID, "start"],
                                                   sd.loc[sd["frac_match"] >= HI_ID, "end"])))
gap_cov = overlap_fraction(start, end, ASM_GAPS)
is_n = (start < SHORT_ARM_GAP_END) | (gap_cov >= 0.5) | (end == L22)      # reference is N (last bin = telomere)
is_sat = (start < P_ARM_END) & ~is_n                                       # satellites / centromere model
is_sd = (segdup_cov >= 0.5) & ~is_n & ~is_sat
mappable = ~(is_n | is_sat | is_sd)

gc_raw = 0.48 + ar1(n, 0.995, 0.055, rng) + ar1(n, 0.8, 0.035, rng)
gc = (0.50 + 0.15 * np.tanh((gc_raw - 0.50) / 0.15)).round(4)    # smooth squeeze into (0.35, 0.65), no pile-up
gc[is_n] = np.nan

# copy number of the patient per bin (bp-weighted, so partial bins get e.g. 1.4)
del_frac = np.clip((np.minimum(end, DEL[1]) - np.maximum(start, DEL[0])) / blen, 0, 1)
cn_patient = 2.0 - del_frac

# ------------------------------------------------------------------------------------------------
# 3. Read depth: patient and a panel of 10 normals
# ------------------------------------------------------------------------------------------------
base = np.ones(n)                                      # what the reference genome allows (shared by all samples)
base[is_sd] = rng.uniform(0.10, 0.60, is_sd.sum())      # MAPQ >= 20 filter removes most reads in near-identical repeats
sat = rng.lognormal(0.0, 1.0, is_sat.sum())
sat[rng.random(is_sat.sum()) < 0.15] = 0.0
base[is_sat] = sat
base[is_n] = 0.0
unreliable = is_sd | is_sat
GC_FILL = np.where(np.isnan(gc), 0.45, gc)


def simulate_sample(mean_depth, cn, r, outliers=False):
    """Mean depth (x) per bin: depth x GC bias x copy number x slow waves x over-dispersed counting noise."""
    centre, width = r.normal(0.40, 0.004), r.normal(0.20, 0.004)   # depth falls at high GC (shared shape)
    gc_bias = np.exp(-0.5 * ((GC_FILL - centre) / width) ** 2)
    wave = np.exp(ar1(n, 0.99, 0.025, r))
    cn_factor = np.where(is_sd, 0.5 + 0.25 * cn, cn / 2.0)  # paralogous reads dilute the signal in repeats
    mu = mean_depth * base * gc_bias * wave * cn_factor
    mu = mu * r.gamma(59.0, 1 / 59.0, n)                    # over-dispersion: CV 0.13 per bin
    extra = np.where(unreliable, r.lognormal(0.0, 0.4, n), 1.0)
    mu = mu * extra
    if outliers:
        hit = mappable & (r.random(n) < 0.003)
        mu = mu * np.where(hit, np.exp(r.normal(0, 0.7, n)), 1.0)
    reads = r.poisson(mu * blen / READ_LEN)
    return reads * READ_LEN / blen


patient_mean = rng.uniform(29.0, 33.0)
depth_patient = simulate_sample(patient_mean, cn_patient, rng, outliers=True)
normal_means = rng.uniform(24.0, 36.0, N_NORMALS)
normals = np.vstack([simulate_sample(m, np.full(n, 2.0), rng) for m in normal_means])
med = np.median(normals[:, mappable], axis=1)
depth_reference = np.median(normals / med[:, None], axis=0) * med.mean()
depth_patient[is_n] = 0.0
depth_reference[is_n] = 0.0

bins = pd.DataFrame({"chrom": CHROM, "start": start, "end": end, "gc": gc, "mappable": mappable,
                     "depth_patient": depth_patient.round(3), "depth_reference": depth_reference.round(3)})
assert (bins.loc[bins["mappable"], ["depth_patient", "depth_reference"]] > 0).all().all()

# ------------------------------------------------------------------------------------------------
# 4. Genome-wide 1-Mb file: the same patient and, for contrast, a colorectal tumour (purity 0.6)
# ------------------------------------------------------------------------------------------------
MB = 1_000_000
ok = bins["mappable"].to_numpy()
p_norm = bins["depth_patient"] / bins.loc[ok, "depth_patient"].median()
r_norm = bins["depth_reference"] / bins.loc[ok, "depth_reference"].median()
log2_10kb = np.log2((p_norm / r_norm).where(ok))
agg = pd.DataFrame({"mb": start // MB, "log2": log2_10kb}).groupby("mb")["log2"].agg(["mean", "count"])
chr22_1mb = agg["mean"].where(agg["count"] >= 10)          # need >= 10 mappable 10-kb bins per 1-Mb bin

# approximate hg38 centromere midpoints (Mb), used only to choose which gene-free gap is the centromere
CEN_MID = {1: 123.4, 2: 93.9, 3: 90.9, 4: 50.0, 5: 48.8, 6: 59.8, 7: 60.1, 8: 45.2, 9: 43.0, 10: 39.8, 11: 53.4,
           12: 35.5, 13: 17.7, 14: 17.2, 15: 19.0, 16: 36.8, 17: 25.1, 18: 18.5, 19: 26.2, 20: 28.1, 21: 12.0,
           22: 15.0}


def centromere_gap(chrom):
    """The gene-free stretch (MANE Select) that contains the approximate centromere midpoint."""
    gg = mane[mane["chrom"] == chrom].sort_values("start0")
    st = np.concatenate([[0], gg["start0"].to_numpy()])
    en = np.concatenate([[0], np.maximum.accumulate(gg["end"].to_numpy())])
    mid = CEN_MID[int(chrom[3:])] * MB
    for a, b in zip(en[:-1], st[1:]):
        if a <= mid < b:
            return int(a), int(b)
    raise ValueError(chrom)


EVENTS = {"chr7": {"p": 3, "q": 3}, "chr8": {"p": 1, "q": 3}, "chr13": {"p": 3, "q": 3},
          "chr17": {"p": 1, "q": 2}, "chr18": {"p": 2, "q": 1}, "chr20": {"p": 2, "q": 4}}
rows, truth_rows = [], []
for k in range(1, 23):
    chrom = f"chr{k}"
    length = int(sizes[chrom])
    s0 = np.arange(0, length, MB, dtype=np.int64)
    e0 = np.minimum(s0 + MB, length)
    g0, g1 = centromere_gap(chrom)
    in_gap = (s0 >= g0) & (e0 <= g1)
    arm = np.where((s0 + e0) / 2 < (g0 + g1) / 2, "p", "q")
    cn_t = np.array([EVENTS.get(chrom, {}).get(a, 2) for a in arm], dtype=float)
    tum = np.log2((PURITY * cn_t + 2 * (1 - PURITY)) / 2)
    tum = tum + ar1(len(s0), 0.6, 0.05, rng) + rng.normal(0, 0.04, len(s0))
    if chrom == CHROM:
        pat = chr22_1mb.reindex(range(len(s0))).to_numpy(dtype=float, copy=True)  # writable copy (pandas 3)
        missing = np.isnan(pat)
    else:
        pat = ar1(len(s0), 0.6, 0.035, rng) + rng.normal(0, 0.025, len(s0))
        missing = in_gap
    pat[missing] = np.nan
    tum[missing] = np.nan
    rows.append(pd.DataFrame({"chrom": chrom, "start": s0, "end": e0, "log2_patient": pat, "log2_tumour": tum}))
    for a, cn_arm in EVENTS.get(chrom, {}).items():
        if cn_arm != 2:
            boundary = int(round((g0 + g1) / 2, -5))
            a0, a1 = (0, boundary) if a == "p" else (boundary, length)
            if chrom in ("chr7", "chr13"):
                if a == "q":
                    continue
                a0, a1, a = 0, length, "whole"
            kind = "GAIN" if cn_arm > 2 else "LOSS"
            truth_rows.append({"sample": "tumour", "event": f"{chrom}{'' if a == 'whole' else a}_{kind.lower()}",
                               "chrom": chrom, "start": a0, "end": a1, "type": kind, "copy_number": cn_arm,
                               "purity": PURITY,
                               "expected_log2": round(float(np.log2((PURITY * cn_arm + 2 * (1 - PURITY)) / 2)), 4),
                               "note": "arm boundary = middle of the gene-free centromere gap (approximate)"})
genome = pd.concat(rows, ignore_index=True)
genome["log2_tumour"] = genome["log2_tumour"] - np.nanmedian(genome["log2_tumour"])   # centre the tumour on CN 2
genome[["log2_patient", "log2_tumour"]] = genome[["log2_patient", "log2_tumour"]].round(4)

truth_rows.insert(0, {"sample": "patient", "event": "22q11.2_LCR22A-D_deletion", "chrom": CHROM,
                      "start": DEL[0], "end": DEL[1], "type": "DEL", "copy_number": 1, "purity": 1.0,
                      "expected_log2": -1.0,
                      "note": (f"heterozygous, constitutional; NAHR between {CHROM}:{int(pair['start']) + 1:,}-"
                               f"{int(pair['end']):,} (LCR22A) and {CHROM}:{int(pair['other_start']) + 1:,}-"
                               f"{int(pair['other_end']):,} (LCR22D), {pair['frac_match'] * 100:.2f}% identity; "
                               f"contains ClinGen ISCA-37446 {CHROM}:{ISCA[0] + 1:,}-{ISCA[1]:,}")})
truth = pd.DataFrame(truth_rows)

# ------------------------------------------------------------------------------------------------
# 5. Genes (real MANE Select), segdups (real), write everything
# ------------------------------------------------------------------------------------------------
genes = mane[(mane["chrom"] == CHROM) & (mane["end"] > 17_000_000) & (mane["start0"] < 23_000_000)]
genes = genes.rename(columns={"start0": "start"})[["chrom", "start", "end", "gene", "strand", "gene_id",
                                                     "transcript", "refseq"]].sort_values(["start", "end", "gene"])

(DATA / "truth").mkdir(parents=True, exist_ok=True)
kw = dict(sep="\t", index=False, na_rep="", lineterminator="\n")
bins.to_csv(DATA / "cnv_chr22_10kb.tsv", float_format="%.6g", **kw)
genome.to_csv(DATA / "cnv_genome_1mb.tsv", float_format="%.6g", **kw)
genes.to_csv(DATA / "cnv_genes_chr22.tsv", **kw)
sd.to_csv(DATA / "cnv_segdups_chr22.tsv", float_format="%.6g", **kw)
truth.to_csv(DATA / "truth" / "cnv_truth.tsv", float_format="%.6g", **kw)

if __name__ == "__main__":
    print(f"chr22: {n:,} bins, mappable {int(mappable.sum()):,} (N {int(is_n.sum())}, satellite {int(is_sat.sum())}, "
          f"segdup >= {HI_ID:.0%}: {int(is_sd.sum())}); assembly gaps {ASM_GAPS}")
    print("LCR22 blocks:", {k: f"{a:,}-{b:,}" for k, (a, b) in LCR.items()})
    print(f"planted deletion {CHROM}:{DEL[0]:,}-{DEL[1]:,} ({(DEL[1] - DEL[0]) / 1e6:.3f} Mb); "
          f"patient mean depth {patient_mean:.1f}x; normals {normal_means.min():.1f}-{normal_means.max():.1f}x")
    inside = ok & (start >= ISCA[0]) & (end <= ISCA[1])
    outside = ok & ((end <= 18_000_000) | (start >= 22_000_000))
    print(f"log2 10 kb: median inside ISCA core {np.nanmedian(log2_10kb[inside]):.3f}, "
          f"outside {np.nanmedian(log2_10kb[outside]):.3f}, SD outside {np.nanstd(log2_10kb[outside]):.3f}")
    print(genome[genome.chrom == CHROM].iloc[14:24].to_string(index=False))
    print(truth.drop(columns="note").to_string(index=False))
