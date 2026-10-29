"""Simulated GIAB-style benchmark of four variant callers (practice B6, seed 20261027 + 7).

One genome, a 6-Mb window of hg38 chr20 (10,000,000-16,000,000, BED style), restricted to simulated
high-confidence (HC) regions, as in a Genome in a Bottle (GIAB) benchmark.

* Every base of the window gets a sequence context: unique, homopolymer/tandem repeat (homopolymer_tr),
  segmental duplication/low mappability (segdup_lowmap) or high GC (high_gc); ~10 % of the window is not HC.
  A variant's context is the context of its position.
* Truth variants: exactly 4,800; contexts multinomial 80/10/5/5 %; indels enriched in homopolymer/TR context
  (overall ~15 % indels).
* Each caller finds a truth variant with a recall that depends on caller, context and type. Misses are
  correlated between callers (Gaussian copula, correlation rho per context, weaker for DeepVariant): a variant
  that is hard for one caller tends to be hard for the others, and some true variants are found by one caller only.
* False positives are planted per context with a fixed pattern of WHICH callers make them (table FP_DESIGN),
  including shared false positives in segmental duplications that 3-4 callers agree on (reads from a
  paralogous copy mapped to the wrong place look like a real heterozygous SNV to every caller).

Alleles are simulated: REF does not match the real hg38 base at that position (no reference sequence is used).
Writes data/variant_callsets.tsv and data/truth/callsets_truth.json (planted design + reference values).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import dviz  # noqa: E402  (only dviz.upset_counts, so the reference counts match the notebook)

SEED = 20261027 + 7
rng = np.random.default_rng(SEED)

CHROM = "chr20"
WIN_START0, WIN_END = 10_000_000, 16_000_000            # BED style: 0-based start, end exclusive
WIN_LEN = WIN_END - WIN_START0
CALLERS = ["GATK", "DeepVariant", "FreeBayes", "Strelka2"]
CONTEXTS = ["unique", "homopolymer_tr", "segdup_lowmap", "high_gc"]
CODE = {c: i + 1 for i, c in enumerate(CONTEXTS)}        # 0 = outside the high-confidence regions
N_TRUTH_MEAN = 4800
CONTEXT_SHARE = [0.80, 0.10, 0.05, 0.05]
P_INDEL = {"unique": 0.105, "homopolymer_tr": 0.52, "segdup_lowmap": 0.09, "high_gc": 0.08}
RHO = {"unique": 0.50, "homopolymer_tr": 0.60, "segdup_lowmap": 0.75, "high_gc": 0.70}
# How strongly each caller follows the shared "difficulty" of a variant (1 = fully). DeepVariant (a neural
# network reading pileup images) makes errors that are less correlated with the three classical callers.
LOADING = np.array([1.0, 0.4, 1.0, 1.0])
MIN_SPACING = 10                                         # bp between any two simulated sites

# Recall (probability of calling a true variant): order GATK, DeepVariant, FreeBayes, Strelka2
RECALL = {
    ("unique", "SNV"): [0.990, 0.995, 0.985, 0.990],
    ("unique", "indel"): [0.960, 0.980, 0.920, 0.945],
    ("homopolymer_tr", "SNV"): [0.960, 0.980, 0.940, 0.960],
    ("homopolymer_tr", "indel"): [0.820, 0.920, 0.700, 0.790],
    ("segdup_lowmap", "SNV"): [0.880, 0.920, 0.860, 0.870],
    ("segdup_lowmap", "indel"): [0.780, 0.840, 0.720, 0.750],
    ("high_gc", "SNV"): [0.930, 0.970, 0.900, 0.930],
    ("high_gc", "indel"): [0.850, 0.910, 0.780, 0.830],
}

# False-positive design: (context, callers that make the same wrong call, expected number, share SNV)
G, D, F, S = CALLERS
FP_DESIGN = [
    ("unique", (G,), 35, 0.80), ("unique", (D,), 4, 0.80), ("unique", (F,), 75, 0.70),
    ("unique", (S,), 35, 0.80), ("unique", (G, F), 5, 0.70), ("unique", (F, S), 6, 0.70),
    ("unique", (G, S), 4, 0.80),
    ("homopolymer_tr", (G,), 20, 0.20), ("homopolymer_tr", (D,), 3, 0.20), ("homopolymer_tr", (F,), 55, 0.15),
    ("homopolymer_tr", (S,), 22, 0.20), ("homopolymer_tr", (G, F), 8, 0.10), ("homopolymer_tr", (F, S), 9, 0.10),
    ("homopolymer_tr", (G, F, S), 5, 0.10), ("homopolymer_tr", (G, D, F, S), 3, 0.10),
    ("segdup_lowmap", (G,), 10, 0.90), ("segdup_lowmap", (D,), 2, 0.90), ("segdup_lowmap", (F,), 18, 0.90),
    ("segdup_lowmap", (S,), 8, 0.90), ("segdup_lowmap", (G, F), 5, 0.90), ("segdup_lowmap", (G, F, S), 10, 0.95),
    ("segdup_lowmap", (G, D, F), 4, 0.95), ("segdup_lowmap", (D, F, S), 3, 0.95), ("segdup_lowmap", (G, D, S), 3, 0.95),
    ("segdup_lowmap", (G, D, F, S), 28, 0.95),
    ("high_gc", (G,), 8, 0.80), ("high_gc", (D,), 2, 0.80), ("high_gc", (F,), 14, 0.80),
    ("high_gc", (S,), 8, 0.80), ("high_gc", (G, F), 3, 0.80), ("high_gc", (G, D, F, S), 3, 0.85),
]
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
GC_SHARE = {"unique": 0.41, "homopolymer_tr": 0.30, "segdup_lowmap": 0.45, "high_gc": 0.70}


# ------------------------------------------------------------------------------------------------
# 1. Context of every base in the window (int8 array; 0 = not high-confidence)
# ------------------------------------------------------------------------------------------------
ctx = np.full(WIN_LEN, CODE["unique"], dtype=np.int8)


def place_blocks(n, lo, hi, code, margin=200):
    """Place n blocks of random length [lo, hi] entirely inside 'unique' sequence (plus a margin)."""
    blocks = []
    while len(blocks) < n:
        ln = int(rng.integers(lo, hi + 1))
        s = int(rng.integers(margin, WIN_LEN - ln - margin))
        if np.all(ctx[s - margin:s + ln + margin] == CODE["unique"]):
            ctx[s:s + ln] = code
            blocks.append((s, s + ln))
    return blocks


segdup_blocks = place_blocks(10, 15_000, 45_000, CODE["segdup_lowmap"], margin=2_000)
gc_blocks = place_blocks(70, 2_000, 7_000, CODE["high_gc"], margin=500)
tr_blocks = place_blocks(5_000, 12, 60, CODE["homopolymer_tr"], margin=20)
gaps = []
for _ in range(110):                                     # regions outside the GIAB-style HC set
    ln = int(np.clip(rng.lognormal(np.log(4_000), 0.9), 300, 60_000))
    s = int(rng.integers(0, WIN_LEN - ln))
    ctx[s:s + ln] = 0
    gaps.append((s, s + ln))
hc_bp = int((ctx > 0).sum())
context_bp = {c: int((ctx == CODE[c]).sum()) for c in CONTEXTS}

# ------------------------------------------------------------------------------------------------
# 2. Sites: positions for truth variants and false-positive calls, at least MIN_SPACING bp apart
# ------------------------------------------------------------------------------------------------
n_truth_total = N_TRUTH_MEAN                              # exactly 4,800 truth variants
n_truth = dict(zip(CONTEXTS, rng.multinomial(n_truth_total, CONTEXT_SHARE).tolist()))
fp_rows = [(c, callers, int(rng.poisson(mu)), p_snv) for c, callers, mu, p_snv in FP_DESIGN]
n_fp_ctx = {c: sum(n for cc, _, n, _ in fp_rows if cc == c) for c in CONTEXTS}

accepted_bins: dict[int, list[int]] = {}


def free(p):
    b = p // MIN_SPACING
    return all(abs(q - p) >= MIN_SPACING for bb in (b - 1, b, b + 1) for q in accepted_bins.get(bb, []))


sites = {}
for c in CONTEXTS:
    need = n_truth[c] + n_fp_ctx[c]
    pool = np.flatnonzero(ctx == CODE[c])
    cand = rng.choice(pool, size=int(need * 1.5) + 50, replace=False)
    got = []
    for p in cand.tolist():
        if free(p):
            accepted_bins.setdefault(p // MIN_SPACING, []).append(p)
            got.append(p)
            if len(got) == need:
                break
    assert len(got) == need, c
    sites[c] = np.array(got, dtype=np.int64)            # offsets in the window, in random order


# ------------------------------------------------------------------------------------------------
# 3. Alleles (VCF style: indels carry the base before the event as anchor; POS = anchor position)
# ------------------------------------------------------------------------------------------------
def random_base(gc):
    return str(rng.choice(["G", "C"])) if rng.random() < gc else str(rng.choice(["A", "T"]))


def snv_alleles(context):
    ref = random_base(GC_SHARE[context])
    if rng.random() < 2 / 3:                               # transitions : transversions ~ 2 : 1
        return ref, TRANSITION[ref]
    return ref, str(rng.choice([b for b in "ACGT" if b not in (ref, TRANSITION[ref])]))


def indel_alleles(context):
    if context == "homopolymer_tr":
        if rng.random() < 0.6:                             # one extra / missing base of a homopolymer run
            base = str(rng.choice(list("ACGT"), p=[0.35, 0.15, 0.15, 0.35]))
            seq = base * int(rng.choice([1, 2, 3], p=[0.75, 0.18, 0.07]))
        else:                                              # one or more units of a short tandem repeat
            unit = str(rng.choice(["CA", "TG", "AT", "AG", "TTC", "AAAT"]))
            seq = unit * int(rng.choice([1, 2, 3], p=[0.7, 0.2, 0.1]))
        anchor = str(rng.choice([b for b in "ACGT" if b != seq[0]]))
    else:
        ln = int(min(rng.geometric(0.45), 12))
        seq = "".join(rng.choice(list("ACGT"), size=ln).tolist())
        anchor = random_base(GC_SHARE[context])
    if rng.random() < 0.55:                                # deletion
        return anchor + seq, anchor
    return anchor, anchor + seq                            # insertion


def make_alleles(context, is_snv):
    return snv_alleles(context) if is_snv else indel_alleles(context)


# ------------------------------------------------------------------------------------------------
# 4. Truth variants and which callers find them (correlated misses)
# ------------------------------------------------------------------------------------------------
rows = []
for c in CONTEXTS:
    offs = sites[c][:n_truth[c]]
    is_indel = rng.random(offs.size) < P_INDEL[c]
    shared = rng.standard_normal(offs.size)
    own = rng.standard_normal((offs.size, len(CALLERS)))
    r = RHO[c] * LOADING                                   # correlation of each caller with the shared factor
    z = np.sqrt(r) * shared[:, None] + np.sqrt(1 - r) * own
    for i, off in enumerate(offs.tolist()):
        t = "indel" if is_indel[i] else "SNV"
        called = (z[i] < norm.ppf(RECALL[(c, t)])).astype(int)
        ref, alt = make_alleles(c, t == "SNV")
        rows.append([WIN_START0 + off + 1, ref, alt, t, c, *called.tolist(), 1])

# ------------------------------------------------------------------------------------------------
# 5. False-positive calls (not in the truth set)
# ------------------------------------------------------------------------------------------------
fp_used = {c: n_truth[c] for c in CONTEXTS}
for c, callers, n, p_snv in fp_rows:
    offs = sites[c][fp_used[c]:fp_used[c] + n]
    fp_used[c] += n
    for off in offs.tolist():
        t = "SNV" if rng.random() < p_snv else "indel"
        ref, alt = make_alleles(c, t == "SNV")
        rows.append([WIN_START0 + off + 1, ref, alt, t, c, *[int(k in callers) for k in CALLERS], 0])

cols = ["pos", "ref", "alt", "type", "context", *CALLERS, "in_truth_set"]
vc = pd.DataFrame(rows, columns=cols).sort_values(["pos", "ref", "alt"]).reset_index(drop=True)
vc.insert(0, "chrom", CHROM)
vc.insert(0, "variant_id", vc["chrom"] + ":" + vc["pos"].astype(str) + ":" + vc["ref"] + ":" + vc["alt"])
assert vc["variant_id"].is_unique and vc["pos"].is_unique
assert not ((vc[CALLERS].sum(axis=1) == 0) & (vc["in_truth_set"] == 0)).any()   # every row is a call or truth
off = vc["pos"].to_numpy() - 1 - WIN_START0
assert (ctx[off] == vc["context"].map(CODE).to_numpy()).all()                   # context follows position
(ROOT / "data" / "truth").mkdir(parents=True, exist_ok=True)
vc.to_csv(ROOT / "data" / "variant_callsets.tsv", sep="\t", index=False, na_rep="")


# ------------------------------------------------------------------------------------------------
# 6. Reference values for the CHECK cells (computed exactly as the notebook does)
# ------------------------------------------------------------------------------------------------
def r6(x):
    return None if x is None or not np.isfinite(x) else round(float(x), 6)


def prf(called, truth):
    called, truth = np.asarray(called, bool), np.asarray(truth, bool)
    tp, fp, fn = int((called & truth).sum()), int((called & ~truth).sum()), int((~called & truth).sum())
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if tp else 0.0
    return {"calls": tp + fp, "TP": tp, "FP": fp, "FN": fn, "precision": r6(prec), "recall": r6(rec), "F1": r6(f1)}


truth = vc["in_truth_set"].astype(bool)
n_callers = vc[CALLERS].sum(axis=1)
per_caller = {k: {"all": prf(vc[k], truth),
                  "SNV": prf(vc.loc[vc.type == "SNV", k], truth[vc.type == "SNV"]),
                  "indel": prf(vc.loc[vc.type == "indel", k], truth[vc.type == "indel"]),
                  "recall_by_context": {c: r6(vc.loc[truth & (vc.context == c), k].mean()) for c in CONTEXTS}}
              for k in CALLERS}
up = dviz.upset_counts(vc, CALLERS)
upset_top = [{"callers": "+".join(k for k in CALLERS if bool(r[k])) if not all(bool(r[k]) for k in CALLERS) else "all four",
              "count": int(r["count"]),
              "true": int(vc.loc[(vc[CALLERS].astype(bool) == r[CALLERS].astype(bool)).all(axis=1), "in_truth_set"].sum())}
             for _, r in up.head(10).iterrows()]
all4 = n_callers == 4
only = {k: (vc[k] == 1) & (n_callers == 1) for k in CALLERS}
consensus = {}
for name, mask in [(">=1", n_callers >= 1), (">=2", n_callers >= 2), (">=3", n_callers >= 3), ("=4", n_callers == 4)]:
    consensus[name] = prf(mask, truth)
shared_fp = (n_callers >= 3) & ~truth
missed_all = truth & (n_callers == 0)
truth_json = {
    "description": "Planted design and reference values for data/variant_callsets.tsv (practice B6). SIMULATED data.",
    "seed": SEED, "genome_build": "hg38", "chrom": CHROM,
    "window_bed": [WIN_START0, WIN_END], "window_1based": f"{CHROM}:{WIN_START0 + 1:,}-{WIN_END:,}",
    "hc_bp": hc_bp, "context_bp": context_bp, "context_bp_share_of_hc": {c: r6(v / hc_bp) for c, v in context_bp.items()},
    "n_blocks": {"segdup_lowmap": len(segdup_blocks), "high_gc": len(gc_blocks), "homopolymer_tr": len(tr_blocks),
                 "non_hc_gaps": len(gaps)},
    "design": {"n_truth": N_TRUTH_MEAN, "context_share": dict(zip(CONTEXTS, CONTEXT_SHARE)), "p_indel": P_INDEL,
               "rho_error_correlation": RHO, "caller_loading": dict(zip(CALLERS, LOADING.tolist())), "recall": {f"{c}|{t}": dict(zip(CALLERS, v)) for (c, t), v in RECALL.items()},
               "fp_design": [{"context": c, "callers": "+".join(k), "expected": mu, "planted": n, "share_SNV": p}
                             for (c, k, mu, p), (_, _, n, _) in zip(FP_DESIGN, fp_rows)]},
    "n_rows": int(len(vc)), "n_truth": int(truth.sum()), "n_false_calls": int((~truth).sum()),
    "n_truth_by_context": {c: int((truth & (vc.context == c)).sum()) for c in CONTEXTS},
    "n_truth_by_type": {t: int((truth & (vc.type == t)).sum()) for t in ["SNV", "indel"]},
    "n_truth_missed_by_all": int(missed_all.sum()),
    "missed_by_all_by_context": {c: int((missed_all & (vc.context == c)).sum()) for c in CONTEXTS},
    "per_caller": per_caller,
    "upset_top10": upset_top,
    "largest_intersection": upset_top[0]["callers"],
    "all_four": {"n": int(all4.sum()), "TP": int((all4 & truth).sum()), "FP": int((all4 & ~truth).sum()),
                 "precision": r6((all4 & truth).sum() / all4.sum()),
                 "FP_by_context": {c: int((all4 & ~truth & (vc.context == c)).sum()) for c in CONTEXTS}},
    "single_caller": {"n": int((n_callers == 1).sum()), "TP": int(((n_callers == 1) & truth).sum()),
                      "precision": r6(((n_callers == 1) & truth).sum() / (n_callers == 1).sum()),
                      "per_caller": {k: {"n": int(m.sum()), "TP": int((m & truth).sum())} for k, m in only.items()}},
    "deepvariant_only": {"n": int(only["DeepVariant"].sum()), "TP": int((only["DeepVariant"] & truth).sum())},
    "consensus": consensus,
    "shared_fp_min_callers": 3,
    "shared_fp_by_context": {c: int((shared_fp & (vc.context == c)).sum()) for c in CONTEXTS},
    "n_shared_fp": int(shared_fp.sum()),
    "shared_fp_top_context_planted": "segdup_lowmap",
}
(ROOT / "data" / "truth" / "callsets_truth.json").write_text(json.dumps(truth_json, indent=1) + "\n")

if __name__ == "__main__":
    print(f"rows {len(vc):,}  truth {int(truth.sum()):,}  false calls {int((~truth).sum()):,}  HC {hc_bp:,} bp")
    print("truth by context:", truth_json["n_truth_by_context"], "by type:", truth_json["n_truth_by_type"])
    for k in CALLERS:
        a = per_caller[k]
        print(f"{k:12s} calls {a['all']['calls']:5d}  P {a['all']['precision']:.4f}  R {a['all']['recall']:.4f}  "
              f"F1 {a['all']['F1']:.4f} | indel P {a['indel']['precision']:.3f} R {a['indel']['recall']:.3f}")
    print("upset:", [(u["callers"], u["count"], u["true"]) for u in upset_top])
    print("all four:", truth_json["all_four"])
    print("DV only:", truth_json["deepvariant_only"], " single-caller:", truth_json["single_caller"]["n"],
          truth_json["single_caller"]["precision"])
    for k, v in consensus.items():
        print(f"consensus {k:4s} calls {v['calls']:5d} P {v['precision']:.4f} R {v['recall']:.4f} F1 {v['F1']:.4f}")
    print("shared FP by context:", truth_json["shared_fp_by_context"], " missed by all:", truth_json["missed_by_all_by_context"])
