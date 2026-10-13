#!/usr/bin/env python3
"""Build data/practice2/ and checkpoints/practice2_*.tsv for Mini practice 2.

Mini practice 2 . Cohort QC: who does not belong?

A synthetic exome cohort of 2,000 sequenced samples, modelled on the UK Biobank exome QC in
Wang et al. 2021, Nature 597:527-532 (doi:10.1038/s41586-021-03855-y). The rules and
thresholds come from that paper and from the KING kinship guide; every sample and every
number is simulated. No UK Biobank data are used. Deterministic (fixed seed), no network.

Usage
  python build/build_practice2.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D2 = ROOT / "data" / "practice2"
CHECK = ROOT / "checkpoints"
SEED = 20261013

CENTRES = ["Centre A", "Centre B", "Centre C", "Centre D", "Centre E", "Centre F"]
CENTRE_W = [0.22, 0.20, 0.18, 0.16, 0.14, 0.10]
BATCHES = ["B1", "B2", "B3", "B4"]
BATCH_W = [0.26, 0.26, 0.24, 0.24]
BATCH_DEPTH = {"B1": 58.0, "B2": 57.0, "B3": 47.0, "B4": 58.0}
BATCH_PCT_OFFSET = {"B1": 0.0, "B2": 0.05, "B3": 0.80, "B4": 0.0}   # B3 = the batch effect
ANC_CENTRE = {"EUR": (0, 0, 0, 0), "AFR": (28, -6, 3, -1), "EAS": (9, 26, -2, 1),
              "SAS": (7, 9, -5, 2)}
ANC_SD = {"EUR": (1, 1, 1, 1), "AFR": (2.0, 1.5, 1.2, 1.2), "EAS": (1.2, 2.0, 1.2, 1.2),
          "SAS": (1.5, 1.5, 1.5, 1.2)}

# thresholds (paper / KING guide / teaching) - the notebook uses the same numbers
FREEMIX_MAX, CCDS10_MIN = 0.04, 94.5             # Wang et al. 2021
ARRAY_MIN = 0.95                                  # teaching threshold
Y_PRESENT, X_TWO = 0.20, 0.30                     # teaching thresholds
KIN_DUP, KIN_REL = 0.354, 0.0884                  # KING: duplicate/MZ; 2nd degree or closer
PR_EUR_MIN, PC_SD = 0.99, 4.0                     # Wang et al. 2021 (peddy + top four PCs)
PCS = ["PC1", "PC2", "PC3", "PC4"]


def draw_eur_pcs(rng, n):
    out = rng.normal(0, 1, (n, 4))
    bad = np.abs(out) > 3.5
    while bad.any():
        out[bad] = rng.normal(0, 1, bad.sum())
        bad = np.abs(out) > 3.5
    return out


def build(rng):
    # ---------------- participants ----------------
    n_part = 1998                                 # + 2 technical duplicate samples = 2,000
    anc = np.array(["EUR"] * 1838 + ["SAS"] * 60 + ["AFR"] * 45 + ["EAS"] * 30 + ["ADMIX"] * 25)
    anc = anc[rng.permutation(n_part)]
    P = pd.DataFrame({"pidx": np.arange(n_part), "ancestry": anc})
    P["sex"] = np.where(rng.random(n_part) < 0.54, "F", "M")
    P["age"] = rng.integers(40, 70, n_part)
    P["centre"] = rng.choice(CENTRES, n_part, p=CENTRE_W)
    P["year"] = rng.integers(2006, 2011, n_part)
    P["role"] = "singleton"
    pools = {a: list(rng.permutation(P.pidx[P.ancestry == a].to_numpy())) for a in ANC_CENTRE}

    def take(k, a="EUR"):
        return [int(pools[a].pop()) for _ in range(k)]

    pairs = []                                    # (pidx_a, pidx_b, kinship, relation)

    def kin(lo, hi):
        return round(float(rng.uniform(lo, hi)), 4)

    def fs():
        return round(float(np.clip(rng.normal(0.25, 0.015), 0.21, 0.29)), 4)

    def set_role(ids, role):
        P.loc[ids, "role"] = role

    # trio: father, mother, child (parents are unrelated to each other)
    f, m, c = take(3)
    P.loc[[f, m, c], "sex"] = ["M", "F", P.at[c, "sex"]]
    P.loc[f, "age"], P.loc[m, "age"], P.loc[c, "age"] = rng.integers(64, 70), rng.integers(62, 70), rng.integers(40, 44)
    pairs += [(f, c, kin(0.240, 0.260), "parent-offspring"), (m, c, kin(0.240, 0.260), "parent-offspring")]
    set_role([f, m, c], "trio")
    # full-sibling pairs (13 EUR + 1 SAS)
    for a in ["EUR"] * 13 + ["SAS"]:
        x, y = take(2, a)
        base = int(rng.integers(42, 68))
        P.loc[x, "age"], P.loc[y, "age"] = base, int(np.clip(base + rng.integers(-6, 7), 40, 69))
        pairs.append((x, y, fs(), "full siblings"))
        set_role([x, y], "sib pair")
    # two sibships of three
    for _ in range(2):
        trio = take(3)
        base = int(rng.integers(45, 64))
        for i, p in enumerate(trio):
            P.loc[p, "age"] = int(np.clip(base + rng.integers(-5, 6), 40, 69))
        for i in range(3):
            for j in range(i + 1, 3):
                pairs.append((trio[i], trio[j], fs(), "full siblings"))
        set_role(trio, "sibship of 3")
    # parent-offspring pairs
    for _ in range(4):
        par, kid = take(2)
        P.loc[par, "age"] = int(rng.integers(62, 70))
        P.loc[kid, "age"] = int(rng.integers(40, P.at[par, "age"] - 19))
        pairs.append((par, kid, kin(0.240, 0.260), "parent-offspring"))
        set_role([par, kid], "parent-offspring")
    # participant X will be sequenced twice; Y is X's full sibling
    X, Y = take(2)
    set_role([X, Y], "duplicate + sibling")
    # second-degree pairs (half-sibs, avuncular, grandparent-grandchild)
    for _ in range(14):
        x, y = take(2)
        pairs.append((x, y, kin(0.105, 0.150), "second degree"))
        set_role([x, y], "second-degree pair")
    # third-degree pairs (first cousins): 19 EUR + 1 AFR
    for a in ["EUR"] * 19 + ["AFR"]:
        x, y = take(2, a)
        pairs.append((x, y, kin(0.050, 0.080), "third degree"))
        set_role([x, y], "third-degree pair")

    planted = {}

    def plant(name, k):
        ids = take(k)
        set_role(ids, name)
        planted[name] = ids
        return ids

    withdrawn = plant("consent withdrawn", 3)
    contam = plant("contaminated", 4)
    lowcov = plant("low coverage", 3)
    border = plant("borderline coverage (passes)", 1)
    swapA, swapB = plant("swapped pair", 2)
    meta = plant("metadata error", 1)[0]
    xxy = plant("XXY-like", 1)[0]
    samesex = plant("same-sex array discordance", 1)[0]
    Z = plant("technical duplicate", 1)[0]
    T1, T2 = plant("twin or unrecorded duplicate", 2)
    pc_out = plant("PC outlier (EUR-like)", 6)
    low_eur = plant("Pr(EUR) just under 0.99", 12)
    P.loc[contam, "sex"] = ["M", "M", "F", "F"]
    P.loc[swapA, "sex"], P.loc[swapB, "sex"] = "F", "M"
    P.loc[meta, "sex"], P.loc[xxy, "sex"], P.loc[samesex, "sex"] = "M", "M", "M"
    P.loc[[T1, T2], "sex"] = "F"
    for col in ["age", "centre", "year"]:
        P.loc[T2, col] = P.at[T1, col]
    P.loc[swapB, "centre"] = P.at[swapA, "centre"]
    pairs.append((T1, T2, 0.4931, "MZ twins or unrecorded duplicate"))

    # ---------------- samples ----------------
    S = pd.DataFrame({"pidx": np.arange(n_part)})
    S["copy_id"] = "a"
    S = pd.concat([S, pd.DataFrame({"pidx": [X, Z], "copy_id": ["b", "b"]})], ignore_index=True)
    n = len(S)
    S["dna_of"] = S.pidx                       # whose DNA is actually in the tube
    S.loc[(S.pidx == swapA), "dna_of"] = swapB
    S.loc[(S.pidx == swapB), "dna_of"] = swapA
    S["batch"] = rng.choice(BATCHES, n, p=BATCH_W)
    S["specimen"] = np.where(rng.random(n) < 0.06, "saliva", "blood")
    for pid, b_a, b_b in [(X, "B1", "B3"), (Z, "B2", "B4")]:
        S.loc[(S.pidx == pid) & (S.copy_id == "a"), "batch"] = b_a
        S.loc[(S.pidx == pid) & (S.copy_id == "b"), "batch"] = b_b
        spec = S.loc[(S.pidx == pid) & (S.copy_id == "a"), "specimen"].iloc[0]
        S.loc[S.pidx == pid, "specimen"] = spec
    S.loc[S.pidx.isin([swapA, swapB]), "batch"] = "B2"
    S.loc[S.pidx.isin([swapA, swapB]), "specimen"] = "blood"
    for pid, b in zip(lowcov, ["B3", "B3", "B1"]):
        S.loc[S.pidx == pid, "batch"] = b
    S.loc[S.pidx == border[0], "batch"] = "B3"
    S["consent_status"] = np.where(S.pidx.isin(withdrawn), "withdrawn", "active")

    # sample IDs: random 5-digit numbers; the swapped pair sat next to each other on a plate
    nums = rng.choice(np.arange(10000, 99990), n, replace=False)
    ia, ib = S.index[S.pidx == swapA][0], S.index[S.pidx == swapB][0]
    taken = set(nums.tolist())
    k = nums[ia] + 1
    while k in taken:
        k += 1
    nums[ib] = k
    S["sample_id"] = ["WES%05d" % v for v in nums]
    assert S.sample_id.is_unique
    pids = rng.choice(np.arange(100000, 999999), n_part, replace=False)
    P["participant_id"] = ["P%06d" % v for v in pids]

    # ---------------- QC metrics ----------------
    dna = P.set_index("pidx").loc[S.dna_of]
    gsex = dna.sex.to_numpy().copy()
    gsex[S.dna_of.to_numpy() == xxy] = "XXY"
    depth = np.array([BATCH_DEPTH[b] for b in S.batch]) + rng.normal(0, 7, n)
    depth -= 4 * (S.specimen == "saliva").to_numpy()
    depth = np.clip(depth, 28, 95)
    pct = 100 - 60 / depth - np.array([BATCH_PCT_OFFSET[b] for b in S.batch]) + rng.normal(0, 0.25, n)
    pct = np.clip(pct, 95.5, 99.9)
    freemix = np.clip(np.exp(rng.normal(np.log(0.004), 0.6, n)), 0.0005, 0.030)
    xhet_f = np.clip(rng.normal(0.60, 0.04, n), 0.45, 0.75)
    xhet_m = np.clip(rng.normal(0.020, 0.008, n), 0.002, 0.06)
    y_m = np.clip(rng.normal(0.48, 0.03, n), 0.38, 0.58)
    y_f = np.clip(rng.normal(0.010, 0.005, n), 0.0, 0.03)
    xhet = np.where(np.isin(gsex, ["F", "XXY"]), xhet_f, xhet_m)
    ydep = np.where(np.isin(gsex, ["M", "XXY"]), y_m, y_f)
    arr = np.clip(rng.normal(0.9975, 0.0008, n), 0.993, 0.9995)

    # ancestry PCs: one draw per participant; samples inherit the PCs of the DNA they carry
    pcs = np.zeros((n_part, 4))
    pr = np.zeros(n_part)
    for a in ANC_CENTRE:
        idx = np.flatnonzero(P.ancestry.to_numpy() == a)
        if a == "EUR":
            pcs[idx] = draw_eur_pcs(rng, len(idx))
            pr[idx] = 1 - np.minimum(rng.exponential(0.0008, len(idx)), 0.004)
        else:
            pcs[idx] = np.array(ANC_CENTRE[a]) + rng.normal(0, 1, (len(idx), 4)) * np.array(ANC_SD[a])
            pr[idx] = rng.uniform(0.001, 0.03, len(idx)) if a == "SAS" else rng.uniform(0, 0.001, len(idx))
    idx = np.flatnonzero(P.ancestry.to_numpy() == "ADMIX")
    other = np.where(np.arange(len(idx)) < 15, "AFR", "SAS")
    frac = rng.uniform(0.25, 0.75, len(idx))
    pcs[idx] = frac[:, None] * np.array([ANC_CENTRE[o] for o in other]) + rng.normal(0, 1, (len(idx), 4))
    pr[idx] = np.clip(1 - frac + rng.normal(0, 0.08, len(idx)), 0.02, 0.95)
    for i, p in enumerate(pc_out):                # off-cluster on PC3/PC4 only
        pc = 2 if i < 4 else 3
        pcs[p, pc] = rng.choice([-1, 1]) * rng.uniform(6.5, 8.5)
        pr[p] = 0.999
    pr[low_eur] = rng.uniform(0.93, 0.985, len(low_eur))
    pcs[T2] = pcs[T1] + rng.normal(0, 0.02, 4)     # MZ twins: identical genomes
    pr[T2] = pr[T1]
    spc = pcs[S.dna_of.to_numpy()] + rng.normal(0, 0.01, (n, 4)) * (S.copy_id == "b").to_numpy()[:, None]
    spr = pr[S.dna_of.to_numpy()]

    # planted metric values
    def at(pid, copy="a"):
        return int(S.index[(S.pidx == pid) & (S.copy_id == copy)][0])

    for pid, fm, ac in zip(contam, [0.052, 0.067, 0.089, 0.118], [0.9851, 0.9812, 0.9760, 0.9690]):
        i = at(pid)
        freemix[i], arr[i] = fm, ac
        if gsex[i] == "M":
            xhet[i] = 0.02 + 1.2 * fm               # contamination inflates X heterozygosity
    for pid, d, pc_ in zip(lowcov, [24.1, 27.6, 29.8], [92.81, 93.94, 94.12]):
        depth[at(pid)], pct[at(pid)] = d, pc_
    depth[at(border[0])], pct[at(border[0])] = 31.2, 94.81
    for pid, ac in [(swapA, 0.6123), (swapB, 0.6047), (samesex, 0.5832)]:
        arr[at(pid)] = ac
    arr[at(meta)], arr[at(xxy)] = 0.9991, 0.9986
    xhet[at(xxy)], ydep[at(xxy)] = 0.552, 0.462
    pct[at(X, "a")], pct[at(X, "b")] = 98.91, 97.62
    pct[at(Y)] = 98.60
    pct[at(Z, "a")], pct[at(Z, "b")] = 99.02, 98.37
    pct[at(T1)], pct[at(T2)] = 98.75, 98.41

    reported = P.set_index("pidx").loc[S.pidx, "sex"].to_numpy().copy()
    reported[S.pidx.to_numpy() == meta] = "F"       # the record says F; the DNA (and array) say XY

    part = P.set_index("pidx").loc[S.pidx]
    manifest = pd.DataFrame({
        "sample_id": S.sample_id, "participant_id": part.participant_id.to_numpy(),
        "reported_sex": reported, "age_at_recruitment": part.age.to_numpy(),
        "recruitment_centre": part.centre.to_numpy(), "specimen": S.specimen,
        "sequencing_batch": S.batch, "consent_status": S.consent_status,
        "collection_year": part.year.to_numpy()})
    metrics = pd.DataFrame({
        "sample_id": S.sample_id, "mean_depth": np.round(depth, 1),
        "pct_ccds_10x": np.round(pct, 2), "freemix": np.round(freemix, 4),
        "chrX_het_rate": np.round(xhet, 3), "chrY_depth_ratio": np.round(ydep, 3),
        "array_concordance": np.round(arr, 4)})
    for j, c in enumerate(PCS):
        metrics[c] = np.round(spc[:, j], 3)
    metrics["pr_eur"] = np.round(spr, 4)

    # kinship pairs (sample level); duplicates of X are both related to Y
    sid = {(int(r.pidx), r.copy_id): r.sample_id for r in S.itertuples()}
    rows = [(sid[(a, "a")], sid[(b, "a")], k, rel) for a, b, k, rel in pairs]
    rows += [(sid[(X, "a")], sid[(X, "b")], 0.4986, "technical duplicate"),
             (sid[(Z, "a")], sid[(Z, "b")], 0.4962, "technical duplicate"),
             (sid[(X, "a")], sid[(Y, "a")], 0.2463, "full siblings"),
             (sid[(X, "b")], sid[(Y, "a")], 0.2481, "full siblings")]
    kin_df = pd.DataFrame([(min(a, b), max(a, b), k, rel) for a, b, k, rel in rows],
                          columns=["sample_id_1", "sample_id_2", "kinship", "relation"])
    kin_df = kin_df.sort_values(["sample_id_1", "sample_id_2"]).reset_index(drop=True)
    assert (kin_df.kinship > 0.0442).all()

    # planted-case answer key
    def sids(pid_list, copy="a"):
        return [sid[(p, copy)] for p in pid_list]

    key = []
    for s_ in sids(withdrawn):
        key.append(("consent withdrawn", s_, "1_consent", "remove first, before any QC"))
    for s_ in sids(contam):
        key.append(("contaminated (FREEMIX > 0.04)", s_, "2_contamination", "exclude"))
    for s_ in sids(lowcov):
        key.append(("low coverage (< 94.5% CCDS at 10x)", s_, "3_coverage", "exclude"))
    key.append(("borderline coverage 94.81% in B3", sid[(border[0], "a")], "kept", "passes; do not move the threshold"))
    key.append(("swapped pair, sample 1 (sex + array discordant)", sid[(swapA, "a")], "4_identity", "exclude both"))
    key.append(("swapped pair, sample 2 (sex + array discordant)", sid[(swapB, "a")], "4_identity", "exclude both"))
    key.append(("metadata error (sex discordant, array 0.999)", sid[(meta, "a")], "kept", "correct reported_sex to M, keep"))
    key.append(("XXY-like (Y present + high X het)", sid[(xxy, "a")], "kept", "flag for review; never call it a swap"))
    key.append(("same-sex array discordance", sid[(samesex, "a")], "4_identity", "exclude"))
    key.append(("technical duplicate (kept copy; has a full sibling)", sid[(X, "a")], "kept",
                "keep better pct_ccds_10x; it then wins the sibling pair"))
    key.append(("technical duplicate (dropped copy)", sid[(X, "b")], "5_duplicates", "drop lower pct_ccds_10x"))
    key.append(("full sibling of the duplicated participant", sid[(Y, "a")], "6_relatedness",
                "counted once, after the duplicate is resolved"))
    key.append(("technical duplicate (kept copy)", sid[(Z, "a")], "kept", "keep better pct_ccds_10x"))
    key.append(("technical duplicate (dropped copy)", sid[(Z, "b")], "5_duplicates", "drop lower pct_ccds_10x"))
    key.append(("twin or unrecorded duplicate (kept)", sid[(T1, "a")], "kept", "review; keep one"))
    key.append(("twin or unrecorded duplicate (dropped)", sid[(T2, "a")], "5_duplicates", "review; keep one"))
    for s_ in sids(pc_out):
        key.append(("PC outlier with Pr(EUR) >= 0.99", s_, "8_ancestry_pc_outlier", "outside the EUR analysis set"))
    for s_ in sids(low_eur):
        key.append(("European cluster, Pr(EUR) < 0.99", s_, "7_ancestry_pr_eur", "outside the EUR analysis set"))
    key_df = pd.DataFrame(key, columns=["case", "sample_id", "expected_step", "expected_handling"])
    return manifest, metrics, kin_df, key_df, P


# ---------------------------------------------------------------------------------------
# reference QC: an independent implementation of the notebook's steps
# ---------------------------------------------------------------------------------------
def genetic_sex(x_het, y_ratio):
    y, xx = y_ratio >= Y_PRESENT, x_het >= X_TWO
    if y and not xx:
        return "XY"
    if xx and not y:
        return "XX"
    return "ambiguous"


def reference_qc(manifest, metrics, kin):
    d = manifest.merge(metrics, on="sample_id", validate="one_to_one").set_index("sample_id").sort_index()
    alive = {s: True for s in d.index}
    excluded, funnel, flags = [], [("0_sequenced", "all delivered exomes", 0, len(d))], []

    def live():
        return d[[alive[s] for s in d.index]]

    def drop(items, step, rule):
        for s, why in items:
            alive[s] = False
            excluded.append((s, step, why))
        funnel.append((step, rule, len(items), sum(alive.values())))

    c = live()
    drop([(s, "consent withdrawn: removed before any QC") for s in c.index[c.consent_status == "withdrawn"]],
         "1_consent", "consent_status == withdrawn")
    c = live()
    drop([(s, f"freemix {v:.3f} > 0.04") for s, v in c.freemix[c.freemix > FREEMIX_MAX].items()],
         "2_contamination", "freemix > 0.04 (VerifyBAMID)")
    c = live()
    drop([(s, f"pct_ccds_10x {v:.2f} < 94.5") for s, v in c.pct_ccds_10x[c.pct_ccds_10x < CCDS10_MIN].items()],
         "3_coverage", "pct_ccds_10x < 94.5")
    c = live()
    items = []
    for s, r in c.iterrows():
        g = genetic_sex(r.chrX_het_rate, r.chrY_depth_ratio)
        expect = {"F": "XX", "M": "XY"}[r.reported_sex]
        status = "ambiguous" if g == "ambiguous" else ("concordant" if g == expect else "discordant")
        if r.array_concordance < ARRAY_MIN:
            if status == "discordant":
                items.append((s, f"array_concordance {r.array_concordance:.3f} < 0.95 and reported "
                                 f"{r.reported_sex} but genetic {g}: likely swap"))
            else:
                items.append((s, f"array_concordance {r.array_concordance:.3f} < 0.95 with {status} "
                                 f"sex: identity not confirmed"))
        elif status == "discordant":
            flags.append((s, "metadata error: reported sex corrected"))
        elif status == "ambiguous":
            flags.append((s, "sex chromosomes ambiguous: review"))
    drop(items, "4_identity", "array_concordance < 0.95")

    pct = d.pct_ccds_10x
    pid = d.participant_id
    items, gone = [], set()
    for r in kin[kin.kinship > KIN_DUP].itertuples():
        a, b = r.sample_id_1, r.sample_id_2
        if not (alive[a] and alive[b]) or a in gone or b in gone:
            continue
        lose, keep = (a, b) if pct[a] < pct[b] else (b, a)
        gone.add(lose)
        if pid[a] == pid[b]:
            items.append((lose, f"duplicate of {keep} (same participant {pid[a]}), lower pct_ccds_10x"))
        else:
            items.append((lose, f"kinship {r.kinship:.3f} with {keep} (different participant): "
                                f"MZ twin or unrecorded duplicate, keep one, review"))
            flags.append((keep, "kept member of a twin/duplicate pair: review"))
    drop(items, "5_duplicates", "kinship > 0.354: keep the better sample")

    rel = kin[kin.kinship > KIN_REL]
    edges = [(a, b) for a, b in zip(rel.sample_id_1, rel.sample_id_2) if alive[a] and alive[b]]
    items = []
    while edges:
        deg = {}
        for a, b in edges:
            deg[a] = deg.get(a, 0) + 1
            deg[b] = deg.get(b, 0) + 1
        top = max(deg.values())
        victim = sorted((s for s, k in deg.items() if k == top), key=lambda s: (pct[s], s))[0]
        items.append((victim, f"kinship > 0.0884 with {top} remaining sample(s)"))
        edges = [e for e in edges if victim not in e]
    drop(items, "6_relatedness", "kinship > 0.0884: greedy, most relatives first")
    funnel.append(("pan_ancestry_set", "QC-passed, unrelated, any ancestry", 0, sum(alive.values())))

    c = live()
    drop([(s, f"pr_eur {v:.3f} < 0.99") for s, v in c.pr_eur[c.pr_eur < PR_EUR_MIN].items()],
         "7_ancestry_pr_eur", "pr_eur < 0.99 (peddy)")
    c = live()
    z = (c[PCS] - c[PCS].mean()) / c[PCS].std(ddof=1)
    items = []
    for s, row in z.iterrows():
        worst = row.abs().idxmax()
        if abs(row[worst]) > PC_SD:
            items.append((s, f"{worst} {row[worst]:+.1f} SD from the mean (limit 4)"))
    drop(items, "8_ancestry_pc_outlier", "outside mean +/- 4 SD on PC1-PC4")
    funnel.append(("eur_analysis_set", "European-ancestry analysis set", 0, sum(alive.values())))

    funnel_df = pd.DataFrame(funnel, columns=["step", "rule", "n_removed", "n_remaining"])
    excl_df = pd.DataFrame(excluded, columns=["sample_id", "step", "reason"])
    excl_df = excl_df.sort_values(["step", "sample_id"]).reset_index(drop=True)
    flags_df = pd.DataFrame(flags, columns=["sample_id", "flag"])
    return funnel_df, excl_df, flags_df


def main():
    D2.mkdir(parents=True, exist_ok=True)
    CHECK.mkdir(exist_ok=True)
    manifest, metrics, kin, key, P = build(np.random.default_rng(SEED))
    manifest.sort_values("sample_id").to_csv(D2 / "cohort_manifest.tsv", sep="\t", index=False)
    metrics.sort_values("sample_id").to_csv(D2 / "sample_qc_metrics.tsv", sep="\t", index=False)
    kin[["sample_id_1", "sample_id_2", "kinship"]].to_csv(D2 / "kinship_pairs.tsv", sep="\t", index=False)

    # re-read as a student would
    manifest = pd.read_csv(D2 / "cohort_manifest.tsv", sep="\t")
    metrics = pd.read_csv(D2 / "sample_qc_metrics.tsv", sep="\t")
    kin = pd.read_csv(D2 / "kinship_pairs.tsv", sep="\t")
    funnel, excl, flags = reference_qc(manifest, metrics, kin)

    # ---- the cohort must contain exactly what the README promises ----
    deg = pd.cut(kin.kinship, [0.0442, 0.0884, 0.177, 0.354, 1], labels=["3rd", "2nd", "1st", "dup"])
    counts = deg.value_counts().to_dict()
    assert counts == {"1st": 28, "3rd": 20, "2nd": 14, "dup": 3}, counts
    assert len(manifest) == 2000 and manifest.participant_id.nunique() == 1998
    assert manifest.age_at_recruitment.between(40, 69).all()
    assert manifest.collection_year.between(2006, 2010).all()
    steps = dict(zip(funnel.step, funnel.n_removed))
    assert steps["1_consent"] == 3 and steps["2_contamination"] == 4 and steps["3_coverage"] == 3
    assert steps["4_identity"] == 3 and steps["5_duplicates"] == 3 and steps["6_relatedness"] == 38
    assert steps["8_ancestry_pc_outlier"] == 6
    fr = funnel.n_remaining.to_numpy()
    assert all(fr[i] == fr[i - 1] - funnel.n_removed[i] for i in range(1, len(fr)))
    assert len(excl) + fr[-1] == len(manifest)
    for r in key.itertuples():                      # every planted case lands where documented
        hit = excl[excl.sample_id == r.sample_id]
        if r.expected_step == "kept":
            assert hit.empty, (r.case, r.sample_id, hit.to_dict("records"))
        else:
            assert len(hit) == 1 and hit.step.iloc[0] == r.expected_step, (r.case, r.sample_id, hit.to_dict("records"))
    flagged = dict(zip(flags.sample_id, flags.flag))
    assert len(flagged) == 3, flagged
    q = manifest.merge(metrics, on="sample_id")
    b = q[q.consent_status == "active"].groupby("sequencing_batch").pct_ccds_10x.median()
    assert b["B3"] < b.drop("B3").min() - 0.5, b     # visible batch effect
    eur_share = (q.pr_eur >= 0.99).mean()
    assert 0.90 <= eur_share <= 0.93, eur_share

    funnel.to_csv(CHECK / "practice2_qc_funnel.tsv", sep="\t", index=False)
    excl.to_csv(CHECK / "practice2_excluded_samples.tsv", sep="\t", index=False)
    key.to_csv(CHECK / "practice2_planted_cases.tsv", sep="\t", index=False)
    print(funnel.to_string(index=False))
    print(flags.to_string(index=False))
    print(b.round(2).to_dict(), "EUR share (pr_eur>=0.99):", round(eur_share, 3))


if __name__ == "__main__":
    main()
