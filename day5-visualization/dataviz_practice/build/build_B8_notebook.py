"""Write electives/B8_tp53_lollipop.ipynb (Practice B8 · TP53 lollipop: where on the protein?).

Numbers quoted in the markdown are read from results/B8.json of the previous execution (two passes, as for B5):
    python build/build_B8_notebook.py && python build/run_notebook.py electives/B8_tp53_lollipop.ipynb electives/B8_tp53_lollipop_executed.ipynb
    python build/build_B8_notebook.py && python build/run_notebook.py electives/B8_tp53_lollipop.ipynb electives/B8_tp53_lollipop_executed.ipynb
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from nbbuild import NB  # noqa: E402

RES_PATH = ROOT / "results" / "B8.json"
RES = json.loads(RES_PATH.read_text()) if RES_PATH.exists() else {}
SN = RES.get("slide_numbers", {})


def v(key, fmt="{:.2f}"):
    """A number from the previous run (CHECK value or slide number), formatted; '…' before the first run."""
    if key in RES and isinstance(RES[key], dict) and "value" in RES[key]:
        x = RES[key]["value"]
    elif key in SN:
        x = SN[key]["value"]
    else:
        return "…"
    out = fmt.format(x)
    return "−" + out[1:] if out.startswith("-") else out


def p_fmt(key):
    """p-value as a × 10^b in plain text (e.g. 3.1 × 10⁻¹⁵²)."""
    if key not in SN and key not in RES:
        return "…"
    x = SN[key]["value"] if key in SN else RES[key]["value"]
    if x >= 0.001:
        return f"{x:.2f}"
    m, e = f"{x:.1e}".split("e")
    sup = str(int(e)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{m} × 10{sup}"


def hot_list(key):
    """'R175H (71), R248Q/W (68), …' from the saved hotspot table."""
    rows = SN.get(key, {}).get("value")
    return ", ".join(f"{r['label']} ({r['n']})" for r in rows) if rows else "…"


nb = NB()
nb.md("""
# Practice B8 · TP53 lollipop: where on the protein?

**Variants → tumours per residue → binomial test and Fisher's exact test → lollipop plot with domains**

*Genomic Data Processing · Day 5 (10/27) · Data visualization in genetics and genomics · about 12 minutes*

> **All data in this practice are simulated (computer-generated). No real patient data are used.**
> The tumour counts are simulated, but the TP53 protein (393 amino acids, UniProt P04637), its domains (IARC TP53
> Database convention) and the hotspot residues are real, so the picture can be compared with cBioPortal or the
> TP53 database.

**How to use this notebook.** Run the cells from top to bottom (Shift + Enter). The steps follow the four beats of
Practice B: **Question → Test → Check → Visualize**. Cells that start with `# CHECK ·` compare your numbers with a
reference value and print ✓ (fine) or ✗ (look again). Each step ends with **What changed?** and a **Your turn** task;
try it first, then click *Answer* to unfold the answer.

| File | What is inside |
|---|---|
| `data/tp53_variants.tsv` | one TP53 variant per tumour: `tumour_id`, `protein_change` (e.g. R175H, R196*, P72fs, X125_splice), `position` (amino acid, 1-based), `variant_class` |
| `data/tp53_domains.tsv` | four protein domains: `aa_start`, `aa_end` (1-based, both ends included), `length_aa` |
| `data/truth/tp53_truth.tsv` | the simulation catalogue with the planted counts – open it only in Step B8.3 |

| Step | Beat | Topic | min |
|---|---|---|---|
| B8.1 | Question | two ways to break a tumour suppressor | 2 |
| B8.2 | Test | tumours per residue; binomial and Fisher tests | 4 |
| B8.3 | Check | against the truth and a hand calculation | 2 |
| B8.4 | Visualize | the lollipop plot | 4 |
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
from scipy import stats
""")

# ------------------------------------------------------------------------------------------------ B8.1
nb.md("""
## Step B8.1 · Question: two ways to break a tumour suppressor  *(2 min)*

### Question

*Do missense and truncating variants fall in different parts of the TP53 protein?*

**Why this step.** *TP53* is a tumour suppressor: its protein, p53, binds DNA and switches on genes that stop damaged
cells from dividing. It is mutated in about half of all cancers. A tumour can disable p53 in two ways. A **missense**
change in the DNA-binding domain leaves a full-length protein that can no longer bind DNA – and the same few residues
are hit again and again (**hotspots**). A **truncating** change (nonsense or frameshift) stops the protein early,
wherever it happens. If missense variants are selected for *what they do* in the DNA-binding domain, they should be
concentrated there; truncating variants only need to destroy the protein, so they should fall roughly in proportion
to length.

**New terms**
* **Residue (amino-acid position)** – the number of an amino acid along the protein, counted from 1 (here 1–393).
* **Domain** – a part of the protein with its own job: transactivation (TAD), proline-rich (PRD), DNA-binding (DBD), tetramerisation (TET).
* **Missense** – one amino acid replaced by another: R175H = arginine (R) at residue 175 → histidine (H).
* **Truncating** – nonsense (a new stop codon, `R196*`) or frameshift (`P72fs`): the protein is cut short.
* **Splice/other** – changes at exon edges (`X125_splice`, position = nearest codon) or small in-frame deletions.
* **Hotspot** – a residue mutated in many unrelated tumours.
""")
nb.code("""
# Step B8.1 · load the variants (one row per tumour) and the protein domains
var = pd.read_csv(DATA / "tp53_variants.tsv", sep="\\t")
dom = pd.read_csv(DATA / "tp53_domains.tsv", sep="\\t")
LENGTH = 393                                      # TP53 protein length in amino acids (UniProt P04637)
print(f"{var['tumour_id'].nunique():,} tumours, {len(var):,} variants (one per tumour)")
print((var["variant_class"].value_counts() / len(var) * 100).round(1).rename("% of tumours").to_string())
var.head()
""")
nb.code("""
# Step B8.1 · the four domains (amino-acid positions are 1-based and both ends are included)
dom[["domain", "short", "aa_start", "aa_end", "length_aa"]]
""")
nb.code("""
# CHECK · the cohort, the spectrum and the length of the DNA-binding domain
dbd = dom.loc[dom["short"] == "DBD"].iloc[0]
dbd_len = int(dbd["aa_end"] - dbd["aa_start"] + 1)          # both ends count: 292 - 94 + 1
p_dbd = dbd_len / LENGTH                                    # the DBD's share of the protein
print(f"DNA-binding domain: residues {dbd['aa_start']}-{dbd['aa_end']} = {dbd_len} of {LENGTH} = {p_dbd:.3f}")
dviz.check("B8.n_tumours", int(var["tumour_id"].nunique()), 1200)
dviz.check("B8.pct_missense", 100 * float((var["variant_class"] == "missense").mean()), between=(70, 76))
dviz.check("B8.dbd_length", dbd_len, 199)
""")
nb.md(f"""
**What changed?** Nothing is computed yet: {v('B8.n_tumours', '{:,}')} tumours, one TP53 variant each –
{v('pct_missense', '{:.0f}')}% missense, {v('pct_truncating', '{:.0f}')}% truncating, {v('pct_splice', '{:.0f}')}%
splice/other. The DNA-binding domain covers {v('B8.dbd_length', '{}')} of 393 residues = {v('dbd_fraction', '{:.3f}')}
of the protein: this is what we expect inside it if position does not matter.

**Your turn.** The domain table says 94–292. Why is the length 199 and not 292 − 94 = 198?

<details><summary>Answer</summary>

Protein positions are **1-based and inclusive**: residue 94 and residue 292 both belong to the domain, so the
length is 292 − 94 + 1 = 199. (A BED interval is different: 0-based start, end excluded, so there end − start is the
length.) An off-by-one here changes the expectation only slightly, but in genomic coordinates the same mistake moves
every breakpoint by one base.
</details>
""")

# ------------------------------------------------------------------------------------------------ B8.2
nb.md("""
## Step B8.2 · Test: tumours per residue, then test the DNA-binding domain  *(4 min)*

### Test

**Why this step.** Counting tumours per residue shows the hotspots. To say "more than expected" we need an
expectation: if variants fell at random positions, the share inside the DBD would equal its share of the length,
199/393 = 0.506. The **binomial test** asks how surprising the observed number inside the DBD is under that
expectation. **Fisher's exact test** compares missense with truncating directly, without assuming 0.506.

**New terms**
* **Binomial test** – *k* "successes" (variants inside the DBD) out of *n*, compared with an expected probability *p₀*.
* **2 × 2 table** – counts split two ways: class (missense / truncating) × place (inside / outside the DBD).
* **Fisher's exact test** – the exact p-value for "the two rows of a 2 × 2 table have the same proportions".
* **Odds ratio** – (inside ÷ outside for missense) divided by (inside ÷ outside for truncating); 1 = no difference.
""")
nb.code("""
# Step B8.2 · count tumours per residue and class; the 10 most frequently mutated residues
per_pos = var.groupby(["position", "variant_class"]).size().rename("n").reset_index()
mine = var.groupby("position").size()                              # tumours per residue (all classes)
top = var.groupby(["position", "protein_change"]).size().rename("n").reset_index()
top = top.sort_values(["position", "n", "protein_change"], ascending=[True, False, True])
top["allele"] = top["protein_change"] + " (" + top["n"].astype(str) + ")"
alleles = top.groupby("position")["allele"].agg(", ".join)           # most frequent allele first
pd.DataFrame({"tumours": mine, "alleles (tumours)": alleles}).sort_values("tumours", ascending=False, kind="stable").head(10)
""")
nb.code("""
# Step B8.2 · binomial test: is the share inside the DNA-binding domain larger than its share of the length?
var["in_dbd"] = var["position"].between(dbd["aa_start"], dbd["aa_end"])     # inclusive on both ends
tests = {}
for cls in ["missense", "truncating"]:
    sub = var.loc[var["variant_class"] == cls]
    k, n = int(sub["in_dbd"].sum()), len(sub)
    res = stats.binomtest(k, n, p_dbd)                                     # two-sided by default
    tests[cls] = {"inside": k, "n": n, "fraction": k / n, "p": float(res.pvalue)}
    print(f"{cls:10s} {k:3d} of {n} inside = {k / n:.3f} (expected {p_dbd:.3f})   binomial p = {res.pvalue:.2g}")
""")
nb.code("""
# Step B8.2 · Fisher's exact test on the 2 x 2 table: class (missense / truncating) x place (inside / outside)
table = pd.DataFrame({"inside DBD": [tests[c]["inside"] for c in tests],
                      "outside DBD": [tests[c]["n"] - tests[c]["inside"] for c in tests]}, index=list(tests))
odds, p_fisher = stats.fisher_exact(table.to_numpy())
print(table.to_string())
print(f"odds ratio {odds:.1f} (missense vs truncating, odds of being inside the DBD);  Fisher p = {p_fisher:.2g}")
""")
nb.md(f"""
**What changed?** {v('B8.n_tumours', '{:,}')} rows → tumours per residue, two binomial tests and one 2 × 2 table.
{v('missense_in_dbd', '{}')} of {v('n_missense', '{}')} missense variants ({v('B8.missense_frac_dbd', '{:.0%}')}) lie in
the DBD against {v('dbd_fraction', '{:.0%}')} expected (binomial p = {p_fmt('missense_binom_p')}); truncating variants
follow the length ({v('B8.truncating_frac_dbd', '{:.0%}')}, p = {p_fmt('truncating_binom_p')}). Fisher's test: odds
ratio {v('B8.fisher_odds_ratio', '{:.1f}')}, p = {p_fmt('fisher_p')}.

**Your turn.** Both tests say the same thing here. Why might you still prefer Fisher's test (missense vs truncating)
to the binomial test against 199/393?

<details><summary>Answer</summary>

The binomial test assumes that, without selection, every residue is **equally likely** to be hit. That is not true:
some DNA sites mutate more often (for example CpG sites, where many hotspot arginines sit), and the exons that code
for the DBD make up much of the coding sequence. Fisher's test compares two classes **within the same protein and the
same tumours**, so biases that affect both classes largely cancel. It still answers the biological question: are
missense variants placed differently from truncating ones?
</details>
""")

# ------------------------------------------------------------------------------------------------ B8.3
nb.md("""
## Step B8.3 · Check: compare with the truth and with a hand calculation  *(2 min)*

### Check

**Open the truth file only now.** It lists every variant of the simulation with its planted number of tumours. Your
counts per residue must match it exactly – a common mistake is to count per `protein_change` (R248Q and R248W
separately) instead of per residue. Then we redo the binomial test by hand with the normal approximation:
z = (observed − expected) ÷ √(n p₀ (1 − p₀)).
""")
nb.code("""
# CHECK · counts per residue equal the planted truth, and both tests point the expected way
truth = pd.read_csv(DATA / "truth" / "tp53_truth.tsv", sep="\\t")
planted = truth.groupby("position")["n_tumours"].sum()
planted = planted[planted > 0]                                # residues with at least one tumour
dviz.check("B8.counts_match_truth", bool(mine.sort_index().equals(planted.sort_index())), True)
dviz.check("B8.top_residue", int(mine.idxmax()), int(planted.idxmax()))
dviz.check("B8.n_residues_30plus", int((mine >= 30).sum()), int((planted >= 30).sum()))
dviz.check("B8.missense_frac_dbd", tests["missense"]["fraction"], between=(0.85, 0.95))
dviz.check("B8.missense_binom_p", tests["missense"]["p"], between=(0, 1e-50))
dviz.check("B8.truncating_frac_dbd", tests["truncating"]["fraction"], between=(0.45, 0.62))
dviz.check("B8.truncating_binom_p", tests["truncating"]["p"], between=(0.05, 1))
dviz.check("B8.fisher_odds_ratio", float(odds), between=(3, 20))
""")
nb.code("""
# CHECK · the binomial test by hand: how many standard errors is the observed count from the expected one?
z = {}
for cls, r in tests.items():
    expected = r["n"] * p_dbd
    z[cls] = (r["inside"] - expected) / np.sqrt(r["n"] * p_dbd * (1 - p_dbd))
    print(f"{cls:10s} expected {expected:5.1f} inside, observed {r['inside']}  ->  z = {z[cls]:+.1f}")
dviz.check("B8.missense_z", float(z["missense"]), between=(10, 40))
dviz.check("B8.truncating_z", float(z["truncating"]), between=(-1.96, 1.96), note="|z| < 1.96 means p > 0.05")
""")
nb.md(f"""
**What changed?** Nothing new was computed – we confirmed it. Your counts per residue equal the planted counts; the
most frequently hit residue is {v('B8.top_residue', '{}')}; {v('B8.n_residues_30plus', '{}')} residues have ≥ 30 tumours.
By hand, missense variants lie z = {v('B8.missense_z', '{:+.1f}')} standard errors above the expectation (far beyond any
threshold) and truncating variants z = {v('B8.truncating_z', '{:+.1f}')} (within ± 1.96, so p > 0.05).

**Your turn.** Count per `protein_change` instead of per `position`
(`var.groupby("protein_change").size().sort_values(ascending=False).head()`). Which residues move up when you count
per residue, and why does the lollipop use residues?

<details><summary>Answer</summary>

Per allele the top five are {hot_list('top_changes')}; per residue they are {hot_list('top_residues')}. Residues
such as R248 and R273 combine two common alleles (Q/W, H/C), so their counts nearly double and R273 moves ahead of R282. The lollipop is a map of *where* the protein
is hit, so one stem per residue is the honest unit; the alleles go into the label.
</details>
""")

# ------------------------------------------------------------------------------------------------ B8.4
nb.md("""
## Step B8.4 · Visualize: the lollipop plot  *(4 min)*

### Visualize

**Why this step.** A lollipop plot shows the answer and the evidence together: position along the protein on the
x-axis (with the domains drawn as boxes, like a ruler), one stem per residue and class with height = tumours,
dot size = tumours, and colour = variant class (Okabe–Ito colours, safe for colour-blind readers). The tested region is
shaded, and the title states n and the test. Only residues with ≥ 30 tumours get a label, so labels never overlap.
""")
nb.code("""
# Step B8.4 · labels for residues with >= 30 tumours: residue + alleles seen in >= 10 tumours (e.g. R248Q/W)
hot_labels = {}
for pos in mine[mine >= 30].index:
    a = top.loc[top["position"] == pos]                          # sorted: most frequent allele first
    keep = a.loc[a["n"] >= 10, "protein_change"].tolist() or a["protein_change"].tolist()[:1]
    hot_labels[pos] = keep[0][:-1] + "/".join(c[-1] for c in keep)
hot_labels
""")
nb.code("""
# Step B8.4 · lollipop: stems, dots sized by tumours, colour = variant class (Okabe-Ito), domains as boxes
COL = {"missense": dviz.OI_BLUE, "truncating": dviz.VERMILION, "splice/other": dviz.GREEN}
fig, ax = plt.subplots(figsize=(7.2, 3.3), layout="constrained")
ax.axvspan(dbd["aa_start"] - 0.5, dbd["aa_end"] + 0.5, color=dviz.HIGHLIGHT, zorder=0)       # the tested region
for _, r in per_pos.sort_values("n", ascending=False).iterrows():                # big dots first, small ones on top
    ax.vlines(r["position"], 0, r["n"], color=dviz.GREY, lw=0.6, zorder=1)
    ax.scatter(r["position"], r["n"], s=8 + 2.5 * r["n"], color=COL[r["variant_class"]], edgecolor="white", lw=0.4, zorder=2)
ax.hlines(-7, 1, LENGTH, color=dviz.GREY, lw=3, zorder=1)                       # the protein backbone, 1-393
for _, d in dom.iterrows():
    ax.add_patch(plt.Rectangle((d["aa_start"] - 0.5, -11), d["length_aa"], 8, color=dviz.SLATE, lw=0, zorder=2))
    ax.text((d["aa_start"] + d["aa_end"]) / 2, -7, d["short"], ha="center", va="center", fontsize=7, color="white", zorder=3)
for pos, label in hot_labels.items():
    ax.annotate(label, (pos, mine[pos]), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=7, color=dviz.INK)
for cls, colour in COL.items():
    ax.scatter([], [], color=colour, s=30, label=f"{cls} (n = {int((var['variant_class'] == cls).sum())})")
ax.legend(loc="upper right", fontsize=7, title="Variant class · dot size = tumours", title_fontsize=7)
ax.set(xlim=(0, LENGTH + 5), ylim=(-14, mine.max() * 1.2), xlabel="Amino-acid position (TP53, 393 aa)", ylabel="Tumours")
ax.set_title(f"TP53 variants in {len(var):,} tumours (simulated) · shaded: DNA-binding domain ({p_dbd:.0%} of the length)\\n"
             f"inside it: {tests['missense']['fraction']:.0%} of missense (binomial p = {tests['missense']['p']:.0e}), "
             f"{tests['truncating']['fraction']:.0%} of truncating (p = {tests['truncating']['p']:.2f})", fontsize=8)
dviz.savefig(fig, "B8_lollipop")
plt.show()
""")
nb.code("""
# CHECK · every residue with >= 30 tumours has a label, and the plot shows every tumour once
dviz.check("B8.n_labels", len(hot_labels), int((planted >= 30).sum()))
dviz.check("B8.tumours_in_plot", int(per_pos["n"].sum()), 1200)
""")
nb.md(f"""
**What changed?** A table of {v('B8.n_tumours', '{:,}')} variants → one picture: the missense hotspots
({hot_list('top_residues')}) stand in the shaded DNA-binding domain, while truncating and splice variants are low
stems spread along the whole protein.

**Your turn.** Change the label threshold from 30 to 15 (in the label cell: `mine[mine >= 15]`) and re-run both cells.
How many labels do you get, and what goes wrong? Then set it back to 30.

<details><summary>Answer</summary>

{v('n_residues_15plus', '{}')} residues have ≥ 15 tumours, so {v('n_residues_15plus', '{}')} labels. Around residues
245–249 and 273–282 the labels start to collide, and the eye is drawn to text instead of the stems. Label only what you
discuss; put the full list in a table or the legend.
</details>
""")

# ------------------------------------------------------------------------------------------------ wrap-up
nb.md(f"""
## Summary

| Step | What you did | Number that proves it |
|---|---|---|
| B8.1 Question | loaded variants and domains | {v('B8.n_tumours', '{:,}')} tumours; {v('pct_missense', '{:.0f}')}% missense; DBD = {v('B8.dbd_length', '{}')}/393 = {v('dbd_fraction', '{:.3f}')} |
| B8.2 Test | binomial test per class; Fisher's exact test | missense {v('B8.missense_frac_dbd', '{:.0%}')} inside (p = {p_fmt('missense_binom_p')}); truncating {v('B8.truncating_frac_dbd', '{:.0%}')} (p = {p_fmt('truncating_binom_p')}); OR {v('B8.fisher_odds_ratio', '{:.1f}')} |
| B8.3 Check | counts vs truth; z by hand | counts match; z = {v('B8.missense_z', '{:+.1f}')} and {v('B8.truncating_z', '{:+.1f}')} |
| B8.4 Visualize | lollipop with domains | {v('B8.n_labels', '{}')} labelled hotspots: {hot_list('top_residues')} |

**Presentation sentence (fill in your numbers).**
*Question* – Do missense and truncating TP53 variants fall in different parts of the protein?
*Takeaway* – {v('B8.missense_frac_dbd', '{:.0%}')} of missense variants lie in the DNA-binding domain, which is
{v('dbd_fraction', '{:.0%}')} of the protein (binomial p = {p_fmt('missense_binom_p')}), with hotspots at R175, R248 and
R273; truncating variants follow the length ({v('B8.truncating_frac_dbd', '{:.0%}')}, p = {p_fmt('truncating_binom_p')}).
*Limitation* – counts are simulated, one variant per tumour; the expectation assumes every residue is equally mutable,
and domain boundaries differ by a few residues between databases.
""")
nb.code("""
# Step B8 · collect the numbers for the slides and the instructor: (value, unit, meaning)
cls_pct = var["variant_class"].value_counts() / len(var) * 100
per_allele = var.groupby("protein_change").size().sort_values(ascending=False, kind="stable")
top5 = mine.sort_values(ascending=False, kind="stable").head(5)
slide = {
    "n_tumours": (len(var), "tumours", "one TP53 variant each"),
    "pct_missense": (round(float(cls_pct["missense"]), 1), "%", "missense share"),
    "pct_truncating": (round(float(cls_pct["truncating"]), 1), "%", "nonsense + frameshift share"),
    "pct_splice": (round(float(cls_pct["splice/other"]), 1), "%", "splice / in-frame share"),
    "dbd_fraction": (round(p_dbd, 4), "fraction", "199/393, expected share inside the DBD"),
    "missense_in_dbd": (tests["missense"]["inside"], "variants", "missense inside the DBD"),
    "n_missense": (tests["missense"]["n"], "variants", "all missense"),
    "missense_binom_p": (tests["missense"]["p"], "p", "binomial test, missense inside DBD vs 199/393"),
    "truncating_in_dbd": (tests["truncating"]["inside"], "variants", "truncating inside the DBD"),
    "n_truncating": (tests["truncating"]["n"], "variants", "all truncating"),
    "truncating_binom_p": (tests["truncating"]["p"], "p", "binomial test, truncating inside DBD vs 199/393"),
    "fisher_p": (float(p_fisher), "p", "Fisher's exact test, missense vs truncating x inside/outside"),
    "n_residues_15plus": (int((mine >= 15).sum()), "residues", "residues with >= 15 tumours"),
}
""")
nb.code("""
# Step B8 · add the hotspot tables, then save everything to results/B8.json
slide["top_residues"] = ([{"label": hot_labels.get(p, str(p)), "position": int(p), "n": int(k)} for p, k in top5.items()],
                         "tumours", "the five most frequently hit residues (all classes)")
slide["top_changes"] = ([{"label": c, "n": int(k)} for c, k in per_allele.head(5).items()],
                        "tumours", "the five most frequent alleles")
extra = {"slide_numbers": {k: {"value": a, "unit": b, "meaning": c} for k, (a, b, c) in slide.items()}}
dviz.save_results("B8", extra=extra)
print(f"Notebook time: {time.time() - T_START:.0f} s")
""")

for c in nb.cells:                                   # notebook conventions: <= 25 lines, first line "# Step" or "# CHECK ·"
    if c.cell_type == "code":
        lines = c.source.splitlines()
        assert len(lines) <= 25, f"{len(lines)} lines: {lines[0]}"
        assert lines[0].startswith(("# Step", "# CHECK ·")), lines[0]
        assert not any(l.lstrip().startswith(("%", "!")) for l in lines), "no magics"

out = ROOT / "electives" / "B8_tp53_lollipop.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
nb.write(str(out))
print("wrote", out.relative_to(ROOT), "| numbers from previous run:", bool(RES))
