"""Build electives/B6_variant_callers.ipynb (practice B6 · Variant callers: agreement is not truth).

Numbers quoted in the markdown are read from data/truth/callsets_truth.json and data/variant_callsets.tsv
(written by build/make_callsets.py), so the text always matches what the CHECK cells verify.
usage:  python build/build_B6_notebook.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
from nbbuild import NB  # noqa: E402
import dviz  # noqa: E402

T = json.loads((ROOT / "data" / "truth" / "callsets_truth.json").read_text())
VC = pd.read_csv(ROOT / "data" / "variant_callsets.tsv", sep="\t")
CALLERS = ["GATK", "DeepVariant", "FreeBayes", "Strelka2"]
PC = T["per_caller"]
n_caller = VC[CALLERS].sum(axis=1)
single_fp = VC[(n_caller == 1) & (VC.in_truth_set == 0)]["context"].value_counts()
single_fp_share = (100 * single_fp / single_fp.sum()).round(0)
truth_share = (100 * VC[VC.in_truth_set == 1]["context"].value_counts() / T["n_truth"]).round(0)
up_all = dviz.upset_counts(VC, CALLERS)
fb = PC["FreeBayes"]["indel"]
best_p = max(CALLERS, key=lambda k: PC[k]["all"]["precision"])
cons = T["consensus"]
a4 = T["all_four"]
sfp = T["shared_fp_by_context"]
seg_share = 100 * sfp["segdup_lowmap"] / T["n_shared_fp"]
seg_truth = 100 * T["n_truth_by_context"]["segdup_lowmap"] / T["n_truth"]
up = {u["callers"]: u for u in T["upset_top10"]}
fbo = T["single_caller"]["per_caller"]["FreeBayes"]
dvo = T["deepvariant_only"]
vote_rules = {">=1": "≥ 1", ">=2": "≥ 2", ">=3": "≥ 3", "=4": "= 4"}
best_rule = max(vote_rules, key=lambda k: cons[k]["F1"])
dv_f1 = PC["DeepVariant"]["all"]["F1"]
missed_unique = T["missed_by_all_by_context"]["unique"]


def p3(x):
    return f"{x:.3f}"


snv_rec = [PC[k]["SNV"]["recall"] for k in CALLERS]
ind_rec = {k: PC[k]["indel"]["recall"] for k in CALLERS}

nb = NB()
# ------------------------------------------------------------------------------------------------ title
nb.md(f"""
# Practice B6 · Variant callers: agreement is not truth

**Truth set → precision, recall and F1 per caller → UpSet of the four call sets → consensus voting → where the shared errors are**

*Genomic Data Processing · Day 5 (10/27) · Data visualization in genetics and genomics · about 18 minutes*

> **All data in this practice are simulated (computer-generated). No real patient data are used.**
> The positions are real hg38 coordinates in a 6-Mb window of chromosome 20 ({T['window_1based']}), styled after a
> Genome in a Bottle (GIAB) benchmark, so the table looks exactly like a real benchmark; the alleles, sequence contexts
> and calls are simulated with a known answer.

**How to use this notebook.** Run the cells from top to bottom (Shift + Enter). Every step follows the same four beats:
**Question → Test → Check → Visualize**. A cell that starts with `# CHECK ·` compares your number with the reference
value and prints ✓ or ✗; if you see ✗, stop and find out why before going on. Each step ends with **What changed?**
and a short **Your turn** task; the answer is folded (click *Answer*). Prompt card 7 (`prompts/prompt_cards.md`) is one
way to ask an AI assistant for the plotting code; the provided cells are the other. Either way, the CHECK cells decide.

| File | What is inside |
|---|---|
| `data/variant_callsets.tsv` | one row per variant that is in the truth set **or** was called by at least one caller ({T['n_rows']:,} rows) |
| `data/truth/callsets_truth.json` | the planted design and reference values – open it only in the CHECK cells |

| Step | Topic | min |
|---|---|---|
| B6.1 | the yardstick: truth set, precision, recall, F1 | 3 |
| B6.2 | each caller against the truth set, SNVs and indels | 4 |
| B6.3 | UpSet plot of the four call sets | 4 |
| B6.4 | agreement is not truth: consensus voting | 4 |
| B6.5 | where the shared errors are | 3 |
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
CALLERS = ["GATK", "DeepVariant", "FreeBayes", "Strelka2"]              # the four call sets (0/1 columns)
CONTEXTS = ["unique", "homopolymer_tr", "segdup_lowmap", "high_gc"]      # sequence context of each position
""")

# ------------------------------------------------------------------------------------------------ B6.1
nb.md("""
## Step B6.1 · The yardstick: truth set, precision, recall and F1  *(3 min)*

### Question

Four programs were run on the same sequenced genome and each produced its own list of variants. **Which list should we
trust, and is a variant more likely to be real when several programs agree?**

**Why this step.** A *variant caller* reads the sequencing reads that were aligned to the reference genome and decides,
position by position, whether this person's DNA differs from the reference. Callers use different statistics and make
different mistakes. To say who is right we need an independent answer key: a **truth set**. The Genome in a Bottle
(GIAB) Consortium, led by NIST, sequenced reference genomes such as HG002 with many technologies and publishes
high-confidence variants together with the **high-confidence regions** where those variants can be trusted. A
benchmark compares a caller's list with the truth set **inside these regions only**.

**New terms**

* **Variant caller** – program that turns aligned reads (BAM) into a list of variants (VCF). GATK HaplotypeCaller,
  DeepVariant, FreeBayes and Strelka2 are four widely used callers.
* **SNV / indel** – single-nucleotide variant (one base changed) / small insertion or deletion of bases.
* **TP, FP, FN** – *true positive*: called and in the truth set; *false positive*: called but not in the truth set;
  *false negative*: in the truth set but not called.
* **Precision** = TP / (TP + FP) – of the variants a caller reports, the fraction that are real.
* **Recall** = TP / (TP + FN) – of the real variants, the fraction the caller found (also called sensitivity).
* **F1** = 2 × precision × recall / (precision + recall) – one number that is high only when both are high.

A worked example with 10 positions and one caller:

| site | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| in the truth set | yes | yes | yes | yes | yes | yes | yes | no | no | no |
| called | yes | yes | yes | yes | yes | no | no | yes | no | no |
| outcome | TP | TP | TP | TP | TP | FN | FN | FP | – | – |

TP = 5, FP = 1, FN = 2, so precision = 5/6 = 0.833, recall = 5/7 = 0.714 and F1 = 0.769. Sites 9 and 10 are true
negatives: they are not used.

### Test
""")
nb.code("""
# Step B6.1 · precision, recall and F1 from two 0/1 columns (this helper is used in every later step)
def prf(called, truth):
    \"\"\"Compare calls with the truth set. Both inputs: one 0/1 (or True/False) value per variant.\"\"\"
    called = np.asarray(called).astype(bool)      # 0/1 -> False/True (do not index a table with 0/1 integers)
    truth = np.asarray(truth).astype(bool)
    tp = int(np.sum(called & truth))              # called and real
    fp = int(np.sum(called & ~truth))             # called but not real
    fn = int(np.sum(~called & truth))             # real but missed
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * precision * recall / (precision + recall)
    return {"TP": tp, "FP": fp, "FN": fn, "precision": precision, "recall": recall, "F1": f1}

toy = pd.DataFrame({"site": range(1, 11),
                    "in_truth_set": [1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
                    "called":       [1, 1, 1, 1, 1, 0, 0, 1, 0, 0]})
toy_result = prf(toy["called"], toy["in_truth_set"])
print({k: round(v, 3) for k, v in toy_result.items()})
""")
nb.md("### Check")
nb.code("""
# CHECK · the worked example: TP 5, FP 1, FN 2
dviz.check("B6.toy_precision", round(toy_result["precision"], 3), 0.833, tol=0.0005)
dviz.check("B6.toy_recall", round(toy_result["recall"], 3), 0.714, tol=0.0005)
dviz.check("B6.toy_F1", round(toy_result["F1"], 3), 0.769, tol=0.0005)
""")
nb.md("""
**What changed?** Nothing in the data yet: you now have a yardstick. The helper `prf()` turns two 0/1 columns into
TP, FP, FN, precision, recall and F1, and it reproduced the hand calculation (0.833, 0.714, 0.769).

**Your turn.** Sites 9 and 10 are true negatives. Why do precision and recall ignore them, and why does nobody
report the *accuracy* of a variant caller?

<details><summary>Answer</summary>

In a genome almost every position is a true negative: a person differs from the reference at roughly one position in
a thousand. Accuracy = (TP + TN) / all positions would be above 99.9% even for a caller that reports nothing at all.
Precision and recall look only at variants (called or real), so they separate good callers from bad ones.
</details>
""")

# ------------------------------------------------------------------------------------------------ B6.2
nb.md(f"""
## Step B6.2 · Each caller against the truth set, SNVs and indels  *(4 min)*

### Question

How precise and how complete is each caller, and are indels harder than SNVs?

**Why this step.** All four callers saw the same reads from the same genome, so any difference between their lists
comes from the callers' methods. Indels are harder than SNVs: a missing base in a run such as `AAAAAAA` gives reads that
fit several alignments, and sequencing errors in such runs look like real indels. We therefore score SNVs and indels
separately.

**New terms**

* **Benchmark table** – one row per variant that is in the truth set *or* was called; the columns `GATK`,
  `DeepVariant`, `FreeBayes`, `Strelka2` are 1 if that caller reported it, `in_truth_set` is 1 if it is real.
* **Truth variants missed by every caller** – rows with `in_truth_set = 1` and four zeros. They must stay in the table,
  otherwise recall is overestimated.

### Test
""")
nb.code("""
# Step B6.2 · load the benchmark (positions are 1-based, hg38) and look at the first rows
vc = pd.read_csv(DATA / "variant_callsets.tsv", sep="\\t")
n_truth = int(vc["in_truth_set"].sum())
print(f"{len(vc):,} rows: {n_truth:,} truth variants and {len(vc) - n_truth:,} calls that are not in the truth set")
print("truth variants by type:", vc.loc[vc["in_truth_set"] == 1, "type"].value_counts().to_dict())
print("truth variants missed by every caller:", int(((vc[CALLERS].sum(axis=1) == 0) & (vc["in_truth_set"] == 1)).sum()))
vc.head()
""")
nb.code("""
# Step B6.2 · precision, recall and F1 for every caller: all variants, SNVs only, indels only
rows = []
for caller in CALLERS:
    for vtype in ["all", "SNV", "indel"]:
        sub = vc if vtype == "all" else vc[vc["type"] == vtype]
        r = prf(sub[caller], sub["in_truth_set"])
        rows.append({"caller": caller, "type": vtype, "calls": r["TP"] + r["FP"], **r})
metrics = pd.DataFrame(rows)
print(metrics[metrics["type"] == "all"].round(4).to_string(index=False))
metrics.pivot(index="caller", columns="type", values="recall").loc[CALLERS, ["all", "SNV", "indel"]].round(3)
""")
nb.md("### Check")
nb.code("""
# CHECK · your table against the reference values in the truth file (open it only in CHECK cells)
ref = json.loads((DATA / "truth" / "callsets_truth.json").read_text())
dviz.check("B6.n_truth", n_truth, ref["n_truth"])
diffs = [abs(r[m] - ref["per_caller"][r["caller"]][r["type"]][m]) for _, r in metrics.iterrows()
         for m in ["precision", "recall", "F1"]]
dviz.check("B6.max_diff_vs_reference", max(diffs), 0, tol=1e-5)
allv = metrics[metrics["type"] == "all"].set_index("caller")
dviz.check("B6.highest_precision", allv["precision"].idxmax(), "DeepVariant")
dviz.check("B6.lowest_precision", allv["precision"].idxmin(), "FreeBayes")
snv = metrics[metrics["type"] == "SNV"].set_index("caller")
ind = metrics[metrics["type"] == "indel"].set_index("caller")
dviz.check("B6.indel_recall_lower_for_every_caller", bool((ind["recall"] < snv["recall"]).all()), True)
""")
nb.md("### Visualize")
nb.code("""
# Step B6.2 · figure: precision, recall and F1 per caller for SNVs and indels (dot plot)
STYLE = {"SNV": dict(color=dviz.OI_BLUE, marker="o"), "indel": dict(color=dviz.ORANGE, marker="D")}
DODGE = {"SNV": 0.13, "indel": -0.13}                              # small vertical offset: dots do not hide each other
order = allv.sort_values("F1").index.tolist()                      # best caller at the top
fig, axs = plt.subplots(1, 3, figsize=(7.4, 3.0), sharey=True, layout="constrained")
for ax, m in zip(axs, ["precision", "recall", "F1"]):
    for vtype, style in STYLE.items():
        sub = metrics[metrics["type"] == vtype].set_index("caller").loc[order]
        n = int(vc.loc[vc["in_truth_set"] == 1, "type"].eq(vtype).sum())
        ax.scatter(sub[m], np.arange(len(order)) + DODGE[vtype], s=40, zorder=2,
                   label=f"{vtype} (n = {n:,} true variants)", **style)
    ax.set_title(m.capitalize() if m != "F1" else "F1")
    ax.set_xlim(0.80, 1.0)
    ax.grid(axis="x", color=dviz.LIGHT, lw=0.8, zorder=0)
axs[0].set_yticks(range(len(order)), order)
axs[1].set_xlabel("Fraction (axis starts at 0.80: dots, not bars)")
fig.legend(*axs[0].get_legend_handles_labels(), loc="outside lower center", ncol=2)
fig.suptitle("Four callers vs a GIAB-style truth set, chr20 10-16 Mb (hg38)", x=0.01, ha="left",
             fontweight="bold", color=dviz.INK)
dviz.savefig(fig, "B6_precision_recall")
plt.show()
""")
nb.md(f"""
**What changed?** The {T['n_rows']:,}-row benchmark became a 12-row table (4 callers × all / SNV / indel).
{best_p} has the highest precision ({p3(PC['DeepVariant']['all']['precision'])}) and recall
({p3(PC['DeepVariant']['all']['recall'])}); FreeBayes the lowest precision ({p3(PC['FreeBayes']['all']['precision'])}).
Every caller is worse on indels: indel recall is {p3(ind_rec['FreeBayes'])} (FreeBayes) to
{p3(ind_rec['DeepVariant'])} (DeepVariant), against {p3(min(snv_rec))}–{p3(max(snv_rec))} for SNVs. Your table matches
the reference to within 1e-5.

**Your turn.** Find FreeBayes' indel row (`metrics[(metrics.caller == "FreeBayes") & (metrics.type == "indel")]`).
How many of its indel calls are false, and what is its indel precision?

<details><summary>Answer</summary>

FreeBayes reported {fb['calls']:,} indels: {fb['TP']:,} real and **{fb['FP']} false**, so indel precision =
{fb['TP']:,} / {fb['calls']:,} = {p3(fb['precision'])}. Most of these false indels sit in homopolymers and short tandem
repeats, where a sequencing error looks like one extra or one missing base.
</details>
""")

# ------------------------------------------------------------------------------------------------ B6.3
nb.md(f"""
## Step B6.3 · UpSet plot of the four call sets  *(4 min)*

### Question

How do the four call sets overlap, and which overlap is the largest?

**Why this step.** With four sets there are 15 possible combinations ("called by GATK and DeepVariant only", "called
by all four", …). A four-set Venn diagram draws them as 15 regions of overlapping ellipses: the regions have
different shapes and sizes, and finding one number means searching the picture. An **UpSet plot** (Lex et al., 2014)
turns every combination into a bar on one common axis, sorted by size; a dot matrix below each bar shows which sets it
contains.

**New terms**

* **Intersection (exclusive)** – the variants called by exactly this combination of callers and by no other caller.
  Every variant belongs to exactly one intersection, so the bars add up to all calls.
* **Set size** – all calls of one caller (the sum of the bars that have a dot in that caller's row).

### Test
""")
nb.code("""
# Step B6.3 · count the variants in every combination of callers (exclusive intersections, largest first)
counts = dviz.upset_counts(vc, CALLERS)
counts["callers"] = [" + ".join(c for c in CALLERS if row[c]) for _, row in counts.iterrows()]
counts.loc[counts[CALLERS].all(axis=1), "callers"] = "all four"
print(f"{len(counts)} intersections, {int(counts['count'].sum()):,} called variants")
print(counts[["callers", "count"]].to_string(index=False))
""")
nb.md("### Check")
nb.code("""
# CHECK · the largest bar is "all four callers"; the counts match the reference
dviz.check("B6.largest_bar", counts.loc[0, "callers"], "all four")
dviz.check("B6.n_all_four", int(counts.loc[0, "count"]), ref["all_four"]["n"])
dv_only = int(counts.loc[counts["callers"] == "DeepVariant", "count"].iloc[0])
dviz.check("B6.n_deepvariant_only", dv_only, ref["deepvariant_only"]["n"])
dviz.check("B6.bars_add_up_to_all_calls", int(counts["count"].sum()), int((vc[CALLERS].sum(axis=1) > 0).sum()))
""")
nb.md("### Visualize")
nb.code("""
# Step B6.3 · UpSet plot of the 10 largest intersections (dviz.upset) and a small layout fix
fig, top = dviz.upset(vc, CALLERS, max_bars=10)
ax_bar, ax_dot, ax_set = fig.axes                                   # bars, dot matrix, set sizes
ax_bar.set_ylabel("Variants called")
fig.suptitle(f"Four call sets, chr20 10-16 Mb (hg38): the 10 largest of {len(counts)} intersections",
             x=0.02, ha="left", fontweight="bold", color=dviz.INK)
dviz.savefig(fig, "B6_upset")
plt.show()
""")
nb.md(f"""
**What changed?** The four 0/1 columns became {len(up_all)} exclusive intersections. The largest bar is **all four
callers** ({a4['n']:,} variants), followed by FreeBayes alone ({up['FreeBayes']['count']}) and GATK + DeepVariant +
Strelka2 ({up['GATK+DeepVariant+Strelka2']['count']}); DeepVariant alone has {dvo['n']} calls. The bars add up to all
{int((n_caller > 0).sum()):,} called variants. On the slide, the same numbers are hard to find in a four-set Venn.

**Your turn.** Without looking at the truth column: which bar would you trust least, *FreeBayes only*
({up['FreeBayes']['count']} calls) or *DeepVariant only* ({dvo['n']} calls)? Write down a guess; Step B6.4 gives the answer.

<details><summary>Answer</summary>

FreeBayes-only calls are almost all false ({fbo['TP']} of {fbo['n']} are real), while DeepVariant-only calls are mostly
real ({dvo['TP']} of {dvo['n']}): DeepVariant finds hard true variants (indels in repeats, variants in duplicated
sequence) that the others miss. The UpSet plot shows *how many* callers agree, not *whether they are right* – that
needs the truth set.
</details>
""")

# ------------------------------------------------------------------------------------------------ B6.4
nb.md(f"""
## Step B6.4 · Agreement is not truth: consensus voting  *(4 min)*

### Question

If we keep only variants that several callers agree on, do we get a better list?

**Why this step.** A popular shortcut is *consensus voting*: keep a variant only if at least *k* of the four callers
reported it. This removes errors that one caller makes alone, but it **cannot** remove errors that all callers share
(they saw the same misleading reads), and it throws away real variants that only one caller found.

**New terms**

* **Consensus rule "≥ k"** – keep a variant if at least *k* callers called it; "≥ 1" is the union of all four lists,
  "= 4" is the intersection.
* **Shared error** – a false call made by several callers at the same position with the same allele.

### Test
""")
nb.code("""
# Step B6.4 · how often is a call real when 1, 2, 3 or 4 callers made it?
vc["n_callers"] = vc[CALLERS].sum(axis=1)                    # 0 = a truth variant that no caller found
by_n = vc[vc["n_callers"] > 0].groupby("n_callers")["in_truth_set"].agg(calls="size", real="sum").reset_index()
by_n["false"] = by_n["calls"] - by_n["real"]
by_n["precision"] = (by_n["real"] / by_n["calls"]).round(3)
print(by_n.to_string(index=False))
single = vc[vc["n_callers"] == 1]                             # calls made by ONE caller only
only = pd.DataFrame({"calls": [int(single[c].sum()) for c in CALLERS],
                     "real": [int(single.loc[single[c] == 1, "in_truth_set"].sum()) for c in CALLERS]},
                    index=CALLERS)
only["precision"] = (only["real"] / only["calls"]).round(3)
print(only.to_string())
""")
nb.code("""
# Step B6.4 · consensus rules: keep a variant if at least k callers called it (DeepVariant alone for comparison)
rules = {"≥ 1 (union)": vc["n_callers"] >= 1, "≥ 2": vc["n_callers"] >= 2, "≥ 3": vc["n_callers"] >= 3,
         "= 4 (all four)": vc["n_callers"] == 4, "DeepVariant alone": vc["DeepVariant"] == 1}
consensus = pd.DataFrame({name: prf(mask, vc["in_truth_set"]) for name, mask in rules.items()}).T
consensus = consensus.astype({"TP": int, "FP": int, "FN": int})
consensus.insert(0, "calls", consensus["TP"] + consensus["FP"])
consensus.round(4)
""")
nb.md("### Check")
nb.code("""
# CHECK · agreement raises precision, but false calls survive in "all four" and recall drops
c4 = consensus.loc["= 4 (all four)"]
dviz.check("B6.precision_all_four", round(float(c4["precision"]), 4), round(ref["all_four"]["precision"], 4), tol=1e-4)
dviz.check("B6.false_calls_all_four", int(c4["FP"]), ref["all_four"]["FP"])
p_single = float(by_n.loc[by_n["n_callers"] == 1, "precision"].iloc[0])
dviz.check("B6.precision_single_caller", p_single, round(ref["single_caller"]["precision"], 3), tol=5e-4)
dviz.check("B6.all_four_more_precise_than_single", bool(c4["precision"] > p_single), True)
dviz.check("B6.all_four_recall_below_DeepVariant",
           bool(c4["recall"] < consensus.loc["DeepVariant alone", "recall"]), True)
""")
nb.md("### Visualize")
nb.code("""
# Step B6.4 · figure: precision and recall of each consensus rule; dotted = DeepVariant alone
names = ["≥ 1 (union)", "≥ 2", "≥ 3", "= 4 (all four)"]
x = np.arange(len(names))
dv = consensus.loc["DeepVariant alone"]
fig, ax = plt.subplots(figsize=(5.4, 3.4))
ax.plot(x, consensus.loc[names, "precision"], "-o", color=dviz.BLUE, label="precision")
ax.plot(x, consensus.loc[names, "recall"], "-o", color=dviz.SLATE, label="recall")
ax.axhline(dv["precision"], color=dviz.BLUE, ls=":", lw=1.2, label=f"DeepVariant alone: precision {dv['precision']:.3f}")
ax.axhline(dv["recall"], color=dviz.SLATE, ls=":", lw=1.2, label=f"DeepVariant alone: recall {dv['recall']:.3f}")
for xi, name in zip(x, names):                                    # higher value labelled above, lower below
    p, r = consensus.loc[name, "precision"], consensus.loc[name, "recall"]
    for v, col, dy in [(p, dviz.BLUE, 7 if p >= r else -16), (r, dviz.SLATE, -16 if p >= r else 7)]:
        ax.annotate(f"{v:.3f}", (xi, v), xytext=(0, dy), textcoords="offset points", ha="center",
                    fontsize=7.5, color=col)
ax.set_xticks(x, ["≥ 1\\n(union)", "≥ 2", "≥ 3", "= 4\\n(all four)"])
ax.set(xlabel="Keep a variant if this many callers called it", ylabel="Fraction (axis starts at 0.88)",
       ylim=(0.88, 1.012), title=f"Consensus voting vs the truth set (n = {n_truth:,} true variants)")
ax.legend(loc="lower center", bbox_to_anchor=(0.45, 0.0), fontsize=7.5)
dviz.savefig(fig, "B6_consensus")
plt.show()
""")
nb.md(f"""
**What changed?** We added one column, `n_callers`, and turned it into precision and recall per rule. Calls made by a
single caller are mostly false (precision {p3(T['single_caller']['precision'])}), calls made by all four are mostly real
(precision {p3(a4['precision'])}) – but the all-four set still contains **{a4['FP']} false calls**, and requiring all
four callers loses {cons['=4']['FN']:,} true variants (recall {p3(cons['=4']['recall'])}). Among the voting rules, {vote_rules[best_rule]}
has the best F1 ({p3(cons[best_rule]['F1'])}) – and DeepVariant alone is {'even slightly better' if dv_f1 > cons[best_rule]['F1'] else 'about as good'} (F1 {p3(dv_f1)}).

**Your turn.** You must choose one rule for (a) a clinical report, where a false variant is worse than a missed one, and
(b) a research screen that should miss as little as possible. Which rules do you pick, and what would you still check?

<details><summary>Answer</summary>

(a) "= 4" or "≥ 3" (precision {p3(cons['=4']['precision'])} / {p3(cons['>=3']['precision'])}), but they still contain
false calls that all callers share, so a clinical variant is confirmed with an independent method (for example Sanger
sequencing or long reads), and we should know *where* those shared errors are (next step). (b) "≥ 1" (recall
{p3(cons['>=1']['recall'])}) costs {cons['>=1']['FP']} false calls; "≥ 2" keeps recall at {p3(cons['>=2']['recall'])}
with far fewer false calls. Agreement is evidence, not proof: only the truth set measures correctness.
</details>
""")

# ------------------------------------------------------------------------------------------------ B6.5
nb.md(f"""
## Step B6.5 · Where are the shared errors?  *(3 min)*

### Question

In which kind of sequence do several callers make the same mistake?

**Why this step.** Callers share an error when they all see the same misleading reads. In a **segmental duplication**
(a long stretch of DNA with a near-identical copy elsewhere in the genome), reads from the other copy can be mapped
here; the small differences between the two copies then look like a real heterozygous SNV to every caller. Knowing
where shared errors concentrate tells us which calls to confirm.

**New terms**

* **Homopolymer / tandem repeat** – a run of one base (`AAAAAAA`) or a short unit repeated (`CACACA`); polymerase
  slippage and sequencing errors create false indels here.
* **Segmental duplication / low mappability** – sequence ≥ 1 kb that is ≥ 90% identical to another place in the
  genome; reads cannot be placed uniquely (low mapping quality).
* **High GC** – G + C-rich sequence; it is covered unevenly by short-read sequencing.
* **Enrichment** – share of the errors in a context divided by its share of the truth variants; 1 = as expected.

### Test
""")
nb.code("""
# Step B6.5 · contexts of the truth variants, of the shared false calls (>= 3 callers) and of the missed variants
truth_v = vc[vc["in_truth_set"] == 1]
shared_fp = vc[(vc["in_truth_set"] == 0) & (vc["n_callers"] >= 3)]
ctx = pd.DataFrame({"truth_variants": truth_v["context"].value_counts(),
                    "shared_false_calls": shared_fp["context"].value_counts(),
                    "missed_by_all": truth_v.loc[truth_v["n_callers"] == 0, "context"].value_counts()})
ctx = ctx.reindex(CONTEXTS).fillna(0).astype(int)
ctx["% of truth"] = (100 * ctx["truth_variants"] / ctx["truth_variants"].sum()).round(1)
ctx["% of shared false"] = (100 * ctx["shared_false_calls"] / ctx["shared_false_calls"].sum()).round(1)
ctx["enrichment"] = (ctx["% of shared false"] / ctx["% of truth"]).round(1)
ctx
""")
nb.md("### Check")
nb.code("""
# CHECK · the shared false calls concentrate in segmental duplications (the planted design)
dviz.check("B6.shared_fp_top_context", ctx["shared_false_calls"].idxmax(), ref["shared_fp_top_context_planted"])
dviz.check("B6.n_shared_false_calls", int(ctx["shared_false_calls"].sum()), ref["n_shared_fp"])
dviz.check("B6.segdup_pct_of_shared_false", float(ctx.loc["segdup_lowmap", "% of shared false"]), between=(50, 100))
dviz.check("B6.segdup_pct_of_truth", float(ctx.loc["segdup_lowmap", "% of truth"]), between=(3, 7))
dviz.check("B6.segdup_false_in_all_four", int(((vc["n_callers"] == 4) & (vc["in_truth_set"] == 0)
                                               & (vc["context"] == "segdup_lowmap")).sum()),
           ref["all_four"]["FP_by_context"]["segdup_lowmap"])
""")
nb.md("### Visualize")
nb.code("""
# Step B6.5 · figure: share of truth variants vs share of shared false calls, per sequence context
LABEL = {"unique": "unique", "homopolymer_tr": "homopolymer / tandem repeat",
         "segdup_lowmap": "segmental duplication", "high_gc": "high GC"}
y = np.arange(len(CONTEXTS))[::-1]                                  # first context at the top
fig, ax = plt.subplots(figsize=(6.2, 2.9))
ax.barh(y + 0.2, ctx["% of truth"], height=0.38, color=dviz.GREY,
        label=f"truth variants (n = {ctx['truth_variants'].sum():,})")
ax.barh(y - 0.2, ctx["% of shared false"], height=0.38, color=dviz.VERMILION,
        label=f"false calls made by ≥ 3 callers (n = {ctx['shared_false_calls'].sum()})")
for yi, a, b in zip(y, ctx["% of truth"], ctx["% of shared false"]):
    ax.text(a + 1, yi + 0.2, f"{a:.0f}%", va="center", fontsize=8, color=dviz.SLATE)
    ax.text(b + 1, yi - 0.2, f"{b:.0f}%", va="center", fontsize=8, color=dviz.VERMILION)
ax.set_yticks(y, [LABEL[c] for c in CONTEXTS])
ax.set(xlabel="% of variants in this context", xlim=(0, 100),
       title="Shared false calls sit in segmental duplications (chr20, hg38)")
ax.legend(loc="lower right", fontsize=7.5)
dviz.savefig(fig, "B6_contexts")
plt.show()
""")
nb.md(f"""
**What changed?** Grouping by sequence context shows where agreement fails: {sfp['segdup_lowmap']} of the
{T['n_shared_fp']} false calls made by ≥ 3 callers ({seg_share:.0f}%) lie in segmental duplications, which hold only
{seg_truth:.0f}% of the truth variants (about {seg_share / seg_truth:.0f}-fold enriched), and {a4['FP_by_context']['segdup_lowmap']}
of the {a4['FP']} false calls in the all-four set are there. {'The ' + str(T['n_truth_missed_by_all']) + ' truth variants that no caller found are also all in hard contexts (none in unique sequence).' if missed_unique == 0 else 'Of the ' + str(T['n_truth_missed_by_all']) + ' truth variants that no caller found, only ' + str(missed_unique) + ' lie in unique sequence.'}

**Your turn.** Change `>= 3` to `== 1` in the first line of the Test cell (false calls made by one caller only) and
run it again. Which context now holds the most errors, and is it still enriched?

<details><summary>Answer</summary>

Single-caller false calls are spread differently: unique sequence holds the most in absolute numbers
({single_fp_share.get('unique', 0):.0f}% of them) but less than its share of truth variants
({truth_share.get('unique', 0):.0f}%), while homopolymers and tandem repeats are enriched
({single_fp_share.get('homopolymer_tr', 0):.0f}% of the errors vs {truth_share.get('homopolymer_tr', 0):.0f}% of the
truth). Private errors are mostly random sequencing noise and repeat indels; shared errors come from mis-mapped reads in
duplicated sequence. (Set the line back to `>= 3` afterwards.)
</details>
""")

# ------------------------------------------------------------------------------------------------ summary
nb.md(f"""
## Summary

| Step | What you did | Number |
|---|---|---|
| B6.1 | precision, recall, F1 on 10 sites by hand and with `prf()` | 0.833 / 0.714 / 0.769 |
| B6.2 | each caller vs the truth set ({T['n_truth']:,} variants), SNVs and indels | precision {p3(PC['FreeBayes']['all']['precision'])} (FreeBayes) to {p3(PC['DeepVariant']['all']['precision'])} (DeepVariant) |
| B6.3 | UpSet plot of {len(up_all)} intersections | largest bar = all four ({a4['n']:,}); DeepVariant only {dvo['n']} |
| B6.4 | consensus voting | all four: precision {p3(a4['precision'])}, recall {p3(cons['=4']['recall'])}; single-caller calls: precision {p3(T['single_caller']['precision'])} |
| B6.5 | context of shared false calls | {sfp['segdup_lowmap']} of {T['n_shared_fp']} ({seg_share:.0f}%) in segmental duplications ({seg_truth:.0f}% of truth variants) |

**Your presentation sentence** (fill in from your outputs):

> **Question.** Which of four variant callers should we trust on chr20 10–16 Mb, and does requiring agreement help?
> **Takeaway.** Against a GIAB-style truth set of {T['n_truth']:,} variants, DeepVariant had the highest precision
> ({p3(PC['DeepVariant']['all']['precision'])}); requiring all four callers raised precision only to {p3(a4['precision'])}
> but lost {100 * (1 - cons['=4']['recall']):.0f}% of the true variants, and {a4['FP_by_context']['segdup_lowmap']} of
> its {a4['FP']} false calls lie in segmental duplications.
> **Limitation.** Simulated benchmark of one genome in one 6-Mb window; real performance depends on coverage, read
> length, the truth-set version and the high-confidence regions.
""")
nb.code("""
# Step B6 · save the numbers used on the slides to results/B6.json (every slide number comes from here)
cons_out = consensus.reset_index(names="rule").round(6).to_dict(orient="records")
slide = {f"calls_{c}": (int(allv.loc[c, "calls"]), "calls", f"set size: variants called by {c}") for c in CALLERS}
for c in CALLERS:
    for m in ["precision", "recall", "F1"]:
        slide[f"{m}_{c}"] = (round(float(allv.loc[c, m]), 4), "fraction", f"{m} of {c} vs the truth set (all variants)")
slide["n_truth"] = (n_truth, "variants", "truth variants in the high-confidence regions")
slide["all_four_count"] = (int(c4["calls"]), "calls", "variants called by all four callers (largest UpSet bar)")
slide["all_four_precision"] = (round(float(c4["precision"]), 4), "fraction", "precision of the all-four intersection")
slide["all_four_false"] = (int(c4["FP"]), "calls", "false calls inside the all-four intersection")
slide["deepvariant_only_count"] = (dv_only, "calls", "variants called by DeepVariant and no other caller")
slide["deepvariant_only_real"] = (int(only.loc["DeepVariant", "real"]), "calls", "DeepVariant-only calls that are real")
slide["single_caller_precision"] = (p_single, "fraction", "precision of calls made by exactly one caller")
for name, key in [("≥ 2", "ge2"), ("≥ 3", "ge3"), ("= 4 (all four)", "eq4"), ("≥ 1 (union)", "ge1")]:
    for m in ["precision", "recall", "F1"]:
        slide[f"consensus_{key}_{m}"] = (round(float(consensus.loc[name, m]), 4), "fraction", f"{m} of rule {name}")
slide["segdup_pct_shared_false"] = (float(ctx.loc["segdup_lowmap", "% of shared false"]), "%",
                                    "share of false calls made by >= 3 callers that lie in segmental duplications")
slide_numbers = {k: {"value": v, "unit": u, "meaning": m} for k, (v, u, m) in slide.items()}
dviz.save_results("B6", extra={"slide_numbers": slide_numbers, "consensus_table": cons_out,
                               "upset_counts": counts[["callers", "count"]].to_dict(orient="records"),
                               "context_table": ctx.reset_index(names="context").to_dict(orient="records")})
print(f"Notebook time: {time.time() - T_START:.0f} s")
""")

out = ROOT / "electives" / "B6_variant_callers.ipynb"
out.parent.mkdir(parents=True, exist_ok=True)
nb.write(str(out))
print("wrote", out.relative_to(ROOT), len(nb.cells), "cells")
