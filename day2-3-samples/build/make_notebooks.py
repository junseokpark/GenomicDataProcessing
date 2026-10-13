#!/usr/bin/env python3
"""Write practice1_sample_to_platform.ipynb and practice2_cohort_qc.ipynb with nbformat.

The notebooks are written without outputs (students run them). Cell ids are fixed so that
rebuilding produces byte-identical files.

Usage
  python build/make_notebooks.py
"""
from pathlib import Path

import nbformat
from nbformat.v4 import new_code_cell, new_markdown_cell, new_notebook

ROOT = Path(__file__).resolve().parents[1]

STYLE = r'''plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                     "legend.frameon": False, "axes.spines.top": False,
                     "axes.spines.right": False})'''

# =====================================================================================
# Mini practice 1
# =====================================================================================
P1 = []
md = lambda s: P1.append(("md", s))
code = lambda s: P1.append(("code", s))

md(r'''# Mini practice 1 · Is this sample fit for this platform?

**Genomic Data Processing · Day 2 (10/13) · about 15 minutes · practice slide 19**

Twelve tubes arrive at a sequencing core with a sample sheet. Before anything goes on a
machine, someone decides, tube by tube: **PASS**, **REVIEW**, **FAIL** or **BLOCKED**, and
writes down why. In Part A that someone is you. In Part B you look at what three platforms
return from the same blood DNA.

One idea runs through both parts:

> A genomics sample is **specimen + molecule + metadata**, and the sample sets the ceiling for
> what any platform can read. **You can clean DNA; you cannot make it longer.**

### What is real and what is synthetic
| Part of the package | Status |
|---|---|
| `hbb_window.hg38.fa`: hg38 chr11:5,218,771–5,278,770 | **Real** UCSC sequence, the same 60-kb β-globin window as Day 1 |
| `genes.hg38.bed`: β-globin genes and exons | **Real** NCBI RefSeq Select models, from Day 1 |
| Input, length and purity numbers in `platform_requirements.tsv` | From the vendor document named in its `source` column, otherwise marked **teaching threshold** |
| The 12 samples, every read, the PAF alignments | **Synthetic**: simulated, no real person |

Nothing in this notebook is a clinical result.

### How to work
* Run the cells from top to bottom (Shift+Enter). Each one finishes in about a second.
* **Predict, then run.** When a cell asks for a prediction, say it or write it down first.
* Each step says **what changes** and **how to check it**. If a cell will not run, the verified
  answer is in `checkpoints/`: load it and keep going.''')

code(r'''import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA  = Path("data/practice1")     # open this notebook from the day2-3-samples/ folder
READS = DATA / "reads"
CHECK = Path("checkpoints")
OUT   = Path("outputs")
OUT.mkdir(exist_ok=True)
pd.set_option("display.max_colwidth", 120)
''' + STYLE + r'''

print("data :", sorted(p.name for p in DATA.iterdir() if p.is_file()))
print("reads:", sorted(p.name for p in READS.iterdir()))''')

md(r'''---
## Part A · From sample sheet to platform decision

*Concept slides: "Measure the molecule before you sequence" (slide 6), "Metadata: the sample's
passport" (slide 7), "Match the sample to the platform" (slide 18).*

### A1 · Read the sample sheet
Each row is one tube, and it carries three kinds of information:

* **specimen**: where the material came from (blood, saliva, an FFPE tumour block, plasma, cells)
* **molecule**: what was extracted (DNA, RNA, or intact cells that have not been extracted yet)
* **metadata**: what the study wants (`requested_assay`, `study_question`) and whether we may
  use the tube at all (`consent_genomic_research`)

`NA` means *not measured*. It never means zero.

**Predict, then run.** One tube must not be processed whatever its quality numbers say.
Which column will tell you?''')

code(r'''sheet = pd.read_csv(DATA / "sample_sheet.tsv", sep="\t")
print(sheet.shape[0], "tubes,", sheet.shape[1], "columns")
sheet[["sample_id", "subject_id", "specimen", "molecule", "requested_assay",
       "study_question", "matched_normal_id", "consent_genomic_research"]]''')

md(r'''It is `consent_genomic_research`: **S12** says `no`. Keep it in mind; we come back to it.

Now the measurements. What each column tells you:

| column | what it measures | can it be fixed? |
|---|---|---|
| `conc_ng_per_ul` × `volume_ul` | how much material is in the tube (ng) | sometimes: re-extract, concentrate, pool |
| `a260_280`, `a260_230` | absorbance ratios. Protein or phenol carry-over lowers 260/280; guanidine salts, phenol or carbohydrates lower 260/230 | **yes**, with a clean-up |
| `din`, `rin` | DNA / RNA integrity number, from 1 (degraded) to 10 (intact) | no |
| `dv200` | % of RNA fragments longer than 200 nucleotides, used when RNA is too degraded for RIN to help | no |
| `modal_fragment_kb` | the most common DNA fragment length | **no**: a read cannot be longer than the molecule |
| `cell_count` | intact cells, for protocols that extract the DNA themselves | recollect |''')

code(r'''sheet[["sample_id", "specimen", "molecule", "conc_ng_per_ul", "volume_ul", "a260_280",
       "a260_230", "din", "rin", "dv200", "modal_fragment_kb", "cell_count"]]''')

md(r'''### A2 · What each platform needs
Read the `source` column before you trust a number. A threshold is only as good as the
document it came from; anything without a vendor source is labelled *teaching threshold*.''')

code(r'''req = pd.read_csv(DATA / "platform_requirements.tsv", sep="\t").set_index("assay")
req[["platform", "kit_or_chemistry", "molecule", "min_input", "input_unit",
     "min_fragment_kb", "integrity_rule", "purity_rule"]]''')

code(r'''for assay, r in req.iterrows():
    print(f"{assay}\n  notes : {r['notes']}\n  source: {r['source']}")''')

md(r'''### A3 · Total mass = concentration × volume
Kits ask for nanograms; a tube label gives ng/µL. Worked example for **S08**:
30 ng/µL × 20 µL = **600 ng**. PacBio asks for 500 ng, so S08 has enough, with little to spare.

**What changes:** the sheet gains one column, `total_ng`.
**How to check:** S08 must show 600, and the cell line S07 must show `NaN`, because cells are
counted, not weighed.''')

code(r'''sheet["total_ng"] = sheet["conc_ng_per_ul"] * sheet["volume_ul"]

s08 = sheet[sheet.sample_id == "S08"].iloc[0]
assert s08.total_ng == 600
sheet[["sample_id", "specimen", "conc_ng_per_ul", "volume_ul", "total_ng", "requested_assay"]]''')

md(r'''### A4 · The decision rules, in the order a core facility applies them
1. **Consent.** No consent: **BLOCKED**. Stop here; do not even look at the quality numbers.
2. **Right molecule, enough cells.** The assay must receive the molecule it expects. For
   cell-input protocols, compare `cell_count` with the kit minimum.
3. **Enough mass.** Below the minimum: **REVIEW** if you have at least half, **FAIL** if less.
4. **Long or intact enough.** Below the integrity rule: **FAIL**. *This cannot be repaired.*
5. **Clean enough.** Purity outside the range: **REVIEW** (clean up, re-measure). *Never FAIL.*
6. **Somatic question on a tumour.** No matched normal from the same person: add a **REVIEW**.

Otherwise the tube is a **PASS**. One tube can collect several reasons; its decision is the
worst of them.

**Predict, then run.** Write down your call for **S03**, **S06**, **S11** and **S12** before
you run the cell after next.''')

code(r'''# purity windows (low, high) copied from the purity_rule column; None = no upper limit
PURITY = {"DNA": {"a260_280": (1.8, 2.0), "a260_230": (2.0, 2.2)},
          "RNA": {"a260_280": (1.9, 2.1), "a260_230": (1.8, None)}}
SUBJECT = dict(zip(sheet.sample_id, sheet.subject_id))      # which person each tube is from

def decide(s, r):
    """Decide one tube. s = one row of the sample sheet, r = one row of the requirements.
    Returns (decision, list of reasons)."""
    # 1. Consent first: without consent we never look at quality
    if s["consent_genomic_research"] != "yes":
        return "BLOCKED", ["BLOCKED: no consent for genomic research (do not process, do not assess quality)"]

    # 2. The assay must receive the molecule it expects (DNA, RNA or intact cells)
    if s["molecule"] != r["molecule"]:
        return "FAIL", [f"FAIL: assay needs {r['molecule']} but the tube holds {s['molecule']}"]
    reasons = []

    # 2-3. Enough material? Cells are counted; DNA and RNA are weighed (ng = conc x volume)
    unit = r["input_unit"]
    have = s["cell_count"] if unit == "cells" else s["conc_ng_per_ul"] * s["volume_ul"]
    need = r["min_input"]
    if have < need:
        level = "REVIEW" if have >= 0.5 * need else "FAIL"
        reasons.append(f"{level}: input {have:,.0f} {unit} is {have / need:.0%} of the {need:,.0f} {unit} minimum")

    # 4. Long or intact enough? A rule reads like "modal_fragment_kb >= 15". Cannot be repaired.
    if pd.notna(r["integrity_rule"]):
        metric, threshold = r["integrity_rule"].split(" >= ")
        threshold = float(threshold)
        value = s[metric]
        if pd.isna(value):
            reasons.append(f"REVIEW: {metric} not measured")
        elif value < threshold:
            reasons.append(f"FAIL: {metric} {value:g} < {threshold:g} required (cannot be repaired)")

    # 5. Clean enough? A clean-up fixes purity, so this is REVIEW, never FAIL
    for ratio, (low, high) in PURITY.get(s["molecule"], {}).items():
        value = s[ratio]
        if pd.isna(value):
            reasons.append(f"REVIEW: {ratio} not measured")
        elif value < low:
            reasons.append(f"REVIEW: {ratio} {value:.2f} below {low:.2f} (clean up, re-measure)")
        elif high is not None and value > high:
            reasons.append(f"REVIEW: {ratio} {value:.2f} above {high:.2f} (clean up, re-measure)")

    # 6. A somatic question on a tumour needs a matched normal from the same person
    if "tumour" in s["specimen"] and "somatic" in s["study_question"].lower():
        normal_id = s["matched_normal_id"]
        if pd.isna(normal_id):
            reasons.append("REVIEW: tumour sample for a somatic question but no matched normal")
        elif SUBJECT.get(normal_id) != s["subject_id"]:
            reasons.append(f"REVIEW: matched normal {normal_id} is not from the same subject")

    levels = [reason.split(":")[0] for reason in reasons]
    decision = "FAIL" if "FAIL" in levels else "REVIEW" if "REVIEW" in levels else "PASS"
    return decision, reasons

print("decide() is ready")''')

md(r'''Look at steps 4 and 5 side by side. Both compare a number with a threshold, but only step 4
can return FAIL. That one difference *is* the lecture's message: a clean-up fixes purity;
nothing fixes length.

Now apply `decide()` to every tube, each against the assay it requested.
**What changes:** a new table with one row per tube. **How to check:** 12 rows, and every
non-PASS row carries at least one reason.''')

code(r'''rows = []
for _, s in sheet.iterrows():
    decision, reasons = decide(s, req.loc[s["requested_assay"]])
    rows.append({"sample_id": s["sample_id"], "requested_assay": s["requested_assay"],
                 "total_ng": s["total_ng"], "decision": decision,
                 "reasons": "; ".join(reasons) if reasons else "-"})
decisions = pd.DataFrame(rows)

assert len(decisions) == 12
assert (decisions.loc[decisions.decision != "PASS", "reasons"] != "-").all()
for r in decisions[decisions.decision != "PASS"].itertuples():
    print(f"{r.sample_id} {r.decision}")
    for reason in r.reasons.split("; "):
        print("    " + reason)
decisions[["sample_id", "requested_assay", "total_ng", "decision"]]''')

code(r'''tally = decisions["decision"].value_counts()
print(tally.to_string())
assert tally.to_dict() == {"PASS": 6, "FAIL": 3, "REVIEW": 2, "BLOCKED": 1}, \
    "different tally: compare your table with checkpoints/practice1_decisions.tsv"''')

md(r'''*If `decide()` will not run, load the verified table and continue:*
`decisions = pd.read_csv(CHECK / "practice1_decisions.tsv", sep="\t")`

### Checkpoint A · read three rows aloud
* **S03** (FFPE tumour DNA → PacBio HiFi) fails twice: 240 ng is under half the input, and the
  fragments are 0.35 kb when the library needs 15–20 kb. More FFPE sections could fix the mass;
  nothing fixes the length. It also lacks a matched normal for a somatic question.
* **S06** (plasma cfDNA → ONT long reads). Cell-free DNA is about 167 bp because it was cut
  between nucleosomes before it left the cell. No platform can read across 50 kb from it, and
  20 ng is 2% of the 1 µg the kit asks for.
* **S12** has excellent numbers and is **BLOCKED**. Consent is checked before quality, and a
  blocked tube's quality data are not used for anything.

> **Q.** S02 and S11 are both REVIEW. Is it the same problem?
> **A.** No. S02 is *dirty* (A260/230 = 1.45, typical of salt or guanidine carry-over): a bead
> clean-up and a re-measurement fix it. S11 is *short of mass* (300 of 500 ng): re-quantify with
> a fluorometer, concentrate, recollect, or move to an assay that needs less.''')

md(r'''### A5 · From a decision to a suggestion
Rules tell you when to stop. Choosing what to do instead needs the study question, so the
suggestions below are written by a person and kept *next to* the rule-based decision, never
mixed into it. S06 shows why: mutation detection and 50-kb phasing lead to different specimens.

**What changes:** a `suggested_alternative` column, then the table is saved.
**How to check:** the printed line says `True` when your file equals the checkpoint.''')

code(r'''ALTERNATIVE = {
    "S02": "Bead clean-up, re-measure A260/230, then run ILMN_PCRFREE_WGS as requested",
    "S03": "ILMN_WES or a panel with an FFPE-compatible prep, plus a matched normal (e.g. blood) from P003",
    "S06": "Mutation detection: Illumina cfDNA targeted panel. Phasing: collect blood and extract long DNA",
    "S10": "RNA_CAPTURE_FFPE: DV200 62% (medium quality) needs 20-40 ng; see the re-check",
    "S11": "Re-quantify (fluorometer) and concentrate; or recollect as blood; or switch to ILMN_PCRFREE_WGS (resolves fewer SVs)",
    "S12": "Do not process; refer back to the study team (consent)",
}
decisions["suggested_alternative"] = [ALTERNATIVE.get(sid, "-") for sid in decisions.sample_id]
decisions.to_csv(OUT / "practice1_decisions.tsv", sep="\t", index=False, na_rep="NA")

mine = pd.read_csv(OUT / "practice1_decisions.tsv", sep="\t")
print("matches checkpoint:", mine.equals(pd.read_csv(CHECK / "practice1_decisions.tsv", sep="\t")))
decisions[["sample_id", "decision", "suggested_alternative"]]''')

md(r'''### A6 · Same tube, better-matched assay
**S10** is FFPE tumour RNA with RIN 2.3. Poly(A) selection grabs each transcript by its 3′ tail,
so on fragmented RNA you keep only the pieces nearest the tail (3′ bias). Capture-based RNA-seq
uses probes along the whole transcript and tolerates fragments. For degraded RNA the useful
number is **DV200**, not RIN.

**Predict, then run.** S10 against `RNA_CAPTURE_FFPE`: PASS or not? And S03 against the exome
(`ILMN_WES`), which accepts short FFPE DNA?''')

code(r'''def dv200_band(dv200):
    """Recommended input for capture RNA-seq by DV200 (Illumina technical note)."""
    if dv200 > 70:
        return "high quality: 20 ng"
    if dv200 >= 50:
        return "medium quality: 20-40 ng"
    if dv200 >= 30:
        return "low quality: 40-100 ng"
    return "too degraded: not recommended"

def tube(sample_id):
    return sheet[sheet.sample_id == sample_id].iloc[0]

for sample_id, assay in [("S10", "RNA_POLYA"), ("S10", "RNA_CAPTURE_FFPE"),
                         ("S03", "PACBIO_HIFI_WGS"), ("S03", "ILMN_WES")]:
    decision, reasons = decide(tube(sample_id), req.loc[assay])
    print(f"{sample_id} on {assay:<16} -> {decision:<6} {'; '.join(reasons)}")

s10 = tube("S10")
print(f"\nS10: DV200 {s10.dv200:g}% -> {dv200_band(s10.dv200)}; the tube holds {s10.total_ng:g} ng")
assert decide(s10, req.loc["RNA_CAPTURE_FFPE"])[0] == "PASS"''')

md(r'''### Checkpoint A6
* **S10:** FAIL → **PASS** by changing the assay, not the tube. DV200 62% is *medium quality*,
  for which 20–40 ng is recommended, and the tube holds 225 ng.
* **S03:** FAIL → **REVIEW** on the exome. The mass and length problems disappear: the exome
  needs only 100 ng, and short-read libraries are built from short fragments anyway. Two flags
  remain: a clean-up for A260/230 = 1.95, and the missing matched normal. A platform cannot fix
  the study design.''')

md(r'''---
## Part B · What each platform returns

*Concept slide: "Side by side: which platform for which question?" (slide 17).*

We simulated **S08**, the blood DNA from Part A (modal fragment 22 kb), sequenced three ways
over the real 60-kb β-globin window: **Illumina 2×150**, **PacBio HiFi** and **ONT ligation**.
Same molecule, same region, roughly the same number of bases. The question: which reads cross
an 8-kb stretch from inside *HBB* to inside *HBD*, the kind of span you need to describe a
structural variant there? This is a teaching target; we make no clinical claim about it.

### B1 · The reference window and the target
The FASTA header `>chr11:5218771-5278770` is a 1-based, closed address (the way a genome
browser writes it). In BED terms the window is `[5218770, 5278770)`, so character `i` of the
sequence string (counting from 0) is chr11 BED position `5218770 + i`.''')

code(r'''def read_fasta(path):
    name, chunks = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:]
        else:
            chunks.append(line)
    return name, "".join(chunks)

window_name, window = read_fasta(DATA / "hbb_window.hg38.fa")
WIN_START0 = 5218770                        # BED start of the window on chr11
print(window_name, f"{len(window):,} bp")
assert len(window) == 60_000

# the 18 bases you met on Day 1: hg38 chr11:5,226,995-5,227,012 (1-based, closed)
day1 = window[5226995 - 1 - WIN_START0 : 5227012 - WIN_START0]
print("chr11:5,226,995-5,227,012 =", day1)
assert day1 == "CTTCTCCTCAGGAGTCAG"''')

code(r'''target = pd.read_csv(DATA / "target_window.bed", sep="\t", header=None,
                     names=["chrom", "start", "end", "name"])
T_START, T_END = int(target["start"][0]), int(target["end"][0])
print(f"target {target['name'][0]}: chr11 [{T_START:,}, {T_END:,}) = {T_END - T_START:,} bp")

genes = pd.read_csv(DATA / "genes.hg38.bed", sep="\t")
inside = genes[(genes.chromStart < T_END) & (genes.chromEnd > T_START)]
inside[["name", "chromStart", "chromEnd", "strand"]]''')

md(r'''The target holds the whole of *HBB* exon 1 and *HBD* exons 2–3, the edge of *HBB* exon 2,
and the *HBB* promoter.

### B2 · One FASTQ record is four lines
```
@read_name        line 1: '@' then the read name
ACGTTGCA...       line 2: the bases
+                 line 3: a separator
??+5@CC>...       line 4: one quality character per base
```
Each quality character encodes a **Phred score**: Q = (ASCII code − 33), and the chance the base
is wrong is p = 10^(−Q/10). So Q10 = 1 error in 10, Q20 = 1 in 100, Q30 = 1 in 1,000,
Q40 = 1 in 10,000. Real read names are longer (PacBio: `movie/hole/ccs`; ONT: a random ID);
ours are short so that you can read them.''')

code(r'''FASTQ = {"illumina_2x150": [READS / "illumina_R1.fastq.gz", READS / "illumina_R2.fastq.gz"],
         "pacbio_hifi":    [READS / "pacbio_hifi.fastq.gz"],
         "ont_ligation":   [READS / "ont_ligation.fastq.gz"]}

def read_fastq(path):
    """Yield (name, bases, quality) for each 4-line FASTQ record."""
    with gzip.open(path, "rt") as fh:
        while True:
            header = fh.readline().rstrip()
            if not header:
                return
            bases = fh.readline().rstrip()
            fh.readline()                           # the '+' separator line
            quality = fh.readline().rstrip()
            assert header.startswith("@") and len(bases) == len(quality)
            yield header[1:], bases, quality

name, bases, quality = next(read_fastq(FASTQ["illumina_2x150"][0]))
print("name   :", name, f"({len(bases)} bases)")
print("bases  :", bases[:40], "...")
print("quality:", quality[:40], "...")
print("Q      :", [ord(c) - 33 for c in quality[:12]], "...")''')

md(r'''### B3 · Read N50: your turn
The **read N50** is the length L such that reads of length L or longer hold at least half of
all the bases. It describes the read a typical *base* sits in, which is what matters for
spanning, rather than the typical read.

Worked example: lengths 2, 3, 4, 5, 6 (total 20, half = 10). Longest first: 6 (running total 6),
then 5 (running total 11, now ≥ 10). **N50 = 5**, while the mean is 4.

Write `read_n50` from the three comments; the two `assert` lines check you. A working version is
filled in so the notebook runs end to end: cover it and try first.''')

code(r'''def read_n50(lengths):
    """Length L such that reads of length >= L hold at least half of all bases."""
    # 1. sort the lengths from longest to shortest
    # 2. walk down the list, keeping a running total of bases
    # 3. return the first length at which the running total reaches half of all bases
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length

assert read_n50([2, 3, 4, 5, 6]) == 5
assert read_n50([150] * 10) == 150         # equal lengths: the N50 is that length
print("read_n50 passes both checks")''')

md(r'''Now read every record of every file and keep, per read, its length and two sums over its bases:
the Q values and the error probabilities. **What changes:** a table with one row per read.
**How to check:** 4,000 Illumina reads (2,000 pairs), 40 HiFi reads, 50 ONT reads.''')

code(r'''records = []
for platform, paths in FASTQ.items():
    for path in paths:
        for name, bases, quality in read_fastq(path):
            q = np.frombuffer(quality.encode(), dtype=np.uint8).astype(int) - 33
            records.append({"platform": platform, "read": name, "length": len(bases),
                            "sum_q": q.sum(), "sum_error": (10 ** (-q / 10)).sum()})
reads = pd.DataFrame(records)

print(reads.groupby("platform", sort=False).size().to_dict())
reads.groupby("platform", sort=False)["length"].describe()[["count", "mean", "min", "max"]]''')

md(r'''### B4 · Average the error probabilities, not the Q values
Q is a logarithm, so averaging Q values hides bad bases. Four bases at Q40, Q40, Q40 and Q10:

* mean Q = 32.5, which sounds like about 1 error in 1,800 bases;
* mean error probability = (0.0001 + 0.0001 + 0.0001 + 0.1) / 4 = 0.025, which is 1 in 40.

The one bad base dominates. The honest summary is the mean of p (convert it back to a Q if you
like). **Predict, then run:** which platform does the naive average flatter most?''')

code(r'''q = np.array([40, 40, 40, 10])
p = 10 ** (-q / 10)
print(f"worked example: mean Q = {q.mean():.1f} -> error {10 ** (-q.mean() / 10):.5f}; "
      f"mean error = {p.mean():.5f} -> Q{-10 * np.log10(p.mean()):.1f}")

totals = reads.groupby("platform", sort=False)[["length", "sum_q", "sum_error"]].sum()
totals["mean_Q_naive"] = totals.sum_q / totals.length
totals["error_if_you_average_Q"] = 10 ** (-totals.mean_Q_naive / 10)
totals["mean_error"] = totals.sum_error / totals.length
totals["honest_Q"] = -10 * np.log10(totals.mean_error)
totals[["mean_Q_naive", "error_if_you_average_Q", "mean_error", "honest_Q"]].round(5)''')

md(r'''**PacBio HiFi** is flattered most: Q43 by the naive average, about Q31 honestly. Most of its bases
carry very high Q values, and its few likely errors sit in homopolymers (runs such as `AAAAAA`),
at bases with much lower Q. Averaging Q values lets the many excellent bases hide the few
doubtful ones; averaging error probabilities does not.

### B5 · Where did each read land? The PAF file
PAF is a plain tab-separated table with one line per alignment. Its first 12 columns:

| # | column | meaning |
|---|---|---|
| 1 | `qname` | read name |
| 2 | `qlen` | read length |
| 3–4 | `qstart`, `qend` | aligned part of the read (0-based, end excluded) |
| 5 | `strand` | `+` the read matches the reference as written; `-` it matches the reverse complement |
| 6 | `tname` | reference sequence name (`chr11`) |
| 7 | `tlen` | reference length (hg38 chr11 = 135,086,622 bp) |
| 8–9 | `tstart`, `tend` | aligned part of the reference: **0-based, end excluded, exactly like BED** |
| 10 | `nmatch` | bases that match the reference |
| 11 | `alnlen` | alignment length: matches, mismatches, insertions and deletions |
| 12 | `mapq` | mapping quality; **255 means not available** |

Our PAF was written from the simulation's truth rather than by an aligner, which is why every
`mapq` is 255. With real reads you would run, for example,
`minimap2 -x map-hifi ref.fa reads.fq.gz > aln.paf` (`map-ont` for nanopore, `sr` for short
reads). The coordinates are real chr11 positions, not positions inside our 60-kb string.''')

code(r'''PAF_COLUMNS = ["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen",
               "tstart", "tend", "nmatch", "alnlen", "mapq"]
paf = pd.read_csv(READS / "alignments.paf", sep="\t", header=None, names=PAF_COLUMNS)
prefix_to_platform = {"ilmn": "illumina_2x150", "hifi": "pacbio_hifi", "ont": "ont_ligation"}
paf["platform"] = paf.qname.str.split("_").str[1].map(prefix_to_platform)

assert paf.platform.notna().all() and len(paf) == len(reads)
assert ((paf.tstart >= WIN_START0) & (paf.tend <= WIN_START0 + len(window))).all()
print(len(paf), "alignments, all inside the window")
paf.groupby("platform", sort=False).head(2)''')

md(r'''A read **spans** the target only if it starts at or before the target's start *and* ends at or
after the target's end: `tstart <= 5,226,000` and `tend >= 5,234,000`.

**Predict, then run.** How many Illumina reads span 8 kb? How many HiFi reads?''')

code(r'''paf["spans_target"] = (paf.tstart <= T_START) & (paf.tend >= T_END)
spanning = paf.groupby("platform", sort=False)["spans_target"].sum()
spanning''')

md(r'''### B6 · One summary row per platform
Depth here is approximate: total bases ÷ window length (60,000 bp).
**What changes:** `outputs/practice1_read_stats.tsv` is written.
**How to check:** it equals the checkpoint; Illumina spans 0 reads, HiFi and ONT at least 3 each.
*If an earlier Part B cell will not run, the verified table is
`pd.read_csv(CHECK / "practice1_read_stats.tsv", sep="\t")`.*''')

code(r'''by_platform = reads.groupby("platform", sort=False)["length"]
stats = pd.DataFrame({
    "reads": by_platform.size(),
    "total_bases": by_platform.sum(),
    "mean_length": by_platform.mean().round(1),
    "max_length": by_platform.max(),
    "read_n50": by_platform.apply(read_n50),
    "mean_error_rate": (totals.sum_error / totals.length).round(6),
    "mean_accuracy_pct": (100 * (1 - totals.sum_error / totals.length)).round(3),
    "approx_depth_x": (by_platform.sum() / len(window)).round(2),
    "reads_spanning_target": spanning,
}).rename_axis("platform").reset_index()
stats.to_csv(OUT / "practice1_read_stats.tsv", sep="\t", index=False)

mine = pd.read_csv(OUT / "practice1_read_stats.tsv", sep="\t")
print("matches checkpoint:", mine.equals(pd.read_csv(CHECK / "practice1_read_stats.tsv", sep="\t")))
span = stats.set_index("platform")["reads_spanning_target"]
assert span["illumina_2x150"] == 0 and span["pacbio_hifi"] >= 3 and span["ont_ligation"] >= 3
stats''')

md(r'''### Checkpoint B · same depth, different answers
All three platforms read this window to roughly 10–12×, yet only long reads cross the 8-kb
target. Depth tells you how much you read; length tells you what you can connect. HiFi reads are
long *and* accurate (about 99.9% per base). The ONT reads here are about 99% accurate per base,
with errors concentrated in homopolymers, and their lengths vary far more.

Two honest limits of this simulation:
* Our ONT reads stop below 58 kb only because the window is 60 kb. Real ultra-long libraries
  (SQK-ULK114, started from cells like S07) reach read N50s of 50–100+ kb and single reads of
  several megabases.
* S08's DNA has a modal fragment of 22 kb, and the ONT read N50 sits just below it. On a
  ligation kit, read length is fragment length: the tube, not the machine, set the ceiling.''')

md(r'''### B7 · Picture it: read lengths on a log scale
A log axis lets 150-bp and 50-kb reads share one plot. The dashed line marks the 8-kb target: a
read to its left cannot span the target, whatever the depth.''')

code(r'''order = ["illumina_2x150", "pacbio_hifi", "ont_ligation"]
bins = np.logspace(np.log10(100), np.log10(60_000), 41)
fig, axes = plt.subplots(3, 1, figsize=(7.5, 6), sharex=True)
for ax, platform in zip(axes, order):
    lengths = reads.loc[reads.platform == platform, "length"]
    ax.hist(lengths, bins=bins, color="#4C72B0")
    ax.axvline(T_END - T_START, color="0.2", linestyle="--", linewidth=1)
    ax.set_ylabel("reads")
    ax.set_title(f"{platform}: {len(lengths):,} reads, read N50 {read_n50(lengths):,} bp, "
                 f"{int(spanning[platform])} span the target", loc="left")
axes[0].annotate("8-kb target", xy=(T_END - T_START, 0.75), xycoords=("data", "axes fraction"),
                 xytext=(6, 0), textcoords="offset points", fontsize=8)
axes[-1].set_xscale("log")
axes[-1].set_xticks([100, 1_000, 10_000, 50_000])
axes[-1].set_xticklabels(["100 bp", "1 kb", "10 kb", "50 kb"])
axes[-1].set_xlabel("read length (log scale)")
fig.suptitle("Only the long reads reach the 8-kb target (S08 blood DNA, same 60-kb window)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "read_length_hist.png", dpi=150)
plt.show()''')

md(r'''### B8 · Optional, if you have two minutes: draw the reads over the target
Each line is one read at its PAF coordinates; the shaded band is the target. Dark lines span it.
For Illumina only the reads near the target are drawn.''')

code(r'''fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharex=True)
for ax, platform in zip(axes, order):
    hits = paf[paf.platform == platform].sort_values("tstart")
    if platform == "illumina_2x150":
        hits = hits[(hits.tend > T_START - 2_000) & (hits.tstart < T_END + 2_000)]
    for y, r in enumerate(hits.itertuples()):
        ax.plot([r.tstart / 1e6, r.tend / 1e6], [y, y], color="0.1" if r.spans_target else "0.65",
                linewidth=1.6 if r.spans_target else 0.7)
    ax.axvspan(T_START / 1e6, T_END / 1e6, color="#4C72B0", alpha=0.15)
    ax.set_title(f"{platform}: {int(hits.spans_target.sum())} of {len(hits):,} drawn span it", loc="left")
    ax.set_xlabel("chr11 position (Mb, hg38)")
    ax.set_yticks([])
axes[0].set_ylabel("reads, sorted by start")
fig.tight_layout()
fig.savefig(OUT / "target_spanning_reads.png", dpi=150)
plt.show()''')

md(r'''---
### What you produced
| file | what it shows |
|---|---|
| `outputs/practice1_decisions.tsv` | one decision per tube, every reason kept, suggestions kept separate |
| `outputs/practice1_read_stats.tsv` | reads, bases, N50, honest error rate, depth and spanning reads per platform |
| `outputs/read_length_hist.png` | why depth alone cannot answer a spanning question |
| `outputs/target_spanning_reads.png` | (optional) the spanning reads themselves |

**Say it in one sentence.** *"I would not put S03 on PacBio, because its fragments (0.35 kb) are
far shorter than the 15–20 kb library; I checked the fragment size, which no clean-up can
change; I would run an FFPE exome with a matched normal instead."* Now say one for a tube of
your choice.''')

# =====================================================================================
# Mini practice 2
# =====================================================================================
P2 = []
md = lambda s: P2.append(("md", s))
code = lambda s: P2.append(("code", s))

md(r'''# Mini practice 2 · Cohort QC: who does not belong?

**Genomic Data Processing · Day 2 (10/13) · about 15 minutes · practice slide 26**

You receive 2,000 exomes from a large cohort. Before anyone runs an association test, every
sample has to earn its place. Were we allowed to use it? Is it clean? Is it the person it claims
to be? Do we already have this person? Which analysis set does it belong to? You will apply the
rules Wang et al. (2021) used for the UK Biobank exomes, and keep a written reason for every
sample you remove.

> **Source of the rules.** Wang Q, Dhindsa RS, Carss K, *et al.* Rare variant contribution to
> human disease in 281,104 UK Biobank exomes. *Nature* 597, 527–532 (2021).
> doi:10.1038/s41586-021-03855-y

### What is real and what is synthetic
| Part | Status |
|---|---|
| QC rules: FREEMIX > 4%, < 94.5% of CCDS bases at ≥ 10×, KING kinship > 0.0884 (greedy), Pr(EUR) ≥ 0.99 and ±4 SD on PC1–PC4 | **Real**, from the Methods of Wang et al. 2021 |
| KING kinship categories (0.354 / 0.177 / 0.0884 / 0.0442) | **Real**, the KING guide |
| Sex-check and array-concordance cut-offs (0.20, 0.30, 0.95) | **Teaching thresholds** |
| All 2,000 samples, every metric, every kinship value | **Synthetic**. No UK Biobank data are used |

### How to work
Run the cells top to bottom; each takes about a second. **Predict, then run** where asked. Every
step says what it removes and how to check it. If a cell will not run, the verified tables are
in `checkpoints/`.''')

code(r'''from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA  = Path("data/practice2")     # open this notebook from the day2-3-samples/ folder
CHECK = Path("checkpoints")
OUT   = Path("outputs")
OUT.mkdir(exist_ok=True)
pd.set_option("display.max_colwidth", 110)
''' + STYLE)

md(r'''---
## 1 · Load and describe the cohort
*Concept slide: "Metadata: the sample's passport" (slide 7).*

Three files. The **manifest** says who each sample is supposed to be. The **QC metrics** say
what the sequencing data look like. The **kinship pairs** list every pair of samples that look
related (KING kinship > 0.0442). A *participant* is a person; a *sample* is one tube of their
DNA, and one person can have two samples.''')

code(r'''manifest = pd.read_csv(DATA / "cohort_manifest.tsv", sep="\t")
metrics  = pd.read_csv(DATA / "sample_qc_metrics.tsv", sep="\t")
kin      = pd.read_csv(DATA / "kinship_pairs.tsv", sep="\t")
print("manifest", manifest.shape, "| metrics", metrics.shape, "| kinship pairs", kin.shape)
print("samples:", manifest.sample_id.nunique(), "| participants:", manifest.participant_id.nunique())
manifest.head()''')

code(r'''print("batches :", manifest.sequencing_batch.value_counts().sort_index().to_dict())
print("centres :", manifest.recruitment_centre.value_counts().sort_index().to_dict())
print("specimen:", manifest.specimen.value_counts().to_dict())
print("sex     :", manifest.reported_sex.value_counts().to_dict())
print("age     :", manifest.age_at_recruitment.min(), "to", manifest.age_at_recruitment.max())
print("years   :", manifest.collection_year.min(), "to", manifest.collection_year.max())
print("consent :", manifest.consent_status.value_counts().to_dict())''')

md(r'''### What each QC metric means
| column | in plain words | typical good value |
|---|---|---|
| `mean_depth` | how many reads cover an average targeted base | tens of reads |
| `pct_ccds_10x` | % of the consensus coding sequence (CCDS) covered by at least 10 reads: can we see every coding base well enough to call a genotype? | ≥ 94.5% |
| `freemix` | VerifyBAMID's estimate of the fraction of DNA that comes from *another person*. A mixture makes heterozygous sites look odd: allele fractions drift away from 0, 50 and 100% | ≤ 0.04 |
| `chrX_het_rate` | fraction of chrX variant sites (outside the pseudo-autosomal regions) that are heterozygous. One X copy cannot be heterozygous, so XY samples sit near 0 | XY ≈ 0, XX ≈ 0.6 |
| `chrY_depth_ratio` | reads on chrY relative to the autosomes. One Y copy against two copies of each autosome gives about 0.5 | XY ≈ 0.5, XX ≈ 0 |
| `array_concordance` | how often the exome genotypes agree with the SNP-array genotypes of the same participant. The same person agrees almost perfectly; two different people agree only by chance (about 0.6 here) | ≥ 0.95 |
| `PC1`–`PC4` | principal components: a few numbers summarising genome-wide genotype patterns. People with similar genetic ancestry sit close together. PCs describe genetic similarity, not identity or race | – |
| `pr_eur` | peddy's probability that the genotypes resemble the European reference group of the 1000 Genomes Project | – |
| `kinship` | KING kinship coefficient: about 0.5 for the same person or identical twins, 0.25 for parent–child or full siblings, 0.125 for second-degree and 0.0625 for third-degree relatives | – |

### The bookkeeping rule
Every removal goes through one function that records **which step** removed the sample and
**why**, and adds a line to the funnel. Nothing leaves the cohort without a written reason.
Run the cells in order; to redo a step, use *Kernel → Restart & Run All*.

*If a step will not run, catch up from the checkpoint and carry on, e.g. for step 4:*
`ck = pd.read_csv(CHECK / "practice2_excluded_samples.tsv", sep="\t")`, then
`remaining -= set(ck.sample_id[ck.step <= "4_identity"])`.''')

code(r'''qc = manifest.merge(metrics, on="sample_id", validate="one_to_one").set_index("sample_id").sort_index()
remaining = set(qc.index)
excluded = []      # one row per removed sample
flags = {}         # sample_id -> why it is kept but must be looked at
funnel = [{"step": "0_sequenced", "rule": "all delivered exomes", "n_removed": 0,
           "n_remaining": len(remaining)}]

def current():
    """The samples still in the analysis set."""
    return qc[qc.index.isin(remaining)]

def remove(reasons, step, rule):
    """reasons = {sample_id: why}. Remove them, record step and reason, add a funnel line."""
    if any(line["step"] == step for line in funnel):
        print(f"{step} was already applied. To redo it, restart the kernel and run all cells.")
        return
    for sample_id, why in reasons.items():
        remaining.discard(sample_id)
        excluded.append({"sample_id": sample_id, "step": step, "reason": why})
    funnel.append({"step": step, "rule": rule, "n_removed": len(reasons),
                   "n_remaining": len(remaining)})
    print(f"{step}: removed {len(reasons)}, {len(remaining)} remain")

def record_set(step, rule):
    """A funnel line that names an analysis set (nothing is removed here)."""
    if not any(line["step"] == step for line in funnel):
        funnel.append({"step": step, "rule": rule, "n_removed": 0, "n_remaining": len(remaining)})

print(len(remaining), "samples to start with")''')

md(r'''---
## 2 · Consent first
Governance comes before statistics. A participant who withdrew consent is removed before anyone
looks at their data, so withdrawn samples appear in no later table, plot or count.

**Predict, then run.** Suppose their QC metrics are excellent. Does that change anything?''')

code(r'''c = current()
withdrawn = c[c.consent_status == "withdrawn"]
remove({s: "consent withdrawn: removed before any QC" for s in withdrawn.index},
       "1_consent", "consent_status == withdrawn")''')

md(r'''It changes nothing. Notice that we did not print their metrics, and we will not.

---
## 3 · Contamination and coverage
Two data-quality rules from the paper: VerifyBAMID **FREEMIX > 4%** (the tube holds a mixture of
two people) and **< 94.5% of CCDS bases at ≥ 10×** (too much of the coding sequence is too thin
to genotype). Contamination goes first for a reason: a mixture inflates heterozygosity, which
would confuse the sex check in step 4.

**What changes:** two funnel lines. **How to check:** each removed sample breaks its rule.''')

code(r'''c = current()
contaminated = c[c.freemix > 0.04]
remove({s: f"freemix {v:.3f} > 0.04" for s, v in contaminated.freemix.items()},
       "2_contamination", "freemix > 0.04 (VerifyBAMID)")

c = current()
low_coverage = c[c.pct_ccds_10x < 94.5]
remove({s: f"pct_ccds_10x {v:.2f} < 94.5" for s, v in low_coverage.pct_ccds_10x.items()},
       "3_coverage", "pct_ccds_10x < 94.5")

qc.loc[list(contaminated.index) + list(low_coverage.index),
       ["reported_sex", "sequencing_batch", "mean_depth", "pct_ccds_10x", "freemix", "chrX_het_rate"]]''')

md(r'''Look at the contaminated **male** samples: their `chrX_het_rate` is several times the usual male
value of about 0.02, because a second person's DNA adds extra alleles at chrX sites. That is
why contamination is checked before the sex check, not after.

Two of the three low-coverage samples come from batch **B3**. Keep that in mind for step 6.

---
## 4 · Is every sample who it claims to be?
*Concept slide: "Cohort QC: is every sample who it claims to be?" (slide 24).*

Two independent checks:
* **Sex check:** genetic sex from the data (is Y present? are there two X copies?) against the
  sex in the manifest.
* **Array concordance:** exome genotypes against the SNP array of the same participant.

Cross the two checks and each cell of the table has its own meaning:

| | array agrees (≥ 0.95) | array disagrees (< 0.95) |
|---|---|---|
| **sex concordant** | pass | identity not confirmed: a same-sex swap the sex check cannot see → **exclude** |
| **sex discordant** | the *record* is wrong, not the DNA → **correct the manifest, keep** | likely sample swap → **exclude** |
| **sex ambiguous** | unusual sex chromosomes (e.g. XXY-like) → **keep, flag for review; never call it a swap** | identity not confirmed → **exclude** |''')

code(r'''Y_PRESENT, X_TWO_COPIES = 0.20, 0.30      # teaching thresholds

def genetic_sex(x_het, y_ratio):
    y_present = y_ratio >= Y_PRESENT
    two_x = x_het >= X_TWO_COPIES
    if y_present and not two_x:
        return "XY"
    if two_x and not y_present:
        return "XX"
    return "ambiguous"

c = current().copy()
c["genetic_sex"] = [genetic_sex(x, y) for x, y in zip(c.chrX_het_rate, c.chrY_depth_ratio)]
expected = c.reported_sex.map({"F": "XX", "M": "XY"})
c["sex_check"] = np.where(c.genetic_sex == "ambiguous", "ambiguous",
                          np.where(c.genetic_sex == expected, "concordant", "discordant"))
c["array_ok"] = c.array_concordance >= 0.95
pd.crosstab(c.sex_check, c.array_ok)''')

md(r'''**Predict, then run.** The next cell lists every sample outside the "concordant + array agrees"
cell. For each one: exclude, correct the record, or flag for review?''')

code(r'''look = c[(c.sex_check != "concordant") | ~c.array_ok]
look[["participant_id", "reported_sex", "chrX_het_rate", "chrY_depth_ratio", "genetic_sex",
      "array_concordance", "sex_check"]]''')

code(r'''reasons = {}
for s, r in look.iterrows():
    if not r.array_ok:
        if r.sex_check == "discordant":
            reasons[s] = (f"array_concordance {r.array_concordance:.3f} < 0.95 and reported "
                          f"{r.reported_sex} but genetic {r.genetic_sex}: likely swap")
        else:
            reasons[s] = (f"array_concordance {r.array_concordance:.3f} < 0.95 with {r.sex_check} "
                          f"sex: identity not confirmed")
    elif r.sex_check == "discordant":
        flags[s] = "metadata error: reported sex corrected"
    else:
        flags[s] = "sex chromosomes ambiguous: review"
remove(reasons, "4_identity", "array_concordance < 0.95")

# correct the record in a copy and tell the study team; never edit the source file silently
manifest_corrected = manifest.set_index("sample_id")
for s, flag in flags.items():
    if flag.startswith("metadata error"):
        manifest_corrected.loc[s, "reported_sex"] = {"XX": "F", "XY": "M"}[c.loc[s, "genetic_sex"]]
pd.Series(flags, name="flag").rename_axis("sample_id").to_frame()''')

md(r'''### Checkpoint 4 · four different stories behind "something is off"
* **The swapped pair** (consecutive sample IDs, as if neighbours on a plate). Each is
  sex-discordant *and* array-discordant: the DNA in each tube belongs to someone else. Exclude
  both and ask the lab to trace the plate. With the full array data you could test whether each
  tube matches the *other* participant's array, which would confirm the swap.
* **The metadata error.** Sex-discordant, yet the exome matches the participant's own array at
  0.999. The DNA is right and the record is wrong: correct the manifest (in a copy, and tell the
  study team) and keep the sample.
* **The XXY-like sample.** Y present *and* two X copies, with array concordance 0.999, so it is
  who it claims to be. This is biology, not a lab error: flag it for review by the study team,
  never call it a swap, and never "fix" it to fit a box.
* **The same-sex discordance.** The sex check passes, the array check fails. A sex check only
  sees swaps between people of different sex; the array check sees them all.''')

md(r'''---
## 5 · Duplicates and relatives
KING kinship guide:

| kinship | relationship |
|---|---|
| > 0.354 | the same person twice (duplicate) or monozygotic (identical) twins |
| 0.177 – 0.354 | first degree: parent–child, full siblings |
| 0.0884 – 0.177 | second degree: half-siblings, grandparent, aunt or uncle |
| 0.0442 – 0.0884 | third degree: first cousins |

Most association tests assume unrelated samples, so related pairs are broken up. The paper kept
no pair with kinship above **0.0884**, which leaves third-degree relatives in.''')

code(r'''def kin_degree(k):
    if k > 0.354:
        return "duplicate / MZ twin"
    if k > 0.177:
        return "1st degree"
    if k > 0.0884:
        return "2nd degree"
    return "3rd degree"

kin["degree"] = kin.kinship.apply(kin_degree)
live = kin[kin.sample_id_1.isin(remaining) & kin.sample_id_2.isin(remaining)]
print("pairs in the file:", len(kin), "| pairs among remaining samples:", len(live))
live.degree.value_counts()''')

md(r'''### 5a · Duplicates first
A pair above 0.354 is the same genome twice. **Same `participant_id`:** a technical duplicate;
keep the copy with the better `pct_ccds_10x`. **Different `participant_id`:** identical twins, or
one person enrolled twice. Genetics cannot tell these apart, so keep one and flag the pair for the
study team (same age? same centre?).''')

code(r'''pct = qc.pct_ccds_10x
dups = live[live.kinship > 0.354]
reasons = {}
for r in dups.itertuples():
    a, b = r.sample_id_1, r.sample_id_2
    if a in reasons or b in reasons:
        continue
    lose, keep = (a, b) if pct[a] < pct[b] else (b, a)
    if qc.participant_id[a] == qc.participant_id[b]:
        reasons[lose] = f"duplicate of {keep} (same participant {qc.participant_id[a]}), lower pct_ccds_10x"
    else:
        reasons[lose] = (f"kinship {r.kinship:.3f} with {keep} (different participant): "
                         f"MZ twin or unrecorded duplicate, keep one, review")
        flags[keep] = "kept member of a twin/duplicate pair: review"
remove(reasons, "5_duplicates", "kinship > 0.354: keep the better sample")

shown = sorted(set(dups.sample_id_1) | set(dups.sample_id_2))
qc.loc[shown, ["participant_id", "reported_sex", "age_at_recruitment", "recruitment_centre",
               "sequencing_batch", "pct_ccds_10x"]]''')

md(r'''### 5b · Relatives: remove as few people as possible
The paper used a greedy rule (the `ukbtools` R package): repeatedly remove the sample with the
**most** remaining close relatives. When two samples tie, we remove the one with the lower
`pct_ccds_10x` (our tie-break). In a trio this removes the child, who is related to both
parents, and keeps both parents: two people kept instead of one.

**Predict, then run.** In a family of three full siblings, how many are removed?''')

code(r'''close = live[(live.kinship > 0.0884) & live.sample_id_1.isin(remaining) & live.sample_id_2.isin(remaining)]
pairs = list(zip(close.sample_id_1, close.sample_id_2))
reasons = {}
while pairs:
    counts = pd.Series([s for pair in pairs for s in pair]).value_counts()
    most = counts[counts == counts.max()].index
    victim = sorted(most, key=lambda s: (pct[s], s))[0]      # tie -> lower pct_ccds_10x
    reasons[victim] = f"kinship > 0.0884 with {counts.max()} remaining sample(s)"
    pairs = [p for p in pairs if victim not in p]
remove(reasons, "6_relatedness", "kinship > 0.0884: greedy, most relatives first")

pan = current()
record_set("pan_ancestry_set", "QC-passed, unrelated, any ancestry")
left = kin[kin.sample_id_1.isin(remaining) & kin.sample_id_2.isin(remaining)]
print("pairs left among the remaining samples:", left.degree.value_counts().to_dict())
assert (left.kinship <= 0.0884).all()''')

md(r'''Two of three siblings are removed: the first has two relatives left, then the remaining pair is
broken. The pairs left are all third degree, which the 0.0884 rule allows. One more detail to
notice: a duplicated participant who also has a sibling was counted **once**, because the
duplicate was resolved before relatedness.

---
## 6 · Ancestry and batch: structure you keep track of
*Concept slide: "What a cohort can and cannot tell you" (slide 25).*

Neither genetic ancestry nor sequencing batch is an error. Both are structure that can create
false associations if ignored, and the paper handled them differently:

* **Ancestry.** The main analyses used a genetically homogeneous European-ancestry set: peddy
  Pr(European) ≥ 0.99, then within ±4 SD of the mean on each of the top four PCs.
* **Batch.** The exomes arrived in four delivery batches. The authors made a "dummy phenotype" for
  each batch to find genes and variants that tracked batch, and excluded those genes and variants.

**What changes:** two funnel lines and the European analysis set. **How to check:** the samples
removed here are *not* QC failures; they remain in the pan-ancestry set.''')

code(r'''c = current()
not_eur = c[c.pr_eur < 0.99]
remove({s: f"pr_eur {v:.3f} < 0.99" for s, v in not_eur.pr_eur.items()},
       "7_ancestry_pr_eur", "pr_eur < 0.99 (peddy)")

eur_like = current()                                      # Pr(EUR) >= 0.99
PCS = ["PC1", "PC2", "PC3", "PC4"]
mean, sd = eur_like[PCS].mean(), eur_like[PCS].std()
z = (eur_like[PCS] - mean) / sd                           # distance from the mean, in SDs
furthest = z.abs().idxmax(axis=1)                         # the PC on which each sample is furthest out
outliers = {s: f"{pc} {z.loc[s, pc]:+.1f} SD from the mean (limit 4)"
            for s, pc in furthest.items() if abs(z.loc[s, pc]) > 4}
remove(outliers, "8_ancestry_pc_outlier", "outside mean +/- 4 SD on PC1-PC4")
record_set("eur_analysis_set", "European-ancestry analysis set")
eur = current()''')

code(r'''groups = [("European analysis set", pan.index.isin(eur.index), "#4C72B0", "o", 5),
          ("Pr(EUR) < 0.99", pan.index.isin(not_eur.index), "#DD8452", "o", 9),
          ("Pr(EUR) >= 0.99 but outside 4 SD", pan.index.isin(list(outliers)), "#C44E52", "x", 28)]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
for ax, (x, y) in zip(axes, [("PC1", "PC2"), ("PC3", "PC4")]):
    for label, mask, colour, marker, size in groups:
        ax.scatter(pan.loc[mask, x], pan.loc[mask, y], s=size, c=colour, marker=marker,
                   linewidths=1, alpha=0.75, label=f"{label} (n={mask.sum():,})")
    box = plt.Rectangle((mean[x] - 4 * sd[x], mean[y] - 4 * sd[y]), 8 * sd[x], 8 * sd[y],
                        fill=False, linestyle="--", edgecolor="0.3", linewidth=0.8)
    ax.add_patch(box)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
axes[0].set_title("PC1-PC2: the ancestry clusters", loc="left")
axes[1].set_title("PC3-PC4: the outliers only show here (dashed box = 4 SD)", loc="left")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=3)
fig.suptitle(f"Pan-ancestry analysis set, n = {len(pan):,}", fontsize=10)
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(OUT / "pca.png", dpi=150)
plt.show()''')

code(r'''batches = sorted(pan.sequencing_batch.unique())
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, column in zip(axes, ["mean_depth", "pct_ccds_10x"]):
    ax.boxplot([pan.loc[pan.sequencing_batch == b, column] for b in batches], widths=0.5)
    ax.set_xticks(range(1, len(batches) + 1))
    ax.set_xticklabels(batches)
    ax.set_xlabel("sequencing batch")
    ax.set_ylabel(column)
lc = qc.loc[list(low_coverage.index)]
axes[1].scatter([batches.index(b) + 1 for b in lc.sequencing_batch], lc.pct_ccds_10x, s=30,
                facecolors="none", edgecolors="#C44E52", label="excluded in step 3 (not in boxes)")
axes[1].axhline(94.5, color="0.3", linestyle="--", linewidth=0.8)
axes[1].annotate("94.5% rule", xy=(0.02, 94.5), xycoords=("axes fraction", "data"),
                 xytext=(0, 4), textcoords="offset points", fontsize=8)
axes[1].legend(loc="lower left")
axes[0].set_title("B3 has lower depth ...", loc="left")
axes[1].set_title("... and lower coverage, yet nearly all of it passes", loc="left")
fig.suptitle(f"Batch check on the pan-ancestry analysis set (n = {len(pan):,})", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "batch_ccds.png", dpi=150)
plt.show()

batch_median = pan.groupby("sequencing_batch")[["mean_depth", "pct_ccds_10x"]].median()
batch_median.round(2)''')

md(r'''**Read the plot.** B3 sits about one percentage point lower on `pct_ccds_10x` and about ten reads
lower on mean depth, yet almost all of its samples pass the 94.5% line. Removing B3 would throw
away a quarter of the cohort for a difference that a model can handle. The usual response: keep
the samples, report the batch effect, include batch as a covariate, and check whether any gene or
variant result tracks batch, which is what the paper's per-batch dummy phenotypes did. Do not
move the 94.5% line to catch B3's weakest passing sample (94.81%): thresholds are set before you
look at the data, not after.''')

md(r'''---
## The QC funnel
Every line of our funnel is one rule, and every number traces back to rows of the
excluded-samples table. For comparison, the paper's funnel (Wang et al. 2021, Methods):

| step | removed | remaining |
|---|---:|---:|
| exomes delivered | | 302,355 |
| sequencing-centre QC: sex discordance, contamination, unresolved duplicates, array discordance | flagged before delivery | |
| VerifyBAMID FREEMIX > 4% | 106 | |
| < 94.5% of CCDS bases at ≥ 10× | 5 | |
| KING kinship > 0.0884 (greedy, most relatives first) | 14,326 | **287,917** any ancestry |
| peddy Pr(European) < 0.99 | 18,212 | 269,706 |
| outside ±4 SD on the top four PCs | 535 | **269,171** European ancestry |''')

code(r'''funnel_df = pd.DataFrame(funnel)
excluded_df = pd.DataFrame(excluded).sort_values(["step", "sample_id"]).reset_index(drop=True)
funnel_df.to_csv(OUT / "practice2_qc_funnel.tsv", sep="\t", index=False)
excluded_df.to_csv(OUT / "practice2_excluded_samples.tsv", sep="\t", index=False)

# the funnel must add up, and every sample must be accounted for exactly once
n = funnel_df.n_remaining.to_numpy()
assert all(n[i] == n[i - 1] - funnel_df.n_removed[i] for i in range(1, len(n)))
assert len(excluded_df) + len(remaining) == len(qc) and excluded_df.sample_id.is_unique

for name in ["practice2_qc_funnel.tsv", "practice2_excluded_samples.tsv"]:
    mine = pd.read_csv(OUT / name, sep="\t")
    print(name, "matches checkpoint:", mine.equals(pd.read_csv(CHECK / name, sep="\t")))
funnel_df''')

md(r'''**Check the paper's arithmetic yourself.** Predict first: should the numbers in the table above
add up exactly?''')

code(r'''print(f"302,355 - 106 - 5 - 14,326 = {302_355 - 106 - 5 - 14_326:,}   reported: 287,917")
print(f"287,917 - 18,212           = {287_917 - 18_212:,}   reported: 269,706")
print(f"269,706 - 535              = {269_706 - 535:,}   reported: 269,171")''')

md(r'''Two of the three subtractions are off by one sample; the third is exact. The paper also says the
sequencing centre had already flagged samples before these steps, so the gap may sit in
bookkeeping we cannot see. The lesson is not that the paper is wrong. It is that a funnel kept by
hand drifts. Ours is produced by the same code that removes the samples, and the `assert` above
stops the notebook if it does not add up.

### Equity note: who is in "the cohort"?
In the paper, 269,171 of the 287,917 QC-passed exomes (93.5%) formed the European-ancestry set.
The authors also analysed 4,744 African, 1,475 East Asian and 5,714 South Asian ancestry
participants (Pr > 0.95 for that group) and ran a pan-ancestry analysis of 281,104 people. That
still leaves 6,813 QC-passed participants in none of the four groups, for example people whose
genotypes do not sit inside a single reference cluster.

Restricting to one ancestry makes a test statistically cleaner, and it leaves every finding least
tested in the people who were left out. In our cohort the `7_ancestry_pr_eur` rows of the
excluded table are not errors: they are participants, and an admixed participant fits no single
group at all. Report who was left out, analyse them where the numbers allow, and say what your
result cannot tell you.''')

code(r'''print(f"pan-ancestry set {len(pan):,} | European analysis set {len(eur):,} "
      f"({len(eur) / len(pan):.1%}) | left out of the European analysis: {len(pan) - len(eur):,}")''')

md(r'''---
## Cohort card: fill it in
One line each, from your outputs above.

* **Delivered:** ___ samples from ___ participants
* **Governance:** ___ removed for withdrawn consent, before any QC
* **Data quality:** ___ contaminated (FREEMIX > 0.04), ___ low coverage (< 94.5% CCDS at 10×)
* **Identity:** ___ excluded (array discordance); ___ record corrected; ___ flagged for review
* **Duplicates and relatives:** ___ duplicates resolved; ___ relatives removed (kinship > 0.0884)
* **Pan-ancestry analysis set:** ___ samples
* **European analysis set:** ___ samples (___% of the pan-ancestry set); left out: ___
* **Batch:** B3 median `pct_ccds_10x` ___ vs ___ elsewhere; what you would do: ___
* **One thing this cohort cannot tell you:** ___

Then run the next cell and compare.''')

code(r'''removed = funnel_df.set_index("step").n_removed
print(f"delivered        {len(qc):,} samples from {qc.participant_id.nunique():,} participants")
print(f"governance       {removed['1_consent']} withdrawn")
print(f"data quality     {removed['2_contamination']} contaminated, {removed['3_coverage']} low coverage")
print(f"identity         {removed['4_identity']} excluded; flags kept for review:")
for s, flag in flags.items():
    print(f"                 {s}  {flag}")
print(f"duplicates       {removed['5_duplicates']} resolved; relatives removed {removed['6_relatedness']}")
print(f"pan-ancestry set {len(pan):,}")
print(f"European set     {len(eur):,} ({len(eur) / len(pan):.1%}); left out {len(pan) - len(eur):,}")
others = batch_median.pct_ccds_10x.drop("B3")
print(f"batch            B3 median pct_ccds_10x {batch_median.pct_ccds_10x['B3']:.2f} "
      f"vs {others.min():.2f}-{others.max():.2f} elsewhere")''')

md(r'''---
### What you produced
| file | what it shows |
|---|---|
| `outputs/practice2_qc_funnel.tsv` | every rule, how many samples it removed, how many remain |
| `outputs/practice2_excluded_samples.tsv` | every removed sample with its step and reason |
| `outputs/pca.png` | the ancestry structure and the European restriction |
| `outputs/batch_ccds.png` | the B3 batch effect, and why it is not an exclusion |''')


# =====================================================================================
def write(cells, path, prefix):
    nb = new_notebook()
    nb.metadata = {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                   "language_info": {"name": "python"}}
    for i, (kind, src) in enumerate(cells):
        cell = new_markdown_cell(src) if kind == "md" else new_code_cell(src)
        cell.id = f"{prefix}-{i:02d}"
        nb.cells.append(cell)
    nbformat.validate(nb)
    with open(path, "w", encoding="utf-8") as fh:
        nbformat.write(nb, fh)
    n_code = sum(k == "code" for k, _ in cells)
    print(f"wrote {path.name}: {len(cells)} cells ({n_code} code)")


if __name__ == "__main__":
    write(P1, ROOT / "practice1_sample_to_platform.ipynb", "p1")
    write(P2, ROOT / "practice2_cohort_qc.ipynb", "p2")
