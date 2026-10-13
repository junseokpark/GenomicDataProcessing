"""Cells of practice3_dna_rna_short_long.ipynb (Mini practice 3 . One person, six datasets).

Imported by make_notebooks.py, which writes the notebook. Kept in its own file because Practice 3
is as long as Practices 1 and 2 together. The CHECK numbers printed by these cells must equal
checkpoints/practice3_checks.tsv, which build_practice3.py computes with independent code.
"""

import csv
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
P3 = []
md = lambda s: P3.append(("md", s))
code = lambda s: P3.append(("code", s))

# Slide titles of the 43-slide Day 2 deck, used for citations ("Title" (slide N)).
SLIDE = {
    2: "What is a genomics sample?",
    4: "Human specimens and what they give you",
    7: "Metadata: the sample's passport",
    9: "RNA sequencing",
    10: "RNA-seq applications",
    11: "One gene, many transcripts",
    13: "Three ways to read a molecule",
    14: "Phred quality: how sure is each base?",
    15: "Illumina: short, accurate, massively parallel",
    16: "Paired-end reads and insert size",
    17: "Duplicates and library complexity",
    18: "PacBio HiFi: long and accurate",
    19: "Oxford Nanopore: the longest reads",
    20: "Read length: mean versus N50",
    21: "Side by side: which platform for which question?",
    22: "Match the sample to the platform",
    24: "Coverage is uneven: the Poisson view",
    33: "Mini practice 3 · One person, six datasets",
    34: "Lightning talks: one figure, two minutes",
    35: "Debrief: what each dataset can and cannot see",
}


def cite(*numbers):
    return " · ".join(f'"{SLIDE[n]}" (slide {n})' for n in numbers)


TITLES = {
    "F01": "Read length: six datasets",
    "F02": "Errors: what kind of mistakes?",
    "F03": "DNA vs RNA: coverage across the locus",
    "F04": "Allele-specific expression",
    "F05": "Isoforms: junctions vs full-length reads",
    "F06": "RNA coverage along the transcript",
    "F07": "Sizing a CAG repeat",
    "F08": "A duplicated block",
    "F09": "Phasing: are two variants on the same copy?",
    "F10": "GC content and coverage",
    "F11": "Methylation from long-read DNA",
    "F12": "Gene expression: do the RNA datasets agree?",
}
SLUG = {
    "F01": "read_length", "F02": "errors", "F03": "dna_vs_rna_coverage", "F04": "allele_specific_expression",
    "F05": "isoforms", "F06": "transcript_coverage", "F07": "cag_repeat", "F08": "duplicated_block",
    "F09": "phasing", "F10": "gc_coverage", "F11": "methylation", "F12": "expression",
}
FIG_SLIDES = {
    "F01": (13, 18, 19, 20), "F02": (14, 15, 18, 19), "F03": (2, 4, 9), "F04": (9, 10), "F05": (10, 11),
    "F06": (9, 19), "F07": (13, 18, 19, 21), "F08": (15, 21), "F09": (16, 18, 19), "F10": (15, 17, 24),
    "F11": (18, 19, 21), "F12": (9, 10, 22),
}

TALK = '''**For your lightning talk** (the pair slide of the template; see **Make your slide** at the end)
* **CLAIM**: one sentence a classmate can check against this figure.
* **1 · Axes and colours**: what do the axes and colours show?
* **2 · Evidence: one number**: point at the figure and quote one number from `CHECK`.
* **3 · So what**: which dataset would you choose for this question, and one thing it cannot do?
* **Q / A**: a question you expect from the floor, and your one-sentence answer.'''


def f01_example_numbers():
    """Numbers for the instructor's F01 example slide, read from the build's own records so the
    text cannot drift from the data: read N50s from checkpoints/practice3_checks.tsv and the
    longest transcript from data/practice3/genes.gtf."""
    with open(PKG / "checkpoints" / "practice3_checks.tsv", newline="") as fh:
        n50 = {r["library"]: int(r["value"]) for r in csv.DictReader(fh, delimiter="\t")
               if r["figure"] == "F01" and r["metric"] == "read_n50"}
    tx = {}
    for line in (PKG / "data" / "practice3" / "genes.gtf").read_text().splitlines():
        f = line.split("\t")
        if f[2] == "exon":
            tid = f[8].split('transcript_id "')[1].split('"')[0]
            tx[tid] = tx.get(tid, 0) + int(f[4]) - int(f[3]) + 1
    return {"ont_kb": f"{n50['DNA-ONT'] / 1000:.1f}", "kinnex_kb": f"{n50['RNA-KINNEX'] / 1000:.1f}",
            "longest_nt": f"{max(tx.values()):,}"}


F01 = f01_example_numbers()
F01_EXAMPLE = {      # each card on the template slide holds two lines, so keep these short
    "claim": "Long DNA reads follow the length of the DNA molecule, but long RNA reads stop at the length of "
             "the transcript.",
    "axes": "Read length on a log scale, one row per dataset. Colour = platform; open circles = RNA; "
            "bar = N50, dotted = median.",
    "evidence": f"Nanopore DNA reads have an N50 of {F01['ont_kb']} kb; Kinnex RNA reads stop at "
                f"{F01['kinnex_kb']} kb, the longest transcript here ({F01['longest_nt']} nt).",
    "so_what": "To cross a 6-kb duplicated block in DNA, choose HiFi or Nanopore DNA. Limit: fewer reads, "
               "higher cost per base than Illumina.",
    "q": "Why are the long RNA reads not 20 kb long, like the long DNA reads?",
    "a": "Each RNA read is one copy of one transcript, so it cannot be longer than the transcript: at most "
         f"{F01['longest_nt']} nt here.",
}


def section(fid, question, biology, cell, try_this):
    """One figure section: heading, question, biology, the figure cell, try this, talk prompts, slides."""
    md(f'''## {fid} · {TITLES[fid]}

**Question.** {question}

**Biology first.** {biology}''')
    code(cell)
    md(f'''**Try this.** {try_this}

{TALK}

*Slides:* {cite(*FIG_SLIDES[fid])}''')


# =====================================================================================
# Title and set-up
# =====================================================================================
md(f'''# Mini practice 3 · One person, six datasets

**Genomic Data Processing · Day 2 (10/13) · 60 minutes in pairs** · {cite(33)} · {cite(34)}

One person gave one blood draw. Their DNA and their RNA were each read three ways: by Illumina
(short reads), by PacBio (long, accurate reads) and by Oxford Nanopore (the longest reads). That
makes six datasets of the same 120,000 bases. Each dataset sees some things clearly and misses
others. Your job is to find out which, and to explain one example in two minutes.

### What you will do in 25 minutes
1. Run the whole notebook (Kernel → Restart & Run All). It takes under 30 seconds and writes 12
   figures to `outputs/`.
2. Go to the figure your pair was given. Read its question, look at the figure and its `CHECK`
   numbers, and agree on an answer.
3. Do its **Try this**: change one parameter, re-run that one cell, and see what moves.
4. Fill in the pair slide of the template (see **Make your slide** at the end), then give a
   two-minute lightning talk.

### How pairs get a figure
Your instructor hands out figures in this order: F03, F05, F07, F04, F08, F02, F09, F06, F10,
F11, F12, F01. F01 is the instructor's worked example. With more than 12 pairs, two pairs share
F05 or F07 and compare their claims. The figures fall into four groups:

| Comparison | Figures |
|---|---|
| DNA vs RNA | F03, F04 |
| Short vs long reads, DNA | F07, F08, F09 |
| Short vs long reads, RNA | F05, F06, F12 |
| Platforms | F01, F02, F10, F11 |

### What is real and what is synthetic
| Part | Status |
|---|---|
| The person `P3-01`, the locus `chrS`, every read and every methylation call | **Synthetic**: simulated for teaching, no real person |
| Read lengths, error rates and biases of the six datasets | Approximate public platform specifications, but **not calibrated to any real run** |
| The file formats (SAM, GTF, FASTA, a bedMethyl-like table) | Real formats, readable by standard tools |

Nothing in this notebook is a clinical result. Because the data are simulated, the truth is
known: your instructor can tell you afterwards exactly what was put into the locus.

### How to work
* Run cells from top to bottom (Shift+Enter). Cells marked *run, no need to edit* hold helper
  code; you can read them, but you do not have to.
* **Predict, then run.** Before you run your figure cell, say what you expect to see.
* Every figure cell prints `CHECK` numbers. Quote one of them as your evidence.''')

code(r'''from pathlib import Path
from collections import Counter
import gzip
import io
import math
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA  = Path("data/practice3")     # open this notebook from the day2-3-samples/ folder
CHECK = Path("checkpoints")
OUT   = Path("outputs")
OUT.mkdir(exist_ok=True)

plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.labelsize": 13,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 12,
                     "legend.frameon": False, "axes.spines.top": False,
                     "axes.spines.right": False})

LIBS = ["DNA-ILMN", "DNA-HIFI", "DNA-ONT", "RNA-ILMN", "RNA-KINNEX", "RNA-ONT"]
DNA_LIBS, RNA_LIBS = LIBS[:3], LIBS[3:]
PLATFORM = {"DNA-ILMN": "Illumina", "RNA-ILMN": "Illumina", "DNA-HIFI": "PacBio",
            "RNA-KINNEX": "PacBio", "DNA-ONT": "ONT", "RNA-ONT": "ONT"}
COLOR = {"Illumina": "#1F5FEA", "PacBio": "#B25E09", "ONT": "#17202A"}


def look(lib):
    """The same style in every figure: colour = platform; DNA solid and filled, RNA dashed and open."""
    color = COLOR[PLATFORM[lib]]
    dna = lib.startswith("DNA")
    return {"color": color, "ls": "-" if dna else "--", "face": color if dna else "white"}


print("set-up done; figures go to", OUT.resolve())''')

# =====================================================================================
# The sample sheet
# =====================================================================================
md(f'''## 1 · The sample sheet: six datasets, one person

A dataset is only useful if you know where it came from. The sample sheet is the passport of each
library ({cite(7)}): which person, which tube, which molecule, which machine and which chemistry.

**Predict.** Which columns tell you that two libraries come from the same person? Which column
tells you whether a library was copied by PCR before sequencing?''')

code(r'''samples = pd.read_csv(DATA / "samples.tsv", sep="\t")
samples[["library_id", "person", "molecule", "specimen_tube", "platform", "chemistry",
         "library_prep", "reads_in_locus"]]''')

md('''**What to notice.** Every row has the same `person`, P3-01, but two different tubes: an EDTA
tube for DNA and a PAXgene tube that keeps RNA from breaking down. So the RNA tells us what blood
cells were doing at the moment of the draw. `reads_in_locus` counts SAM records: an Illumina
fragment gives two records (read 1 and read 2), a long read gives one.''')

# =====================================================================================
# A SAM file is a table
# =====================================================================================
md('''## 2 · A SAM file is a table

An aligner takes each read and records where it sits on the reference sequence. The result is a
SAM file: a few header lines starting with `@`, then **one line per read** with 11 fixed columns
separated by tabs. Here are the first lines of the Illumina DNA file:''')

code(r'''with gzip.open(DATA / "DNA-ILMN.sam.gz", "rt") as fh:
    for i, line in enumerate(fh):
        print(line.rstrip("\n")[:110] + (" ..." if len(line) > 110 else ""))
        if i == 7:
            break''')

md('''The columns we use, in plain words:

| Column | What it says | Example |
|---|---|---|
| `FLAG` | yes/no facts about the read, added up as a number: 1 = paired, 4 = **unmapped** (no place found), 16 = read on the reverse strand, 64 = read 1, 128 = read 2 | 99 = 1 + 2 + 32 + 64 |
| `POS` | the reference position of the first aligned base (counting from 1) | 4 |
| `MAPQ` | how sure the aligner is about the position: 60 = sure, **0 = fits two or more places equally well** | 60 |
| `CIGAR` | how the read lines up with the reference, as counts and letters (below) | `150M` |
| `SEQ` | the bases of the read, written in reference orientation | `CAGGCC...` |

CIGAR letters: **M** read base aligned to a reference base (same or different), **I** extra base
in the read (insertion), **D** reference base missing from the read (deletion), **N** reference
skipped by splicing (an intron, RNA only), **S** read bases at the end that are not aligned
(soft clip). The tag `NM` at the end of each line counts the differences from the reference:
mismatches + inserted bases + deleted bases.

### Walking a CIGAR: `5M2I3M`
A read with CIGAR `5M2I3M` at `POS` 100 and `SEQ` = `ACGTAGGTCA` (10 bases):
* `5M`: read bases 0-4 (`ACGTA`) sit on reference positions 100-104.
* `2I`: read bases 5-6 (`GG`) are extra; they have **no** reference position.
* `3M`: read bases 7-9 (`TCA`) sit on reference positions 105-107.

The helper functions below do this walk for every read. *Run, no need to edit.*''')

code(r'''# ---- helpers: run, no need to edit ----------------------------------------------------
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cigar):
    """'5M2I3M' -> [(5, 'M'), (2, 'I'), (3, 'M')]"""
    return [(int(n), op) for n, op in CIGAR_OP.findall(cigar)]


def aligned_blocks(pos, cigar):
    """Where a read sits on the reference: one (ref_start, ref_end, read_start) per M run.
    M moves along both; I and S move along the read only; D and N along the reference only."""
    blocks, ref, read = [], pos, 0
    for n, op in parse_cigar(cigar):
        if op == "M":
            blocks.append((ref, ref + n - 1, read))
            ref, read = ref + n, read + n
        elif op in "DN":
            ref += n
        elif op in "IS":
            read += n
    return blocks


def read_index_at(blocks, ref_pos):
    """Index (from 0) of the read base aligned to reference position ref_pos, or None."""
    for start, end, read in blocks:
        if start <= ref_pos <= end:
            return read + (ref_pos - start)
    return None


def base_at(seq, blocks, ref_pos):
    """The read's base at a reference position, or None if the read has no base there."""
    i = read_index_at(blocks, ref_pos)
    return None if i is None else seq[i]


def repeat_length(blocks, left_anchor, right_anchor):
    """Number of read bases between two reference positions (both must be aligned)."""
    a, b = read_index_at(blocks, left_anchor), read_index_at(blocks, right_anchor)
    return None if a is None or b is None else b - a - 1


example = aligned_blocks(100, "5M2I3M")
print("blocks:", example)
for ref_pos in (103, 104, 105, 106):
    print("reference", ref_pos, "-> read index", read_index_at(example, ref_pos),
          "base", base_at("ACGTAGGTCA", example, ref_pos))
print("read bases between reference 104 and 105:", repeat_length(example, 104, 105))''')

md('''**Check.** Reference 104 is read base 4 (`A`) and reference 105 is read base 7 (`T`): the two
inserted bases (`GG`) sit between them, so `repeat_length` says 2. The same function will count
the bases inside a repeat in F07.

Now load the six SAM files and the annotation. *Run, no need to edit* (about 3 seconds).''')

code(r'''# ---- load the six datasets and the annotation: run, no need to edit ---------------------
SAM_COLUMNS = ["qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen",
               "seq", "qual", "tags"]


def read_sam(path):
    """One row per SAM record, plus columns we compute: mapped, blocks, end, length, n_M, n_I, n_D."""
    with gzip.open(path, "rt") as fh:
        body = "".join(line for line in fh if not line.startswith("@"))
    df = pd.read_csv(io.StringIO(body), sep="\t", header=None, names=SAM_COLUMNS,
                     quoting=3, keep_default_na=False, dtype=str)
    for col in ("flag", "pos", "mapq", "pnext", "tlen"):
        df[col] = df[col].astype(int)
    df["nm"] = pd.to_numeric(df["tags"].str[5:], errors="coerce")     # "NM:i:3" -> 3
    df["mapped"] = (df["flag"] & 4) == 0
    df["length"] = df["seq"].str.len()
    df["blocks"] = [aligned_blocks(p, c) if m else [] for p, c, m in
                    zip(df["pos"], df["cigar"], df["mapped"])]
    df["end"] = [b[-1][1] if b else 0 for b in df["blocks"]]
    totals = [Counter() for _ in range(len(df))]
    for t, cigar in zip(totals, df["cigar"]):
        for n, op in parse_cigar(cigar):
            t[op] += n
    for op in "MID":
        df["n_" + op] = [t[op] for t in totals]
    return df


def depth(df, min_mapq=0):
    """Read depth at every position: depth[p] for p = 1 .. 120,000 (depth[0] is unused)."""
    use = df[df["mapped"] & (df["mapq"] >= min_mapq)]
    starts = [s for blocks in use["blocks"] for s, e, r in blocks]
    ends = [e for blocks in use["blocks"] for s, e, r in blocks]
    change = np.zeros(LOCUS_LEN + 2)
    np.add.at(change, starts, 1)
    np.add.at(change, np.array(ends, dtype=int) + 1, -1)
    return np.cumsum(change)[:LOCUS_LEN + 1]


def n50(lengths):
    """Half of all sequenced bases are in reads at least this long."""
    s = sorted(lengths, reverse=True)
    half, run = sum(s) / 2, 0
    for x in s:
        run += x
        if run >= half:
            return x


def wilson(k, n, z=1.96):
    """95% interval for a fraction k / n (Wilson score interval)."""
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def kb(x):
    return f"{x / 1000:.1f} kb" if x >= 1000 else f"{int(x)} bp"


# the reference sequence, genes, special regions and the two variants
fasta = (DATA / "chrS.fa").read_text().splitlines()
CHRS = "".join(line.strip() for line in fasta[1:])
LOCUS_LEN = len(CHRS)
gtf = pd.read_csv(DATA / "genes.gtf", sep="\t", header=None,
                  names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"])
gtf["gene"] = gtf["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
gtf["transcript"] = gtf["attributes"].str.extract(r'transcript_id "([^"]+)"', expand=False)
gtf["exon"] = gtf["attributes"].str.extract(r'exon_id "[^"-]+-([^"]+)"', expand=False)
exons = gtf[gtf["feature"] == "exon"]
GENE_ORDER = ["ASE1", "ALTX", "LOW1", "OFF1"]
GENE_SPAN = {g: (int(x["start"].min()), int(x["end"].max())) for g, x in exons.groupby("gene")}
GENE_EXONS = {g: sorted(set(zip(x["start"], x["end"]))) for g, x in exons.groupby("gene")}
EXON_POS = {g: np.unique(np.concatenate([np.arange(s, e + 1) for s, e in GENE_EXONS[g]]))
            for g in GENE_ORDER}
TX_EXONS = {t: list(zip(x["start"], x["end"])) for t, x in exons.groupby("transcript")}
FIRST_TX = {g: sorted(exons.loc[exons["gene"] == g, "transcript"].unique())[0] for g in GENE_ORDER}
TX_LENGTH = {g: sum(e - s + 1 for s, e in TX_EXONS[FIRST_TX[g]]) for g in GENE_ORDER}
regions = pd.read_csv(DATA / "regions.tsv", sep="\t").set_index("name")
sites = pd.read_csv(DATA / "sites.tsv", sep="\t").set_index("name")

sam = {lib: read_sam(DATA / f"{lib}.sam.gz") for lib in LIBS}
DEPTH = {lib: depth(sam[lib]) for lib in LIBS}

overview = pd.DataFrame({
    "records": [len(sam[lib]) for lib in LIBS],
    "unmapped": [int((~sam[lib]["mapped"]).sum()) for lib in LIBS],
    "MAPQ_0": [int((sam[lib]["mapped"] & (sam[lib]["mapq"] == 0)).sum()) for lib in LIBS],
    "mean_depth": [round(float(DEPTH[lib][1:].mean()), 1) for lib in LIBS]}, index=LIBS)
assert (overview["records"] == samples.set_index("library_id").loc[LIBS, "reads_in_locus"]).all()
print(f"chrS: {LOCUS_LEN:,} bases; genes: {', '.join(GENE_ORDER)}")
print(overview.to_string())

def draw_genes(ax, xlim):
    """Gene models: a line for each gene, a box for each exon, the gene name above."""
    for g in GENE_ORDER:
        s, e = GENE_SPAN[g]
        if e < xlim[0] or s > xlim[1]:
            continue
        ax.plot([max(s, xlim[0]), min(e, xlim[1])], [0, 0], color="0.3", lw=1)
        for es, ee in GENE_EXONS[g]:
            ax.add_patch(plt.Rectangle((es, -0.3), ee - es + 1, 0.6, color="0.3", lw=0))
        ax.text((max(s, xlim[0]) + min(e, xlim[1])) / 2, 0.45, g, ha="center", va="bottom",
                fontsize=12, style="italic")
    ax.set_ylim(-0.6, 1.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


# ---- CHECK numbers: every figure stores its numbers here and prints them -------------------
CHECKS = {}


def check(fig, lib, metric, value):
    CHECKS[(fig, lib, metric)] = value


def show_checks(fig):
    """Print this figure's CHECK numbers as a table and save all CHECK numbers so far."""
    rows = [(lib, metric, value) for (f, lib, metric), value in CHECKS.items() if f == fig]
    metrics = list(dict.fromkeys(m for _, m, _ in rows))
    libs = [lib for lib in LIBS if any(l == lib for l, _, _ in rows)]
    table = pd.DataFrame("", index=metrics, columns=libs, dtype=object)
    for lib, metric, value in rows:
        table.loc[metric, lib] = value
    print(f"CHECK {fig}")
    print(table.to_string())
    pd.DataFrame([(f, l, m, v) for (f, l, m), v in CHECKS.items()],
                 columns=["figure", "library", "metric", "value"]).to_csv(
        OUT / "practice3_checks.tsv", sep="\t", index=False)


def title(fig, fid, text):
    fig.suptitle(f"{fid} · {text}", x=0.01, ha="left", fontsize=13)


def save(fig, fid, slug):
    path = OUT / f"practice3_{fid}_{slug}.png"
    fig.savefig(path, dpi=200)
    print("saved", path)''')

md('''**Check.** The table above should list 24,000 records for DNA-ILMN (12,000 read pairs) and
20,000 for RNA-ILMN (10,000 fragments). In DNA-ILMN, 31 reads are unmapped and about 2,300 have
MAPQ 0: F07 and F08 explain why. Mean depth is about 30 for Illumina DNA and about 20 for the
long-read DNA libraries.''')

# =====================================================================================
# F01-F12
# =====================================================================================
section(
    "F01",
    "How long are the reads, and why are long RNA reads shorter than long DNA reads?",
    "Illumina reads a fixed number of bases from each end of a fragment: the run is set to 150 or "
    "100 cycles. PacBio and Nanopore read one molecule from end to end, so the read is as long as "
    "the molecule that went into the machine. Extracted DNA can stay tens of kilobases long. What "
    "sets the length of an RNA molecule?",
    r'''LOG_X = True        # Try this: False puts the same reads on a linear axis

fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
rng = np.random.default_rng(1)                 # fixed seed: the same dots every time
labels = []
for i, lib in enumerate(LIBS):
    lengths = sam[lib]["length"].to_numpy()
    median, read_n50 = float(np.median(lengths)), n50(lengths)
    check("F01", lib, "n_reads", len(lengths))
    check("F01", lib, "median_length", round(median, 1))
    check("F01", lib, "read_n50", int(read_n50))
    st, row = look(lib), len(LIBS) - 1 - i    # first library on top
    shown = rng.choice(lengths, size=min(len(lengths), 800), replace=False)
    ax.scatter(shown, row + rng.uniform(-0.3, 0.3, len(shown)), s=12, facecolors=st["face"],
               edgecolors=st["color"], linewidths=0.7, alpha=0.4)
    ax.plot([median, median], [row - 0.42, row + 0.42], color=st["color"], lw=1.5, ls=":")
    ax.plot([read_n50, read_n50], [row - 0.42, row + 0.42], color=st["color"], lw=3.5)
    labels.append(f"{lib}\nN50 {kb(read_n50)}")

longest = max(TX_LENGTH.values())              # the longest transcript in genes.gtf
ax.axvline(longest, color="0.55", lw=1.2, ls="--", zorder=0)
ax.text(longest * (1.08 if LOG_X else 1) + (0 if LOG_X else 1_000), len(LIBS) - 0.55,
        f"longest transcript\n{longest:,} nt", color="0.35", fontsize=12, va="top")
ax.plot([], [], color="0.3", lw=3.5, label="N50")
ax.plot([], [], color="0.3", lw=1.5, ls=":", label="median")
ax.legend(loc="lower right")
ax.set_yticks(range(len(LIBS)))
ax.set_yticklabels(labels[::-1])
ax.set_ylim(-0.6, len(LIBS) - 0.4)
if LOG_X:
    ax.set_xscale("log")
    ax.set_xlim(60, 200_000)
    ax.set_xticks([100, 1_000, 10_000, 100_000])
    ax.set_xticklabels(["100 bp", "1 kb", "10 kb", "100 kb"])
    ax.xaxis.set_minor_formatter(plt.NullFormatter())
else:
    ax.set_xlim(-2_000, 105_000)
    ax.set_xticks([0, 25_000, 50_000, 75_000, 100_000])
    ax.set_xticklabels(["0", "25 kb", "50 kb", "75 kb", "100 kb"])
ax.set_xlabel("read length (one dot per read, up to 800 per dataset)")
title(fig, "F01", "Read length: six datasets")
save(fig, "F01", "read_length")
plt.show()
show_checks("F01")''',
    "Set `LOG_X = False` and re-run. What happens to the two Illumina rows, and why do we use a "
    "log axis for read length?",
)

section(
    "F02",
    "Which datasets make mostly insertion/deletion errors, and where do they show up?",
    "A sequencing error can be a wrong base (mismatch), an extra base (insertion) or a missing base "
    "(deletion). Runs of one base, called homopolymers, are hard to count when a machine reads a "
    "signal as the molecule moves past. At HP1 the person really has ten A's on both copies, so a "
    "read showing 9 or 11 is wrong. Reads that touch the person's real variants (REP1, SNV1, SNV2) "
    "are left out of the error count.",
    r'''SHOW = LIBS        # Try this: SHOW = ["DNA-ONT", "RNA-ONT"]

rep, hp1 = regions.loc["REP1"], regions.loc["HP1"]
KINDS = [("mismatch", "o"), ("insertion", "^"), ("deletion", "s")]
fig, (ax, ax_hp) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True,
                                gridspec_kw={"width_ratios": [1.5, 1]})
ylabels = []
for i, lib in enumerate(SHOW):
    df = sam[lib]
    touches = (df["pos"] <= rep["end"]) & (df["end"] >= rep["start"])
    for name in ("SNV1", "SNV2"):
        touches |= (df["pos"] <= sites.loc[name, "pos"]) & (df["end"] >= sites.loc[name, "pos"])
    use = df[df["mapped"] & ~touches]
    aligned, ins, dels = use["n_M"].sum(), use["n_I"].sum(), use["n_D"].sum()
    mism = (use["nm"] - use["n_I"] - use["n_D"]).sum()      # NM counts all three kinds
    per_kb = {"mismatch": 1000 * mism / aligned, "insertion": 1000 * ins / aligned,
              "deletion": 1000 * dels / aligned}
    accuracy = round(float(100 * (1 - (mism + ins + dels) / aligned)), 2)   # rounded once: label = CHECK = scorecard
    check("F02", lib, "reads_used", len(use))
    for kind, _ in KINDS:
        check("F02", lib, f"{kind}_per_kb", round(float(per_kb[kind]), 3))
    check("F02", lib, "accuracy_pct", accuracy)
    st, row = look(lib), len(SHOW) - 1 - i
    for k, (kind, marker) in enumerate(KINDS):
        ax.plot(per_kb[kind], row + (1 - k) * 0.22, marker=marker, ms=9, ls="none",
                color=st["color"], mfc=st["face"], mew=1.5)
    ylabels.append(f"{lib}\n{accuracy:.2f}% correct")
for kind, marker in KINDS:
    ax.plot([], [], marker=marker, ls="none", ms=8, color="0.35", mfc="0.35", label=kind)
ax.legend(loc="upper right", handletextpad=0.3)
ax.set_yticks(range(len(SHOW)))
ax.set_yticklabels(ylabels[::-1])
ax.set_ylim(-0.6, len(SHOW) - 0.4)
ax.set_xlim(left=-0.6)
ax.set_xlabel("errors per 1,000 aligned bases")

for lib in [l for l in SHOW if l in DNA_LIBS]:              # HP1 lies between genes: DNA only
    df = sam[lib][sam[lib]["mapped"]]
    n_a = [repeat_length(b, hp1["start"] - 1, hp1["end"] + 1) for b in df["blocks"]]
    n_a = np.array([x for x in n_a if x is not None])
    check("F02", lib, "hp1_reads", len(n_a))
    for k in (9, 10, 11):
        check("F02", lib, f"hp1_pct_{k}A", round(100 * float((n_a == k).mean()), 1))
    xs = np.arange(7, 14)
    st = look(lib)
    ax_hp.plot(xs, [100 * float((n_a == x).mean()) for x in xs], marker="o", ms=7, lw=1.8,
               color=st["color"], mfc=st["face"], ls=st["ls"], label=lib)
ax_hp.set_xticks(range(7, 14))
ax_hp.set_xlabel("A's read at HP1 (truth: 10)")
ax_hp.set_ylabel("% of reads")
ax_hp.set_ylim(-3, 140)                                     # head room for the legend
ax_hp.set_yticks(range(0, 101, 20))
ax_hp.legend(loc="upper right", fontsize=12, handlelength=1.5)
title(fig, "F02", "Errors: what kind of mistakes?")
save(fig, "F02", "errors")
plt.show()
show_checks("F02")''',
    "Set `SHOW = [\"DNA-ONT\", \"RNA-ONT\"]`. Same platform, different molecule: which kind of "
    "error grows most when Nanopore reads RNA directly?",
)

section(
    "F03",
    "OFF1 has no RNA reads. Is the gene missing from this person's genome?",
    "Every cell carries the same genome, so DNA reads should land evenly along the locus: on exons, "
    "introns and the DNA between genes. RNA reads come from transcripts, so they pile up on the "
    "exons of genes that the sampled cells (here, blood cells) are using. The RNA axis is "
    "logarithmic because RNA depth spans from 0 to more than 1,000.",
    r'''XLIM = (0, 120_000)     # Try this: XLIM = (37_000, 49_500) zooms in on ASE1

BIN = 500 if XLIM[1] - XLIM[0] > 20_000 else 50
centres = np.arange(LOCUS_LEN // BIN) * BIN + BIN / 2
all_exons = np.concatenate([EXON_POS[g] for g in GENE_ORDER])
fig, (ax_dna, ax_rna, ax_gene) = plt.subplots(3, 1, figsize=(8, 4.5), sharex=True,
                                              constrained_layout=True,
                                              gridspec_kw={"height_ratios": [2.2, 2.2, 1.1]})
for lib in LIBS:
    d = DEPTH[lib]
    for g in GENE_ORDER:
        check("F03", lib, f"{g}_exon_depth", round(float(d[EXON_POS[g]].mean()), 2))
    check("F03", lib, "pct_bases_in_exons", round(100 * float(d[all_exons].sum() / d[1:].sum()), 1))
    st = look(lib)
    ax = ax_dna if lib in DNA_LIBS else ax_rna
    ax.plot(centres, d[1:].reshape(-1, BIN).mean(axis=1), color=st["color"], ls=st["ls"],
            lw=1.5, label=lib)
ax_dna.set_ylim(0, 45)
ax_dna.set_ylabel("DNA\ndepth")
ax_rna.set_yscale("symlog", linthresh=1)
ax_rna.set_ylim(0, 4_000)
ax_rna.set_yticks([0, 1, 10, 100, 1000])
ax_rna.set_yticklabels(["0", "1", "10", "100", "1000"])
ax_rna.set_ylabel("RNA\ndepth")
for ax in (ax_dna, ax_rna):
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1), ncol=3, borderaxespad=0.1, handlelength=2.5)
draw_genes(ax_gene, XLIM)
ax_gene.set_xlim(*XLIM)
ax_gene.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1000:g}"))
ax_gene.set_xlabel("position on chrS (kb)")
title(fig, "F03", "DNA vs RNA: coverage across the locus")
save(fig, "F03", "dna_vs_rna_coverage")
plt.show()
show_checks("F03")''',
    "Set `XLIM = (37_000, 49_500)` to zoom in on ASE1. Where do the RNA reads sit, and where do the "
    "DNA reads sit?",
)

section(
    "F04",
    "Why do DNA and RNA give different allele fractions at SNV1?",
    "At a heterozygous site one copy of the chromosome carries the reference base (G at SNV1) and "
    "the other copy carries the alternative base (A). DNA reads sample the two copies of the genome. "
    "RNA reads sample transcripts instead: how many come from each copy depends on how much each "
    "copy is transcribed. SNV1 sits in an exon of ASE1, so both DNA and RNA reads cover it.",
    r'''MIN_MAPQ = 0      # Try this: 20, then 60 (keep only reads the aligner placed confidently)

snv = sites.loc["SNV1"]
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
for i, lib in enumerate(LIBS):
    df = sam[lib][sam[lib]["mapped"] & (sam[lib]["mapq"] >= MIN_MAPQ)]
    bases = [base_at(s, b, snv["pos"]) for s, b in zip(df["seq"], df["blocks"])]
    n_ref, n_alt = bases.count(snv["ref"]), bases.count(snv["alt"])
    frac = n_alt / (n_ref + n_alt)
    lo, hi = wilson(n_alt, n_ref + n_alt)
    check("F04", lib, "n_ref", n_ref)
    check("F04", lib, "n_alt", n_alt)
    check("F04", lib, "alt_fraction", round(frac, 3))
    check("F04", lib, "ci_low", round(lo, 3))
    check("F04", lib, "ci_high", round(hi, 3))
    st, row = look(lib), len(LIBS) - 1 - i
    bars = ax.errorbar(frac, row, xerr=[[frac - lo], [hi - frac]], fmt="o", ms=10, color=st["color"],
                       mfc=st["face"], mew=1.8, elinewidth=2, capsize=5)
    bars[2][0].set_linestyle(st["ls"])
    ax.annotate(f"{n_alt} of {n_ref + n_alt} reads", (lo, row), xytext=(-12, 0),
                textcoords="offset points", ha="right", va="center", fontsize=12, color=st["color"])
ax.axvline(0.5, color="0.55", lw=1.2, ls=":", zorder=0)
ax.text(0.51, len(LIBS) - 0.45, "0.5: as many reads from each copy", ha="left", fontsize=12,
        color="0.35")
ax.set_xlim(0, 1)
ax.set_ylim(-0.6, len(LIBS))
ax.set_yticks(range(len(LIBS)))
ax.set_yticklabels(LIBS[::-1])
ax.set_xlabel("fraction of reads with the alternative base A at SNV1 (95% interval)")
title(fig, "F04", "Allele-specific expression")
save(fig, "F04", "allele_specific_expression")
plt.show()
show_checks("F04")''',
    "Set `MIN_MAPQ = 20`, then `60`. Does the answer depend on the mapping-quality filter here? "
    "(Compare with F08.)",
)

section(
    "F05",
    "From short reads alone, how many ALTX isoforms are possible? What do long reads add?",
    "Splicing joins exons and removes introns; one gene can join different exons into different "
    "mRNAs (isoforms). ALTX has two alternative exon pairs: E2a or E2b, and E4a or E4b. A short read "
    "that crosses a junction proves that two exons were joined, but it covers one junction at a "
    "time. A full-length read shows the whole chain of exons of one molecule.",
    r'''MIN_OVERHANG = 8     # Try this: 30 (count a junction only with 30 aligned bases on each side)

altx = exons[exons["gene"] == "ALTX"].drop_duplicates("exon").set_index("exon")
EX = {name: (int(r["start"]), int(r["end"])) for name, r in altx.iterrows()}
JUNCTIONS = [("E1", "E2a"), ("E1", "E2b"), ("E2a", "E3"), ("E2b", "E3"),
             ("E3", "E4a"), ("E3", "E4b"), ("E4a", "E5"), ("E4b", "E5")]
lookup = {(EX[a][1], EX[b][0]): (a, b) for a, b in JUNCTIONS}

# short reads: count each N (intron) that joins two annotated ALTX exons
junction_reads = Counter()
df = sam["RNA-ILMN"][sam["RNA-ILMN"]["mapped"]]
for pos, cigar in zip(df["pos"], df["cigar"]):
    ops, ref = parse_cigar(cigar), pos
    for k, (n, op) in enumerate(ops):
        if op == "N" and 0 < k < len(ops) - 1:
            key = (ref - 1, ref + n)              # last base before the intron, first base after
            left, right = ops[k - 1], ops[k + 1]
            if key in lookup and left[1] == "M" and right[1] == "M" and min(left[0], right[0]) >= MIN_OVERHANG:
                junction_reads[lookup[key]] += 1
        if op in "MDN":
            ref += n
for a, b in JUNCTIONS:
    check("F05", "RNA-ILMN", f"junction_{a}_{b}", junction_reads[(a, b)])

# every RNA dataset: which alternative exons does each read contain (at least 20 aligned bases)?
CHAINS = ["E2a+E4a", "E2b+E4b", "E2a+E4b", "E2b+E4a"]
chains = {}
for lib in RNA_LIBS:
    df = sam[lib]
    df = df[df["mapped"] & (df["end"] >= GENE_SPAN["ALTX"][0]) & (df["pos"] <= GENE_SPAN["ALTX"][1])]
    c = Counter()
    for blocks in df["blocks"]:
        has = {x: sum(max(0, min(e, EX[x][1]) - max(s, EX[x][0]) + 1) for s, e, _ in blocks) >= 20
               for x in ("E2a", "E2b", "E4a", "E4b")}
        two = [x for x in ("E2a", "E2b") if has[x]]
        four = [x for x in ("E4a", "E4b") if has[x]]
        c[f"{two[0]}+{four[0]}" if len(two) == 1 and len(four) == 1 else "incomplete"] += 1
    for name in CHAINS + ["incomplete"]:
        check("F05", lib, "chain_" + name.replace("+", "_"), c[name])
    chains[lib] = c

fig, (ax_j, ax_c) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [1.45, 1]})
st = look("RNA-ILMN")
top = max(junction_reads.values()) if junction_reads else 1
for a, b in JUNCTIONS:                 # junctions of E2a / E4a above the exons, of E2b / E4b below
    up = 1 if ("E2a" in (a, b) or "E4a" in (a, b)) else -1
    x0, x1 = EX[a][1], EX[b][0]
    xs = np.linspace(x0, x1, 60)
    n = junction_reads[(a, b)]
    ax_j.plot(xs, up * (0.15 + 0.9 * np.sin(np.pi * (xs - x0) / (x1 - x0))), color=st["color"],
              ls=st["ls"], lw=1.2 + 3 * n / top)
    ax_j.text((x0 + x1) / 2, up * 1.12, str(n), ha="center", va="bottom" if up > 0 else "top",
              fontsize=12, color=st["color"])
for name, (s, e) in EX.items():
    ax_j.add_patch(plt.Rectangle((s, -0.15), e - s + 1, 0.3, facecolor="0.8", edgecolor="0.3", lw=1))
ax_j.set_xticks([(s + e) / 2 for s, e in EX.values()])
ax_j.set_xticklabels(list(EX))
ax_j.set_xlim(EX["E1"][0] - 300, EX["E5"][1] + 300)
ax_j.set_ylim(-1.55, 1.55)
ax_j.set_yticks([])
ax_j.spines["left"].set_visible(False)
ax_j.set_xlabel("ALTX exons (a-exon joins above, b below)")
ax_j.set_title("RNA-ILMN: reads across each junction", loc="left", fontsize=12)

base = np.arange(len(CHAINS))[::-1]
for j, lib in enumerate(RNA_LIBS):
    st = look(lib)
    ys = base + (1 - j) * 0.26
    vals = [chains[lib][name] for name in CHAINS]
    ax_c.barh(ys, vals, height=0.24, facecolor="white", edgecolor=st["color"], lw=1.6, ls=st["ls"],
              label=lib)
    for y, v in zip(ys, vals):
        ax_c.text(v + 5, y, str(v), va="center", fontsize=11, color=st["color"])
ax_c.set_yticks(base)
ax_c.set_yticklabels(CHAINS)
ax_c.set_xlim(0, max(max(c[n] for n in CHAINS) for c in chains.values()) * 1.25 + 10)
ax_c.set_xlabel("reads with this exon chain")
ax_c.legend(loc="lower right", fontsize=12, handlelength=1.5)
ax_c.set_title("full chains of alternative exons", loc="left", fontsize=12)
title(fig, "F05", "Isoforms: junctions vs full-length reads")
save(fig, "F05", "isoforms")
plt.show()
show_checks("F05")''',
    "Set `MIN_OVERHANG = 30` and re-run. Junction counts drop; do the 60/40-type proportions "
    "between the alternative exons change?",
)

section(
    "F06",
    "Why does one RNA dataset have more reads near the 3' end?",
    "If every molecule were read from end to end, coverage along a transcript would be flat. Library "
    "chemistry decides where reads start and stop: Illumina reads cDNA fragments, Kinnex reads "
    "full-length cDNA copies, and direct RNA sequencing feeds the native molecule into the pore tail "
    "first. An mRNA's tail is its poly(A) end, the 3' end.",
    r'''GENE = "ASE1"      # Try this: GENE = "ALTX"

tx = FIRST_TX[GENE]
pos = np.concatenate([np.arange(s, e + 1) for s, e in TX_EXONS[tx]])     # 5' to 3' (+ strand)
index = {p: i for i, p in enumerate(pos)}
T, k = len(pos), int(len(pos) * 0.2)                                      # k = 20% of the length
x = (np.arange(T) + 0.5) / T * 100
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
ax.axvspan(0, 20, color="0.93", zorder=0)
ax.axvspan(80, 100, color="0.93", zorder=0)
for lib in RNA_LIBS:
    cov = DEPTH[lib][pos]
    check("F06", lib, "ratio_3to5", round(float(cov[-k:].mean() / cov[:k].mean()), 2))
    if lib == "RNA-ILMN":
        check("F06", lib, "pct_reads_reach_5prime", "NA")      # fragments: start anywhere by design
    else:
        firsts = []
        for blocks in sam[lib].loc[sam[lib]["mapped"], "blocks"]:
            hits = [max(s, es) for s, e, _ in blocks for es, ee in TX_EXONS[tx] if s <= ee and e >= es]
            if hits:
                firsts.append(index[min(hits)])                   # most 5' base of the read
        reach = int(np.sum(np.array(firsts) < 50))               # starts within 50 nt of the 5' end
        check("F06", lib, "pct_reads_reach_5prime", round(100 * reach / len(firsts), 1))
    st = look(lib)
    smooth = np.convolve(cov / cov.mean(), np.ones(25) / 25, mode="same")   # 25-nt running mean
    ax.plot(x[12:-12], smooth[12:-12], color=st["color"], ls=st["ls"], lw=2.2, label=lib)
ax.text(10, 0.03, "first 20%", ha="center", fontsize=12, color="0.35")
ax.text(90, 0.03, "last 20%", ha="center", fontsize=12, color="0.35")
ax.set_xlim(0, 100)
ax.set_ylim(0, None)
ax.legend(loc="upper left")
ax.set_xlabel(f"position along the {tx} transcript (% of its length, 5' to 3')")
ax.set_ylabel("coverage ÷ mean coverage")
title(fig, "F06", "RNA coverage along the transcript")
save(fig, "F06", "transcript_coverage")
plt.show()
show_checks("F06")''',
    "Set `GENE = \"ALTX\"`. Does the same dataset show the same 3' bias on a different gene?",
)

section(
    "F07",
    "Is one allele expanded, and which dataset can measure it?",
    "A short tandem repeat such as (CAG)n can grow from one generation to the next; long CAG "
    "expansions cause diseases such as Huntington's disease. To size a repeat, a read must cover "
    "the whole repeat plus unique DNA on both sides (the flanks). The reference genome has 20 CAG "
    "units at REP1; the person has two alleles, which may differ.",
    r'''FLANK = 10      # Try this: 30 (a read must also cover 30 unique bases on each side)


def cag_identity(seq):
    """Share of the read that matches a pure CAG repeat (or CTG: the same repeat, other strand)."""
    best = 0.0
    for motif in ("CAG", "CTG"):
        if seq.count(motif) < 0.25 * len(seq):        # quick skip: clearly not a repeat
            continue
        for shift in range(3):
            pattern = (motif * (len(seq) // 3 + 3))[shift:shift + len(seq)]
            best = max(best, sum(a == b for a, b in zip(seq, pattern)) / len(seq))
    return best


rep = regions.loc["REP1"]
ref_units = (rep["end"] - rep["start"] + 1) // 3
rng = np.random.default_rng(7)
fig, (ax_s, ax_l) = plt.subplots(1, 2, figsize=(8, 4.5), sharey=True, constrained_layout=True)
for i, lib in enumerate(DNA_LIBS):
    units = []
    for blocks in sam[lib].loc[sam[lib]["mapped"], "blocks"]:
        if read_index_at(blocks, rep["start"] - FLANK) is None or read_index_at(blocks, rep["end"] + FLANK) is None:
            continue
        n = repeat_length(blocks, rep["start"] - 1, rep["end"] + 1)
        if n is not None:
            units.append(n / 3)
    units = np.array(units)
    all_cag = sum(len(s) >= 30 and cag_identity(s) >= 0.95 for s in sam[lib]["seq"])
    check("F07", lib, "spanning_reads", len(units))
    for name, sel in (("short", units < 100), ("long", units >= 100)):
        u = units[sel]
        check("F07", lib, f"{name}_n", len(u))
        if len(u):
            values, counts = np.unique(np.floor(u + 0.5).astype(int), return_counts=True)
            check("F07", lib, f"{name}_mode", int(values[np.argmax(counts)]))
            check("F07", lib, f"{name}_sd", round(float(u.std()), 2))
        else:
            check("F07", lib, f"{name}_mode", "NA")
            check("F07", lib, f"{name}_sd", "NA")
    check("F07", lib, "all_cag_reads", all_cag)
    st, row = look(lib), len(DNA_LIBS) - 1 - i
    y = row + rng.uniform(-0.25, 0.25, len(units))
    for ax in (ax_s, ax_l):
        ax.scatter(units, y, s=45, color=st["color"], edgecolors="white", linewidths=0.6, zorder=3)
    if (units >= 100).sum() == 0:
        ax_l.text(180, row, f"no {lib} read spans\nthis allele; {all_cag} reads\nare all CAG (unmapped)",
                  ha="center", va="center", fontsize=12, color=st["color"])
ax_s.axvline(ref_units, color="0.55", lw=1.2, ls="--", zorder=0)
ax_s.text(ref_units - 0.4, len(DNA_LIBS) - 0.55, f"reference\ngenome: {ref_units}", ha="right",
          va="top", fontsize=12, color="0.35")
ax_s.set_xlim(14, 36)
ax_l.set_xlim(166, 194)
ax_s.set_yticks(range(len(DNA_LIBS)))
ax_s.set_yticklabels(DNA_LIBS[::-1])
ax_s.set_ylim(-0.6, len(DNA_LIBS) - 0.4)
ax_s.set_title("14-36 units", loc="left", fontsize=12)
ax_l.set_title("166-194 units", loc="left", fontsize=12)
fig.supxlabel("CAG units in each spanning read (bases between the flanks ÷ 3)", fontsize=13)
title(fig, "F07", "Sizing a CAG repeat")
save(fig, "F07", "cag_repeat")
plt.show()
show_checks("F07")''',
    "Set `FLANK = 30`. Which dataset loses spanning reads, and why does a longer flank requirement "
    "cost short reads more than long reads?",
)

section(
    "F08",
    "Why do short reads disappear here when you filter on mapping quality?",
    "Some DNA blocks exist twice in the genome as near-identical copies (segmental duplications). "
    "Here DUP_A and DUP_B are 6 kb long and differ at only 3 bases. An aligner reports how sure it "
    "is about a read's position as MAPQ, and most analyses drop reads with low MAPQ.",
    r'''MIN_MAPQ = 20       # Try this: 0 (keep reads that fit both copies equally well)

dup_a, dup_b = regions.loc["DUP_A"], regions.loc["DUP_B"]
a0, a1 = int(dup_a["start"]), int(dup_a["end"])
inside = np.arange(a0, a1 + 1)
flanks = np.r_[a0 - 3_000:a0, a1 + 1:a1 + 3_001]           # 3 kb of unique DNA on each side
BIN, XLIM = 500, (50_000, 95_000)
centres = np.arange(LOCUS_LEN // BIN) * BIN + BIN / 2
fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
for name, r in (("DUP_A", dup_a), ("DUP_B (copy)", dup_b)):
    ax.axvspan(r["start"], r["end"], color="0.9", zorder=0)
    ax.text((r["start"] + r["end"]) / 2, 46, name, ha="center", va="top", fontsize=12, color="0.3")
for lib in DNA_LIBS:
    d = depth(sam[lib], MIN_MAPQ)
    a, f = float(d[inside].mean()), float(d[flanks].mean())
    check("F08", lib, "dupA_depth", round(a, 2))
    check("F08", lib, "flank_depth", round(f, 2))
    check("F08", lib, "ratio_inside_flank", round(a / f, 3))
    st = look(lib)
    ax.plot(centres, d[1:].reshape(-1, BIN).mean(axis=1), color=st["color"], ls=st["ls"], lw=2,
            label=lib)
ax.set_xlim(*XLIM)
ax.set_ylim(0, 48)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1000:g}"))
ax.set_xlabel("position on chrS (kb)")
ax.set_ylabel(f"depth, reads with MAPQ ≥ {MIN_MAPQ}")
ax.legend(loc="lower left", bbox_to_anchor=(0, 1), ncol=3, borderaxespad=0.1)
title(fig, "F08", "A duplicated block")
save(fig, "F08", "duplicated_block")
plt.show()
show_checks("F08")''',
    "Set `MIN_MAPQ = 0`. The Illumina depth comes back: are those reads now placed correctly? "
    "(Look at how MAPQ 0 is defined in section 2, and at "
    + cite(15) + ".)",
)

section(
    "F09",
    "Are SNV1 and SNV2 on the same copy of the chromosome?",
    "Each person has two copies of each chromosome. When a gene carries two heterozygous variants, "
    "it matters whether they sit on the same copy (cis) or on different copies (trans): two "
    "damaging variants in trans leave no intact copy of the gene. Only a read, or a read pair, that "
    "covers both sites shows which bases travel together. SNV2 and SNV1 are 3,000 bases apart.",
    r'''COUNT_PAIRS = False     # Try this: True (count Illumina read pairs with one mate on each SNV)

s1, s2 = sites.loc["SNV1"], sites.loc["SNV2"]
COMBOS = [("alt", "ref"), ("ref", "alt"), ("ref", "ref"), ("alt", "alt")]
lo_x, hi_x = s2["pos"] - 1_500, s1["pos"] + 1_500
fig, (ax_r, ax_c) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True,
                                 gridspec_kw={"width_ratios": [1.2, 1]})
row, centres = 0, []
for lib in DNA_LIBS:
    df = sam[lib][sam[lib]["mapped"]]
    combos = Counter()
    for seq, blocks in zip(df["seq"], df["blocks"]):
        b1, b2 = base_at(seq, blocks, s1["pos"]), base_at(seq, blocks, s2["pos"])
        if b1 is None or b2 is None:
            continue
        a1 = "alt" if b1 == s1["alt"] else "ref" if b1 == s1["ref"] else None
        a2 = "alt" if b2 == s2["alt"] else "ref" if b2 == s2["ref"] else None
        combos[f"{a1}_{a2}" if a1 and a2 else "other"] += 1
    check("F09", lib, "reads_covering_both", sum(combos.values()))
    for a1, a2 in COMBOS:
        check("F09", lib, f"{a1}_{a2}", combos[f"{a1}_{a2}"])
    check("F09", lib, "other", combos["other"])
    # left panel: up to 14 reads that touch either site, drawn as lines
    near = df[(df["pos"] <= hi_x) & (df["end"] >= lo_x)]
    touch = [read_index_at(b, s1["pos"]) is not None or read_index_at(b, s2["pos"]) is not None
             for b in near["blocks"]]
    near = near[touch].sort_values("pos")
    pick = near.iloc[np.unique(np.linspace(0, len(near) - 1, min(14, len(near))).astype(int))]
    st = look(lib)
    start = row
    for p, e in zip(pick["pos"], pick["end"]):
        ax_r.plot([max(p, lo_x), min(e, hi_x)], [row, row], color=st["color"], ls=st["ls"], lw=2.2,
                  solid_capstyle="butt")
        row += 1
    centres.append((start + row - 1) / 2)
    row += 2
    # right panel: how many reads show each combination of bases
    base = np.arange(len(COMBOS))
    j = DNA_LIBS.index(lib)
    vals = [combos[f"{a1}_{a2}"] for a1, a2 in COMBOS]
    ax_c.bar(base + (j - 1) * 0.27, vals, width=0.25, color=st["face"], edgecolor=st["color"], lw=1.2,
             label=lib)
    for xb, v in zip(base + (j - 1) * 0.27, vals):
        ax_c.text(xb, v + 0.3, str(v), ha="center", va="bottom", fontsize=11, color=st["color"])
for snv in (s2, s1):
    ax_r.axvline(snv["pos"], color="0.55", lw=1.2, ls=":", zorder=0)
    ax_r.text(snv["pos"], -1.2, snv.name, ha="center", va="bottom", fontsize=12, color="0.3")
ax_r.set_xlim(lo_x, hi_x)
ax_r.set_ylim(-3.2, row - 1)
ax_r.set_yticks(centres)
ax_r.set_yticklabels(DNA_LIBS)
ax_r.invert_yaxis()
ax_r.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1000:g}"))
ax_r.set_xlabel("position on chrS (kb)")
ax_r.set_title("reads touching SNV2 or SNV1", loc="left", fontsize=12)
ax_c.set_xticks(base)
ax_c.set_xticklabels([f"{s1[a1]} + {s2[a2]}\n{a1}/{a2}" for a1, a2 in COMBOS])
ax_c.set_xlabel("base at SNV1 + base at SNV2")
ax_c.set_ylabel("reads covering both sites")
ax_c.set_ylim(0, max(1, ax_c.get_ylim()[1]) * 1.15)
ax_c.legend(loc="upper right", fontsize=12, handlelength=1.2)
title(fig, "F09", "Phasing: are two variants on the same copy?")
save(fig, "F09", "phasing")
plt.show()
if COUNT_PAIRS:
    df = sam["DNA-ILMN"][sam["DNA-ILMN"]["mapped"]]
    seen = {}
    for q, blocks in zip(df["qname"], df["blocks"]):
        c = seen.setdefault(q, [False, False])
        c[0] |= read_index_at(blocks, s1["pos"]) is not None
        c[1] |= read_index_at(blocks, s2["pos"]) is not None
    check("F09", "DNA-ILMN", "pairs_one_mate_each", sum(a and b for a, b in seen.values()))
show_checks("F09")''',
    "Set `COUNT_PAIRS = True`. Could an Illumina read pair, with an insert of about 350 bp ("
    + cite(16) + "), ever link two sites 3,000 bases apart?",
)

section(
    "F10",
    "Why does one dataset lose coverage in GC-rich DNA?",
    "The GC content of a window is the fraction of its bases that are G or C; it varies along every "
    "chromosome. Some library preparations copy DNA by PCR before sequencing, and PCR does not copy "
    "every sequence equally well. The `library_prep` column of the sample sheet says which "
    "libraries were amplified.",
    r'''WINDOW = 1_000      # Try this: 500 or 2_000 (must divide 120,000)

gc = np.array([(CHRS[i:i + WINDOW].count("G") + CHRS[i:i + WINDOW].count("C")) / WINDOW
               for i in range(0, LOCUS_LEN, WINDOW)])
high, mid = gc >= 0.60, (gc >= 0.40) & (gc < 0.50)
cgi = regions.loc["CGI_ALTX"]
cgi_window = int(((cgi["start"] + cgi["end"]) // 2 - 1) // WINDOW)     # window holding the island
fig, axes = plt.subplots(1, 3, figsize=(8, 4.5), sharey=True, constrained_layout=True)
norms = {}
for ax, lib in zip(axes, DNA_LIBS):
    w = DEPTH[lib][1:].reshape(-1, WINDOW).mean(axis=1)
    norm = norms[lib] = w / np.median(w)
    check("F10", lib, "n_windows_gc60", int(high.sum()))
    check("F10", lib, "n_windows_gc40_50", int(mid.sum()))
    check("F10", lib, "norm_depth_gc60", round(float(norm[high].mean()), 3))
    check("F10", lib, "norm_depth_gc40_50", round(float(norm[mid].mean()), 3))
    check("F10", lib, "ratio", round(float(norm[high].mean() / norm[mid].mean()), 3))
    check("F10", lib, "sd_gc40_50", round(float(norm[mid].std()), 3))   # window-to-window scatter
    st = look(lib)
    ax.axvspan(40, 50, color="0.92", zorder=0)
    ax.axvspan(60, 100 * gc.max() + 1, color="0.92", zorder=0)
    ax.scatter(100 * gc, norm, s=24, color=st["color"], alpha=0.75, edgecolors="none")
    ax.scatter(100 * gc[cgi_window], norm[cgi_window], s=130, facecolors="none",
               edgecolors=st["color"], lw=1.6)
    ax.set_title(lib, loc="left", fontsize=13, color=st["color"])
axes[0].annotate("ALTX promoter\n(CpG island)", (100 * gc[cgi_window], norms["DNA-ILMN"][cgi_window]),
                 xytext=(52, 0.2), fontsize=12, color="0.3", arrowprops={"arrowstyle": "-", "color": "0.5"})
for ax in axes:
    ax.text(45, 1.68, "40-50%", ha="center", fontsize=11, color="0.4")
    ax.text(62 + (100 * gc.max() - 60) / 2, 1.68, "≥60%", ha="center", fontsize=11, color="0.4")
axes[0].set_ylim(0, 1.8)
axes[0].set_ylabel(f"depth ÷ median depth ({kb(WINDOW)} windows)")
fig.supxlabel("GC content of the window (%)", fontsize=13)
title(fig, "F10", "GC content and coverage")
save(fig, "F10", "gc_coverage")
plt.show()
show_checks("F10")''',
    "Set `WINDOW = 500`, then `2_000`. Which numbers move? Compare the drop at ≥60% GC with "
    "`sd_gc40_50`, the window-to-window scatter: is the drop bigger than the noise?",
)

section(
    "F11",
    "Which promoter is methylated, and does that fit the RNA data?",
    "A cytosine followed by a guanine (a CpG site) can carry a methyl group: 5-methylcytosine, or "
    "5mC. A methylated promoter usually goes with a silent gene; an unmethylated CpG-island promoter "
    "usually goes with an active one. PacBio and Nanopore read 5mC from the native DNA molecule; "
    "there is no methylation table for DNA-ILMN (check `library_prep` in the sample sheet).",
    r'''MIN_READS = 5      # Try this: 15 (use only CpGs read by at least 15 molecules)

meth = {lib: pd.read_csv(DATA / f"methylation_{lib}.tsv", sep="\t") for lib in ("DNA-HIFI", "DNA-ONT")}
cgi, off1p = regions.loc["CGI_ALTX"], regions.loc["OFF1_PROMOTER"]
REGION_ORDER = ["CGI_ALTX", "OFF1_PROMOTER", "ALTX_body", "other"]


def region_of(pos):
    name = np.full(len(pos), "other", dtype=object)
    s, e = GENE_SPAN["ALTX"]
    name[(pos >= s) & (pos <= e)] = "ALTX_body"
    name[(pos >= off1p["start"]) & (pos <= off1p["end"])] = "OFF1_PROMOTER"
    name[(pos >= cgi["start"]) & (pos <= cgi["end"])] = "CGI_ALTX"
    return name


rna = sam["RNA-ILMN"]
read1 = rna.loc[rna["mapped"] & ((rna["flag"] & 64) > 0), "blocks"]     # one read per fragment
fragments = {g: sum(any(s <= ee and e >= es for s, e, _ in b for es, ee in GENE_EXONS[g]) for b in read1)
             for g in ("ALTX", "OFF1")}
WINDOWS = {"ALTX": (17_000, 32_000), "OFF1": (63_500, 77_500)}
fig, axes = plt.subplots(2, 2, figsize=(8, 4.5), sharex="col", sharey="row", constrained_layout=True,
                         gridspec_kw={"height_ratios": [4, 1.3]})
for lib, tab in meth.items():
    tab = tab[tab["n_reads"] >= MIN_READS]
    region = region_of(tab["pos"].to_numpy())
    for name in REGION_ORDER:
        sel = region == name
        check("F11", lib, f"n_cpg_{name}", int(sel.sum()))
        check("F11", lib, f"mean_{name}", round(float(tab["fraction"].to_numpy()[sel].mean()), 3))
    st = look(lib)
    for col, (gene, (x0, x1)) in enumerate(WINDOWS.items()):
        part = tab[(tab["pos"] >= x0) & (tab["pos"] <= x1)]
        axes[0, col].scatter(part["pos"], part["fraction"], s=9, color=st["color"], alpha=0.3,
                             edgecolors="none")
        smooth = part["fraction"].rolling(15, center=True, min_periods=5).mean()
        axes[0, col].plot(part["pos"], smooth, color=st["color"], ls=st["ls"], lw=2, label=lib)
for col, (gene, (x0, x1)) in enumerate(WINDOWS.items()):
    r = cgi if gene == "ALTX" else off1p
    axes[0, col].axvspan(r["start"], r["end"], color="0.88", zorder=0)
    axes[0, col].set_title(f"{gene}: {fragments[gene]:,} RNA-ILMN fragments", loc="left", fontsize=12)
    draw_genes(axes[1, col], (x0, x1))
    axes[1, col].set_xlim(x0, x1)
    axes[1, col].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1000:g}"))
axes[0, 0].text(cgi["end"] + 300, 0.03, "promoter\n(CpG island)", fontsize=11, color="0.3")
axes[0, 1].text(off1p["start"] - 300, 0.03, "promoter", ha="right", fontsize=11, color="0.3")
axes[0, 0].set_ylim(-0.02, 1.05)
axes[0, 0].set_ylabel("methylated\nfraction per CpG")
axes[0, 1].legend(loc="lower right", fontsize=12)
fig.supxlabel("position on chrS (kb)", fontsize=13)
title(fig, "F11", "Methylation from long-read DNA")
save(fig, "F11", "methylation")
plt.show()
show_checks("F11")''',
    "Set `MIN_READS = 15`. Fewer CpGs pass: do the region means change?",
)

section(
    "F12",
    "Do the three RNA datasets agree on which genes are on, and how sure is each?",
    "The share of a gene's reads among all reads estimates how much of the RNA comes from that gene. "
    "An Illumina fragment is a piece of a molecule, so longer transcripts give more fragments and "
    "the counts are divided by transcript length; a long read is one whole molecule. With more "
    "reads, the 95% interval gets narrower.",
    r'''SUBSAMPLE = 1.0       # Try this: 0.25 (keep a random quarter of the long reads)

fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
for j, lib in enumerate(RNA_LIBS):
    df = sam[lib][sam[lib]["mapped"]]
    if lib == "RNA-ILMN":
        df = df[(df["flag"] & 64) > 0]                       # read 1 only: one count per fragment
    elif SUBSAMPLE < 1:
        df = df.sample(frac=SUBSAMPLE, random_state=1)
    counts = {g: 0 for g in GENE_ORDER}
    for blocks in df["blocks"]:
        for g in GENE_ORDER:
            if any(s <= ee and e >= es for s, e, _ in blocks for es, ee in GENE_EXONS[g]):
                counts[g] += 1
                break
    per = {g: counts[g] / (TX_LENGTH[g] if lib == "RNA-ILMN" else 1) for g in GENE_ORDER}
    total = sum(per.values())
    st = look(lib)
    for gi, g in enumerate(GENE_ORDER):
        share, c = per[g] / total, counts[g]
        if c:
            lo, hi = share * max(0.0, 1 - 1.96 / math.sqrt(c)), share * (1 + 1.96 / math.sqrt(c))
        else:                                                # zero reads: "rule of three" upper bound
            lo, hi = 0.0, (3 / (TX_LENGTH[g] if lib == "RNA-ILMN" else 1)) / total
        check("F12", lib, f"{g}_reads", c)
        check("F12", lib, f"{g}_share", round(share, 4))
        check("F12", lib, f"{g}_share_lo", round(lo, 4))
        check("F12", lib, f"{g}_share_hi", round(hi, 4))
        bars = ax.errorbar(gi + (j - 1) * 0.22, 100 * share, yerr=[[100 * (share - lo)], [100 * (hi - share)]],
                           fmt="o", ms=9, color=st["color"], mfc=st["face"], mew=1.8, elinewidth=1.8,
                           capsize=4, label=lib if gi == 0 else None)
        bars[2][0].set_linestyle(st["ls"])
ax.set_xticks(range(len(GENE_ORDER)))
ax.set_xticklabels(GENE_ORDER, fontstyle="italic")
ax.set_xlim(-0.5, len(GENE_ORDER) - 0.5)
ax.set_ylim(-2, None)
ax.set_ylabel("share of this locus's RNA (%)")
ax.legend(loc="upper right")
ax.text(len(GENE_ORDER) - 1, 9, "Illumina counts are divided\nby transcript length", ha="center",
        fontsize=11, color="0.35")
title(fig, "F12", "Gene expression: do the RNA datasets agree?")
save(fig, "F12", "expression")
plt.show()
show_checks("F12")''',
    "Set `SUBSAMPLE = 0.25`. How much wider do the long-read intervals get? Short reads give the "
    "most reads for the money; what can only the long reads tell you (see F05)?",
)

# =====================================================================================
# Scorecard, answer check, slide
# =====================================================================================
md(f'''## Scorecard: what each dataset can and cannot see

The scorecard turns the `CHECK` numbers into one answer per question per dataset. The rules are
in the cell, so you can argue with them. It is the starting point for the debrief
({cite(35)}).''')

code(r'''def get(fig, lib, metric):
    return CHECKS[(fig, lib, metric)]


rows = []
for lib in LIBS:
    row = {"library": lib, "molecule": lib[:3], "platform": PLATFORM[lib],
           "read_n50": get("F01", lib, "read_n50"),
           "accuracy_pct": round(float(get("F02", lib, "accuracy_pct")), 2)}
    if lib in DNA_LIBS:
        row["sizes_expanded_repeat"] = "yes" if get("F07", lib, "long_n") >= 3 else "no"
        row["resolves_duplication"] = "yes" if get("F08", lib, "ratio_inside_flank") >= 0.5 else "no"
        row["phases_snv1_snv2"] = "yes" if get("F09", lib, "alt_ref") + get("F09", lib, "ref_alt") >= 3 else "no"
        row["methylation_from_reads"] = "yes" if (DATA / f"methylation_{lib}.tsv").exists() else "no"
        for k in ("allele_specific_expression", "full_length_isoforms", "gene_expression"):
            row[k] = "n/a"
    else:
        lo, hi = get("F04", lib, "ci_low"), get("F04", lib, "ci_high")
        row["allele_specific_expression"] = "yes" if lo > 0.5 or hi < 0.5 else "no"
        full = sum(get("F05", lib, f"chain_{k}") for k in ("E2a_E4a", "E2b_E4b", "E2a_E4b", "E2b_E4a"))
        row["full_length_isoforms"] = "yes" if full >= 10 else "no"
        share = [get("F12", lib, f"{g}_share") for g in ("ASE1", "ALTX", "LOW1")]
        row["gene_expression"] = ("yes" if share[0] > share[1] > share[2] and get("F12", lib, "OFF1_reads") == 0
                                  else "no")
        for k in ("sizes_expanded_repeat", "resolves_duplication", "phases_snv1_snv2", "methylation_from_reads"):
            row[k] = "n/a"
    rows.append(row)
scorecard = pd.DataFrame(rows)[["library", "molecule", "platform", "read_n50", "accuracy_pct",
                                "sizes_expanded_repeat", "resolves_duplication", "phases_snv1_snv2",
                                "methylation_from_reads", "allele_specific_expression",
                                "full_length_isoforms", "gene_expression"]]
scorecard.to_csv(OUT / "practice3_scorecard.tsv", sep="\t", index=False)
print("rules: repeat sized = at least 3 reads span the long allele (F07); duplication resolved = "
      "MAPQ>=20 depth inside DUP_A at least half the flank depth (F08); phased = at least 3 reads link "
      "the two SNVs (F09); allele-specific expression = the SNV1 interval excludes 0.5 (F04); "
      "full-length isoforms = at least 10 reads with a complete exon chain (F05); expression = "
      "ASE1 > ALTX > LOW1 and no OFF1 reads (F12)")
scorecard.set_index("library").T''')

md('''## Compare with the answer file

`checkpoints/practice3_checks.tsv` holds every `CHECK` number from an independent calculation.
If you ran the notebook unchanged, all of them should match. If you changed a **Try this**
parameter, the numbers of that figure will differ, which is the point.''')

code(r'''mine = pd.read_csv(OUT / "practice3_checks.tsv", sep="\t", dtype=str, keep_default_na=False)
answer = pd.read_csv(CHECK / "practice3_checks.tsv", sep="\t", dtype=str, keep_default_na=False)
both = answer.merge(mine, on=["figure", "library", "metric"], how="left", suffixes=("_answer", "_yours"))
same = both["value_answer"] == both["value_yours"]
print(f"{int(same.sum())} of {len(both)} CHECK numbers match checkpoints/practice3_checks.tsv")
if not same.all():
    print("These differ (did you change a Try this parameter?):")
    print(both.loc[~same].head(12).to_string(index=False))''')

md(f'''## Make your slide

Your instructor provides `presentation/practice3_lightning_talk_template.pptx`. It has two slides:
the **example slide**, the instructor's worked example on F01 (below), and the **pair slide**,
which your pair fills in for your figure. Each talk has two minutes ({cite(34)}).

| On the pair slide | What to write |
|---|---|
| Title | your figure ID and title, for example *F03 · DNA vs RNA: coverage across the locus* |
| CLAIM | one sentence a classmate can check against the figure |
| The figure | `outputs/practice3_F0X_<name>.png`; the file name is printed under your figure |
| 1 · Axes and colours | what the axes, rows and colours show |
| 2 · Evidence: one number | point at the figure and quote one number from `CHECK` |
| 3 · So what: choice and one limit | which dataset you would choose for this question, and one thing it cannot do |
| Q / A bar | a question you expect from the floor, and your one-sentence answer |

### The example slide: F01 · {TITLES['F01']}

* **CLAIM.** {F01_EXAMPLE['claim']}
* **1 · Axes and colours.** {F01_EXAMPLE['axes']}
* **2 · Evidence: one number.** {F01_EXAMPLE['evidence']}
* **3 · So what: choice and one limit.** {F01_EXAMPLE['so_what']}
* **Q / A.** *{F01_EXAMPLE['q']}* {F01_EXAMPLE['a']}

After the talks, the debrief ({cite(35)}) puts the 12 claims side by side. Keep the scorecard open.''')
