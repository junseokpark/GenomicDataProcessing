# Practice package — Genomic Data Processing, Day 5 (10/27)

The practice for Part 6 of the Day 5 lecture *Data visualization in genetics and genomics*: about 90 minutes, in pairs.
Students first make simple figures from small tables they can read by eye. Then they work with realistic simulated
data: they run the test or method a real analysis would use, check the result against the planted answer, and only
then draw the figure. Every Practice B exercise has the same four beats, **Question → Test → Check → Visualize**, and a
`# CHECK ·` cell prints ✓ or ✗ for every number that matters.

| Notebook (in `dataviz_practice/`) | Lecture slides | What students do |
|---|---|---|
| `practice.ipynb` (everyone) | 83–89, 97–101 | Setup · Exercise 1, spot the misleader (M1–M8) · Practice A, simple data and simple figures (A1–A5) · B1, two groups · B2, bulk RNA-seq from counts to a volcano plot · D1–D3, defend one figure |
| `electives/B3_gwas.ipynb` | 91 | GWAS: real signal or ancestry? |
| `electives/B4_single_cell.ipynb` | 92 | Single cell: neighbours, not distances |
| `electives/B5_cnv_tracks.ipynb` | 93 | Copy number: find the 22q11.2 deletion |
| `electives/B6_variant_callers.ipynb` | 94 | Variant callers: agreement is not truth |
| `electives/B7_classifier.ipynb` | 95 | Classifiers: rare positives need PR |
| `electives/B8_tp53_lollipop.ipynb` | 96 | TP53 lollipop: where do mutations fall? |

Each pair works through `practice.ipynb` and **one** elective. The session ends with each pair presenting one figure,
reviewed with the 16-item `checklist.md`.

| Part | Minutes |
|---|---|
| Setup and file check | 5 |
| Exercise 1 · diagnose eight misleading figures | 8 |
| Practice A · simple data (A1–A5) | 20 |
| Practice B · B1 two groups, B2 bulk RNA-seq | 27 |
| One elective per pair (B3–B8) | 20 |
| Defend one figure: verify, accessibility, optional AI remake | 10 |

AI is optional. The prompt cards (`prompts/prompt_cards.md`) are one way to write the plotting code; the provided cells
are the other. Whatever writes the code, the CHECK cells decide whether the figure is right.

## How to run

With **conda** (Miniforge is recommended) on macOS, Linux or WSL2:

```bash
cd day5-visualization
conda env create -f environment.yml      # once
conda activate dataviz-practice
cd dataviz_practice && jupyter lab       # open practice.ipynb, then your elective in electives/
```

With **pip**: `python -m pip install -r dataviz_practice/requirements.txt`, then `jupyter lab` in `dataviz_practice/`.
For **Google Colab**, follow the quick start in [`dataviz_practice/README.md`](dataviz_practice/README.md).

Run each notebook from top to bottom. The core notebook runs in about 10 seconds, each elective in under a minute, and
no internet connection is needed. The elective notebooks work whether Jupyter starts in `electives/` or in
`dataviz_practice/`.

## Where the answers are

| File | What it holds |
|---|---|
| `solutions/executed/` | every notebook with all of its outputs |
| `solutions/ANSWER_KEY_bad_figures.md` | Exercise 1: what is wrong with M1–M8 and how to fix it |
| `solutions/figures.py` | reference code for the prompt-card figures (`python solutions/figures.py 1` … `8`) |
| `prompts/prompt_cards.md` | the eight prompt cards from the slides, with the values to check against |
| `INSTRUCTOR_NOTES.md` | timing, expected CHECK values, common mistakes, *Your turn* answers, simplifications in the data |
| `docs/DATA_DICTIONARY.md` | every data file and column, with units and coordinate conventions |

## The data

All data are simulated: no real people, patients or samples. Each dataset sits on real anchors so that names and
coordinates look familiar: hg38 chromosome lengths and MANE Select genes, glucocorticoid-response genes after Himes et
al. (2014), AMD loci after Fritsche et al. (2016), the 22q11.2 region with real UCSC segmental duplications and the
ClinGen region ISCA-37446, the TP53 protein (UniProt P04637), and Anscombe's published quartet. Every count, genotype,
cell, score and variant call is simulated with a fixed seed and a planted answer. **No output of this practice is a
clinical result.**

## Checked

From a fresh copy, all seven notebooks run without errors and all 180 pass/fail CHECKs print ✓ (practice 39, B3 22,
B4 36, B5 21, B6 22, B7 25, B8 15). They were tested with two sets of library versions; see "Tested with" in
`dataviz_practice/README.md`.

## Rebuilding the data, notebooks and lecture figures

The generators, notebook builders and lecture-figure scripts are in
[`dataviz_practice/build/`](dataviz_practice/build); see its README. `bash build/make_data.sh` (run in
`dataviz_practice/`) rebuilds `data/` byte for byte.
