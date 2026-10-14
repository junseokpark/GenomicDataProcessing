# Practice 2 · CNVs from reads to interpretation

**Call on hg19 → lift to hg38 → test expression in matched RNA → score with TESSERACT-CNV → read the tracks**

*Genomic Data Processing · Day 4 (10/15) · Understanding genetic variants · about 75 minutes in class*

> **Everything about the people in this practice is synthetic.** Three fictional patients (P1–P3) and five fictional
> controls (C1–C5) were simulated by computer. No real patient data are used. The DNA reads were simulated from **real
> hg19 sequence**, and the copy-number variants (CNVs) were planted at **real loci** taken from public ClinVar and LOVD
> records, so you can check every result against public annotation. The RNA counts are synthetic too ("fibroblast-like"
> baseline, not GTEx).

## What you will do

| Step | Question | Main tool |
|---|---|---|
| 0 | What are we looking for, and why two genome builds? | – |
| 1 | Where in hg19 did each read come from? | `minimap2`, `samtools` |
| 2 | How many DNA fragments per 5-kb bin, compared with controls? | `pysam`, panel of normals |
| 3 | Where does the copy number change? | a small hidden Markov model (HMM) |
| 4 | What are the same CNVs in hg38 coordinates – and did the conversion work? | `CrossMap` + UCSC chain file |
| 5 | Do genes inside the CNVs change their RNA level? | `pandas` |
| 6 | What would TESSERACT-CNV need as input? | input table, 60-kb segment grid, adapter cell |
| 7 | Which public annotation tracks overlap each CNV? | UCSC tracks (pre-downloaded) |
| 8 | All evidence in one table per CNV | – |
| 9 | Wrap-up and exit questions | – |

## Why this design?

An earlier version of the slides called CNVs on a toy "hg19-like" contig and then switched to unrelated real CNVs for
TESSERACT-CNV, because TESSERACT reads **real hg38 annotation tracks**. Here the synthetic patients are simulated from real
hg19 sequence with CNVs at real loci, so **the same CNV flows through every step**: DNA (hg19) → liftover → RNA (hg38) →
TESSERACT input (hg38) → track discussion.

| Sample | Planted CNV (hg19) | Public source |
|---|---|---|
| P1 | 22q11.2 proximal deletion (LCR22A–D), CN 1, chr22:18,916,842–21,465,662 | ClinVar VCV000155158 |
| P2 | PMP22 whole-gene duplication, CN 3, chr17:15,133,094–15,168,674 | LOVD PMP22_000101 |
| P2 | UGT2B17 deletion, CN 1, chr4:69,375,169–69,483,248 (common, benign) | ClinVar VCV000147015 |
| P3 | distal 22q11.2 deletion, CN 1, chr22:21,029,655–22,481,498 (the lecture's "case 3") | ClinVar VCV000816205 (Pathogenic) and VCV000443513 (Uncertain significance) – same coordinates |
| P3 | UGT2B17 deletion on both copies, CN 0 | as above |
| C2 | UGT2B17 deletion, CN 1 (healthy people carry CNVs too) | as above |

## Setup

You need **conda** (Miniforge recommended) or **mamba**.

**macOS or Linux**

```bash
cd practice2_cnv_tesseract
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate d4-practice
jupyter lab                              # then open practice2_cnv_tesseract_student.ipynb
```

> **One environment for both Day 4 practices.** If you already created `d4-practice` for Practice 1, skip `env create` and just run `conda activate d4-practice` (the two `environment.yml` files are identical).

**Windows:** use **WSL2** (Windows Subsystem for Linux) with Ubuntu: open the Ubuntu terminal, install Miniforge for
Linux inside WSL, and follow the macOS/Linux commands. `samtools`, `minimap2` and `CrossMap` are not available for
native Windows.

**Run without opening Jupyter (optional):**

```bash
jupyter nbconvert --to notebook --execute practice2_cnv_tesseract_student.ipynb --output my_run.ipynb
```

On a laptop the whole notebook needs about 1–3 minutes of computing. No internet connection is needed: all public data are
included. (Only the optional live track refresh in Step 7 downloads data.)

## Files

| Path | Content |
|---|---|
| `practice2_cnv_tesseract_student.ipynb` | the notebook you work in (no outputs) |
| `practice2_cnv_tesseract_executed.ipynb` | the same notebook after a complete run (to compare with) |
| `data/ref/hg19_cnv_windows.fa.gz` (+ `.fai`, `.gzi`) | practice reference: real hg19 sequence for chr22:18,000,001–23,000,000, chr17:14,600,001–15,700,000 and chr4:69,000,001–69,800,000 inside full-length chr22/chr17/chr4 that are otherwise `N` |
| `data/ref/windows_hg19.bed` | the three real-sequence windows |
| `data/ref/hg19ToHg38.over.chain.gz` | UCSC liftover chain file |
| `data/fastq/{P1,P2,P3,C1}_R{1,2}.fastq.gz` | synthetic paired-end reads, 2 × 100 bp, ≈ 3× depth |
| `data/pon/pon_bincounts_hg19_5kb.tsv.gz` | panel of normals: fragment counts per 5-kb bin for C1–C5 (+ GC and N fraction) |
| `data/truth/planted_cnvs_hg19.{bed,tsv}` | the planted CNVs (look only in Step 3) |
| `data/truth/rna_truth_cnv_genes.tsv` | how the synthetic RNA effect was simulated (instructors) |
| `data/annotation/genes_hg19_refseqselect_windows.tsv` | hg19 protein-coding genes in the windows (UCSC `ncbiRefSeqSelect`, NM_ only) |
| `data/annotation/genes_hg38_mane_select.tsv.gz` | hg38 MANE Select protein-coding genes, genome-wide (MANE v1.5) |
| `data/rna/counts_hg38_mane.tsv.gz`, `data/rna/samples.tsv` | synthetic RNA-seq gene counts for all 8 samples, sample sheet |
| `data/tracks_hg38/*.tsv.gz` | UCSC hg38 tracks around each CNV (+ `provenance_tracks_hg38.json`) |
| `data/tracks_hg19/*.tsv.gz` | hg19 segmental duplications and ClinVar CNVs in the windows |
| `data/provenance_public_data.json` | URL, retrieval time and data version of every public download |
| `scripts/` | data generator (`make_all.sh`), UCSC downloader with live track refresh (`fetch_public_data.py`), plotting helpers (`p2_figures.py`) |
| `results/p2_results.json`, `results/figures/` | results and slide figures of the reference run (your run overwrites them with the same values) |
| `fonts/` | SUIT font used by the lecture slides (SIL Open Font License 1.1); figures fall back to DejaVu Sans if absent |

The notebook writes alignments to `data/bam/` and intermediate tables (CNV calls, liftover files, TESSERACT input) to
`data/results/`.

## TESSERACT-CNV

TESSERACT-CNV (Buralkin I, Botas J, Chang K-L, Deng Y, Papastathopoulos-Katsaros A, Liu Z, Park J. *Vision-Based Genomic
Model for Copy Number Variant Pathogenicity Prediction.* bioRxiv 2026, doi:10.64898/2026.05.21.726953) is **not public**.
The notebook prepares its input and calls the lab's installation through **one adapter cell** in Step 6
(`TESSERACT_CMD`). If the command is not set, the notebook prints "TESSERACT-CNV not run" and continues. No score is ever
simulated or filled in.

## Rebuilding the data (instructors)

`bash scripts/make_all.sh` downloads the public inputs again (UCSC API and hgdownload) and rebuilds every synthetic file
with the fixed seed **20261015**. See `INSTRUCTOR_NOTES.md`.

## Sources and terms

Sequence, chain file and annotation tracks: UCSC Genome Browser (https://genome.ucsc.edu; data are subject to the terms of
the original providers: NCBI RefSeq/MANE, ClinGen, ClinVar, gnomAD, ENCODE, DGV, UCSC). ClinVar/LOVD record identifiers
are given for teaching; check the current records online before using them for anything else. This material is for
teaching only and is not for clinical use.
