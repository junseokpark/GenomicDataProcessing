# Practice packages — Genomic Data Processing, Day 4 (10/15)

Two 75-minute practices for the Day 4 lecture *Understanding genetic variants*. Both start from the
raw sequencing reads of synthetic people at real genomic loci and end with evidence you can defend.
Every step ends with a **Check** cell: count what changed before you move on.

Open each notebook from **its own folder**; every path inside a notebook is relative to that folder.

| Folder | Notebook | Lecture slides | The question students answer |
|---|---|---|---|
| [`practice1_variant_calling/`](practice1_variant_calling) | `practice1_variant_calling_student.ipynb` | 40–52 · Practice 1 · From reads to variants | A child has sickle cell disease and both parents are healthy. Which variant explains it, and can we trust the call? |
| [`practice2_cnv_tesseract/`](practice2_cnv_tesseract) | `practice2_cnv_tesseract_student.ipynb` | 53–69 · Practice 2 · CNV case study | Does each of three referred patients carry a copy-number variant that could explain the features, and how strong is the evidence? |

The step numbers on the slides match the step numbers in the notebooks.

## How to run

You need **conda** (Miniforge is recommended) on macOS or Linux. On Windows, use WSL2 with Ubuntu. Both
practices use the same environment, so create it once:

```bash
cd day4-VariantCalling
conda env create -f environment.yml      # once (or: mamba env create -f environment.yml)
conda activate d4-practice

cd practice1_variant_calling && jupyter lab        # Practice 1
cd ../practice2_cnv_tesseract && jupyter lab       # Practice 2
```

Run the student notebook from top to bottom (Shift + Enter). Each notebook finishes within a few minutes
on a laptop. If you get stuck, the `*_executed.ipynb` copy shows every expected output. The practice
READMEs cover Apple-silicon Macs, WSL2 and common error messages.

No internet connection is needed. The exceptions are optional: the Ensembl VEP query in Practice 1,
Step 11, and the live track refresh in Practice 2, Step 7.

## Practice 1 · From reads to variants (GRCh38 trio)

Reads were simulated from the real GRCh38 β-globin window (chr11:5,200,001–5,300,000) for a child (S01),
father (S02) and mother (S03): 2 × 150 bp at about 30×. The reference is full-length chr11 with real
sequence only in the window, so every position is a real GRCh38 position.

| Step | What students do | Main tool |
|---|---|---|
| 1 | Look at the reference, the reads and their qualities | `samtools faidx`, `fastp` |
| 2 | Align, sort and index; label each sample with a read group | `minimap2`, `samtools` |
| 3 | Mark PCR duplicates | `samtools markdup` |
| 4 | Check the coverage along the window | `samtools depth` |
| 5 | Read the pileup at HBB and work out genotype likelihoods by hand | `samtools mpileup`, Python |
| 6 | Call variants jointly in the trio | `bcftools mpileup` / `call` |
| 7 | Normalise how each variant is written | `bcftools norm` |
| 8–9 | Filter, then benchmark against the planted truth | `bcftools filter`, Python |
| 10–11 | Check inheritance; name the variant in HGVS | `bcftools +mendelian2`, Ensembl VEP (optional) |

The answer is NM_000518.5(HBB):c.20A>T p.(Glu7Val). Three planted traps each cost one call when the
matching safeguard is skipped: a pile of PCR duplicates, systematic low-quality bases, and a
low-coverage patch.

## Practice 2 · CNVs from reads to interpretation (hg19 → hg38)

Reads were simulated from real hg19 sequence (2 × 100 bp, about 3×). CNVs were planted at loci taken from
public ClinVar and LOVD records:

| Sample | Planted CNV |
|---|---|
| P1 | proximal 22q11.2 deletion |
| P2 | PMP22 duplication, plus a common UGT2B17 deletion |
| P3 | distal 22q11.2 deletion, plus loss of both copies of UGT2B17 |
| C1–C5 | five controls, used as the panel of normals |

| Step | What students do | Main tool |
|---|---|---|
| 0 | The clinical question, and why two genome builds | — |
| 1 | Map the DNA reads to hg19 | `minimap2`, `samtools` |
| 2 | Turn read depth into a copy-number signal against a panel of normals | `pysam`, `pandas` |
| 3 | Segment and call CNVs with a small hidden Markov model | Python |
| 4 | Lift the calls to hg38, and check the conversion | `CrossMap` + UCSC chain file |
| 5 | Test expression in matched RNA (synthetic) on hg38 | `pandas` |
| 6 | Prepare the TESSERACT-CNV input; run it through the adapter cell | input table, 60-kb segment grid |
| 7 | Read the public annotation tracks around each CNV | UCSC tracks (pre-downloaded) |
| 8 | Put all evidence for each CNV into one table | — |

**TESSERACT-CNV** (Buralkin et al., bioRxiv 2026, doi:10.64898/2026.05.21.726953) is not publicly
available. Step 6 writes its input on hg38 and calls a lab installation through one adapter cell
(`TESSERACT_CMD`). If the command is not set, the notebook prints "TESSERACT-CNV not run" and continues.
No score is ever simulated.

## Answers and instructor notes

Each practice folder contains `INSTRUCTOR_NOTES.md`, with expected values, timings, common errors and
answers. Each also has an answer key under `data/truth/`:

- **Practice 1:** students open the key only in the benchmark step (Step 9).
- **Practice 2:** students open the key only after calling CNVs (Step 3).

## Real data vs teaching data

**Real:**
- the GRCh38 and hg19 reference sequence;
- the gene models (MANE Select, RefSeq Select) and the UCSC hg19 → hg38 chain file;
- the ClinVar and LOVD loci;
- the UCSC annotation tracks (ClinGen dosage and gene–disease, ClinVar, gnomAD, DGV, ENCODE, conservation).

Each real source is recorded with its URL and retrieval date in the package's provenance files.

**Synthetic:** every person, read, genotype and RNA count, generated with the fixed seed 20261015. The
phenotypes are fictional.

The evidence tables are teaching exercises, not classifications. **No output of these materials is a
clinical result.**

## Reproducing the data

```bash
cd practice1_variant_calling
python scripts/generate_data.py      # rebuilds data/ from the saved UCSC downloads; --download re-fetches them

cd ../practice2_cnv_tesseract
bash scripts/make_all.sh             # re-downloads the public inputs and rebuilds every synthetic file
```

Tested in October 2026 with:
- Python 3.12, pysam 0.24.1
- samtools, bcftools and htslib 1.24
- minimap2 2.31, fastp 1.4.0
- CrossMap (conda package 0.7.4)

`environment.yml` pins these versions.
