# Practice package — Genomic Data Processing, Day 2 (10/13)

Two mini practice sessions of about 15 minutes each, one notebook per session. Open the
notebooks from **this folder** (`day2-3-samples/`): every path inside them is relative to it.

| Notebook | Practice slide | The question students answer |
|---|---|---|
| `practice1_sample_to_platform.ipynb` | 25 · Mini practice 1 · Is this sample fit for this platform? | Which tubes can go on which platform, and what does each platform return from the same DNA? |
| `practice2_cohort_qc.ipynb` | 32 · Mini practice 2 · Cohort QC: who does not belong? | Which of 2,000 exomes must leave the analysis, at which step, and why? |

## How to run

```bash
cd day2-3-samples
python -m pip install -r requirements.txt      # pandas, numpy, matplotlib
jupyter lab                                    # or any notebook front-end
```

Open `practice1_sample_to_platform.ipynb`, run it top to bottom (Shift+Enter), then
`practice2_cohort_qc.ipynb`. Everything the notebooks write lands in `outputs/`.

## Requirements

| Needed | Why |
|---|---|
| Python 3.9+ with `pandas` ≥ 1.3, `numpy` ≥ 1.21, `matplotlib` ≥ 3.4 | both notebooks |
| JupyterLab or Jupyter Notebook *(optional)* | to open the notebooks |
| `nbformat`, `nbclient`, `ipykernel` *(maintainers only)* | `build/make_notebooks.py` and `build/verify_package.py` |

No network access, no reference download, no aligner, no `pysam`, no index build. Every cell
finishes in about a second, and each notebook runs top to bottom in a few seconds.

Tested headless with Python 3.12 / pandas 3.0.6 / numpy 2.5.4 / matplotlib 3.11.2, and with the
oldest supported stack, Python 3.9 / pandas 1.3.5 / numpy 1.21 / matplotlib 3.4. Rebuilding the
package on either stack gives byte-identical data, checkpoints and notebooks.

## What is real and what is synthetic

**Real**
* `data/practice1/hbb_window.hg38.fa` is UCSC hg38 chr11:5,218,771–5,278,770 (60,000 bp), the same
  β-globin window as Day 1. The build checks it by SHA-256 and against the Day 1 checkpoint bases
  (hg38 chr11:5,226,995–5,227,012 = `CTTCTCCTCAGGAGTCAG`).
* `data/practice1/genes.hg38.bed` holds the NCBI RefSeq Select β-globin models from Day 1.
* The input, length and integrity numbers in `platform_requirements.tsv` come from the vendor
  documents named in its `source` column (table below).
* The cohort QC rules come from Wang *et al.* 2021 (*Nature* 597, 527–532,
  doi:10.1038/s41586-021-03855-y) and the kinship categories from the KING guide.

**Synthetic**: the 12 tubes (S01–S12) and their subjects (P001–P012), every read and every
alignment, and the whole 2,000-sample cohort. No real person's data and no UK Biobank data are
used, and no result has clinical meaning. The 8-kb target `HBB_HBD_window` is a teaching target
("imagine a structural variant here"); nothing is claimed about the region.

**Teaching thresholds** are labelled as such in the files: the Illumina `min_fragment_kb`
values, the PacBio rule "the DNA must be at least as long as the 15–20 kb library", poly(A)
RNA-seq RIN ≥ 7 and 100 ng, the RNA purity window, the sex-check cut-offs (chrY ratio 0.20,
chrX heterozygosity 0.30) and the array-concordance cut-off (0.95).

### Platform numbers and where they come from

| assay | kit / chemistry | input | length / integrity rule | source |
|---|---|---|---|---|
| `ILMN_PCRFREE_WGS` | Illumina DNA PCR-Free Prep | 25 ng minimum (25–300 ng recommended; protocol accepts 25–2,000 ng) | modal fragment ≥ 1 kb *(teaching)* | Illumina product page; reference guide 1000000086922 |
| `ILMN_WES` | Illumina DNA Prep with Enrichment | 100 ng (100–500 ng recommended for human DNA) | ≥ 0.15 kb *(teaching; FFPE workflow)* | Illumina DNA Prep reference guide 1000000025416 |
| `PACBIO_HIFI_WGS` | PacBio Revio, SPRQ-Nx chemistry | 500 ng | ≥ 15 kb *(teaching: human WGS libraries are 15–20 kb)* | PacBio Revio product page / specification sheet |
| `ONT_LIGATION_LONG` | Ligation Sequencing Kit V14 (SQK-LSK114), R10.4.1 | 1 µg | ≥ 10 kb (long-fragment libraries) | ONT store page / chemistry technical document |
| `ONT_ULTRALONG` | Ultra-Long DNA Sequencing Kit V14 (SQK-ULK114) | 6 million cells | none: the protocol extracts the DNA itself (read N50 50–100+ kb, reads to 4+ Mb) | ONT store page |
| `RNA_POLYA` | poly(A) selection mRNA-seq | 100 ng *(teaching)* | RIN ≥ 7 *(teaching)* | common core-facility practice |
| `RNA_CAPTURE_FFPE` | capture RNA-seq for degraded RNA (e.g. Illumina RNA exome, TruSeq RNA Access) | by DV200: > 70% 20 ng; 50–70% 20–40 ng; 30–50% 40–100 ng | DV200 ≥ 30% | Illumina technical note 470-2014-001 |

Purity: DNA A260/280 1.8–2.0 and A260/230 2.0–2.2 (Illumina DNA Prep guidance); RNA, teaching
window A260/280 1.9–2.1 and A260/230 ≥ 1.8. Purity outside the window is always REVIEW, never
FAIL, because a clean-up can fix it.

## Layout

```
day2-3-samples/
├── README.md
├── requirements.txt
├── .gitignore                          ignores outputs/* except outputs/.gitkeep
├── practice1_sample_to_platform.ipynb  Mini practice 1
├── practice2_cohort_qc.ipynb           Mini practice 2
├── data/
│   ├── practice1/
│   │   ├── sample_sheet.tsv            12 tubes: specimen + molecule + metadata
│   │   ├── platform_requirements.tsv   7 assays: input, length/integrity and purity rules, sources
│   │   ├── hbb_window.hg38.fa (+ .fai) real hg38 chr11:5,218,771-5,278,770
│   │   ├── genes.hg38.bed              real RefSeq Select beta-globin models (from Day 1)
│   │   ├── target_window.bed           chr11 5226000 5234000 HBB_HBD_window (8 kb)
│   │   └── reads/
│   │       ├── illumina_R1.fastq.gz    2,000 read pairs, 2 x 150 bp (mate 1)
│   │       ├── illumina_R2.fastq.gz    (mate 2)
│   │       ├── pacbio_hifi.fastq.gz    40 HiFi reads
│   │       ├── ont_ligation.fastq.gz   50 nanopore reads
│   │       └── alignments.paf          4,090 alignments, written from the simulation truth
│   └── practice2/
│       ├── cohort_manifest.tsv         2,000 samples from 1,998 participants
│       ├── sample_qc_metrics.tsv       depth, CCDS 10x, FREEMIX, sex metrics, array concordance, PCs, Pr(EUR)
│       └── kinship_pairs.tsv           65 pairs with KING kinship > 0.0442
├── checkpoints/                        verified outputs, for when a cell will not run
│   ├── practice1_decisions.tsv
│   ├── practice1_read_stats.tsv
│   ├── practice2_qc_funnel.tsv
│   ├── practice2_excluded_samples.tsv
│   └── practice2_planted_cases.tsv     answer key: every planted case and the step it must reach
├── outputs/                            what the notebooks write (only .gitkeep is tracked)
└── build/
    ├── build_practice1.py              sample sheet, requirements, reads, PAF, Practice 1 checkpoints
    ├── build_practice2.py              cohort, metrics, kinship, Practice 2 checkpoints + answer key
    ├── make_notebooks.py               writes both notebooks with nbformat (no outputs stored)
    └── verify_package.py               runs both notebooks headless from a clean copy and compares
```

The checkpoints are produced by **independent re-implementations** inside the build scripts,
not by running the notebooks, so `verify_package.py` compares two separate codings of the same
rules.

## Mini practice 1 · Is this sample fit for this platform?

Lecture message: a genomics sample is **specimen + molecule + metadata**, and the sample sets
the ceiling for what any platform can read. *You can clean DNA; you cannot make it longer.*

**Part A · sample sheet → platform decision**

1. A1 Read the sheet: specimen, molecule and metadata side by side; `NA` means not measured.
2. A2 Read the requirements table, source column first.
3. A3 Total mass = concentration × volume (worked example: S08, 30 ng/µL × 20 µL = 600 ng).
4. A4 `decide()` applies six rules in order: consent → BLOCKED; molecule and cell count; mass
   (REVIEW at ≥ 50% of the minimum, else FAIL); length or integrity → FAIL (cannot be repaired);
   purity → REVIEW (clean-up); tumour without a matched normal for a somatic question → REVIEW.
   The decision is the worst of the reasons collected. Tally checked: 6 / 2 / 3 / 1.
5. A5 Suggestions are written by a person and kept next to the rule-based decision; the table
   is saved to `outputs/practice1_decisions.tsv` and compared with the checkpoint.
6. A6 Same tube, better assay: S10 on `RNA_CAPTURE_FFPE` → PASS (DV200 62%: 20–40 ng, tube holds
   225 ng); S03 on `ILMN_WES` → REVIEW (purity flag and the missing matched normal remain).

**Part B · what each platform returns** (S08 blood DNA, simulated over the real window)

1. B1 Window and target: FASTA header vs BED coordinates; the Day 1 bases re-checked.
2. B2 The four-line FASTQ record and Phred decoding (Q = ASCII − 33, p = 10^(−Q/10)).
3. B3 Students write `read_n50()` from three comments; two asserts check it.
4. B4 Average error probabilities, not Q values (worked example Q40, Q40, Q40, Q10).
5. B5 The 12 PAF columns; spanning = `tstart ≤ 5,226,000` and `tend ≥ 5,234,000`.
6. B6 One summary row per platform → `outputs/practice1_read_stats.tsv`.
7. B7 Log-scale read-length histogram → `outputs/read_length_hist.png`.
8. B8 *(optional)* Reads drawn over the target → `outputs/target_spanning_reads.png`.

### Expected Part A result

| tube | requested assay | decision | reasons (as the notebook writes them) |
|---|---|---|---|
| S01 | PACBIO_HIFI_WGS | **PASS** | - |
| S02 | ILMN_PCRFREE_WGS | **REVIEW** | REVIEW: a260_230 1.45 below 2.00 (clean up, re-measure) |
| S03 | PACBIO_HIFI_WGS | **FAIL** | FAIL: input 240 ng is 48% of the 500 ng minimum; FAIL: modal_fragment_kb 0.35 < 15 required (cannot be repaired); REVIEW: a260_230 1.95 below 2.00 (clean up, re-measure); REVIEW: tumour sample for a somatic question but no matched normal |
| S04 | ILMN_WES | **PASS** | - |
| S05 | ILMN_WES | **PASS** | - |
| S06 | ONT_LIGATION_LONG | **FAIL** | FAIL: input 20 ng is 2% of the 1,000 ng minimum; FAIL: modal_fragment_kb 0.167 < 10 required (cannot be repaired); REVIEW: a260_280 1.75 below 1.80 (clean up, re-measure); REVIEW: a260_230 1.60 below 2.00 (clean up, re-measure) |
| S07 | ONT_ULTRALONG | **PASS** | - |
| S08 | PACBIO_HIFI_WGS | **PASS** | - |
| S09 | RNA_POLYA | **PASS** | - |
| S10 | RNA_POLYA | **FAIL** | FAIL: rin 2.3 < 7 required (cannot be repaired) |
| S11 | PACBIO_HIFI_WGS | **REVIEW** | REVIEW: input 300 ng is 60% of the 500 ng minimum |
| S12 | ONT_LIGATION_LONG | **BLOCKED** | BLOCKED: no consent for genomic research (do not process, do not assess quality) |

Expected tally: PASS 6 (S01 S04 S05 S07 S08 S09), REVIEW 2 (S02 S11), FAIL 3 (S03 S06 S10),
BLOCKED 1 (S12).

### Deliberate cases in Part A (do not "fix" them)

| tube | what | correct response |
|---|---|---|
| S12 | perfect blood DNA, `consent_genomic_research = no` | **BLOCKED** before any quality check; its quality data are not used |
| S03 | FFPE tumour DNA: 240 ng, 0.35-kb fragments, PacBio requested, no matched normal | **FAIL**: mass might be fixed, length cannot; suggest an FFPE exome or panel plus a matched normal |
| S06 | plasma cfDNA (~167 bp, 20 ng) for 50-kb phasing on ONT | **FAIL**: no platform makes long reads from cfDNA; the alternative depends on the question |
| S10 | FFPE RNA, RIN 2.3, DV200 62%, poly(A) requested | **FAIL** on poly(A) (3′ bias); **PASS** on capture RNA-seq |
| S02 | saliva DNA, A260/230 = 1.45 | **REVIEW**: clean up and re-measure, never FAIL |
| S11 | buccal DNA, 300 ng for PacBio (60% of 500 ng) | **REVIEW**: re-quantify, concentrate, recollect, or switch assay |
| S04 / S05 | tumour and its matched normal, same subject P004 | **PASS**: the matched-normal check looks up the subject, not just the ID |
| S07 | 6.5 million cells for ultra-long ONT | **PASS**: cell input is counted, not weighed |

### Part B simulation and expected numbers

| platform | how the reads were simulated |
|---|---|
| Illumina 2×150 | 2,000 pairs; insert N(400, 25) bp clipped to 350–450; per-cycle Q falls from about Q38 to Q30 (mate 1) and Q36 to Q27 (mate 2), noise SD 2.5, 0.4% low-Q dips; substitutions drawn with p = 10^(−Q/10) |
| PacBio HiFi | 40 reads, length N(16 kb, 3 kb) clipped to 8–25 kb; Q ≈ N(45, 6) (30–60) outside homopolymers and N(24, 4) inside runs ≥ 4 bp; errors mostly homopolymer indels |
| ONT ligation (R10.4.1-like) | 50 reads, log-normal lengths (σ 0.8) tuned to a read N50 of about 20 kb, 1–58 kb; per-read mean Q drawn from 22–27, per-base SD 5, homopolymers 7 Q lower; indels enriched in homopolymers |

All reads come from both strands at uniform positions in the window. `alignments.paf` is written
from the simulation truth: full-length alignments, real chr11 coordinates (0-based, half-open,
target length 135,086,622), `nmatch` and `alnlen` from the true edit operations, MAPQ 255
("not available"). The ONT maximum stays below 58 kb only because the window is 60 kb; real
ultra-long libraries reach 100 kb to 4+ Mb.

| platform | reads | bases | mean length | max length | read N50 | mean accuracy | depth | reads spanning the 8-kb target |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Illumina 2×150 | 4,000 | 600,000 | 150.0 | 150 | 150 | 99.919% | 10.00× | 0 |
| PacBio HiFi | 40 | 646,752 | 16,168.8 | 21,467 | 16,522 | 99.929% | 10.78× | 6 |
| ONT ligation | 50 | 736,144 | 14,722.9 | 56,177 | 20,212 | 99.077% | 12.27× | 5 |

## Mini practice 2 · Cohort QC: who does not belong?

Modelled on the UK Biobank exome QC of Wang *et al.* 2021: 302,355 exomes; sequencing-centre
checks for sex discordance, contamination, unresolved duplicates and array discordance; then
FREEMIX > 4% (106 removed), < 94.5% of CCDS bases at ≥ 10× (5), KING kinship > 0.0884 with a
greedy most-relatives-first rule (14,326) → 287,917 of any ancestry; peddy Pr(European) < 0.99
(18,212) and ±4 SD on the top four PCs (535) → 269,171 European-ancestry participants; four
delivery batches checked for batch effects with per-batch dummy phenotypes.

The synthetic cohort: 2,000 exomes from 1,998 participants aged 40–69, six recruitment centres,
mostly blood with some saliva, batches B1–B4, collection years 2006–2010. About 92% sit in the
European PC cluster, with South Asian, African, East Asian and admixed clusters.
Kinship file: 3 duplicate/MZ pairs, 28 first-degree pairs (including a trio, two sibships of
three and a duplicated participant's sibling), 14 second-degree and 20 third-degree pairs.

Steps in the notebook:

1. Load and describe the cohort; what each metric means biologically.
2. Consent first: withdrawn samples leave before any QC, plot or count.
3. Contamination (FREEMIX > 0.04) and coverage (pct_ccds_10x < 94.5).
4. Sex check × array concordance: swap, metadata error, ambiguous (review), same-sex discordance.
5. Duplicates (kinship > 0.354; keep the better pct_ccds_10x) then relatives (greedy above 0.0884).
6. Ancestry (Pr(EUR) ≥ 0.99 and ±4 SD on PC1–PC4 → `outputs/pca.png`) and batch
   (`outputs/batch_ccds.png`); the QC funnel, the paper's funnel and its arithmetic, an equity
   note, and a cohort card students fill in.

### Expected funnel (`outputs/practice2_qc_funnel.tsv`)

| step | rule | removed | remaining |
|---|---|---:|---:|
| `0_sequenced` | all delivered exomes |  | 2,000 |
| `1_consent` | consent_status == withdrawn | 3 | 1,997 |
| `2_contamination` | freemix > 0.04 (VerifyBAMID) | 4 | 1,993 |
| `3_coverage` | pct_ccds_10x < 94.5 | 3 | 1,990 |
| `4_identity` | array_concordance < 0.95 | 3 | 1,987 |
| `5_duplicates` | kinship > 0.354: keep the better sample | 3 | 1,984 |
| `6_relatedness` | kinship > 0.0884: greedy, most relatives first | 38 | 1,946 |
| `pan_ancestry_set` | QC-passed, unrelated, any ancestry |  | 1,946 |
| `7_ancestry_pr_eur` | pr_eur < 0.99 (peddy) | 171 | 1,775 |
| `8_ancestry_pc_outlier` | outside mean +/- 4 SD on PC1-PC4 | 6 | 1,769 |
| `eur_analysis_set` | European-ancestry analysis set |  | 1,769 |

The pan-ancestry analysis set is 1,946; the European analysis set is 1,769 (90.9% of
the pan-ancestry set). `outputs/practice2_excluded_samples.tsv` lists every one of the other
samples exactly once with its step and reason; the 171 `7_ancestry_pr_eur` rows are
participants left out of the European analysis, not QC failures.

### Deliberate cases in Practice 2

| case | sample IDs | evidence in the data | where it lands | correct response |
|---|---|---|---|---|
| Consent withdrawn (3) | WES54841, WES54793, WES25052 | consent_status = withdrawn; metrics unremarkable | `1_consent` | remove first, before any QC; never plot or count them |
| Contaminated (4) | WES28515, WES52423, WES94009, WES92427 | FREEMIX 0.052, 0.067, 0.089, 0.118; the two male samples show chrX_het_rate 0.082, 0.100 (male baseline ≈ 0.02) | `2_contamination` | exclude; contamination is checked before the sex check because it inflates heterozygosity |
| Low coverage (3) | WES47545, WES24796, WES94990 | pct_ccds_10x 92.81% (B3), 93.94% (B3), 94.12% (B1) | `3_coverage` | exclude |
| Borderline coverage (1) | WES13035 | pct_ccds_10x 94.81% in B3 | kept | passes; do not move the threshold after looking |
| Swapped pair (2, consecutive IDs) | WES90876, WES90877 | WES90876: reported F, chrX_het 0.026, chrY 0.487, array 0.6123; WES90877: reported M, chrX_het 0.602, chrY 0.008, array 0.6047 | `4_identity` | exclude both; ask the lab to trace the plate |
| Metadata error (1) | WES76523 | reported F, genetic XY (chrX_het 0.025, chrY 0.474), array 0.9991 | kept (flag) | correct reported_sex to M in a copy, tell the study team, keep |
| XXY-like (1) | WES83617 | reported M, chrY 0.462 and chrX_het 0.552, array 0.9986 | kept (flag) | flag for review; never call it a swap, never force it into XX/XY |
| Same-sex array discordance (1) | WES31055 | reported M, genetic XY, array 0.5832 | `4_identity` | exclude: the sex check cannot see a same-sex swap |
| Technical duplicate with a sibling (2 + 1) | WES22141, WES74622, WES60985 | WES22141 and WES74622: same participant P891448, kinship 0.4986, pct_ccds_10x 98.91 vs 97.62; both copies ≈ 0.25 with sibling WES60985 (pct 98.60) | WES74622: `5_duplicates`; WES60985: `6_relatedness` | drop the worse copy first, so the sibling pair is counted once; the greedy step then removes the sibling (lower pct) |
| Technical duplicate (2) | WES27401, WES74705 | same participant P987958, kinship 0.4962, pct_ccds_10x 99.02 vs 98.37 | WES74705: `5_duplicates` | keep the better pct_ccds_10x |
| Different IDs, kinship ≈ 0.5 (2) | WES86470, WES33032 | participants P917121 / P838067, kinship 0.4931, same sex, age (56) and centre | WES33032: `5_duplicates`; WES86470 kept (flag) | MZ twins or one person enrolled twice: keep one (better pct), flag the pair for review |
| PC outliers with Pr(EUR) ≥ 0.99 (6) | WES69127, WES58717, WES81758, WES56035, WES58405, WES68025 | 4 off-cluster on PC3, 2 on PC4 (6.5-8.5 units; the European cluster has SD ≈ 1); invisible on PC1-PC2 | `8_ancestry_pc_outlier` | outside the European analysis set; still in the pan-ancestry set |
| European cluster, Pr(EUR) < 0.99 (12) | WES87840, WES70194, WES69974, WES89241, WES81247, WES28216, WES99566, WES66426, WES30620, WES73196, WES73904, WES90832 | Pr(EUR) 0.930-0.983 | `7_ancestry_pr_eur` | outside the European analysis set; not a QC failure |

Families: one trio (the greedy rule removes the child and keeps both parents), 13 European and 1
South Asian full-sibling pairs, two sibships of three (two removed from each), four
parent–offspring pairs and 14 second-degree pairs; 38 samples are removed for relatedness.
The 20 third-degree pairs (one of them African) stay, as the 0.0884 rule allows. B3 has lower
depth and lower pct_ccds_10x than B1, B2 and B4, yet almost all of it passes: keep it, report it,
adjust for it.

## How the notebooks map onto the slides

Slide numbers refer to the 40-slide Day 2 deck (*Molecular tools for genomics*, 10/13); titles are
given too, so a citation still resolves if slides move. `build/check_slide_refs.py` checks every
citation below and in both notebooks against the deck.

| Notebook part | Practice slide | Slides cited in the notebook |
|---|---|---|
| Practice 1, Part A (A1–A6) | 25 · Mini practice 1 · Is this sample fit for this platform? | 6 "Measure the molecule before you sequence", 7 "Metadata: the sample's passport", 22 "Match the sample to the platform" |
| Practice 1, Part B (B1–B8) | 25 · Mini practice 1 · Is this sample fit for this platform? | 21 "Side by side: which platform for which question?", 14 "Phred quality: how sure is each base?", 20 "Read length: mean versus N50", 23 "Depth: how much sequencing do I need?", 24 "Coverage is uneven: the Poisson view" |
| Practice 2, steps 1–3 | 32 · Mini practice 2 · Cohort QC: who does not belong? | 7 "Metadata: the sample's passport", 24 "Coverage is uneven: the Poisson view" |
| Practice 2, steps 4–5 | 32 · Mini practice 2 · Cohort QC: who does not belong? | 30 "Cohort QC: is every sample who it claims to be?" |
| Practice 2, step 6, funnel, equity note | 32 · Mini practice 2 · Cohort QC: who does not belong? | 31 "What a cohort can and cannot tell you" |

The reads-and-coverage theory behind Part B sits next to the platform slides: slide 14 (Phred
quality), slide 16 (insert size), slide 17 (duplicates), slide 20 (N50), slide 23 (depth) and
slide 24 (coverage). Platform background: slide 15 (Illumina), slide 18 (PacBio HiFi) and slide 19
(Oxford Nanopore). The one message of the day is slide 2 ("What is a genomics sample?"); the
take-home is slide 35.

## Instructor notes

* **Timing.** Practice 1: Part A about 8 minutes, Part B about 7. Practice 2: steps 1–3 about
  4 minutes, step 4 about 4, step 5 about 3, step 6 and the funnel about 4. If time is short,
  skip B8 in Practice 1 and the paper-arithmetic cell in Practice 2; neither feeds a later cell.
* **Computing is not the bottleneck.** Each notebook runs top to bottom in a few seconds; the
  15 minutes are for reading outputs and deciding.
* **Answers to say out loud.** Practice 1: tally 6 / 2 / 3 / 1; S10 passes on capture RNA-seq;
  0 Illumina reads span the 8-kb target versus 6 HiFi and 5 ONT at similar depth (10–12×); HiFi is
  the platform most flattered by averaging Q values (naive Q43, honest Q31). Practice 2: funnel
  2,000 → 1,946 → 1,769; three samples are kept with a review flag (metadata error, XXY-like,
  kept twin); the paper's own subtractions are off by one sample in two places.
* **Where students get stuck.**
  * Reading several reasons for one tube. The rule is "the worst reason wins", and S03 has four.
  * DIN versus fragment length: integrity is a score, length is what limits a read.
  * `NA` versus 0: S07 has no mass because cells are counted, not weighed.
  * N50 versus mean: ask "in what length of read does a typical *base* sit?"
  * PAF coordinates are 0-based, half-open and on real chr11, not positions in the 60-kb string.
  * Treating ancestry rows as QC failures. They are participants left out of one analysis.
  * Wanting to "fix" the XXY-like sample or to move the 94.5% line to catch B3's 94.81% sample.
  * Running cells out of order in Practice 2. Removal steps refuse to run twice and say so;
    *Kernel → Restart & Run All* resets everything.
* **Good questions for the room.** Which REVIEW is a clean-up and which is a mass problem
  (S02 vs S11)? Why does S06's alternative depend on the study question? What can the sex check
  not see that the array check can? Who is missing from a European-only analysis?

## Rebuild and verify

```bash
python build/build_practice1.py && python build/build_practice2.py && \
python build/make_notebooks.py && python build/verify_package.py
```

* Everything is deterministic (fixed seeds; gzip written with mtime 0), so a rebuild reproduces
  the same data, checkpoints and notebooks.
* `build_practice1.py` reads the real window from `data/practice1/` and refuses to run if its
  checksum changed. To bootstrap from the Day 1 package instead:
  `python build/build_practice1.py --day1-seq <hg38_5218771_5278770.txt> --day1-genes <genes.hg38.bed>`.
* `verify_package.py` copies the package to a temporary folder with an empty `outputs/`, runs
  both notebooks headless with `nbclient` (kernel `python3`), and checks: no errors and nothing on
  stderr, under 60 s in total, every output table equal to its checkpoint, the teaching facts
  (tally, spanning reads, funnel arithmetic, every planted case at its documented step), all
  figures written, no stored outputs, and no network calls. It exits non-zero on any failure.
  On machines that forbid Jupyter kernel sockets it says so and runs the same cells through
  IPython's `InteractiveShell` in a fresh interpreter instead (`--executor` chooses explicitly).
