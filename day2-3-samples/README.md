# Practice package — Genomic Data Processing, Day 2 (10/13)

Three mini practice sessions, one notebook per session. Practices 1 and 2 take about 15 minutes
each; Practice 3 is a 60-minute session in pairs that ends in two-minute lightning talks. Open the
notebooks from **this folder** (`day2-3-samples/`): every path inside them is relative to it.

| Notebook | Practice slide | The question students answer |
|---|---|---|
| `practice1_sample_to_platform.ipynb` | 25 · Mini practice 1 · Is this sample fit for this platform? | Which tubes can go on which platform, and what does each platform return from the same DNA? |
| `practice2_cohort_qc.ipynb` | 32 · Mini practice 2 · Cohort QC: who does not belong? | Which of 2,000 exomes must leave the analysis, at which step, and why? |
| `practice3_dna_rna_short_long.ipynb` | 33 · Mini practice 3 · One person, six datasets | One person's DNA and RNA, each read by Illumina, PacBio and Nanopore: what can each dataset see, and what does it miss? |

## How to run

```bash
cd day2-3-samples
python -m pip install -r requirements.txt      # pandas, numpy, matplotlib
jupyter lab                                    # or any notebook front-end
```

Open `practice1_sample_to_platform.ipynb`, run it top to bottom (Shift+Enter), then
`practice2_cohort_qc.ipynb`, then `practice3_dna_rna_short_long.ipynb`. Everything the notebooks
write lands in `outputs/`.

## Requirements

| Needed | Why |
|---|---|
| Python 3.9+ with `pandas` ≥ 1.3, `numpy` ≥ 1.21, `matplotlib` ≥ 3.4 | all three notebooks |
| JupyterLab or Jupyter Notebook *(optional)* | to open the notebooks |
| `nbformat`, `nbclient`, `ipykernel` *(maintainers only)* | `build/make_notebooks.py` and `build/verify_package.py` |
| `python-pptx` *(maintainers only)* | `build/check_slide_refs.py` |

No network access, no reference download, no aligner, no `pysam`, no index build: Practice 3
reads its SAM files with pandas. Every cell finishes in a few seconds at most, and each notebook
runs top to bottom in under 10 seconds.

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
("imagine a structural variant here"); nothing is claimed about the region. **All of Practice 3
is synthetic**: the person P3-01, the 120-kb locus `chrS` and its four genes, every read, every
alignment and every methylation call. Its read models approximate public platform specifications
but are **not calibrated to any real run**.

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
├── practice3_dna_rna_short_long.ipynb  Mini practice 3
├── presentation/
│   └── practice3_lightning_talk_template.pptx   the lightning-talk template: a worked example, then the blank slide
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
│   ├── practice2/
│   │   ├── cohort_manifest.tsv         2,000 samples from 1,998 participants
│   │   ├── sample_qc_metrics.tsv       depth, CCDS 10x, FREEMIX, sex metrics, array concordance, PCs, Pr(EUR)
│   │   └── kinship_pairs.tsv           65 pairs with KING kinship > 0.0442
│   └── practice3/                      synthetic person P3-01, synthetic 120-kb locus chrS (2.9 MB)
│       ├── samples.tsv                 the six libraries: person, tube, molecule, platform, chemistry, prep
│       ├── chrS.fa (+ .fai)            the locus sequence
│       ├── genes.gtf                   ALTX (two isoforms), ASE1, OFF1, LOW1: transcript and exon lines
│       ├── regions.tsv                 CGI_ALTX, REP1, HP1, DUP_A, DUP_B, OFF1_PROMOTER, in plain words
│       ├── sites.tsv                   SNV1 and SNV2
│       ├── DNA-ILMN.sam.gz             12,000 Illumina read pairs, 2 x 150 bp (~30x)
│       ├── DNA-HIFI.sam.gz             169 PacBio HiFi reads (~20x)
│       ├── DNA-ONT.sam.gz              183 Nanopore reads (~20x)
│       ├── RNA-ILMN.sam.gz             10,000 stranded poly(A) RNA-seq fragments, 2 x 100 bp
│       ├── RNA-KINNEX.sam.gz           1,500 PacBio Kinnex full-length cDNA reads
│       ├── RNA-ONT.sam.gz              2,000 Nanopore direct RNA reads
│       └── methylation_DNA-{HIFI,ONT}.tsv   5mC per CpG: n_reads, n_methylated, fraction
├── checkpoints/                        verified outputs, for when a cell will not run
│   ├── practice1_decisions.tsv
│   ├── practice1_read_stats.tsv
│   ├── practice2_qc_funnel.tsv
│   ├── practice2_excluded_samples.tsv
│   ├── practice2_planted_cases.tsv     answer key: every planted case and the step it must reach
│   ├── practice3_checks.tsv            every Practice 3 CHECK number (282)
│   ├── practice3_scorecard.tsv         what each dataset can and cannot see
│   └── practice3_truth.json            answer key: what was built into chrS, the read models, the facts
├── outputs/                            what the notebooks write (only .gitkeep is tracked)
└── build/
    ├── build_practice1.py              sample sheet, requirements, reads, PAF, Practice 1 checkpoints
    ├── build_practice2.py              cohort, metrics, kinship, Practice 2 checkpoints + answer key
    ├── build_practice3.py              locus, six libraries, Practice 3 checkpoints + answer key
    ├── make_notebooks.py               writes the three notebooks with nbformat (no outputs stored)
    ├── practice3_cells.py              the cells of the Practice 3 notebook (used by make_notebooks.py)
    ├── verify_package.py               runs the notebooks headless from a clean copy and compares
    └── check_slide_refs.py             checks every slide citation against the Day 2 deck
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

## Mini practice 3 · One person, six datasets

**Purpose.** One synthetic person, P3-01, gave one blood draw: DNA from an EDTA tube and RNA from
a PAXgene tube. Each molecule was read three ways (Illumina, PacBio, Oxford Nanopore), giving six
datasets of the same 120-kb locus. DNA tells you what the genome contains; RNA tells you what the
blood cells are doing now. Long reads span repeats and duplications, link variants and read whole
transcripts; short reads are accurate, numerous and cheap per read. Each pair studies one of 12
figures and presents it as a claim backed by one number, in two minutes.

### The six datasets

| library | molecule | platform · chemistry (teaching label) | layout / prep | in the locus | read model (approximate, not calibrated) |
|---|---|---|---|---|---|
| `DNA-ILMN` | DNA | Illumina NovaSeq X · 2 x 150 bp | PCR-amplified library | 12,000 pairs, ~30x | insert N(350, 50); substitutions, 0.1% at base 1 rising to 0.8% at base 150; GC bias exp(-0.5((gc-0.45)/0.18)^2) |
| `DNA-HIFI` | DNA | PacBio Revio · HiFi reads, 15-20 kb library | PCR-free; 5mC from polymerase kinetics | 169 reads, ~20x | length N(16, 3) kb, 8-30 kb; ~Q33, mostly 1-bp indels, some homopolymer slips |
| `DNA-ONT` | DNA | Oxford Nanopore PromethION · R10.4.1, Ligation Sequencing Kit V14 | PCR-free; 5mC from the current signal | 183 reads, ~20x | log-normal length, N50 ~22 kb, 1-100 kb; ~Q20 (97-99.7%), 45% deletions, homopolymer length errors |
| `RNA-ILMN` | RNA | Illumina NovaSeq X · 2 x 100 bp, stranded poly(A) mRNA | cDNA fragments ~250 bp; read 1 antisense | 10,000 fragments | 3% unspliced pre-mRNA; mild 3' bias (1.2x) |
| `RNA-KINNEX` | RNA | PacBio Revio · Kinnex full-length RNA | full-length cDNA, poly(A) trimmed | 1,500 reads | 85% full length, 15% 5'-truncated; HiFi accuracy |
| `RNA-ONT` | RNA | Oxford Nanopore PromethION · Direct RNA, SQK-RNA004 | native RNA read 3' to 5', written as DNA letters | 2,000 reads | starts at the poly(A) end, ~30% reach the 5' end; ~97% accurate, deletions most common |

### What pairs do (60 minutes)

| time | activity | deck |
|---|---|---|
| 0-5 | the six datasets, the one-slide recipe, instructor demo of a finished slide on F01 | 33 "Mini practice 3 · One person, six datasets", 34 "Lightning talks: one figure, two minutes" |
| 5-30 | pairs run the notebook (Run All takes about 4 s), study their figure, do its **Try this**, fill one slide | notebook |
| 30-55 | lightning talks: 2 minutes per pair plus about 15 s changeover (11-12 pairs) | student slides |
| 55-60 | debrief: what each dataset can and cannot see, with the scorecard | 35 "Debrief: what each dataset can and cannot see" |

Figures are handed out in this order: F03, F05, F07, F04, F08, F02, F09, F06, F10, F11, F12, F01
(F01 is the instructor's worked example). With more than 12 pairs, two pairs share F05 or F07 and
compare their claims.

**The lightning-talk rule: one slide, two minutes.** The slide holds the figure, a **claim** in one
sentence with **one number** from `CHECK`, and a **so what** with **one limit**: which dataset you
would choose for this question, and one thing it cannot do. The one-slide template is
`presentation/practice3_lightning_talk_template.pptx`: its first page is a worked example on F01, its second the page each pair fills.

### The twelve figures

Each section of the notebook has the question, two or three sentences of biology, one code cell
that draws the figure, saves it and prints its `CHECK` numbers, a **Try this** (one parameter to
change), three prompts for the talk, and the slides to look back at. The notebook asks; it never
prints the answer.

| figure | question | CHECK numbers | slides |
|---|---|---|---|
| F01 · Read length: six datasets | How long are the reads, and why are long RNA reads shorter than long DNA reads? | read count, median and N50 per dataset | 13 "Three ways to read a molecule", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads", 20 "Read length: mean versus N50" |
| F02 · Errors: what kind of mistakes? | Which datasets make mostly insertion/deletion errors, and where do they show up? | mismatches, insertions, deletions per 1,000 aligned bases; % of reads with 9, 10, 11 A's at HP1 | 14 "Phred quality: how sure is each base?", 15 "Illumina: short, accurate, massively parallel", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads" |
| F03 · DNA vs RNA: coverage across the locus | OFF1 has no RNA reads. Is the gene missing from this person's genome? | mean exon depth per gene; % of aligned bases in exons | 2 "What is a genomics sample?", 4 "Human specimens and what they give you", 9 "RNA sequencing" |
| F04 · Allele-specific expression | Why do DNA and RNA give different allele fractions at SNV1? | reads with G and A at SNV1, fraction, 95% interval | 9 "RNA sequencing", 10 "RNA-seq applications" |
| F05 · Isoforms: junctions vs full-length reads | From short reads alone, how many ALTX isoforms are possible? What do long reads add? | reads per junction; reads per exon chain | 10 "RNA-seq applications", 11 "One gene, many transcripts" |
| F06 · RNA coverage along the transcript | Why does one RNA dataset have more reads near the 3' end? | 3'/5' coverage ratio; % of reads reaching the 5' end | 9 "RNA sequencing", 19 "Oxford Nanopore: the longest reads" |
| F07 · Sizing a CAG repeat | Is one allele expanded, and which dataset can measure it? | spanning reads, modal size and spread per allele; Illumina all-CAG reads | 13 "Three ways to read a molecule", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads", 21 "Side by side: which platform for which question?" |
| F08 · A duplicated block | Why do short reads disappear here when you filter on mapping quality? | MAPQ ≥ 20 depth inside DUP_A and in its flanks, and their ratio | 15 "Illumina: short, accurate, massively parallel", 21 "Side by side: which platform for which question?" |
| F09 · Phasing: are two variants on the same copy? | Are SNV1 and SNV2 on the same copy of the chromosome? | reads covering both sites; the four combinations of bases | 16 "Paired-end reads and insert size", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads" |
| F10 · GC content and coverage | Why does one dataset lose coverage in GC-rich DNA? | normalized depth in ≥ 60% GC and 40-50% GC windows, ratio, window-to-window scatter | 15 "Illumina: short, accurate, massively parallel", 17 "Duplicates and library complexity", 24 "Coverage is uneven: the Poisson view" |
| F11 · Methylation from long-read DNA | Which promoter is methylated, and does that fit the RNA data? | CpGs and mean methylation per region | 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads", 21 "Side by side: which platform for which question?" |
| F12 · Gene expression: do the RNA datasets agree? | Do the three RNA datasets agree on which genes are on, and how sure is each? | reads per gene, share of the locus's RNA, 95% interval | 9 "RNA sequencing", 10 "RNA-seq applications", 22 "Match the sample to the platform" |

Comparison groups: DNA vs RNA = F03, F04; short vs long reads, DNA = F07, F08, F09; short vs long
reads, RNA = F05, F06, F12; platforms = F01, F02, F10, F11.

### Files the notebook writes

* `outputs/practice3_F01_read_length.png`, `practice3_F02_errors.png`,
  `practice3_F03_dna_vs_rna_coverage.png`, `practice3_F04_allele_specific_expression.png`,
  `practice3_F05_isoforms.png`, `practice3_F06_transcript_coverage.png`,
  `practice3_F07_cag_repeat.png`, `practice3_F08_duplicated_block.png`, `practice3_F09_phasing.png`,
  `practice3_F10_gc_coverage.png`, `practice3_F11_methylation.png`, `practice3_F12_expression.png`:
  8 x 4.5 in at 200 dpi, text 11-13 pt, ready for a slide. Colour = platform (Illumina `#1F5FEA`,
  PacBio `#B25E09`, ONT `#17202A`); DNA solid lines and filled markers, RNA dashed and open.
* `outputs/practice3_checks.tsv`: every CHECK number (282), compared at the end of the notebook
  with `checkpoints/practice3_checks.tsv`.
* `outputs/practice3_scorecard.tsv`: one yes / no / n/a per question per dataset.

### How the data were made, and the honest limits

* `build/build_practice3.py` simulates the two haplotypes, draws reads from them with the models
  above, and writes each read's alignment by projecting it onto `chrS`; there is no aligner. The
  180-unit allele therefore appears as an insertion inside REP1, an insertion at a read end becomes
  a soft clip, and a short read lying entirely inside the expanded repeat is unmapped (FLAG 4,
  placed at its mate). A short read inside one copy of the duplication that avoids the three
  differing bases gets MAPQ 0 and is placed in either copy at random, as an aligner would.
* `samtools` 1.24 reads all six files, sorts and indexes them, and `samtools calmd` recomputes
  every `NM` tag without a difference; every intron in the RNA files is an annotated one.
* **Base qualities.** Illumina reads carry a per-cycle quality that matches the error model.
  Long reads carry one quality value per read (the Phred value of that read's simulated error
  rate) instead of per-base values: none of the figures uses base qualities, and a uniform string
  keeps the files small (2.9 MB in total) without pretending to model base-by-base confidence.
  Real HiFi and nanopore files vary base by base.
* Illumina reads carry substitutions only. The per-cycle model from the design (0.1% rising to
  0.8%) averages 99.56% correct per base, a little below what a good NovaSeq X run achieves.
* Long DNA reads start on an even, jittered grid and alternate between the two haplotypes. A real
  run scatters reads at random; with only ~20 reads over any base that noise alone could hide a
  lesson (a heterozygous site reading 75:25, say), so the teaching data spread reads evenly.
  Nanopore depth still wanders (15-23x) because its read lengths vary so much; F10 reports the
  window-to-window scatter so students can tell noise from bias.
* Methylation tables summarise simulated per-molecule calls with 3% (HiFi) and 5% (ONT) call error.
* All four genes are on the + strand. The build stops if any of the 16 teaching facts is no longer
  visible in the data, and `verify_package.py` tests the same facts on the notebook's output.

### For instructors

The answer key for the twelve figures, the design truth of `chrS` and notes on where students get stuck are in a
separate instructor guide, kept outside this public repository because in this practice students form their own
claims. `checkpoints/practice3_truth.json` is the build record that `verify_package.py` uses; open it after the
talks, not before.

## How the notebooks map onto the slides

Slide numbers refer to the 43-slide Day 2 deck (*Molecular tools for genomics*, 10/13), in which
slides 33-35 hold Mini practice 3; titles are given too, so a citation still resolves if slides
move. `build/check_slide_refs.py` checks every citation below and in the three notebooks against
the deck.

| Notebook part | Practice slide | Slides cited in the notebook |
|---|---|---|
| Practice 1, Part A (A1–A6) | 25 · Mini practice 1 · Is this sample fit for this platform? | 6 "Measure the molecule before you sequence", 7 "Metadata: the sample's passport", 22 "Match the sample to the platform" |
| Practice 1, Part B (B1–B8) | 25 · Mini practice 1 · Is this sample fit for this platform? | 21 "Side by side: which platform for which question?", 14 "Phred quality: how sure is each base?", 20 "Read length: mean versus N50", 23 "Depth: how much sequencing do I need?", 24 "Coverage is uneven: the Poisson view" |
| Practice 2, steps 1–3 | 32 · Mini practice 2 · Cohort QC: who does not belong? | 7 "Metadata: the sample's passport", 24 "Coverage is uneven: the Poisson view" |
| Practice 2, steps 4–5 | 32 · Mini practice 2 · Cohort QC: who does not belong? | 30 "Cohort QC: is every sample who it claims to be?" |
| Practice 2, step 6, funnel, equity note | 32 · Mini practice 2 · Cohort QC: who does not belong? | 31 "What a cohort can and cannot tell you" |
| Practice 3, set-up and sample sheet | 33 · Mini practice 3 · One person, six datasets | 34 "Lightning talks: one figure, two minutes", 7 "Metadata: the sample's passport" |
| Practice 3, DNA vs RNA (F03, F04) | 33 · Mini practice 3 · One person, six datasets | 2 "What is a genomics sample?", 4 "Human specimens and what they give you", 9 "RNA sequencing", 10 "RNA-seq applications" |
| Practice 3, short vs long reads, DNA (F07, F08, F09) | 33 · Mini practice 3 · One person, six datasets | 13 "Three ways to read a molecule", 15 "Illumina: short, accurate, massively parallel", 16 "Paired-end reads and insert size", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads", 21 "Side by side: which platform for which question?" |
| Practice 3, short vs long reads, RNA (F05, F06, F12) | 33 · Mini practice 3 · One person, six datasets | 9 "RNA sequencing", 10 "RNA-seq applications", 11 "One gene, many transcripts", 19 "Oxford Nanopore: the longest reads", 22 "Match the sample to the platform" |
| Practice 3, platforms (F01, F02, F10, F11) | 33 · Mini practice 3 · One person, six datasets | 13 "Three ways to read a molecule", 14 "Phred quality: how sure is each base?", 15 "Illumina: short, accurate, massively parallel", 17 "Duplicates and library complexity", 18 "PacBio HiFi: long and accurate", 19 "Oxford Nanopore: the longest reads", 20 "Read length: mean versus N50", 21 "Side by side: which platform for which question?", 24 "Coverage is uneven: the Poisson view" |
| Practice 3, scorecard and debrief | 33 · Mini practice 3 · One person, six datasets | 35 "Debrief: what each dataset can and cannot see" |

The reads-and-coverage theory behind Part B sits next to the platform slides: slide 14 (Phred
quality), slide 16 (insert size), slide 17 (duplicates), slide 20 (N50), slide 23 (depth) and
slide 24 (coverage). Platform background: slide 15 (Illumina), slide 18 (PacBio HiFi) and slide 19
(Oxford Nanopore). The one message of the day is slide 2 ("What is a genomics sample?"); the
take-home is slide 38.

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
* **Practice 3, timing.** 5 minutes of set-up with the F01 demo, 25 minutes in pairs, 25 minutes
  of talks (keep a visible two-minute timer), 5 minutes of debrief with the scorecard. The
  notebook runs in about 4 seconds; the 25 minutes are for reading one figure well.
## Rebuild and verify

```bash
python build/build_practice1.py && python build/build_practice2.py && \
python build/build_practice3.py && python build/make_notebooks.py && \
python build/verify_package.py
python build/check_slide_refs.py "<path>/Day 2 deck.pptx"        # needs python-pptx
```

* Everything is deterministic (fixed seeds; gzip written with mtime 0), so a rebuild reproduces
  the same data, checkpoints and notebooks.
* `build_practice3.py` takes about 7 seconds. It computes every CHECK number a second time with
  its own code (from the SAM files it has just written), writes `checkpoints/practice3_*`, and
  stops if any of the 16 teaching facts is no longer visible.
* `build_practice1.py` reads the real window from `data/practice1/` and refuses to run if its
  checksum changed. To bootstrap from the Day 1 package instead:
  `python build/build_practice1.py --day1-seq <hg38_5218771_5278770.txt> --day1-genes <genes.hg38.bed>`.
* `verify_package.py` copies the package to a temporary folder with an empty `outputs/`, runs
  the three notebooks headless with `nbclient` (kernel `python3`), and checks: no errors and
  nothing on stderr; Practices 1 and 2 under 60 s together and Practice 3 under 30 s; every output
  table equal to its checkpoint; the teaching facts (tally, spanning reads, funnel arithmetic,
  every planted case at its documented step, and Practice 3's 16 facts tested on the notebook's
  own CHECK numbers); all 16 figures written; no stored outputs; only the allowed imports; and no
  network calls. It exits non-zero on any failure.
  On machines that forbid Jupyter kernel sockets it says so and runs the same cells through
  IPython's `InteractiveShell` in a fresh interpreter instead (`--executor` chooses explicitly).
