# Practice 2 · Instructor notes

*CNVs from reads to interpretation – Day 4 (10/15), "Understanding genetic variants". All samples are synthetic.*

## 1. At a glance

* **Run time.** The reference run of the whole notebook took **24 s** of computing (12-core
  laptop; mapping four samples ≈ 18 s). Expect 1–3 min on a student laptop. Classroom time ≈ 75 min:
  Step 0 5 · 1 9 · 2 12 · 3 9 · 4 10 · 5 9 · 6 6 · 7 9 · 8 4 · 9 2 min.
* **Headline results.** All five planted CNVs in the patients are recovered with the right copy number, no false calls;
  leave-one-out on the five controls finds only C2's UGT2B17 deletion. The 22q11.2 deletion of P1 lifts to
  chr22:18,927,488–21,110,711 (2.18 Mb) by breakpoint lift; the whole-interval CrossMap lift gives 117 pieces.
  RNA: expressed genes in the 22q11.2 deletions have median log2FC -0.85 (P1) and
  -0.8 (P3); PMP22 +0.55; UGT2B17 in P3 has 0 reads.
* **TESSERACT-CNV is not run** in the shipped notebook (no public code/weights). The adapter cell prints
  "TESSERACT-CNV not run …" and every downstream cell continues. No score exists anywhere in the package.

## 2. Expected CHECK values, step by step

### Step 0
Chromosome lengths chr4 191,154,276 · chr17 81,195,210 · chr22 51,304,566; 6.9 Mb of real sequence out of 324 Mb.

### Step 1 · mapping (hg19)
| sample | read pairs | % mapped | % properly paired | % duplicates | mean depth (x) |
|---|---|---|---|---|---|
| P1 | 93857 | 99.97 | 99.93 | 2.78 | 2.64 |
| P2 | 99683 | 99.96 | 99.92 | 2.78 | 2.81 |
| P3 | 104911 | 99.96 | 99.91 | 2.74 | 2.96 |
| C1 | 106631 | 99.96 | 99.92 | 2.79 | 3 |

Mean depth is lower in P1 (one copy of 2.5 Mb is missing) and differs between samples on purpose (per-sample depth factors
0.88–1.15, so library-size normalisation matters). MAPQ, control C1:
| where | reads | % MAPQ < 20 | median MAPQ |
|---|---|---|---|
| inside segmental duplications | 56234 | 37.5 | 50 |
| unique sequence | 156946 | 0.2 | 60 |

`samtools view -c` demo (C1): LCR22A chr22:18700001-18800000: 3,187 reads, 1,025 with MAPQ ≥ 20; unique chr22:19700001-19800000: 2,385 reads, 2,383 with MAPQ ≥ 20.

### Step 2 · depth → copy-number signal
* 1,380 bins of 5 kb. Fragments counted: P1 81,747 · P2 87,019 · P3 92,799 · C1 93,049.
  **Your C1 = panel C1 (identical vectors).** samtools read counts in chr17:15,145,001–15,150,000:
  P2 163, C1 131
  (fragments: P2 79, C1 63).
* Masked bins: gap 21 · low panel depth 166 · variable controls
  6 · **total 193**. 97% of the low-depth bins lie in segmental duplications.
* Centring offsets (log2): P1 +0.269 · P2 +0.007 ·
  P3 +0.166 · C1 -0.011.
  (The deletion removes ~18 % of P1's fragments in this 6.9-Mb mini-genome, which inflates normal bins by about log2(1/0.82) ≈ 0.29; the measured median offset is 0.27.)
* GC effect in C1 (range of 10 GC-group medians): raw depth 39% → ratio to panel
  5%. Spread of C1 log2 ratios: SD 0.19 (σ for CN 2 = 0.1946).
* Worked examples:
| sample | bin (hg19, 1-based) | what | raw | per million | panel | ratio | log2 | log2 centred | expected |
|---|---|---|---|---|---|---|---|---|---|
| P2 | chr17:15,150,001-15,155,000 | PMP22 duplication (CN 3) | 122 | 1402 | 841.5 | 1.67 | 0.74 | 0.73 | 0.58 |
| P1 | chr22:19,755,001-19,760,000 | TBX1, 22q11.2 deletion (CN 1) | 43 | 526 | 758.8 | 0.69 | -0.53 | -0.8 | -1 |
| P3 | chr4:69,415,001-69,420,000 | UGT2B17, both copies deleted (CN 0) | 0 | 0 | 696.6 | 0 | -9.45 | -3 | −∞ (−3) |
| C1 | chr17:15,150,001-15,155,000 | the same PMP22 bin in control C1 (CN 2) | 69 | 741.5 | 841.5 | 0.88 | -0.18 | -0.17 | 0 |

  Single bins scatter by ±0.2–0.4; that is the point of Step 3.

### Step 3 · HMM calls (hg19, BED 0-based)
| cnv_id | chrom | start | end | cn | n_bins | log2 |
|---|---|---|---|---|---|---|
| P1_22q11_DEL | chr22 | 18915000 | 21465000 | 1 | 433 | -0.988 |
| P2_PMP22_DUP | chr17 | 15130000 | 15170000 | 3 | 8 | 0.549 |
| P2_UGT2B17_DEL | chr4 | 69375000 | 69485000 | 1 | 19 | -1.076 |
| P3_22q11_DEL | chr22 | 21030000 | 22480000 | 1 | 211 | -0.993 |
| P3_UGT2B17_DEL | chr4 | 69375000 | 69485000 | 0 | 19 | -2.933 |

Comparison with the planted truth (reciprocal overlap ≥ 50 %, same type):
| planted | expected CN | called CN | bins | mean log2 | recip. overlap | start error (kb) | end error (kb) |
|---|---|---|---|---|---|---|---|
| P1_22q11_LCR22A-D_DEL | 1 | 1 | 433 | -0.988 | 0.999 | -1.8 | -0.7 |
| P2_PMP22_DUP | 3 | 3 | 8 | 0.549 | 0.89 | -3.1 | 1.3 |
| P2_UGT2B17_DEL | 1 | 1 | 19 | -1.076 | 0.983 | -0.2 | 1.8 |
| P3_22q11_distal_DEL | 1 | 1 | 211 | -0.993 | 0.999 | 0.3 | -1.5 |
| P3_UGT2B17_DEL_HOM | 0 | 0 | 19 | -2.933 | 0.983 | -0.2 | 1.8 |

False calls in patients: 0; calls in C1: 0; leave-one-out on all controls:
only C2 chr4:69,375,001–69,485,000 (CN 1).

### Step 4 · liftover
Called CNVs (rule used downstream = **breakpoint lift**):
| cnv_id | hg19 | (a) pieces | (a) pieces off main chrom | (a) lifted bp (kb) | (b) breakpoint lift | (b) size (kb) | (c) outer span | (c) size (kb) | genes hg19 | genes hg38 (b) | genes hg38 (c) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P1_22q11_DEL | chr22:18,915,001-21,465,000 | 117 | chr22_KI270731v1_random:120986-121187 | 2449.8 | chr22:18,927,488-21,110,711 | 2183.2 | chr22:18,339,130-21,343,827 | 3004.7 | 48 | 43 | 48 |
| P2_PMP22_DUP | chr17:15,130,001-15,170,000 | 1 | - | 40 | chr17:15,226,684-15,266,683 | 40 | chr17:15,226,684-15,266,683 | 40 | 1 | 1 | 1 |
| P2_UGT2B17_DEL | chr4:69,375,001-69,485,000 | 1 | - | 110 | chr4:68,509,283-68,619,282 | 110 | chr4:68,509,283-68,619,282 | 110 | 1 | 1 | 1 |
| P3_22q11_DEL | chr22:21,030,001-22,480,000 | 97 | chr1:226114942-226114970 | 1449.6 | chr22:20,675,713-22,125,606 | 1449.9 | chr22:20,668,232-22,125,606 | 1457.4 | 24 | 24 | 24 |
| P3_UGT2B17_DEL | chr4:69,375,001-69,485,000 | 1 | - | 110 | chr4:68,509,283-68,619,282 | 110 | chr4:68,509,283-68,619,282 | 110 | 1 | 1 | 1 |

`.unmap` file: **empty** for every CNV (the key teaching point: pieces in chain gaps are dropped silently).
Gene content, P1: only in hg19 → FAM246B, FAM246C, GGTLC3, RIMBP3, TMEM191B (all 100 % inside segmental duplications; in GRCh38 their
MANE copies sit at 18.5–18.6 Mb in LCR22A, FAM246C has no MANE transcript); extra in the outer span → DGCR6 (outside the
deletion in hg19) plus the same LCR genes. P3, PMP22, UGT2B17: identical gene lists.

Public records lifted the same way:
| cnv_id | hg19 | (a) pieces | (b) breakpoint lift | (b) size (kb) | (c) outer span | (c) size (kb) | public GRCh38 |
|---|---|---|---|---|---|---|---|
| P1_22q11_LCR22A-D_DEL | chr22:18,916,842-21,465,662 | 117 | chr22:18,929,329-21,111,373 | 2182 | chr22:18,339,130-21,343,827 | 3004.7 | ClinVar VCV000155158: 18,339,130-21,111,373 |
| P2_PMP22_DUP | chr17:15,133,094-15,168,674 | 1 | chr17:15,229,777-15,265,357 | 35.6 | chr17:15,229,777-15,265,357 | 35.6 | LOVD PMP22_000101: 15,229,777-15,265,357 |
| P2_UGT2B17_DEL | chr4:69,375,169-69,483,248 | 1 | chr4:68,509,451-68,617,530 | 108.1 | chr4:68,509,451-68,617,530 | 108.1 | ClinVar VCV000147015: 68,509,451-68,617,530 |
| P3_22q11_distal_DEL | chr22:21,029,655-22,481,498 | 97 | chr22:20,675,367-22,127,104 | 1451.7 | chr22:20,668,232-22,127,104 | 1458.9 | - |

ClinVar's GRCh38 start for VCV000155158 (18,339,130) is exactly the first base of the most proximal CrossMap piece, and its
end (21,111,373) is exactly the lifted end breakpoint.

### Step 5 · RNA (synthetic)
| cnv_id | genes | expressed | median log2FC (expressed) | fits copy number |
|---|---|---|---|---|
| P1_22q11_DEL | 43 | 32 | -0.85 | 26 |
| P2_PMP22_DUP | 1 | 1 | 0.55 | 1 |
| P2_UGT2B17_DEL | 1 | 1 | -0.37 | 0 |
| P3_22q11_DEL | 24 | 19 | -0.8 | 14 |
| P3_UGT2B17_DEL | 1 | 1 | -3.75 | 1 |

Key genes:
| cnv_id | gene | expressed | control CPM | patient CPM | log2FC | z | reading |
|---|---|---|---|---|---|---|---|
| P1_22q11_DEL | HIRA | True | 25.2 | 19 | -0.39 | -5 | weaker/other (buffered? noise?) |
| P1_22q11_DEL | TBX1 | True | 20 | 10.7 | -0.87 | -6.7 | fits the copy number |
| P1_22q11_DEL | COMT | True | 141.2 | 73.4 | -0.93 | -6.4 | fits the copy number |
| P1_22q11_DEL | CRKL | True | 62.6 | 34.6 | -0.84 | -7.8 | fits the copy number |
| P2_PMP22_DUP | PMP22 | True | 396.7 | 579.9 | 0.55 | 4.2 | fits the copy number |
| P2_UGT2B17_DEL | UGT2B17 | True | 6.5 | 4.7 | -0.37 | -0.9 | weaker/other (buffered? noise?) |
| P3_22q11_DEL | CRKL | True | 62.6 | 35.6 | -0.8 | -7.4 | fits the copy number |
| P3_UGT2B17_DEL | UGT2B17 | True | 6.5 | 0 | -3.75 | -8.7 | no RNA, as expected |

z-scores from only five controls are unstable (the control SD itself is estimated from five values) – discuss log2FC first.
P2's UGT2B17 change is weaker than expected (-0.37): the control mean includes C2, who also carries one copy, plus
sampling noise. UGT2B17 CPM by sample: P1 7.06 (CN 2) · P2 4.71 (CN 1) · P3 0 (CN 0) · C1 6.99 (CN 2) · C2 3.52 (CN 1) · C3 6.93 (CN 2) · C4 8.19 (CN 2) · C5 6.64 (CN 2).

Wrong-build demo (hg19 numbers intersected with hg38 MANE genes):
| cnv_id | correct (lifted) | wrong build | missed | wrongly added |
|---|---|---|---|---|
| P1_22q11_DEL | 43 | 46 | – | FAM246A, HIC2, RIMBP3B |
| P2_PMP22_DUP | 1 | 0 | PMP22 | – |
| P2_UGT2B17_DEL | 1 | 1 | UGT2B17 | UGT2B4 |
| P3_22q11_DEL | 24 | 17 | AIFM3, CRKL, LZTR1, P2RX6, PI4KA, SERPIND1, SNAP29, THAP7 | VPREB1 |
| P3_UGT2B17_DEL | 1 | 1 | UGT2B17 | UGT2B4 |


### Step 6 · TESSERACT-CNV input
| cnv_id | cnv_type | size_bp | in design range (50 bp-1 Mbp) | windows (500 bp) | segments (60 kb) | inference |
|---|---|---|---|---|---|---|
| P1_22q11_DEL | DEL | 2183224 | False | 4367 | 37 | segmented (score = max over segments) |
| P2_PMP22_DUP | DUP | 40000 | True | 80 | 1 | standard |
| P2_UGT2B17_DEL | DEL | 110000 | True | 220 | 2 | standard |
| P3_22q11_DEL | DEL | 1449894 | False | 2900 | 25 | segmented (score = max over segments) |
| P3_UGT2B17_DEL | DEL | 110000 | True | 220 | 2 | standard |

Segment template rows: 67 (seg_score column empty). Adapter message:
"TESSERACT-CNV not run: TESSERACT_CMD is not configured …".

### Step 7 · tracks (hg38)
| CNV | hg38 | ClinGen items (score, % of item inside the CNV) | ClinVar same-type RO≥50% (P/LP · VUS · benign) | P/LP SNVs inside | max pop. freq. (gnomAD SV AF / DGV Gold) | segdup fraction |
|---|---|---|---|---|---|---|
| P1_22q11_DEL | chr22:18,927,488-21,110,711 | HI 3: ISCA-37446 (22q11.2 recurrent (DGS) region (proximal, A-B, A-C, or A-D) (includes TBX1)) (covered 99.8%); HI 1: TBX1 (covered 100.0%); HI 2: ISCA-37516 (22q11.2 recurrent region (central, B-D or C-D) (includes CRKL)) (covered 99.9%); HI 1: CRKL (covered 100.0%); TS 3: ISCA-37446 (22q11.2 recurrent (DGS) region (proximal, A-B, A-C, or A-D) (includes TBX1)) (covered 99.8%); TS 0: TBX1 (covered 100.0%); TS 1: ISCA-37516 (22q11.2 recurrent region (central, B-D or C-D) (includes CRKL)) (covered 99.9%); TS 0: CRKL (covered 100.0%) | 111 · 5 · 0 | 1115 | 7.2e-05 / 0 | 0.129 |
| P2_PMP22_DUP | chr17:15,226,684-15,266,683 | HI 3: ISCA-37436 (17p12 recurrent (HNPP/CMT1A) region (includes PMP22)) (covered 3.0%); HI 3: PMP22 (covered 100.0%); TS 3: ISCA-37436 (17p12 recurrent (HNPP/CMT1A) region (includes PMP22)) (covered 3.0%); TS 0: PMP22 (covered 100.0%) | 5 · 0 · 0 | 91 | 0 / 0 | 0 |
| P2_UGT2B17_DEL | chr4:68,509,283-68,619,282 | none | 0 · 0 · 5 | 0 | 0.35 / 0.161 | 0.553 |
| P3_22q11_DEL | chr22:20,675,713-22,125,606 | HI 3: ISCA-37446 (22q11.2 recurrent (DGS) region (proximal, A-B, A-C, or A-D) (includes TBX1)) (covered 19.9%); HI 2: ISCA-37516 (22q11.2 recurrent region (central, B-D or C-D) (includes CRKL)) (covered 59.4%); HI 1: CRKL (covered 100.0%); HI 3: ISCA-37397 (22q11.2 recurrent region (distal type I, D-E or D-F)) (covered 53.2%); TS 3: ISCA-37446 (22q11.2 recurrent (DGS) region (proximal, A-B, A-C, or A-D) (includes TBX1)) (covered 19.9%); TS 1: ISCA-37516 (22q11.2 recurrent region (central, B-D or C-D) (includes CRKL)) (covered 59.4%); TS 0: CRKL (covered 100.0%); TS 1: ISCA-37397 (22q11.2 recurrent region (distal type I, D-E or D-F)) (covered 53.2%) | 4 · 1 · 0 | 904 | 0 / 0 | 0.355 |
| P3_UGT2B17_DEL | chr4:68,509,283-68,619,282 | none | 0 · 0 · 5 | 0 | 0.35 / 0.161 | 0.553 |

Note the PMP22 levels: gene HI 3, gene TS **0**, region ISCA-37436 (17p12 HNPP/CMT1A) HI 3 / TS 3. ClinGen's gene page
explains the TS 0: no duplication of PMP22 alone had been published when the gene was curated (2012).

### Step 8 · evidence directions (rule of thumb, not a classification)
| CNV | genes | dosage (ClinGen) | RNA | population | ClinVar |
|---|---|---|---|---|---|
| P1_22q11_DEL | towards pathogenic (≥ 25 genes) | towards pathogenic (score-3 HI item ≥ 95% inside) | supports a dosage effect (not pathogenicity) | neutral (rare or absent) | mostly towards pathogenic (some VUS/benign records) |
| P2_PMP22_DUP | neutral (gene number) | towards pathogenic (contains PMP22, the gene named in a score-3 TS region; but PMP22's own gene-level TS score is 0 - read the curation) | supports a dosage effect (not pathogenicity) | neutral (rare or absent) | towards pathogenic (similar P/LP records) |
| P2_UGT2B17_DEL | neutral (gene number) | neutral (no score-3 dosage item) | weak/unclear dosage effect | towards benign (common, ≥ 1%) | towards benign (similar benign records) |
| P3_22q11_DEL | neutral (gene number) | uncertain (partial overlap with a score-3 HI region) | supports a dosage effect (not pathogenicity) | neutral (rare or absent) | mostly towards pathogenic (some VUS/benign records) |
| P3_UGT2B17_DEL | neutral (gene number) | neutral (no score-3 dosage item) | supports the call (no RNA) | towards benign (common, ≥ 1%) | towards benign (similar benign records) |


## 3. Design notes and parameter choices

* **Reference.** UCSC `getData/sequence` (hg19) for chr22:18,000,001–23,000,000, chr17:14,600,001–15,700,000 and
  chr4:69,000,001–69,800,000, embedded in full-length N-masked chromosomes (real coordinates, works directly with the UCSC
  chain file). bgzip + `samtools faidx`.
* **Reads** (`scripts/simulate_dna.py`, seed `default_rng([20261015, sample_index])`): diploid haplotypes; deletions removed
  from haplotype 1 (both for CN 0); duplication = tandem copy on haplotype 1. Fragments N(350, 50) truncated to 200–600 bp,
  2 × 100 bp, GC acceptance exp(−½((g−0.45)/0.15)²) identical for all samples, binned qualities 2/12/23/37 with errors drawn
  at 10^(−Q/10) (Q2 → N), 3 % duplicate fragments, shuffled order. Base depth 3.0× × per-sample factor
  (P1 1.08, P2 0.94, P3 1.12, C1 1.00, C2 0.90, C3 1.15, C4 0.88, C5 1.05).
* **Why 3×?** Replicate test (`16` independent panels and P2 samples per depth; chr17 + chr4 windows; full pipeline:
  simulate → minimap2 → markdup → counting → caller; PMP22 duplication counted as recovered with ≥ 50 % overlap):
| base depth (x) | P2 depth (x) | P_SWITCH 1e-4 | P_SWITCH 1e-3 | P_SWITCH 1e-2 (used) |
|---|---|---|---|---|
| 1.5 | 1.41 | 3/16 | 8/16 | 12/16 |
| 2.0 | 1.88 | 4/16 | 8/16 | 15/16 |
| 2.5 | 2.35 | 10/16 | 16/16 | 16/16 |
| 3.0 | 2.82 | 13/16 | 15/16 | 16/16 |
| 3.5 | 3.29 | 16/16 | 16/16 | 16/16 |

  No false calls occurred in any setting (P2 outside the planted CNVs, and four leave-one-out controls per replicate).
  With the shipped caller (P_SWITCH = 0.01) recovery is reliable from 2.5× (16/16); 2.0× missed 1/16. We ship 3.0×, the
  lowest depth of the requested 3–4× range, which keeps a margin. With a 100-fold stricter switch penalty (1e-4) the
  caller would need 3.5×: sensitivity depends on depth **and** on caller settings – a good discussion point.
* **Panel of normals** (`scripts/build_pon.py`): C1 from the shipped FASTQ, C2–C5 simulated to a temporary folder,
  identical commands and counting code (`scripts/cnvlib.py`).
* **Counting:** one count per properly paired fragment at its midpoint (left read, TLEN > 0), MAPQ ≥ 20, no duplicates,
  secondary, supplementary or QC-fail reads.
* **Normalisation and masks:** per-million scaling; panel = median of five controls; masks = >10 % N, panel depth < 50 % of
  the typical bin, or robust control SD > 0.60 (≈ 3 × usual noise; robust so that C2's UGT2B17 deletion alone does not mask
  the locus); centring on the copy-neutral peak (smoothed-histogram mode, then median within ±0.3, three times).
* **HMM:** states CN 0–4, expected log2 −3 (clipped floor), −1, 0, +0.585, +1; σ = σ0·√(2/CN) (σ0 = robust SD of
  control-vs-control log2 ratios), σ(CN 0) = 0.6; P_SWITCH = 0.01; Viterbi per window on unmasked bins; segments ≥ 3 bins.
* **RNA** (`scripts/simulate_rna.py`, seed `default_rng([20261015, 99])`): 19,229 MANE Select genes; log10 CPM ~ N(1, 0.75);
  27.5 % not expressed; NB dispersion 0.02·e^N(0,0.4) (teaching genes PMP22, COMT, CRKL, TBX1, UGT2B17, HIRA, DGCR8, MAPK1,
  SNAP29, TXNRD2, GAPDH, ACTB fixed baselines and dispersion 0.01); fold (CN/2)^β with β ~ U(0.85, 1.0) or, for 15 % of
  genes, U(0.2, 0.5); partial duplications unchanged; chrY genes expressed only in samples marked M. Truth in
  `data/truth/rna_truth_cnv_genes.tsv`. UGT2B17 was given a low baseline (8 CPM) so that CN 0 is visible; in real skin
  fibroblasts UGT2B17 is barely expressed – say so when you show it.
* **Tracks:** UCSC hg38 `mane`, `clinGenHaplo`, `clinGenTriplo`, `clinGenGeneDisease`, `pliByGene` (gnomAD v2.1.1),
  `clinvarCnv`, `clinvarMain` (kept P/LP only), `phastConsElements100way`, `encodeCcreCombined`, `genomicSuperDups`,
  `gnomadStructuralVariants` (v4.1), `dgvGold`; names verified with `/list/tracks?genome=hg38`. DECIPHER skipped
  (protectedData). Regions = breakpoint-lifted planted CNV (P1: union with the ClinVar GRCh38 span) ± max(20 %, 50 kb).
  ClinGen dosage map version: 28 Apr 2026; MANE v1.5. Live refresh: `refresh_tracks_hg38()` in
  `scripts/fetch_public_data.py` (or `REFRESH_TRACKS = True` in Step 7).

## 4. Common errors and fixes

| Symptom | Cause | Fix |
|---|---|---|
| `minimap2: command not found` | environment not active | `conda activate d4-practice`, restart the kernel |
| `AssertionError: C1 differs from the panel` | mapping options or `MIN_MAPQ` changed, or Step 1 interrupted | re-run Step 1 completely, keep `MIN_MAPQ = 20` |
| `NameError` (e.g. `bins` not defined) | cells run out of order | *Kernel → Restart and run all* |
| Jupyter kernel never starts / `PermissionError` when binding a port | locked-down machine | `python scripts/run_notebook_inprocess.py practice2_cnv_tesseract_student.ipynb out.ipynb` (runs all cells in one in-process IPython shell) |
| CrossMap prints `[INFO]` lines | normal (stderr) | – |
| Different number of CrossMap pieces | other chain file or CrossMap version | use the shipped chain + CrossMap 0.7.3 |
| `Glyph ... missing from font SUIT` warnings | SUIT lacks the minus sign; fallback is DejaVu Sans | harmless |
| Windows: tools not found | native Windows | use WSL2 (README) |
| Track refresh fails | no internet / UCSC API down | keep `REFRESH_TRACKS = False` (pre-downloaded data are complete) |

## 5. Discussion – suggested answers

* **Step 7 · which tracks drive a score?** P1: the HI-3 22q11.2 DGS region (ISCA-37446, 99.8 % inside the call) and
  segments around TBX1 and CRKL; P3: partial overlap with the HI-3 distal region (53 %) and the central region (HI 2), CRKL
  inside. PMP22: one 40-kb segment with HI 3 / TS-3 region, Definitive gene–disease validity, 91 P/LP SNVs.
* **Circularity.** ClinGen dosage scores and ClinVar are model inputs and also ACMG/ClinGen evidence; do not add a high
  model score on top of the same ClinGen points.
* **What the model does not see.** Population SV frequency (excluded on purpose), inheritance, phenotype, RNA.
* **UGT2B17.** Whole-gene loss, yet gnomAD SV AF ≈ 0.35 for this deletion and DGV Gold losses up to 16 %; ClinVar benign;
  homozygous deletion in a healthy-looking person (P3's other CNV, and common in the population).
* **Breakpoints in LCRs / > 1 Mbp.** Exact breakpoints cannot be resolved by depth; gene content inside LCRs differs between
  builds; > 1 Mbp → segmented inference, maximum segment score.
* **Step 8.** P1 and P2: the dosage-curation line (with care for PMP22's gene-level TS 0) plus similar ClinVar P/LP records;
  P3: inheritance and phenotype are missing – the most informative next step is parental testing; P2's UGT2B17 deletion is
  reported (if at all) as a benign, common CNV, separately from the PMP22 duplication.

## 6. Caveats to mention in class

1. **Practice reference ≠ whole genome.** Segmental duplications whose partner copy is outside the three windows look unique
   here (e.g. around UGT2B17), so MAPQ is higher than in a real whole-genome analysis. The LCR22 copies are all inside the
   chr22 window and behave realistically.
2. **Mini-genome normalisation.** A 2.5-Mb deletion is 18 % of P1's analysed sequence, hence the centring step.
3. **Synthetic RNA.** Fibroblast-like baseline invented for teaching (not GTEx); named teaching genes have fixed baselines
   and lower noise.
4. **Track-to-group mapping** (Step 7) is our interpretation of the preprint's Table 1 groups.
5. **The executed notebook** shipped in the package was produced with `scripts/run_notebook_inprocess.py` (the build
   machine did not allow Jupyter kernels to open sockets); it runs exactly the student notebook's cells in order. A normal
   `jupyter nbconvert --execute` is expected to give the same numbers (every step is deterministic), but this could not be
   tested on the build machine – please run the student notebook once on a classroom laptop before the session.

## 7. "Your turn" and exit questions – short answers

(The full answers are in the notebook, inside the collapsible *Answer* blocks.)

| Where | Question (short) | Answer (short) |
|---|---|---|
| Step 0 | why full-length N-masked chromosomes? | real hg19 coordinates; chain file and browsers work directly; reads cannot map to N |
| Step 1 | which SAM columns hold chromosome, position, MAPQ? | 3, 4 (1-based), 5; FLAG in column 2; `RG:Z:` tag at the end |
| Step 2 | why the median of the panel, not the mean? | robust to one unusual control (C2's UGT2B17 deletion); controls are not CNV-free |
| Step 3 | why is the PMP22 call 40 kb instead of 35.6 kb? | breakpoints are rounded to 5-kb bin edges; exact junctions need split/discordant reads or long reads |
| Step 4 | why is an empty `.unmap` file not proof? | CrossMap drops parts in chain gaps silently; check pieces, sizes and gene content |
| Step 5 | why does HIRA change less than −1? | buffering/dosage compensation; noise of one sample vs five controls; tissue effects |
| Step 6 | windows and segments for P1's deletion? | 4,367 windows and 37 segments; the maximum is used because one dosage-sensitive segment is enough |
| Step 7 | what would help most for P3? | parental testing (inheritance), phenotype match, exact gene content (CRKL), population data |
| Step 8 | why not count TESSERACT and ClinGen HI separately? | ClinGen/ClinVar are model inputs (and shaped its labels): double counting |
| Exit 1 | why +0.58 and not +1? | 3/2 = 1.5, log2 1.5 ≈ 0.58; +1 would be 4 copies |
| Exit 2 | why were LCR22 bins masked? | MAPQ < 20 → almost no usable depth; breakpoints inside repeats cannot be placed by depth |
| Exit 3 | two liftover checks; why is ClinVar's interval longer? | breakpoint lift + gene content; NCBI remapping reports an outer span starting in LCR22A |
| Exit 4 | evidence TESSERACT does not see? | population SV frequency, inheritance, phenotype, RNA |
| Exit 5 | why is a missing RNA change not evidence against a deletion? | gene not expressed in the tissue, buffering, or noise – absence of signal is uninformative |
