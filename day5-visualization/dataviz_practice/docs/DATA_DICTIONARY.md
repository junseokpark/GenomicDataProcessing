# Data dictionary · dataviz_practice

Every file in `data/` and `data/truth/`: columns, units, coordinate conventions and how it was simulated. The short
table is in `README.md`. Scripts named `build/…` are the generators; they are in the instructor archive, and each
reruns byte-identically.

**Conventions.** hg38; intervals BED-style (0-based, half-open `start`, `end`); single positions 1-based; tab-separated;
missing values are empty fields; truth files are opened only in CHECK steps.

## Practice A and B1 · simple tables (`build/make_simple.py`, seed 20261027 + 1)

### `data/anscombe.tsv` · 44 × 3 · 468 B
* `dataset` (I–IV), `x`, `y`: Anscombe’s quartet, the published values (Anscombe 1973, *Am Stat* 27:17–21). All four sets
  share means, variances, the regression line and r = 0.82; only a plot shows that they differ.

### `data/mapping_rates.tsv` · 12 × 3 · 148 B
* `pipeline` (A, B, C), `library` (4 per pipeline), `pct_mapped` (% of reads that mapped). Means 96.1, 97.4 and 98.0%:
  the true ratio of C to A is 1.02, but bars on an axis from 95.5% make C look 4.2 times taller (A3).

### `data/expression_by_group.tsv` · 88 × 4 · 2 kB
* `dataset` (`small_n`: 4 per group; `large_n`: 40 per group), `group` (Control, Case), `sample`, `expression` (TPM of one
  gene, GENE_X). In `large_n` the Case group is bimodal: 13 of 40 cases are responders above 90 TPM
  (`truth/expression_by_group_truth.tsv`). Small study: medians 60.8 vs 89.8 TPM.

### `data/timecourse.tsv` · 54 × 4 · 873 B
* `gene` (FKBP5, TSC22D3, PER1), `time_h` (0, 1, 2, 4, 8, 24 h after dexamethasone), `replicate` (wells 1–3),
  `log2fc_vs_0h`. Shapes follow known glucocorticoid responses: PER1 transient, TSC22D3 fast and sustained, FKBP5
  slow and still rising at 24 h (3.66).

### `data/celltype_composition.tsv` · 32 × 4 · 629 B
* `donor` (D1–D4), `cell_type` (8 blood cell types), `n_cells`, `percent` (of that donor’s cells). The same 2,940 cells as
  B4 (`build/make_cell_labels.py`): CD14 Mono 19.0%, Memory CD4 T 16.3% overall.

### `data/methylation_cpg.tsv.gz` · 60,000 × 6 · 527 kB
* `chrom` (chr1–22), `pos` (1-based CpG position, hg38), `beta_wgbs` and `beta_nanopore` (methylated fraction, 0–1),
  `depth_wgbs` and `depth_nanopore` (reads covering the CpG; medians 22 and 28). Simulated per-CpG true β, mostly near
  0 or 1, read by both methods with sampling noise. Pearson r = 0.984, Spearman ρ = 0.897; the
  methods agree less at intermediate β (median |difference| 0.10 for β 0.3–0.7, 0.034 for β < 0.1 or > 0.9).

## B2 · bulk RNA-seq (`build/make_rnaseq.py`, seed 20261027 + 3)

12 libraries of airway smooth-muscle-like cells, 6 control and 6 dexamethasone-treated, sequenced in two batches
(A and B, 3 + 3 in each, so batch is balanced with treatment).

### `data/rnaseq_counts.tsv.gz` · 8,000 × 15 · 245 kB
* `gene` (HGNC symbol), `gene_id` (Ensembl), `chrom`, then 12 integer count columns `ctrl_A1` … `trt_B3`. Genes: 8,000 real
  protein-coding genes (MANE Select), in genomic order. Counts are negative binomial with a mean–dispersion trend;
  library sizes 11.5–26.1 million reads (size factors 0.65-1.47). The 15 known glucocorticoid targets of Himes et
  al. (2014) are among the induced genes (FKBP5, TSC22D3, ZBTB16, KLF15, PER1, DUSP1, CRISPLD2, KLF9, ERRFI1, CEBPD,
  SGK1, MT2A, ALOX15B, CDKN1C, PDK4); every other changed gene, including every down-regulated one, is random.
  Batch shifts some genes in both groups.

### `data/rnaseq_metadata.tsv` · 12 × 4 · 285 B
* `sample`, `condition` (control, treated), `batch` (A, B), `RIN` (RNA integrity number, 8.4–9.8).

### `data/rnaseq_log2norm.tsv.gz` · 6,940 × 13 · 273 kB
* `gene` + 12 samples: log2(normalised count + 1), median-of-ratios size factors (genes expressed in every sample), kept
  if the mean normalised count is ≥ 10: 6,940 genes. The same matrix as notebook Steps B2.2–B2.3 and prompt cards 3–4.

### `data/de_results.tsv` · 8,000 × 8 · 562 kB
* PyDESeq2, design `~ batch + condition`, treated vs control. `baseMean` (mean normalised count), `log2FC` (maximum
  likelihood), `log2FC_shrunk` (apeGLM shrinkage), `lfcSE`, `pvalue` (Wald test), `padj` (Benjamini–Hochberg).
  Empty `pvalue`/`padj` = not tested: 608 genes, 606 of them with zero counts in all 12 samples. With FDR < 5% and
  |log2FC_shrunk| > 1: 177 up, 106 down.

### `data/de_results_no_batch.tsv` · 8,000 × 7 · 498 kB
* The same analysis with design `~ condition` (no shrunken column): 367 genes at FDR < 5%, against
  396 with batch in the design.

### Truth files
* `truth/expression_by_group_truth.tsv` (40 × 2): `sample`, `responder` (1 = upper subgroup of the large-study cases).
* `truth/rnaseq_truth.tsv` (8,000 × 6): `gene`, `is_de` (1 = planted change), `true_log2fc`, `batch_log2_shift`,
  `baseline_mean`, `known_gr_target` (1 = one of the 15 known targets).


## B3 · GWAS


#### `data/gwas_sumstats.tsv.gz` (8,257,493 bytes)

Genome-wide summary statistics: one row per SNP, 299,185 rows, sorted by `CHR`, then `BP`. Tab-separated, gzip with mtime 0; floats with 6 significant digits; no missing values.

| column | type | unit / convention | meaning |
|---|---|---|---|
| SNP | text | `chr:pos:ref:alt`; pos 1-based, hg38 | SNP identifier (ref/alt bases are random) |
| CHR | integer | 1–22, no `chr` prefix | chromosome (autosomes only) |
| BP | integer | 1-based position, hg38 | position of the SNP |
| A1 | text | A/C/G/T | effect allele = **alt** allele; genotypes count copies of A1; BETA and OR refer to A1 |
| A2 | text | A/C/G/T | other allele = ref allele |
| MAF | float | fraction, 0.01–0.5 | minor-allele frequency in all 3,000 people |
| BETA | float | log odds ratio per copy of A1 | approximate log OR from the 4-PC trend test (definition below) |
| SE | float | same scale as BETA | standard error; BETA / SE = z and z² = χ² |
| P | float | probability | Armitage trend test, EIGENSTRAT form, **4 ancestry PCs** |
| P_no_ancestry_PCs | float | probability | the same test **without** PCs (K = 0) |

*Test.* With covariates C = [1, PC1…PCK], genotype g and case status y (1/0) are replaced by their residuals g̃, ỹ after least-squares regression on C; r = corr(g̃, ỹ), χ² = (N − K − 1)·r², 1 df. *Effect size.* BETA = (Σg̃ỹ / Σg̃²) / (φ(1 − φ)) with φ = case fraction = 0.5 (the linear-to-logistic approximation of the log odds ratio); SE = √(Σỹ² / (N − K − 1) / Σg̃²) / (φ(1 − φ)). exp(BETA) slightly underestimates large ORs (planted CFH OR 2.40; estimated 2.12).

#### `data/gwas_genotypes_subset.npz` (7,289,121 bytes)

`np.load` gives three arrays: `genotypes` (int8, 3,000 people × 10,097 SNPs; 0/1/2 = copies of A1), `iid` (row order = `gwas_samples.tsv`) and `snp` (column order = `gwas_snps_subset.tsv`). Columns: the 10,011 LD-thinned genome-wide SNPs used for the PCA plus every SNP within ±500 kb of the CFH lead SNP (90 SNPs; 4 belong to both sets).

#### `data/gwas_snps_subset.tsv` (439,959 bytes)

One row per column of the genotype matrix, sorted by CHR, BP. Columns `SNP, CHR, BP, A1, A2, MAF` as in the summary statistics; `pca_snp` (1 = LD-thinned SNP: the common SNP, MAF ≥ 0.05, closest to the centre of each 279.5-kb window) and `cfh_window` (1 = within ±500 kb of the CFH lead SNP).

#### `data/gwas_samples.tsv` (51,022 bytes)

One row per person (3,000), sorted by `IID`; row order is random with respect to ancestry and case status. `IID` (ID0001…ID3000), `case` (1 = case, 0 = control), `site` (S1/S2: recruitment site), `sex` (F/M), `age` (years at recruitment, integer; normal with mean 76 for cases and 71 for controls, SD 7, clipped to 55–95; age does not depend on genotype).

#### `data/gwas_genes_hg38.tsv` (592,594 bytes)

18,352 autosomal MANE Select protein-coding genes copied from `build/ref/genes_hg38_mane_select.tsv.gz`: `chrom` (chr1…chr22), `start0` (0-based start, BED), `end` (BED end = 1-based last base), `gene`, `strand`. Used for the nearest-gene labels and the gene track. **Convention:** a 1-based SNP position BP lies inside a gene when start0 < BP ≤ end.

#### `data/truth/gwas_truth.tsv` (448 bytes) — open only in CHECK cells

Planted disease loci: `locus` (literature label), `gene` (MANE gene whose midpoint holds the causal SNP), `CHR`, `BP` (1-based hg38), `SNP`, `A1` (risk-or-protective allele = alt), `A1_FREQ` (in the 3,000 people), `OR` (target per-copy odds ratio of A1, defined as the OR a single-SNP logistic regression within ancestry groups returns).

| locus | gene | CHR | BP | SNP | A1 | A1_FREQ | OR |
|---|---|---|---|---|---|---|---|
| CFH | CFH | 1 | 196,699,774 | 1:196699774:G:A | A | 0.365 | 2.4 |
| ARMS2/HTRA1 | HTRA1 | 10 | 122,488,230 | 10:122488230:C:G | G | 0.2348 | 2.6 |
| C2/CFB | CFB | 6 | 31,949,090 | 6:31949090:G:A | A | 0.0995 | 0.55 |
| C3 | C3 | 19 | 6,699,177 | 19:6699177:A:T | T | 0.2943 | 1.45 |
| APOE | APOE | 19 | 44,907,595 | 19:44907595:T:C | C | 0.1243 | 0.75 |
| LIPC | LIPC | 15 | 58,500,918 | 15:58500918:A:T | T | 0.3102 | 0.85 |
| CETP | CETP | 16 | 56,972,898 | 16:56972898:T:C | C | 0.3542 | 1.2 |
| TIMP3 | TIMP3 | 22 | 32,832,373 | 22:32832373:G:C | C | 0.34 | 1.3 |

#### `data/truth/gwas_ancestry.tsv` (27,019 bytes) — open only in CHECK cells

`IID`, `ancestry_group` (A: 1,800 people, B: 1,200).

#### `data/truth/gwas_truth.json` (6,757 bytes)

All design parameters (`design`) and realised values (`realized`: λ, loci with and without PCs, the ancestry-informative SNPs, PC variance, Hudson Fst). The notebook reads reference values from it in CHECK cells.

#### How the data were simulated (`build/make_gwas.py`)

* **People.** 1,800 people of genetic-ancestry group A and 1,200 of group B; cases are the 756 (A, 42%) and 744 (B, 62%) people with the highest disease liability in their group, so 1,500 cases and 1,500 controls. Site S2 holds 360 people of A and 900 of B (S1 the rest): site follows ancestry imperfectly; case fraction S1 45.3%, S2 56.5%.
* **SNP positions.** 312,000 autosomal positions, numbers per chromosome proportional to the hg38 length, uniform between 1 Mb before the first and 1 Mb after the last MANE gene of the chromosome (start clipped at 1; this leaves out the unassembled short arms of chromosomes 13, 14, 15, 21, 22). Random ref base, random different alt base. QC: SNPs with sample MAF < 0.01 removed (12,815), 299,185 remain.
* **Allele frequencies.** Ancestral alt frequency from a MAF spectrum of 80% Uniform(0.05, 0.5) + 20% log-uniform(0.01, 0.05), alt allele minor or major at random; with probability 0.5 a SNP copies the frequency of its left neighbour (± 0.02), like SNPs that tag the same haplotype. Group frequencies: Balding–Nichols, Beta(p(1−F)/F, (1−p)(1−F)/F) with F = 0.0028 (realised Hudson Fst between A and B = 0.002914).
* **LD.** One latent Gaussian AR(1) chain per haplotype along each chromosome: z_j = ρ_j z_(j−1) + √(1 − ρ_j²) e_j, ρ_j = exp(−d_j / 3 Mb − h_j), d_j = distance to the previous SNP, h_j = summed strength of the recombination hotspots between the two SNPs (one per 50 kb on average, strength ~ Exponential(mean 0.3)). Allele = alt when z_j < Φ⁻¹(p_group); genotype = sum of the two haplotypes (Hardy–Weinberg within groups). Around each planted SNP the 12 SNPs on either side copy its frequency, so the causal SNP sits in a block of proxies and hits form towers.
* **Disease model.** Liability = Σ planted effects + polygenic background (5,000 random SNPs > 2 Mb from the planted loci, total variance 2.0 on the logit scale) + logistic noise; the top-liability people of each group are the cases. Because the other risk factors attenuate single-SNP odds ratios, the planted (conditional) ORs were calibrated in a 1,000,000-person pilot so that the marginal OR equals the target: conditional OR CFH 3.857, ARMS2/HTRA1 4.405, C2/CFB 0.3749, C3 1.843, APOE 0.6216, LIPC 0.7609, CETP 1.349, TIMP3 1.549; pilot marginal OR CFH 2.399, ARMS2/HTRA1 2.598, C2/CFB 0.5504, C3 1.449, APOE 0.7502, LIPC 0.8503, CETP 1.2, TIMP3 1.299.
* **Ancestry-informative regions (no disease effect).** LCT (chr2:135,812,517): A1 frequency 0.2767 in A vs 0.7571 in B, P = 0.38 with PCs, 5.1 × 10⁻¹² without; HERC2 (chr15:28,216,610): A1 frequency 0.2867 in A vs 0.7704 in B, P = 0.96 with PCs, 9.4 × 10⁻¹⁰ without. They mimic lactase-persistence (LCT/MCM6) and eye-colour (HERC2/OCA2) alleles, the textbook stratification artefacts.
* **PCs and tests.** PCA: z-scored genotypes of the 10,011 thinned SNPs, scikit-learn PCA(10, randomized, random_state 0); PC1 0.3107%, PC2 0.07806%. Trend test for all SNPs with K = 0 and K = 4 PCs (λ 1.3285 and 1.0256).

**Real anchors.** hg38 chromosome lengths (`build/ref/hg38_chrom_sizes.tsv`); MANE Select gene coordinates for the causal-SNP positions (gene midpoints), the nearest-gene labels and the gene track; the choice of loci and approximate ORs/frequencies follows the AMD GWAS of Fritsche et al. (2016, *Nat Genet* 48:134–143; CFH, ARMS2/HTRA1, C2/CFB, C3, APOE, LIPC, CETP, TIMP3); LCT and HERC2 positions for the two ancestry-informative regions. Everything else is random.


## B4 · single cell


General: no file of this dataset has genomic interval or position columns. MANE Select coordinates (`build/ref/genes_hg38_mane_select.tsv.gz`, hg38, `start0` 0-based, BED style) are used only to put the genes in genomic order (chromosome 1–22, X; then `start0`), like a 10x `features.tsv`. Generator: `build/make_scrna.py`, seed 20261032 (= 20261027 + 5); a second run gives byte-identical files.

#### `data/scrna_counts.npz` · 863,109 bytes
* UMI counts, **2,940 cells (rows) × 2,000 genes (columns)**; `scipy.sparse` CSR matrix, dtype `int32`; 510,083 non-zero entries (91.3% zeros). Load with `scipy.sparse.load_npz`. Row *i* = row *i* of `scrna_cells.tsv`; column *j* = row *j* of `scrna_genes.tsv`.
* Unit: UMIs (unique molecules) per gene and cell. Median library size 2,188 UMIs (range 544–7,412); median 171 detected genes per cell (2,000-gene panel).
* How simulated: counts ~ Poisson(μ × Gamma(1/φ, φ)) (gamma-Poisson = negative binomial). μ = library size × cell-type expression profile × donor effect × cell state. Library size: log-normal, median 2,200 UMIs, sdlog 0.36, clipped to 550–7,600, times a cell-type factor (Naive CD4 T 0.9, Memory CD4 T 1.0, CD8 T 1.0, NK 1.0, B 1.05, CD14 Mono 1.35, FCGR3A Mono 1.3, DC 1.4). Overdispersion φ by gene class {'hk': 0.04, 'ribo': 0.04, 'marker': 0.3, 'background': 0.2} × log-normal(0, 0.3) per gene. Donor effect: log-normal sd 0.08 per gene and donor.
* Cell states: naive ↔ memory CD4 T and CD14 ↔ FCGR3A monocytes are continua (per-cell Beta score mixing the two profiles in log space; Beta parameters {'Naive CD4 T': [1.3, 10.0], 'Memory CD4 T': [10.0, 1.3], 'CD14 Mono': [1.3, 7.0], 'FCGR3A Mono': [7.0, 1.3]}); CD8 T cells carry an effector score (score s ~ Beta(2, 2); cytotoxic genes x exp(2.4 (s - 0.5)); naive genes and LTB, IL7R, CD27, CD28 x exp(-1.8 (s - 0.5))). NKG7 is exactly 0 in B cells and CD14 monocytes (needed for slide 21).
* Gene classes: 41 housekeeping, 40 ribosomal-protein, 202 marker/programme genes (levels in `data/truth/scrna_truth.json`), 1717 randomly drawn background MANE genes (log-normal level, sdlog 1.4, 1.5% of the UMIs; most are barely expressed in blood, which makes the matrix sparse like a full 10x matrix).
* Real anchors: real gene symbols and Ensembl IDs (MANE Select, hg38); canonical PBMC markers (Seurat PBMC 3k tutorial: IL7R, CCR7, S100A4, CD14, LYZ, MS4A1, CD8A, FCGR3A, MS4A7, GNLY, NKG7, FCER1A, CST3, plus LEF1, TCF7, SELL, GPR183, S100A8/9, VCAN, FCN1, CDKN1C, LST1, CD79A/B, BANK1, CD8B, GZMB, PRF1, KLRD1, CLEC10A) with levels set by hand to resemble 10x PBMC data (approximate UMIs per cell); housekeeping ACTB, B2M, GAPDH, ribosomal genes.

#### `data/scrna_genes.tsv` · 50,312 bytes
* 2,000 rows, one per matrix column, genomic order. Columns: `gene` (HGNC symbol, MANE Select), `gene_id` (Ensembl gene ID with version, as in the MANE file, e.g. `ENSG00000188290.13`). No coordinates.
* Selection: 243 hand-picked genes (41 housekeeping + 202 markers and cell-type programme genes, table in `build/make_scrna.py`), 40 of the 76 ribosomal-protein genes (RPL/RPS, not chrY) and 1717 background genes, both drawn at random (seed 20261032) from MANE Select (chrY excluded).

#### `data/scrna_cells.tsv` · 135,150 bytes
* 2,940 rows, one per matrix row, **sorted by `cell_id`** (not in truth-file order: students must join by barcode). Columns: `cell_id` (10x-style barcode, 16 nt + `-1`); `donor` (D1–D4, 4 donors); `n_counts` (UMIs per cell = row sum of the matrix); `n_genes` (genes with ≥ 1 UMI); `umap1`, `umap2` (fallback UMAP, unitless, 4 decimals; computed by the generator with the notebook's settings: genes in ≥ 3 cells, counts per 10,000 + log1p, UMAP of 30 PCs (1,000 HVGs) · n_neighbors 15 · min_dist 0.3 · seed 0 — identical to the notebook's UMAP on this machine). No cell types.

#### `data/scrna_markers.tsv` · 283 bytes
* 8 rows. Columns: `cell_type`; `markers` (comma-separated gene symbols). Prior knowledge from the literature (Seurat PBMC 3k tutorial and standard PBMC markers), **not** the truth: Naive CD4 T: CCR7, LEF1, TCF7, SELL, IL7R; Memory CD4 T: S100A4, GPR183, ANXA1, KLRB1, IL7R; CD8 T: CD8A, CD8B, CCL5, GZMK; NK: GNLY, NKG7, GZMB, PRF1, KLRD1; B: MS4A1, CD79A, CD79B, BANK1; CD14 Mono: CD14, LYZ, S100A8, S100A9, VCAN, FCN1; FCGR3A Mono: FCGR3A, MS4A7, CDKN1C, LST1; DC: FCER1A, CST3, CLEC10A.

#### `data/truth/scrna_truth.tsv` · 90,961 bytes
* Made by `build/make_cell_labels.py` (md5 `d133e1e248c1470a54d66c259fa19201`). Columns `cell_id`, `donor`, `cell_type`; order differs from `scrna_cells.tsv`. Open only in the Check beats.

#### `data/truth/scrna_truth.json` · 47,964 bytes
* Simulation parameters: seed, model, library-size model, dispersions, continua, the marker table (approximate UMIs per cell for each marker and cell type) and the ribosomal genes used. For instructors.


## B5 · copy number


All intervals are BED style (0-based start, end excluded) on **hg38**; add 1 to `start` for 1-based display.
Generator: `build/make_cnv.py`, seed `np.random.default_rng(20261027 + 6)`; running it twice gives byte-identical files.

#### `data/cnv_chr22_10kb.tsv`

Read depth of one patient and of a reference panel in 10-kb bins over all of chr22.  
**5,082 rows · 232,566 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `chrom` | text | – | chromosome (`chr22`) |
| `start` | integer | bp | bin start, 0-based |
| `end` | integer | bp | bin end, exclusive (last bin 50,810,000–50,818,468 is shorter) |
| `gc` | float | fraction 0–1 | GC fraction of the bin; empty where the reference is N (assembly gap) |
| `mappable` | True/False | – | False = do not use: short arm, centromere and satellites before 15.41 Mb, assembly gaps, telomere bin, bins ≥ 50% covered by segmental duplications with ≥ 98% identity |
| `depth_patient` | float | x (mean depth) | patient, constitutional DNA, ~30x WGS (reads × 150 bp ÷ bin length) |
| `depth_reference` | float | x (mean depth) | median of 10 simulated normals, each divided by its own median, rescaled to the panel's mean depth |

*Coordinates.* `start` 0-based, `end` exclusive, hg38. Bins are 10,000 bp from position 0.

*How it was made.* For every sample and bin: expected depth = mean depth × reference-genome factor (1 in unique sequence; 0.10–0.60 in near-identical repeats because the MAPQ filter removes most reads; log-normal(0, 1) with 15% dropouts in satellites; 0 in N) × GC bias exp(−½((gc − c)/w)²) with c ≈ 0.40, w ≈ 0.20 per sample × slow coverage waves (AR(1), SD 0.025 in log) × copy number (CN/2; in repeats 0.5 + 0.25·CN, diluted by paralogous reads) × gamma over-dispersion (CV 0.13) × extra log-normal noise (SD 0.4) in repeat and satellite bins; 0.3% outlier bins in the patient; reads ~ Poisson. GC = 0.48 + two AR(1) series, squeezed smoothly into (0.35, 0.65). The patient has CN 1 inside the planted deletion (bp-weighted at its two edge bins).

*Real anchors.* hg38 chr22 length 50,818,468 bp (`build/ref/hg38_chrom_sizes.tsv`); UCSC `genomicSuperDups` for chr22 17.98–22.41 Mb (downloaded for Day 4, `build/ref/d4_tracks_hg38/`) for the repeat mask and LCR22 blocks; the three assembly gaps inside LCR22A (holes of exactly 50 or 100 kb between segdup clusters with no phastCons element, see the B5 simplifications in INSTRUCTOR_NOTES.md); ClinGen ISCA-37446.

#### `data/cnv_genes_chr22.tsv`

Real MANE Select protein-coding genes overlapping chr22 17–23 Mb (gene track).  
**82 rows · 6,816 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `chrom` | text | – | chromosome |
| `start` | integer | bp | gene start, 0-based |
| `end` | integer | bp | gene end, exclusive |
| `gene` | text | – | HGNC symbol |
| `strand` | +/− | – | transcription strand |
| `gene_id` | text | – | Ensembl gene ID (versioned) |
| `transcript` | text | – | Ensembl MANE Select transcript |
| `refseq` | text | – | matching RefSeq transcript |

*Coordinates.* `start` 0-based, `end` exclusive, hg38.

*How it was made.* Not simulated: rows of `build/ref/genes_hg38_mane_select.tsv.gz` with chrom chr22, end > 17,000,000 and start < 23,000,000.

*Real anchors.* MANE Select (19,229 protein-coding genes, hg38).

#### `data/cnv_segdups_chr22.tsv`

Real segmental duplications on chr22 17.98–22.41 Mb (subset of UCSC genomicSuperDups).  
**590 rows · 36,532 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `chrom` | text | – | chromosome of this copy |
| `start` | integer | bp | copy start, 0-based |
| `end` | integer | bp | copy end, exclusive |
| `other_chrom` | text | – | chromosome of the paralogous copy (may be an unplaced contig) |
| `other_start` | integer | bp | start of the other copy, 0-based |
| `other_end` | integer | bp | end of the other copy |
| `strand` | +/− | – | + = same orientation (can cause deletions/duplications by NAHR), − = inverted |
| `frac_match` | float | fraction 0–1 | sequence identity of the two copies |
| `lcr22` | A/B/C/D or empty | – | LCR22 block that contains this segdup |

*Coordinates.* `start` 0-based, `end` exclusive, hg38.

*How it was made.* Not simulated. Rows of the Day-4 file for chr22, de-duplicated (the two Day-4 windows overlap, so 239 rows were duplicates). LCR22 blocks = segdups with ≥ 98% identity merged across holes < 20 kb or across an assembly gap; result: LCR22A chr22:18,159,724-18,939,750; LCR22B chr22:20,324,593-20,377,695; LCR22C chr22:20,678,348-20,738,049; LCR22D chr22:21,111,385-21,562,827.

*Real anchors.* UCSC hg38 genomicSuperDups (retrieved 2026-10-14, provenance JSON in `build/ref/d4_tracks_hg38/`).

#### `data/cnv_genome_1mb.tsv`

Genome-wide 1-Mb log2 ratios for chr1–22: the same patient and, for contrast, a colorectal tumour.  
**2,887 rows · 109,286 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `chrom` | text | – | chr1 … chr22 (autosomes only) |
| `start` | integer | bp | bin start, 0-based |
| `end` | integer | bp | bin end, exclusive (last bin of a chromosome is shorter) |
| `log2_patient` | float | log2 ratio | patient vs reference; empty = no unique sequence |
| `log2_tumour` | float | log2 ratio | tumour (purity 0.6) vs its own median; empty = no unique sequence |

*Coordinates.* `start` 0-based, `end` exclusive, hg38, 1,000,000-bp bins from position 0.

*How it was made.* Patient chr22 = mean of the 10-kb log2 ratios (median-normalised, ratio to the panel, masked) in each 1-Mb bin with ≥ 10 mappable 10-kb bins – identical to the notebook's `rebin()`. Patient chr1–21: AR(1) noise (ρ 0.6, SD 0.035) + white noise (SD 0.025). Tumour: log2((0.6·CN + 0.4·2)/2) with CN from the truth file, + AR(1) (SD 0.05) + white noise (SD 0.04), then median-centred. Empty bins: inside the gene-free centromere gap of each chromosome (and chr22 where the 10-kb data have no mappable bin).

*Real anchors.* hg38 chromosome lengths; centromere gaps = the gene-free stretch in MANE Select that contains an approximate hg38 centromere midpoint; arm boundary = middle of that gap (chr8 45.2, chr17 24.5, chr18 17.9, chr20 28.5 Mb). The colorectal pattern (gains 7, 8q, 13, 20q; losses 8p, 17p, 18q) is the classic chromosomal-instability pattern.

#### `data/truth/cnv_truth.tsv`

The planted events (students open it only in the CHECK step of B5.3).  
**8 rows · 1,209 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `sample` | text | – | `patient` or `tumour` |
| `event` | text | – | short name |
| `chrom` | text | – | chromosome |
| `start` | integer | bp | event start, 0-based |
| `end` | integer | bp | event end, exclusive |
| `type` | DEL/GAIN/LOSS | – | kind of event |
| `copy_number` | integer | copies | copy number in the affected cells |
| `purity` | float | fraction | share of cells carrying the event (1 = constitutional) |
| `expected_log2` | float | log2 ratio | log2((purity·CN + (1 − purity)·2)/2), before centring |
| `note` | text | – | how the event was placed |

*Coordinates.* `start` 0-based, `end` exclusive, hg38.

*How it was made.* Patient deletion chr22:18,717,065-21,342,616 (2.626 Mb), placed by non-allelic homologous recombination between two directly oriented copies of one 99.67%-identical segment (chr22:18,709,565-18,725,027 in LCR22A and chr22:21,335,351-21,350,332 in LCR22D), crossover at the same relative position (48.5%) of both copies. Tumour: 7 arm-level events, purity 0.6.

*Real anchors.* ClinGen ISCA-37446 (22q11.2 recurrent DGS region A–D, chr22:18,924,717-21,111,383) lies completely inside the planted deletion; UCSC segdup coordinates for the NAHR pair.



## B6 · variant callers


#### `data/variant_callsets.tsv` (320,392 bytes; 5,210 rows; tab-separated, sorted by `pos`)

A simulated GIAB-style benchmark of **one genome** in the window chr20:10,000,001-16,000,000 (hg38; BED
10,000,000–16,000,000). One row per variant that is in the truth set **or** was called by at
least one caller, inside simulated high-confidence regions (5,332,659 bp of the 6,000,000 bp window). Truth variants
missed by every caller are included (8 rows with four zeros).

| column | type / unit | meaning |
|---|---|---|
| `variant_id` | text | `chrom:pos:ref:alt`, e.g. `chr20:10002015:A:G` (unique) |
| `chrom` | text | always `chr20` (hg38) |
| `pos` | integer, **1-based** | VCF-style position; for indels the position of the anchor base (the base before the event) |
| `ref`, `alt` | text (A/C/G/T) | VCF-style alleles; indels carry the anchor base (`TA` → `T` is a 1-bp deletion). **Simulated: REF does not match the real hg38 base.** |
| `type` | `SNV` / `indel` | SNV = one base changed; indel = insertion or deletion of 1–12 bases |
| `context` | `unique`, `homopolymer_tr`, `segdup_lowmap`, `high_gc` | sequence context of the position (homopolymer/tandem repeat; segmental duplication/low mappability; high GC) |
| `GATK`, `DeepVariant`, `FreeBayes`, `Strelka2` | 0/1 | 1 = this caller reported the variant (stored as integers – convert with `== 1` or `.astype(bool)` before filtering) |
| `in_truth_set` | 0/1 | 1 = the variant is in the truth set (a real variant of this genome) |

**How it was simulated.** Every base of the window gets a context: 10 segmental-duplication
blocks (15–45 kb; 240,862 bp in HC), 70 high-GC blocks (2–7 kb;
268,985 bp), 5,000 short homopolymer/tandem-repeat intervals (12–60 bp;
157,769 bp), the rest unique (4,665,043 bp); 110 gaps are not
high-confidence. A variant's context is the context of its position. Exactly 4,800 truth variants (contexts
multinomial 80/10/5/5%; P(indel) 10.5% unique, 52% homopolymer/TR, 9% segdup, 8% high GC → 4,089
SNVs and 711 indels). Each caller finds a truth variant with a recall that depends on caller ×
context × type (table `design.recall` in the truth JSON); misses are correlated through a Gaussian copula (correlation
0.50 unique, 0.60 homopolymer/TR, 0.75 segdup, 0.70 high GC; DeepVariant loads 0.4 on the shared factor, the other
callers 1.0), so hard variants are hard for everyone and some real variants are found by one caller only. False calls
are planted per context with a fixed pattern of which callers make them (`design.fp_design`: 31 patterns,
410 false calls; by context unique 154, homopolymer/TR
123, segdup 92, high GC 41),
including shared false SNVs in segmental duplications that 3–4 callers agree on (mis-mapped reads from a paralogous copy).
Sites are ≥ 10 bp apart. Deterministic: two runs give byte-identical files.

**Real anchors.** hg38 coordinates and chromosome length of chr20 (`build/ref/hg38_chrom_sizes.tsv`); the window contains
real MANE Select genes (e.g. SNAP25, MKKS, JAG1, FLRT3; `build/ref/genes_hg38_mane_select.tsv.gz`); caller names and the GIAB
benchmarking logic (high-confidence regions, TP/FP/FN). Contexts, alleles and calls are simulated (no chr20 repeat or
segdup track is shipped; the Day 4 tracks cover chr22/chr17 only).

#### `data/truth/callsets_truth.json` (11,213 bytes)

Planted design (`design`: recall table, error correlation, caller loading, false-positive patterns with expected and
planted counts), interval statistics (`hc_bp`, `context_bp`, `n_blocks`) and the reference values used by the CHECK cells
(`per_caller` precision/recall/F1 for all/SNV/indel and recall by context, `upset_top10`, `all_four`, `single_caller`,
`deepvariant_only`, `consensus`, `shared_fp_by_context`, `missed_by_all_by_context`). Students open it only in CHECK cells.

#### Files written by the notebook

* `results/B6.json` – every CHECK value plus `slide_numbers`, `consensus_table`, `upset_counts`, `context_table`.
* `results/figures/B6_consensus.png` (1453 × 1059 px, 300 dpi)
* `results/figures/B6_contexts.png` (1973 × 902 px, 300 dpi)
* `results/figures/B6_precision_recall.png` (2224 × 900 px, 300 dpi)
* `results/figures/B6_upset.png` (1746 × 1097 px, 300 dpi)


## B7 · classifier


#### `data/classifier_scores.tsv` (660,712 bytes; 16,000 rows; tab-separated, sorted by `variant_id`)

| column | type / unit | meaning |
|---|---|---|
| `variant_id` | text | `mv00001` … `mv16000`, assigned in random order (carries no information; not an HGVS name) |
| `gene` | text | real MANE Select gene symbol (hg38), sampled uniformly (10,816 distinct genes; chrY, chrM and LOC symbols excluded) |
| `split` | `train` / `calibration` / `test` | 60/20/20%, stratified by label (9,600 / 3,200 / 3,200) |
| `label_pathogenic` | 0/1 | 1 = pathogenic; 560 of 16,000 = 3.5% in every split |
| `score_modelA` | 0–1, rank | model A: a single conservation-like score, reported as a percentile (a rank, **not a probability**) |
| `score_modelB` | 0–1 | model B: ensemble score, **raw and over-confident** (raw 0.6 means a true probability of 0.25) |

No genomic coordinates (variants are identified by `variant_id` and `gene`). No missing values.

**How it was simulated.** One latent "damage" value d per variant: benign d ~ N(0, 1), with
4.5% hard negatives d ~ N(2.3, 1); pathogenic d ~ N(4.2, 1),
with 5% hard positives d ~ N(0.8, 1). The true probability p_true(d) follows from
Bayes' rule at 3.5% prevalence (strictly increasing in d). Model B's raw score is mis-calibrated on purpose:
logit(p_true) = 2.0 × logit(score_B) − 1.9095; on the **train** split model B "memorised" the labels
(d shifted by ±1.0). Model A = percentile of 1.0 × d + 1.3 × noise.
Labels and splits are exact (560 pathogenic: 336 / 112 / 112). Deterministic (byte-identical over two runs).

**Real anchors.** MANE Select gene symbols; the 3.5% prevalence and the ROC/PR/calibration methods (Saito & Rehmsmeier 2015;
Niculescu-Mizil & Caruana 2005). Variants, labels and scores are simulated.

#### `data/truth/classifier_truth.tsv` (387,819 bytes) and `classifier_truth.json` (1,980 bytes)

TSV: `variant_id`, `p_true` (planted probability of being pathogenic, 0–1), `latent_d` (the latent damage value).
JSON: design parameters, the planted calibration map (`true = expit(K·logit(score_B) + C)` at raw 0.1 … 0.9), reference
metrics per split (n, prevalence, AUROC and AP of A and B), test-split Brier scores (raw, isotonic, true probability), mean
absolute error to p_true, and the reference cut-off for precision ≥ 0.90. CHECK cells only.

#### Files written by the notebook

* `results/B7.json` – every CHECK value plus `slide_numbers`, `per_split`, `calibration_raw`, `calibration_isotonic`.
* `results/figures/B7_calibration.png` (1055 × 1087 px, 300 dpi)
* `results/figures/B7_pr.png` (1031 × 1087 px, 300 dpi)
* `results/figures/B7_roc.png` (1039 × 1087 px, 300 dpi)
* `results/figures/B7_scores.png` (1449 × 948 px, 300 dpi)
* `results/figures/B7_threshold.png` (1216 × 1087 px, 300 dpi)


## B8 · TP53 lollipop


Amino-acid positions are **1-based**; domain intervals are **1-based and inclusive** (length = end − start + 1), the
convention of every protein database (a deliberate exception to the BED rule for genomic intervals).
Generator: `build/make_tp53.py`, seed `np.random.default_rng(20261027 + 9)`; byte-identical on re-run.

#### `data/tp53_variants.tsv`

One somatic TP53 variant per tumour (simulated cohort).  
**1,200 rows · 31,084 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `tumour_id` | text | – | T0001 … T1200 |
| `protein_change` | text | – | cBioPortal-style short notation: missense `R175H`, nonsense `R196*`, frameshift `P72fs`, splice `X125_splice`, in-frame deletion `V157del` |
| `position` | integer | amino acid (1-based) | residue of the change; for splice variants the codon next to the splice site |
| `variant_class` | text | – | `missense`, `truncating` (nonsense + frameshift) or `splice/other` (splice site + in-frame deletion) |

*Coordinates.* `position` 1-based on the 393-aa TP53 protein (UniProt P04637).

*How it was made.* Catalogue of variants with weights (expected tumours): 59 named missense alleles at real hotspots (R175H 62, R248Q 38, R273H 34, R248W 32, R282W 32, R273C 30, G245S 22, Y220C 22, R249S 18, V157F 12, G245D 10, …), 6 recurrent nonsense alleles (R213*, R196*, R306*, R342*, Q192*, E298*), 10 splice sites, and a long tail (one random reachable missense allele per residue: 274 expected tumours spread over the DBD, 91 outside; nonsense at every stop-reachable residue, 54; frameshifts at every residue, 132; 10 in-frame deletions). Class totals fixed to 876 / 264 / 60; tumours drawn within each class by a multinomial; order shuffled. Alleles are checked against the embedded sequence (reference residue) and the genetic code (alternative residue reachable by one nucleotide change).

*Real anchors.* TP53 protein sequence (UniProt P04637, 393 aa, embedded; the reference residue of all 65 named missense/nonsense alleles matches); real hotspot residues and recurrent nonsense positions.

#### `data/tp53_domains.tsv`

Four TP53 protein domains (domain boxes and the expectation of the binomial test).  
**4 rows · 450 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `domain` | text | – | Transactivation, Proline-rich, DNA-binding, Tetramerisation |
| `short` | text | – | TAD, PRD, DBD, TET |
| `aa_start` | integer | amino acid | first residue (1-based, inclusive) |
| `aa_end` | integer | amino acid | last residue (inclusive) |
| `length_aa` | integer | amino acids | aa_end − aa_start + 1 |
| `source` | text | – | where the boundaries come from |

*Coordinates.* 1-based, both ends included.

*How it was made.* Not simulated: TAD 1–61, PRD 64–92, DBD 94–292, TET 325–356.

*Real anchors.* IARC TP53 Database convention (now hosted as the NCI TP53 Database; Bouaoun et al., Hum Mutat 2016). UniProt and Pfam place some edges a few residues differently; the nuclear localisation signal and the C-terminal regulatory region are not drawn.

#### `data/truth/tp53_truth.tsv`

The simulation catalogue with the planted counts (open only in Step B8.3).  
**1,065 rows · 29,857 bytes · tab-separated · empty field = missing**

| column | type | unit | meaning |
|---|---|---|---|
| `protein_change` | text | – | allele |
| `position` | integer | amino acid | residue |
| `variant_class` | text | – | class |
| `expected_n` | float | tumours | expected count from the weights |
| `n_tumours` | integer | tumours | realised count in `tp53_variants.tsv` (0 = allele in the catalogue but not drawn) |

*Coordinates.* `position` 1-based.

*How it was made.* 1,065 catalogue rows; the realised counts sum to 1,200.

*Real anchors.* As for `tp53_variants.tsv`.


