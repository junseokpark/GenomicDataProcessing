# Instructor notes · D5 practice (Genomic Data Processing, 10/27)

For the instructor and teaching assistants. Students do not need this file. The answer key for Exercise 1 is
`solutions/ANSWER_KEY_bad_figures.md`; the executed notebooks are in `solutions/executed/`. You may hold back
`solutions/` until the end of the session.

## Learning goals

1. Recognise how a figure misleads (truncated axes, dynamite plots, rainbow colour maps, draw order, pies, red–green
   colour, distance claims on a UMAP, dual axes) and name the fix.
2. Choose the figure that matches the data and the question (Practice A).
3. For realistic data, run the method a real analysis would use, **check it against a known answer**, and only then
   draw the figure (Practice B: Question → Test → Check → Visualize).
4. Defend a figure: its numbers match the printed output, it reads in greyscale and for colour-vision deficiency,
   and it states n, the test, thresholds, units and the genome build.

## Timing (about 90 minutes, in pairs)

| Part | Min | What to watch | If short on time |
|---|---|---|---|
| Setup (Step 0) | 5 | `S0.files_missing = 0`; Colab needs `pip install umap-learn adjustText` | Pre-install in a shared image |
| Exercise 1 · M1–M8 | 8 | One problem and one fix per figure, with the checklist | Discuss M1, M4, M6 only |
| Practice A · A1–A5 | 20 | Students predict before running each cell | Skip A4; A5 is the overplotting example |
| B1 · two groups | 10 | n = 4: rank test vs t-test; error-bar types | Keep; it sets up the CHECK habit |
| B2 · bulk RNA-seq | 17 | Untested genes, batch in the design, both thresholds | Skip the optional DESeq2 re-run (B2.5b) |
| Elective (one per pair) | 20 | Each pair reports one CHECK number | B8 takes about 12 minutes |
| Defend · D1–D3 | 10 | D1 numbers match the CHECK cells; D2 colour-vision check | D3 (AI) is optional |
| Present | class | One sentence question, axes first, one number | Two or three pairs |

## Core reference values (all ✓ in `solutions/executed/practice_executed.ipynb`)

| Step | CHECK | Value |
|---|---|---|
| A1 | correlation of each Anscombe set | r = 0.82 in all four |
| A2 | medians, large study | Control 61.65, Case 68.85 TPM; 13 cases above 90 TPM |
| A3 | apparent vs true bar ratio | looks 4.2×, true 1.02× |
| A4 | composition | CD14 Mono 19.0%, Memory CD4 T 16.3%; FKBP5 at 24 h 3.66 log2 FC |
| A5 | WGBS vs nanopore | r = 0.984, ρ = 0.897 |
| B1 | n = 4 per group | medians 60.8 vs 89.8; Mann–Whitney p = 0.0286, Welch p = 0.107 |
| B1 | n = 40 per group | Mann–Whitney p = 8.0e-05, Welch p = 2.2e-05 |
| B1 | error bars | 95% CI is 3.18× the SEM bar (t for 3 df); Case SD 34.7 TPM |
| B2 | design | 12 samples, batch balanced; size factors 0.65-1.47 |
| B2 | PCA | PC1 38.9% (treatment), PC2 23.0% (batch) |
| B2 | hand-made test (6,940 genes) | 282 genes at FDR 5% with batch, 271 without; false discovery proportion 0.039 |
| B2 | DESeq2 table | 7,392 tested, 608 not tested; 396 at FDR 5% (5.8% not planted) |
| B2 | volcano (FDR 5% and \|log2FC_shrunk\| > 1) | 177 up, 106 down; 1 gene not planted; all 15 known targets up |

Two counts are easy to confuse. The hand-made test in Step B2.4 uses the 6,940 genes with mean normalised count ≥ 10, so
"5% false positives without adjustment" is 347 genes there; slide 35 uses the 7,392 genes tested in
`de_results.tsv`, so about 370. Both are right for their own table.

## Common mistakes and how to spot them

* **B1: p = 0.029 looks like evidence.** With 4 vs 4, the smallest two-sided Mann–Whitney p is 2/70 = 0.029: every
  Case value is above every Control value, and the test cannot say more. The Welch t-test gives 0.107 because one Case
  value is 150 TPM. Ask: what would p be if you had 40 per group? (Both tests below 10⁻⁴.)
* **B2: raw instead of shrunken fold change.** padj < 0.05 with raw log2FC gives 192 up and 119 down, not
  177 and 106; low-count genes with large raw fold changes get in (see the MA plot).
* **B2: adjusted p alone.** FDR < 5% without a fold-change threshold gives 242 up and 154 down, many with tiny effects.
* **B2: raw p.** 758 genes have p < 0.05 in `de_results.tsv`; if no gene changed at all, about 370 would pass by chance alone (7,392 × 0.05).
* **B2: untested genes.** Filling empty padj with 1 or 0 puts 608 genes on the volcano; they must be removed and counted.
* **Colour.** Red–green volcanoes (M6) collapse for deuteranopia; check with `dviz.show_cvd()`.
* **Joins.** Truth tables are not in the same row order as the data (B4 on purpose): always join on the ID column.

## Electives

Each elective has the same four beats. For each one, the notes below give a minute-by-minute plan, the expected
CHECK values, common mistakes and the answers to the *Your turn* tasks, then the simplifications that were made on
purpose in the simulation (say them if students ask whether the data are real).

All figures use SUIT, with DejaVu Sans as the fallback for glyphs that SUIT lacks (such as −, ≥, ×, β, λ, ρ and χ²).
`dviz.use_style()` sets this, so no "Glyph missing" warning appears.


## B3 · GWAS


#### At a glance

* **Run time.** The whole notebook runs in about 5 s on the test machine (printed as "Notebook time: 5 s"); the slowest cells are the Manhattan plot and the set-up imports. Expect a few seconds more on a student laptop.
* **Classroom time (≈ 20 min).** B3.1 question and stratification idea 3 min · B3.2 genotype PCA 4 min · B3.3 trend test with/without PCs, λ, QQ 5 min · B3.4 Manhattan, loci, truth check 5 min · B3.5 regional plot at CFH 3 min. Summary and presentation sentence fit into the B3.5 slot or the plenary.
* **Headline.** λ 1.33 without PCs → 1.03 with 4 PCs; 53 SNPs pass 5 × 10⁻⁸ in 5 loci (CFH p = 4.2 × 10⁻⁴⁶, OR 2.12, CFB p = 3.7 × 10⁻¹⁷, OR 0.49, HTRA1 p = 3.2 × 10⁻⁵⁰, OR 2.45, C3 p = 5.0 × 10⁻¹³, OR 1.49, TIMP3 p = 5.8 × 10⁻¹⁰, OR 1.39), all planted, none false. Without PCs 7 loci pass, 2 of them false.

#### Expected CHECK values (all must show ✓; `·` = recorded value for a Your-turn cell)

| CHECK | value | reference | |
|---|---|---|---|
| `B3.n_people` | 3,000 | 3000 | ✓ |
| `B3.n_cases` | 1,500 | 1500 | ✓ |
| `B3.n_controls` | 1,500 | 1500 | ✓ |
| `B3.n_pca_snps` | 10,011 | 10011 | ✓ |
| `B3.pc1_ancestry_agreement` | 1 | in [0.99, 1.0] | ✓ |
| `B3.pc1_var_pct` | 0.3107 | 0.3107 ± 0.001 | ✓ |
| `B3.pc1_over_pc2` | 3.981 | in [2, 100] | ✓ |
| `B3.yt_agreement_1000_snps` | 0.9177 | recorded | · |
| `B3.lambda_subset_noPC` | 1.323 | in [1.15, 1.5] | ✓ |
| `B3.lambda_subset_PC` | 1.011 | in [0.95, 1.1] | ✓ |
| `B3.max_abs_log10p_difference` | 2.113e-06 | in [0, 0.001] | ✓ |
| `B3.yt_lambda_site` | 1.179 | recorded | · |
| `B3.n_snps` | 299,185 | 299185 | ✓ |
| `B3.lambda_noPC` | 1.328 | in [1.2, 1.4] | ✓ |
| `B3.lambda_PC` | 1.026 | in [1.0, 1.08] | ✓ |
| `B3.n_snps_p_lt_5e-8` | 53 | 53 | ✓ |
| `B3.n_loci` | 5 | in [4, 6] | ✓ |
| `B3.false_loci` | 0 | 0 | ✓ |
| `B3.strong_loci_found` | 4 | 4 | ✓ |
| `B3.or_direction_agrees` | True | True | ✓ |
| `B3.yt_n_loci_noPC` | 7 | recorded | · |
| `B3.yt_false_loci_noPC` | 2 | recorded | · |
| `B3.cfh_lead_snp` | 1:196699774:G:A | 1:196699774:G:A | ✓ |
| `B3.cfh_window_snps` | 90 | 90 | ✓ |
| `B3.cfh_lead_r2_with_itself` | 1 | 1.0 ± 1e-09 | ✓ |
| `B3.cfh_causal_r2_with_lead` | 1 | in [0.5, 1.0] | ✓ |
| `B3.cfh_n_r2_ge_0.6` | 10 | recorded | · |
| `B3.yt_second_best_r2` | 0.8725 | recorded | · |

#### Your turn – answers

* **B3.1 (toy calculation).** Cases 0.350, controls 0.330: a difference of +0.020 for a SNP with no effect; over ~300,000 SNPs this lifts every statistic a little, and strongly differentiated SNPs can pass 5 × 10⁻⁸.
* **B3.2 (PCA with 1,000 SNPs).** Agreement with the ancestry labels drops from 1.000 to 0.918: each SNP carries little ancestry information (Fst ≈ 0.002914), so PCA needs thousands of SNPs.
* **B3.3 (site instead of PCs).** λ = 1.18, between 1.32 (no adjustment) and 1.01 (4 PCs): site is only a proxy for ancestry.
* **B3.4 (column without PCs).** 7 loci, 2 false: chr2:135,739,644 near UBXN4 (p = 2.4 × 10⁻¹⁴ without PCs, 0.10 with PCs); chr15:28,095,696 near OCA2 (p = 1.0 × 10⁻¹² without PCs, 0.12 with PCs) – the LCT-like and HERC2-like ancestry regions.
* **B3.5 (is the lead causal?).** Here the lead is the planted causal SNP, but the runner-up 1:196714024:G:A (r² = 0.87, 14 kb away) has p = 1.1 × 10⁻⁴⁰; with other noise it could have led. Fine-mapping (SuSiE, FINEMAP) assigns probabilities.

#### Three common student mistakes and how to spot them

1. **Using the wrong p-value column** (`P_no_ancestry_PCs` instead of `P`). Spot: λ = 1.33 instead of 1.03, `B3.lambda_PC` and `B3.false_loci` turn ✗, and two extra towers appear on chromosomes 2 and 15 (nearest genes UBXN4 and OCA2). The same mistake is what the assertion 0.95 < λ(P) < 1.1 on prompt card 5 catches.
2. **Mixing coordinate conventions between tables.** The GWAS table uses `CHR` = 1…22 and 1-based `BP`; the gene table uses `chr1` and 0-based BED starts. Comparing `genes.chrom == 1` gives an empty table (IndexError in `nearest_gene`), and forgetting the +1 mislabels a SNP that sits on a gene's first base. Spot: the CFH lead SNP must lie *inside* CFH (distance 0) and the gene track must not be empty.
3. **Running the PCA on the wrong SNPs or without standardising** (all subset columns including the dense CFH window, or raw 0/1/2 values). Spot: `B3.pc1_var_pct` turns ✗, and `B3.max_abs_log10p_difference` is no longer ≈ 0 because the PCs differ from the ones used for the summary statistics. Also watch for row-order errors between `gwas_samples.tsv` and the genotype rows – the `assert` in B3.2 stops the notebook.

Also seen: forgetting to sort by CHR and BP before merging hits into loci (n_loci ≠ 5); computing λ from p or −log10 p instead of χ² (use `dviz.genomic_lambda`); reading the triangle markers of a capped QQ plot as real values.

#### Talking points and caveats

* PC1 explains only 0.31% of the genotype variance but separates the groups perfectly: in genotype data most variation is between individuals; % variance is not a measure of how well ancestry is captured.
* λ with PCs is 1.03, not exactly 1: part is polygenic background (planted on purpose), part sampling noise. LD-score regression (Bulik-Sullivan et al., 2015) separates polygenicity from confounding.
* The 3 missed planted loci (APOE OR 0.75, LIPC OR 0.85, CETP OR 1.20) are a power lesson: 3,000 people detect only strong effects; the real AMD GWAS needed ~34,000 people for 34 loci.
* Nearest gene ≠ causal gene: the planted loci are labelled CFB and HTRA1 here (literature: C2/CFB, ARMS2/HTRA1); in gene-dense regions such as the MHC class III region on 6p21 several genes are equally close.
* The simulation's LD is a block-like latent-Gaussian model, not real haplotypes; r² patterns are realistic in shape but not in detail.


### Simplifications in the GWAS data

1. Allele-frequency differentiation: Balding–Nichols F = 0.0028 (realised Hudson Fst 0.002914) instead of ≈ 0.01. With N = 3,000 and case fractions 0.42 vs 0.62 the expected inflation is λ ≈ 1 + 110.6·F, so F = 0.01 would give λ ≈ 2.11, far outside the 1.2–1.4 target; the λ target and the exact case fractions were kept.
2. Two ancestry-informative regions without disease effect were added (LCT-like on chr2, HERC2-like on chr15; A1 frequency 0.2767 vs 0.7571 and 0.2867 vs 0.7704) so that stratification produces visible false loci without PCs (2 false of 7); they vanish with PCs.
3. Added `data/gwas_genes_hg38.tsv` (18,352 autosomal MANE Select genes copied from build/ref) so the notebook needs nothing outside data/ for nearest-gene labels and the gene track.
4. SNP count: 312,000 positions simulated and filtered at MAF ≥ 0.01, leaving 299,185 (about 300,000). Positions are drawn between 1 Mb before the first and 1 Mb after the last MANE gene of each chromosome; counts per chromosome stay proportional to the full hg38 length.
5. OR definition: truth `OR` is the target marginal per-allele OR within ancestry groups; the liability model uses larger conditional ORs (calibrated in a 1,000,000-person pilot) because the polygenic background and the other loci attenuate single-SNP ORs. Estimated lead ORs scatter around the targets and the linear-to-logistic BETA slightly underestimates large ORs.
6. LD model: the AR(1) correlation also drops at random recombination hotspots (block-like LD), and the 12 SNPs on either side of each planted/ancestry SNP copy its allele frequency (a tagging haplotype) so that towers form; plain AR(1) with random frequencies gave towers of 1 SNP.
7. Optional polygenic background included (5,000 SNPs, liability variance 2.0).
8. Colours: Manhattan hits are vermilion (as in the original deck image) because the chromosomes already use two alternating blues. Gene labels are upright: SUIT has no italic face.
9. Slide QQ plots are capped (image-37-1 at −log10 p = 20, practice_B3 at 15) with triangles for larger values, so that the λ difference stays visible; the cap and the maximum are stated in the legend.
10. practice_B3.png colours the PCA by case status only (no site markers) to keep the small panel readable (the notebook PCA uses colour = case and marker = site), and its QQ panel uses the genome-wide summary statistics (both P columns) rather than the notebook's subset.
11. 'λ on the subset' (notebook B3.3) is computed on the 10,011 LD-thinned SNPs only; the 90 CFH-window SNPs (one associated locus in strong LD) are left out of λ and of the subset QQ plot.
12. `dviz.qq_points` returns inf for p = 0; this does not happen here (smallest p 3.2 × 10⁻⁵⁰).

## B4 · single cell

> **Version note.** t-SNE and UMAP layouts depend on library versions. With Python 3.11, scikit-learn 1.4 and umap-learn 0.5.5 every CHECK still passed, but the 15-NN-kept shares moved (UMAP 15.8% → 15.2%, t-SNE 29.0% → 29.3%), as did the draw-order visibility gain (3.56 → 4.62). The slides quote the values from `results/B4.json` of the reference run.


#### At a glance
* Run time of the executed notebook: **33 s** in total (this build machine; UMAP 7.9 s including numba compilation, t-SNE 2.7 s). Without `umap-learn` the notebook uses the stored UMAP and varies the t-SNE perplexity in B4.6.
* Headline: k-means (k = 8) on 30 PCs recovers the planted types with ARI **0.962**; marker-based names are correct for **98.5%** of cells (44 wrong: naive ↔ memory CD4 T 26, CD14 ↔ FCGR3A Mono 12, CD8 T named as CD4 T 6). t-SNE/UMAP keep 12.4–29.0% of each cell's 15 nearest neighbours (chance 0.51%; PC1–PC2 6.34%), while the order of the 28 cluster-centre distances agrees with the 30-PC space only partly (ρ 0.36–0.67; PC1–PC2 0.96).

#### Minutes per step (suggested plan, about 20 min)
| Step | What students do | min |
|---|---|---|
| B4.1 | question; load; UMIs and genes per cell; QC histograms | 3 |
| B4.2 | counts per 10,000 + log1p; 1,000 HVGs; PCA (30 PCs), scree plot | 3 |
| B4.3 | UMAP (15 neighbours, min_dist 0.3, seed 0) and t-SNE (perplexity 30); settings and run time | 3 |
| B4.4 | k-means k = 8 on 30 PCs; open truth; ARI; truth × cluster heatmap | 3 |
| B4.5 | dot plot of 35 markers; automatic names; % named correctly; labelled UMAP | 3 |
| B4.6 | UMAP 5 / 200 neighbours, t-SNE, PC1–PC2: 15-NN kept vs centre-distance ρ | 4 |
| B4.7 | NKG7 draw order; summary; presentation sentence | 1 |

#### Expected CHECK values (all ✓ in the executed notebook)
| CHECK | value | reference | meaning |
|---|---|---|---|
| `B4.n_cells` | 2,940 | 2940 |  |
| `B4.n_genes` | 2,000 | 2000 |  |
| `B4.median_umis` | 2187.5 | in [1800, 2600] | UMIs per cell |
| `B4.umis_min` | 544 | in [300, 1000] | smallest library |
| `B4.umis_max` | 7,412 | in [5000, 10000] | largest library |
| `B4.median_genes` | 171.0 | in [100, 400] | of a 2,000-gene panel |
| `B4.pct_zeros` | 91.33 | in [90, 99] | % of matrix entries that are 0 |
| `B4.genes_kept` | 1,840 | in [1700, 2000] | genes seen in >= 3 cells |
| `B4.cells_below_200_genes` | 2,426 | in [1470, 2940] |  |
| `B4.cells_table_matches` | True | True |  |
| `B4.rowsum_max_error` | 1.82e-12 | in [0, 0.01] | largest \|row sum - 10,000\| after scaling |
| `B4.hvg_n` | 1,000 | 1000 |  |
| `B4.markers_in_hvg` | 31 | in [28, 35] | of 35 marker genes |
| `B4.pc1_3_var_pct` | 56.1 | in [30, 90] | % variance, PC1-PC3 |
| `B4.pc8_30_mean_var_pct` | 0.4 | in [0, 1] | % per PC: noise level |
| `B4.umap_ok` | True | True |  |
| `B4.tsne_ok` | True | True |  |
| `B4.umap_runtime_s` | 7.9 | in [0, 60] | seconds (first call also compiles code) |
| `B4.tsne_runtime_s` | 2.7 | in [0, 60] | seconds |
| `B4.ari` | 0.962 | in [0.9, 1.0] | adjusted Rand index, k-means (k = 8) vs truth |
| `B4.labels_unique` | 8 | 8 | 8 clusters -> 8 different names |
| `B4.pct_correct` | 98.5 | in [95, 100] | % of cells named correctly |
| `B4.n_wrong` | 44 | in [0, 147] | cells with a wrong name (5% = 147) |
| `B4.pct_cd14_mono` | 19.0 | 19.0 ± 0.05 | % of cells |
| `B4.pct_memory_cd4` | 16.3 | 16.3 ± 0.05 | % of cells |
| `B4.knn_tsne30_pct` | 29.0 | in [5.1, 100] | t-SNE, perplexity 30: >= 10 x chance |
| `B4.knn_umap15_pct` | 15.8 | in [5.1, 100] | UMAP, 15 neighbours: >= 10 x chance |
| `B4.knn_umap5_pct` | 16.8 | in [5.1, 100] | UMAP, 5 neighbours: >= 10 x chance |
| `B4.knn_umap200_pct` | 12.4 | in [5.1, 100] | UMAP, 200 neighbours: >= 10 x chance |
| `B4.knn_pca2_pct` | 6.3 | in [0, 12.4] | a linear 2-D view keeps fewer |
| `B4.same_label_min_pct` | 98.0 | in [95, 100] | t-SNE / UMAP: neighbours have the same name |
| `B4.rho_range` | 0.6 | in [0.3, 2] | rho from 0.36 to 0.96 |
| `B4.n_isolated_answers` | 3 | in [2, 6] | different answers to 'which cluster is most isolated?' |
| `B4.nkg7_umis_in_B_and_CD14` | 0 | 0 | all-zero groups |
| `B4.nkg7_pct_nk` | 99.6 | in [90, 100] | % of NK cells with NKG7 |
| `B4.visible_gain` | 3.56 | in [1.5, 100] | sorted / zeros on top |

#### Three common mistakes and how to spot them
1. **Truth aligned by row order instead of joined by barcode.** `scrna_cells.tsv` is sorted by `cell_id`, the truth file is not. Aligning rows gives ARI **0.001** (correct join: 0.962). Spot it: CHECK `B4.ari` ✗ and a flat truth × cluster table. Fix: `cells[['cell_id']].merge(truth, on='cell_id', validate='one_to_one')`.
2. **Skipping log1p** (PCA on counts per 10,000). A few very high genes (LYZ, S100A8/9, ribosomal genes) dominate the PCs and k-means gives ARI **0.476**. Spot it: CHECK `B4.ari` ✗, monocyte genes dominate every PC, the dot plot names are mixed up.
3. **Copying QC thresholds or reading the UMAP like a map.** The tutorial rule "≥ 200 genes per cell" would remove **2,426** of 2,940 cells, because this is a 2,000-gene panel (median 171 genes). And on the default UMAP the most isolated island is **B**, whereas in the 30-PC space it is **NK**. Spot it: students who filter cells lose most of the data (downstream CHECKs fail); students who argue from gaps cannot explain the B4.6 table, where the "most isolated cluster" changes with every setting.

#### Answers to the Your-turn tasks
* **B4.1** (≥ 200-gene rule): removes 2,426 of 2,940 cells; the rule is for full-transcriptome matrices (background: 20,000–30,000 genes), this panel shows a median of 171 genes per cell. Derive thresholds from your own distributions.
* **B4.2** (scree flattens after ~PC 7): 8 cell types span at most 7 independent directions; PC 8–30 explain 0.4% each on average (noise). 30 PCs are safe.
* **B4.3** (UMAP seed 1): islands move, rotate or mirror; only 27.8% of each cell's 15 nearest neighbours are the same on both maps, but 98.7% of a cell's neighbours stay in its own cluster (seed 0: 98.8%) and 15-NN kept is 15.6% vs 15.8%. k-means on the PCs does not change.
* **B4.4** (k-means on the UMAP): ARI 0.954 (30 PCs: 0.962). Similar here, but the map depends on n_neighbors, min_dist and seed and distorts distances, so clusters drawn on it inherit arbitrary choices. (Extra: scaling every gene before PCA gives ARI 0.945.)
* **B4.5** (naive vs memory markers): CCR7 in 87% of naive-cluster vs 53% of memory-cluster cells; S100A4 96% vs 67%; also LEF1, TCF7, SELL (naive) and ANXA1, KLRB1, GPR183 (memory); IL7R in both.
* **B4.6** ("B cells most different"): no. Most isolated cluster per space: 30 PCs (reference): NK; PCA (PC1, PC2): NK; t-SNE, perplexity 30: CD8 T; UMAP, 15 neighbours: B; UMAP, 5 neighbours: NK; UMAP, 200 neighbours: NK.
* **B4.7** (Seurat default order and violins): data order shows NKG7 in 27.9% of occupied squares, between zeros-on-top (14.3%) and sorted (51.0%). NKG7 is in 0.0% of B cells and 0.0% of CD14 monocytes: a violin of these groups is a flat line; VlnPlot's added noise draws an invented shape (slide 21).


### Simplifications in the single cell data

* Simulation strength: the naive vs memory CD4 T difference is stronger than in real PBMC data (the expression table gives the naive state 23 and the memory state 29 programme genes, including extra real ones such as NELL2, BACH2, LGALS1 and KLF6, with larger fold changes) so that k-means on 30 PCs separates them (CHECK `B4.ari`: in [0.9, 1.0]). CD8 T and NK do not touch on the UMAP (background: in many real PBMC UMAPs they do).
* Panel matrix: 2,000 genes with a median of 2,188 UMIs per cell and 91.3% zeros means the UMIs concentrate in ribosomal and housekeeping genes; detected genes per cell (median 171) and log-normalised values are not comparable with a full-transcriptome 10x matrix (background: several hundred to a few thousand detected genes per cell, depending on the chemistry).
* PCA is centred but not scaled to unit variance (background: Seurat's ScaleData scales by default); scaling gave ARI 0.945 vs 0.962 on this sparse panel. run_pca(scale=True) is offered.
* HVGs: Scanpy 'seurat'-flavour dispersion z-score on counts per 10,000, but with 20 equal-size bins (pd.qcut) instead of equal-width bins (background: Scanpy uses pd.cut).
* Colours: cell types use Okabe-Ito colours with OI yellow replaced by dark gold #8C6D1F (FCGR3A Mono) for visibility on white; DC = OI black #000000. Sequential maps are viridis (the course default) — slides 17 and 42 originally used grey→blue; slide 42 gains a colour bar that the original lacked.
* Slide 21: the right panel adds Gaussian noise with sd 0.03 (log units) so the invented violin shape is visible at slide size. Seurat's VlnPlot adds far smaller noise (background: rnorm()/100,000 in Seurat's source, not run here) — same effect, much thinner shape.
* Slide 19: the bad pie keeps all 8 slices, shadow, explosion and legend but uses the colour-blind-safe cell-type colours (not matplotlib's default tab10); the sorted bars highlight the two compared types (CD14 Mono, Memory CD4 T) in blue instead of the largest type.
* Slides 18, 41, 42 and practice use the marker-based names of the k-means clusters (98.5% equal to the truth); slides 19 and 21 use the truth file (exact composition; exactly-zero NKG7 groups). Slide 42 shows 14 genes: the original slide's 12 (read from the original image) plus S100A4 and MS4A7; 13 of the 14 are among the 35 genes of data/scrna_markers.tsv, CD3D (pan-T-cell gene) is not.
* Labels: adjustText runs with a fixed iteration limit (iter_lim=25) because without a limit adjustText 1.4.0 stops after 1 s (time_lim), which makes label positions depend on machine speed; slides 41 and practice (UMAP) and the PCA panel of slide 18 use hand-placed labels with leader lines.
* Log axes in the notebook use plain number labels, because matplotlib's math-text tick labels do not use the DejaVu Sans fallback.
* Supporting files: build/scrna_helpers.py (helpers copied into the notebook with inspect.getsource), build/answers_B4.py (writes results/B4_yourturn.json, the numbers behind the folded answers and the mistakes section) and data/truth/scrna_truth.json (simulation parameters).
* The normalisation CHECK reports the largest |row sum − 10,000| rather than the row sums, because `dviz.check` prints floats with four significant digits (10,000.0 would show as 1e+04).
* MANE Select lacks LILRA3 and OAZ1, which were therefore not used as anchors.
* Run times are machine-dependent and not deterministic (UMAP includes numba compilation on the first call); all other values in results/B4.json are identical across executions.

## B5 · copy number


#### Timing (about 20 min)

| step | beat | minutes | what students do |
|---|---|---|---|
| B5.1 | Question | 3 | load the bins; worked example TBX1 vs MAPK1 bin |
| B5.2 | Test | 5 | median-normalise, divide by the panel, log2, mask; GC picture |
| B5.3 | Test + Check | 5 | segment; deletion coordinates, genes; compare with the truth |
| B5.4 | Visualize | 4 | four tracks on one x-axis |
| B5.5 | Visualize + Check | 3 | genome-wide 1-Mb view, patient vs tumour; Your turn re-binning |

Run time of the whole notebook with `build/run_notebook.py`: 1.0 s.

#### Expected CHECK values (all pass)

| CHECK | value | reference | pass |
|---|---|---|---|
| `B5.n_bins` | 5,082 | 5082 | ✓ |
| `B5.chr22_end` | 50,818,468 | 50818468 | ✓ |
| `B5.n_mappable` | 3,410 | in [3000, 4000] | ✓ |
| `B5.first_mappable_Mb` | 15.41 | in [15.0, 16.0] | ✓ |
| `B5.median_log2_outside` | 0.04948 | in [-0.1, 0.1] | ✓ |
| `B5.median_log2_core` | −0.9615 | in [-1.1, -0.9] | ✓ |
| `B5.sd_log2` | 0.2187 | in [0.15, 0.25] | ✓ |
| `B5.gc_effect_raw` | 0.2628 | in [0.1, 0.5] | ✓ |
| `B5.gc_effect_ratio` | 0.03945 | in [0.0, 0.08] | ✓ |
| `B5.start_in_LCR22A` | True | True | ✓ |
| `B5.end_in_LCR22D` | True | True | ✓ |
| `B5.truth_in_intervals` | True | True | ✓ |
| `B5.mean_log2_del` | −0.9647 | in [-1.1, -0.9] | ✓ |
| `B5.copy_number` | 1.025 | 1.0 ± 0.1 | ✓ |
| `B5.n_segments` | 3 | 3 | ✓ |
| `B5.key_genes_inside` | 5 | 5 | ✓ |
| `B5.MAPK1_outside` | True | True | ✓ |
| `B5.bins_1mb_low` | 2 | 2 | ✓ |
| `B5.bins_1mb_overlap` | 4 | 4 | ✓ |
| `B5.tumour_8q_log2` | 0.3821 | 0.379 ± 0.08 | ✓ |
| `B5.rebin_equals_1mb_file` | True | True | ✓ |

Key printed numbers: called deletion chr22:18,940,001-21,110,000 (BED
`chr22 18940000 21110000`), 2.17 Mb,
208 mappable bins, median log2 inside −0.962,
42 MANE genes; planted chr22:18,717,065-21,342,616
(2.626 Mb); uncertainty intervals 780 kb (start, masked LCR22A) and
330 kb (end, masked LCR22D). Genes in the planted deletion but outside the call: DGCR6, PRODH
(printed in the B5.3 CHECK cell).

#### Three common mistakes and how to spot them

1. **Using all bins instead of mappable bins** (normalising by the median or mean of every bin, or segmenting
   without the mask). Zeros from N gaps and satellite pile-ups (hundreds of x before 15.4 Mb) enter the
   calculation. *Spot it:* `B5.median_log2_outside` moves away from 0, `B5.n_segments` is larger than 3, or log2 values
   of ±3 appear before 15.4 Mb. *Fix:* `.where(ok)` and `cn.loc[ok, ...].median()`.
2. **Off-by-one and end-of-bin errors.** Printing the BED start as a 1-based position (18,940,000 instead of
   18,940,001), or taking the end from the *start* of the last bin (10 kb too short). *Spot it:* 1-based starts of
   bin-based calls always end in …,001; the length must be 2.17 Mb, not 10 kb less.
3. **Expecting base-exact breakpoints.** Students think the call is wrong because it is 2.17 Mb while
   the truth is 2.626 Mb. *Spot it:* they compare `t["start"]` with `d_start` and stop.
   *Explain:* both breakpoints lie in masked LCR22 sequence; read depth gives an interval (780 and
   330 kb), and `B5.truth_in_intervals` is ✓. Exact junctions need long reads, optical mapping or
   breakpoint PCR.

#### Your-turn answers (short; full answers are folded in the notebook)

| step | task | answer |
|---|---|---|
| B5.1 | why do masked rows jump around? | LCR22 reads map to several copies → low MAPQ, filtered; depth low and erratic in every sample and depends on each person's LCR22 structure; from 18.94 Mb the sequence is unique and the patient is at half depth |
| B5.2 | log2 for 3 copies and 0 copies? | log2(3/2) ≈ +0.58 (not +1); 0 copies → no reads, log2(0) = −∞ (clipped at −3 in plots) |
| B5.3 | is the mean consistent with one copy? | SE = SD/√n = 0.219/√208 ≈ 0.0152; relative to the copy-neutral level (+0.049) the deletion sits at −1.014 = one copy lost |
| B5.4 | why do all tracks zoom? | `sharex=True` shares one x-axis; *TBX1* (19.761 Mb), *COMT* (19.942 Mb), *DGCR8* (20.08 Mb) at half depth |
| B5.5 | re-bin to 100 kb and 1 Mb | 10 kb: 208 bins fully inside, SD 0.219; 100 kb: 21 bins fully inside, SD 0.075; 1 Mb: 2 bins fully inside, SD 0.04 – the SD falls by less than √(bins averaged) because of correlated coverage waves |


### Simplifications in the copy number data

1. **Repeat mask only where real segdup data exist** (chr22 17.98–22.41 Mb, downloaded for Day 4). Elsewhere on chr22 only
   the short arm/centromere (< 15.41 Mb) and the telomere bin are masked; other chr22 repeats (e.g. LCR22E–H, IGL) are
   simulated as unique sequence.
2. **Assembly gaps inside LCR22A are inferred**, not read from the UCSC gap table (it is not among the downloaded tracks): holes of exactly
   50,000 or 100,000 bp between segdup clusters with no phastCons element (also no cCRE, gnomAD-SV, ClinVar or MANE
   feature). The 10.51-Mb short-arm gap end and the 15.41-Mb start of unique sequence are approximate hg38 landmarks.
3. **Planted deletion 2.626 Mb**, slightly above the intended 2.5–2.6 Mb, because it was placed by NAHR between a
   real pair of directly oriented 99.67%-identical LCR22A/LCR22D copies (direct-orientation LCR22A/LCR22D pairs with
   ≥ 99% identity imply deletions of 2.436 Mb (99.34%), 2.497 Mb (99.42%), 2.626 Mb (99.67%), 2.939 Mb (99.47%), 3.123 Mb (99.07%); none lies strictly
   inside 2.5–2.6 Mb, and we used the most similar pair, 2.626 Mb at 99.67% identity).
   It contains ISCA-37446 completely.
4. **The called deletion (2.17 Mb) is shorter than the planted one**, as in real short-read data: both breakpoints lie in
   masked LCR22 sequence. Slides report the called interval and the planted interval separately.
5. **The normal level is +0.049, not 0**, because the notebook normalises by the median of chr22 only and 6% of the mappable bins are deleted; the notebook explains it instead of re-centring.
6. **1-Mb file:** autosomes only; centromere gaps from the gene-free MANE stretch around approximate (rounded, typed)
   hg38 centromere midpoints; patient chr1–21 and the whole tumour are simulated directly at 1 Mb.
7. The notebook gives the same `results/B5.json` whether it is run from `electives/` or from the package root.

## B6 · variant callers


**Timing (18 min):** B6.1 yardstick 3 · B6.2 callers vs truth 4 · B6.3 UpSet 4 · B6.4 consensus 4 · B6.5 contexts 3.
Whole-notebook run time on the build machine: 0.9 s with `build/run_notebook.py`
(in-process runner, 18 code cells); expect well under a minute on a student laptop.

**Expected CHECK values** (all ✓ in the executed notebook):

| CHECK | value | reference | ok |
|---|---|---|---|
| `B6.toy_precision` | 0.833 | 0.833 ± 0.0005 | ✓ |
| `B6.toy_recall` | 0.714 | 0.714 ± 0.0005 | ✓ |
| `B6.toy_F1` | 0.769 | 0.769 ± 0.0005 | ✓ |
| `B6.n_truth` | 4,800 | 4800 | ✓ |
| `B6.max_diff_vs_reference` | 4.94253e-07 | 0 ± 1e-05 | ✓ |
| `B6.highest_precision` | DeepVariant | DeepVariant | ✓ |
| `B6.lowest_precision` | FreeBayes | FreeBayes | ✓ |
| `B6.indel_recall_lower_for_every_caller` | True | True | ✓ |
| `B6.largest_bar` | all four | all four | ✓ |
| `B6.n_all_four` | 4,399 | 4399 | ✓ |
| `B6.n_deepvariant_only` | 39 | 39 | ✓ |
| `B6.bars_add_up_to_all_calls` | 5,202 | 5202 | ✓ |
| `B6.precision_all_four` | 0.9927 | 0.9927 ± 0.0001 | ✓ |
| `B6.false_calls_all_four` | 32 | 32 | ✓ |
| `B6.precision_single_caller` | 0.11 | 0.11 ± 0.0005 | ✓ |
| `B6.all_four_more_precise_than_single` | True | True | ✓ |
| `B6.all_four_recall_below_DeepVariant` | True | True | ✓ |
| `B6.shared_fp_top_context` | segdup_lowmap | segdup_lowmap | ✓ |
| `B6.n_shared_false_calls` | 52 | 52 | ✓ |
| `B6.segdup_pct_of_shared_false` | 76.9 | in [50, 100] | ✓ |
| `B6.segdup_pct_of_truth` | 5 | in [3, 7] | ✓ |
| `B6.segdup_false_in_all_four` | 24 | 24 | ✓ |

**Three common mistakes and how to spot them**

1. **Filtering with 0/1 integers.** `vc.loc[vc["GATK"], "in_truth_set"]` (as on demo slide 59) treats 0/1 as *row labels*
   and returns rows 0 and 1, so every caller gets the same meaningless "precision" (measured on this file:
   1.0 for all four callers). Spot: four identical precisions, or
   `B6.max_diff_vs_reference` ✗. Fix: `vc[vc["GATK"] == 1]` or `.astype(bool)` (the `prf()` helper does this).
2. **Recall without the missed variants.** Dividing by the number of truth variants *that somebody called*, or dropping the
   8 rows with four zeros, inflates recall (the union then has recall 1.000 instead of
   0.998). Spot: `B6.max_diff_vs_reference` ✗; recall of "≥ 1" equals 1.
3. **Inclusive vs exclusive intersections.** Reading "GATK + DeepVariant" as all variants called by both (thousands) instead
   of the exclusive UpSet bar (38), or reading the set size
   "DeepVariant (4,787)" as the DeepVariant-only bar (39). Spot:
   `B6.bars_add_up_to_all_calls` – exclusive bars add up to 5,202.

**Answers to the Your-turn tasks** (full text in the folded answers of the notebook)

* B6.1 – true negatives are almost every base of the genome; accuracy would be > 99.9% even for a caller that reports
  nothing, so precision/recall/F1 are used.
* B6.2 – FreeBayes indels: 681 calls, 88 false,
  precision 0.871 (mostly homopolymer/tandem-repeat indels).
* B6.3 – FreeBayes-only calls: 2 of 169
  real; DeepVariant-only: 32 of 39 real. Agreement counts do not show correctness.
* B6.4 – clinical: "= 4" / "≥ 3" (precision 0.993 / 0.989) plus
  orthogonal confirmation, because shared errors survive; screen: "≥ 2" (recall 0.990) or "≥ 1"
  (410 false calls). DeepVariant alone (F1 0.989) beats every voting rule.
* B6.5 – false calls made by one caller only: most are in unique sequence (in absolute numbers) but homopolymers/tandem
  repeats are enriched; shared errors (≥ 3 callers) are 76.9% in segmental duplications.

**Discussion points.** (1) The all-four bar is the largest (4,399) and still holds 32
false calls – 24 in segmental duplications. (2) Calls made by one caller are
mostly false (precision 0.110), except DeepVariant's (32 of
39 real). (3) Real benchmarks use hap.py / vcfeval, which also match differently represented
indels and score genotype errors; here every row is one site-allele, genotypes are not modelled.


### Simplifications in the variant callers data

* **Determinism verified:** each generator was run twice; MD5 checksums of every data file were identical.
* **`dviz.upset` hides the set names:** they are drawn as y tick labels of the dot matrix, but the gridspec gap is
  `wspace=0.05`, so the opaque set-size panel covers them. The notebook fixes this in a visible cell
  (`ax_dot.get_subplotspec().get_gridspec().update(wspace=0.6)`); the slide figures use their own UpSet drawing (same
  counts from `dviz.upset_counts`) that matches the original slide layout (no set-size panel).
* **Caller columns are 0/1, not booleans.** On this file, `vc.loc[vc[c], 'in_truth_set'].mean()` silently returns 1.0
  for every caller, and `upsetplot.from_indicators(callers, data=vc)` raises `ValueError: The indicators must all be
  boolean`. One line after `read_csv` fixes both: `vc[callers] = vc[callers].astype(bool)`. Slide 59, prompt card 7
  and `solutions/figures.py` include it.
* **Truth total fixed at exactly 4,800** (not Poisson around 4,800), contexts multinomial; the FP patterns and
  the per-caller recall table were tuned so that the DeepVariant-only bar is inside the top 10 (rank 9), which slide 20 needs.
* **Sequence contexts are simulated intervals**, not real chr20 annotation (no repeat/segdup/GC track for chr20 is in
  `build/ref/`); alleles are random and do not match hg38.
* **Colours.** Okabe–Ito categories (SNV blue, indel orange; Venn callers orange/sky/green/purple) and vermilion for false
  calls; the original Venn used a red–green palette with mismatched label colours – not reproduced.

## B7 · classifier


**Timing (18 min):** B7.1 splits 3 · B7.2 ROC 3 · B7.3 PR 3 · B7.4 calibration 5 · B7.5 cut-off 4.
Whole-notebook run time on the build machine: 1.9 s with `build/run_notebook.py`
(in-process runner, 18 code cells).

**Expected CHECK values** (all ✓ in the executed notebook):

| CHECK | value | reference | ok |
|---|---|---|---|
| `B7.prevalence_test` | 0.035 | 0.035 ± 0.0001 | ✓ |
| `B7.n_pathogenic_test` | 112 | 112 | ✓ |
| `B7.same_prevalence_all_splits` | True | True | ✓ |
| `B7.train_AUROC_B_inflated` | True | True | ✓ |
| `B7.AUROC_A` | 0.9369 | 0.9369 ± 0.0001 | ✓ |
| `B7.AUROC_B` | 0.9666 | 0.9666 ± 0.0001 | ✓ |
| `B7.AUROC_A_in_design_range` | 0.936928 | in [0.92, 0.94] | ✓ |
| `B7.AUROC_B_in_design_range` | 0.966622 | in [0.96, 0.98] | ✓ |
| `B7.AP_A` | 0.5506 | 0.5506 ± 0.0001 | ✓ |
| `B7.AP_B` | 0.8329 | 0.8329 ± 0.0001 | ✓ |
| `B7.AP_A_above_prevalence` | True | True | ✓ |
| `B7.AP_B_above_prevalence` | True | True | ✓ |
| `B7.AP_B_above_AP_A` | True | True | ✓ |
| `B7.brier_raw_B` | 0.0274 | 0.0274 ± 0.0001 | ✓ |
| `B7.brier_calibrated_B` | 0.0127 | 0.0127 ± 0.0001 | ✓ |
| `B7.brier_improved` | True | True | ✓ |
| `B7.isotonic_never_reverses_order` | True | True | ✓ |
| `B7.AUROC_change_after_calibration` | -0.006 | 0 ± 0.01 | ✓ |
| `B7.mean_error_vs_planted_raw` | 0.1053 | 0.1053 ± 0.0001 | ✓ |
| `B7.mean_error_vs_planted_calibrated` | 0.0049 | 0.0049 ± 0.0001 | ✓ |
| `B7.cut_score_B` | 0.779526 | 0.779526 ± 1e-06 | ✓ |
| `B7.precision_calibration_at_cut` | 0.9011 | 0.9011 ± 0.0001 | ✓ |
| `B7.precision_test_at_cut` | 0.8452 | 0.8452 ± 0.0001 | ✓ |
| `B7.recall_test_at_cut` | 0.6339 | 0.6339 ± 0.0001 | ✓ |
| `B7.test_precision_near_target` | 0.845238 | in [0.8, 1.0] | ✓ |

**Three common mistakes and how to spot them**

1. **Evaluating on all variants or on the train split.** Model B then looks almost perfect (train AUROC
   0.997, AP 0.974). Spot: `B7.AUROC_B` ✗ (reference 0.9666).
   Fix: `test = cs[cs["split"] == "test"]` before any metric.
2. **Fitting the isotonic regression on the test split** (or on all data). The calibrated Brier score then looks better
   than it will be on new data (leakage). Spot: `B7.brier_calibrated_B` ✗ (reference 0.0127).
   Fix: `iso.fit(calib[...])`, `iso.predict(test[...])`.
3. **Swapping the PR axes or forgetting the baseline.** `precision_recall_curve` returns *precision first*;
   `plot(precision, recall)` draws a mirrored curve, and an AP compared with 0.5 instead of the prevalence
   (0.035) looks "bad". Spot: curve starting at the bottom right; no dashed line at 3.5%.
   (Related: choosing the cut-off on the test split gives `B7.cut_score_B` ✗.)

**Answers to the Your-turn tasks**

* B7.1 – calibration vs test AUROC (0.992 vs 0.967): sampling noise with only
  112 pathogenic variants per split (the notebook gives a bootstrap 95% interval for the test AUROC).
* B7.2 – adding benign variants leaves ROC/AUROC unchanged (FPR is a rate) but multiplies the false alarms.
* B7.3 – AP must be compared with its baseline (prevalence): AP 0.20 is ~6× better than random at 3.5%, worse than random at 50%.
* B7.4 – the top isotonic step is 1.00 because a handful of calibration variants above a raw score were all pathogenic;
  extremes are uncertain, calibration is valid only for similar data.
* B7.5 – the calibrated probability at the cut-off (0.76) is that of the *lowest* flagged
  variant; precision (0.845 on test) averages over the whole flagged list.

**Discussion points.** (1) ROC looks excellent for both models (A 0.937, B 0.967); at 3.5%
prevalence PR separates them (AP 0.551 vs 0.833). (2) Isotonic regression is monotone but creates ties:
AUROC moved by -0.0060; Platt (logistic) scaling would keep AUROC exactly – and here the planted
mis-calibration is exactly logit-linear, so Platt scaling would be the better model (good extension task). (3) The cut-off chosen
on the calibration split (precision 0.901) gives precision 0.845 on
the test split: optimism. (4) Real benchmarks split by gene (not by variant) and check that test variants were not used to
train the predictor (circularity).


### Simplifications in the classifier data

* **Determinism verified:** each generator was run twice; MD5 checksums of every data file were identical.
* **`score_modelB_calibrated` is computed, not read.** It is not a column of `classifier_scores.tsv`; students create it
  in Step B7.4 with isotonic regression fitted on the calibration split, as prompt card 8 asks.
* **"AUROC unchanged" is checked as |ΔAUROC| ≤ 0.01 plus a no-reversal check.** Isotonic regression merges scores into
  10 steps (ties), so AUROC moves by -0.0060; the order of the variants is never reversed.
* **Model A colour:** slate (the course colour for a reference group) instead of the orange of the original slide images; model B blue,
  calibrated model B Okabe–Ito green, as before.
* **Extra planted truth:** `data/truth/classifier_truth.tsv` (true probability per variant) is an addition that lets B7.4
  check the calibrated scores against the planted probabilities (mean absolute error 0.105 raw →
  0.005 calibrated).
* **Parameters tuned on the seeded data** so that the test split lands in the intended ranges (A AUROC 0.92–0.94, AP 0.5–0.6;
  B AUROC 0.96–0.98, AP 0.80–0.87) and every raw-score calibration bin has ≥ 20 test variants; the realised values are
  sample values with 112 pathogenic test variants (AUROC B bootstrap 95% interval in the notebook).

## B8 · TP53 lollipop


#### Timing (about 12 min)

| step | beat | minutes | what students do |
|---|---|---|---|
| B8.1 | Question | 2 | load variants and domains; class shares; DBD length 199 |
| B8.2 | Test | 4 | tumours per residue; binomial test per class; Fisher's exact test |
| B8.3 | Check | 2 | counts vs truth; z by hand |
| B8.4 | Visualize | 4 | lollipop with domains; labels for residues with ≥ 30 tumours |

Run time of the whole notebook with `build/run_notebook.py`: 1.9 s.

#### Expected CHECK values (all pass)

| CHECK | value | reference | pass |
|---|---|---|---|
| `B8.n_tumours` | 1,200 | 1200 | ✓ |
| `B8.pct_missense` | 73 | in [70, 76] | ✓ |
| `B8.dbd_length` | 199 | 199 | ✓ |
| `B8.counts_match_truth` | True | True | ✓ |
| `B8.top_residue` | 175 | 175 | ✓ |
| `B8.n_residues_30plus` | 5 | 5 | ✓ |
| `B8.missense_frac_dbd` | 0.8995 | in [0.85, 0.95] | ✓ |
| `B8.missense_binom_p` | 1.3 × 10^-137 | in [0, 1e-50] | ✓ |
| `B8.truncating_frac_dbd` | 0.553 | in [0.45, 0.62] | ✓ |
| `B8.truncating_binom_p` | 0.1396 | in [0.05, 1] | ✓ |
| `B8.fisher_odds_ratio` | 7.237 | in [3, 20] | ✓ |
| `B8.missense_z` | 23.28 | in [10, 40] | ✓ |
| `B8.truncating_z` | 1.517 | in [-1.96, 1.96] | ✓ |
| `B8.n_labels` | 5 | 5 | ✓ |
| `B8.tumours_in_plot` | 1,200 | 1200 | ✓ |

Printed test results: missense 788/876 inside the DBD, binomial
p = 1.3 × 10^-137; truncating 146/264,
p = 0.140; Fisher odds ratio 7.24, p = 7.8 × 10^-33.

#### Three common mistakes and how to spot them

1. **Counting per allele instead of per residue** (`groupby("protein_change")`). R248Q and R248W become two stems, R273
   drops behind R282. *Spot it:* `B8.counts_match_truth` is ✗ and the labels differ from R175H, R248Q/W, R273H/C, R282W, G245S.
2. **Off-by-one domain length or exclusive comparisons** (`292 - 94 = 198`, or `(pos > 94) & (pos < 292)`). The
   expectation becomes 198/393 or residues 94 and 292 drop out. *Spot it:* `B8.dbd_length` ✗, or the inside counts differ
   from 788 / 146. *Fix:* `Series.between(94, 292)` is inclusive.
3. **Over-reading the p-value and over-decorating the plot.** p = 1.3 × 10^-137 is not an effect size (report
   90% vs 51%); colouring the domains with the same colours as the
   variant classes, or labelling every residue, hides the message. *Spot it:* domain boxes in class colours; more than
   5 labels.

#### Your-turn answers (short; full answers are folded in the notebook)

| step | task | answer |
|---|---|---|
| B8.1 | why 199 and not 198? | 1-based inclusive: 292 − 94 + 1 = 199 |
| B8.2 | why prefer Fisher's test? | the binomial test assumes every residue is equally mutable (not true: CpG sites, exon composition); Fisher compares the two classes within the same protein and tumours, so shared biases cancel |
| B8.3 | count per allele: which residues move up? | per allele R175H 66, R248Q 42, R282W 35, R273H 28, G245S 25; per residue R175H 71, R248Q/W 68, R273H/C 60, R282W 40, G245S 38; R248 and R273 combine two alleles |
| B8.4 | label threshold 15? | 14 labels; collisions around 245–249 and 273–282; label only what you discuss |


### Simplifications in the TP53 lollipop data

1. **Class totals are fixed** (876 / 264 / 60) and only the counts within each class are random; a
   plain multinomial over all classes drifted visibly from the intended ~73 / 22 / 5% with this seed (splice/other
   came out well above 5%).
2. **One variant per tumour**; real tumours sometimes carry two TP53 variants. Splice variants are placed at the codon next to
   the splice site (cBioPortal convention); in-frame deletions are counted as splice/other.
3. **Domains:** four domains in the IARC TP53 Database convention; boundaries differ by a few residues in UniProt and
   Pfam; NLS and C-terminal regulatory region not drawn.
4. **TP53 sequence** embedded from UniProt P04637 (no network); checked internally against the
   65 named alleles and the genetic code, not re-downloaded.

