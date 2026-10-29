# Prompt cards 1–8 (lecture Part 5, slides 67–74)

Each card is a prompt you can paste into an AI coding assistant to remake one figure type from the practice data. The text below is the same as on the slides. Using an assistant is optional: the notebooks teach every figure without one.

**How to use a card**

1. Start the conversation with the context line below, then paste the card's prompt.
2. Run the code it returns from the package folder (the folder that contains `data/`).
3. Compare the figure and the printed numbers with your notebook and with **Check against** below. If they differ, find out why before you change anything.
4. Ask for one change at a time (for example “labels overlap; use adjustText”), then look again.
5. Record the tool, the prompt, your edits and your checks in `ai_use_log.md`.

**Context line (paste first):**

```text
Write one self-contained Python script (pandas, numpy, matplotlib, scipy, scikit-learn; adjustText only where asked) that reads the file(s) below from data/ and saves the figure as PNG (300 dpi) and SVG in results/figures/. Print every number that the figure shows. Do not invent columns: use only those listed.
```

Paste only code, column names and printed numbers. Never paste patient data or unpublished results into a public tool.

## Card 1 · Distribution  (slide 67 · notebook step B1)

```text
Data: expression_by_group.tsv (dataset, group, sample, expression in TPM); two panels: small_n (n = 4) and large_n (n = 40).
Question: is GENE_X higher in Case than Control?
Plot a box (median, IQR, 1.5 × IQR whiskers) with every sample as jittered points; n under each label. No bar of means.
Mann–Whitney U per panel: print p and medians; no stars.
Control grey, Case #0072B2; shared y-axis from 0; 120 × 55 mm; PNG 300 dpi + SVG.
```

**Principles enforced:** show every point · n visible · test named · no dynamite · shared axis

**Check against:** small_n: medians 60.8 (Control) vs 89.8 (Case), Mann–Whitney p = 0.029; large_n: p = 8.0 × 10⁻⁵.

**Reference code:** `python solutions/figures.py 1`

*Q. Why say “no bar of means” explicitly in the prompt?* Bar ± error is a common default in generated code; naming what you don’t want prevents it.

## Card 2 · Volcano + MA  (slide 68 · notebook step B2.6)

```text
Data: de_results.tsv; padj is BH; NA padj = not tested: exclude and report how many.
Volcano: x = log2FC_shrunk, y = −log10(padj); dashed lines at padj 0.05 and |log2FC| 1; up #D55E00, down #0072B2, rest grey.
Symmetric x-axis; label the 8 most significant genes with adjustText; n up / n down in the legend.
MA: x = log2(baseMean + 1); raw log2FC grey, shrunken blue.
Print tested/untested/up/down; assert no gene is both; 183 × 70 mm.
```

**Principles enforced:** adjusted p · shrinkage · thresholds · CVD-safe · counts

**Check against:** 7,392 tested, 608 not tested (NA padj); 177 up and 106 down at padj < 0.05 and |log2FC_shrunk| > 1.

**Reference code:** `python solutions/figures.py 2`

*Q. Why tell the model that NA padj means “not tested”?* Otherwise it may zero-fill or colour them; untested genes must be excluded and counted, not shown as significant.

## Card 3 · Clustered heatmap  (slide 69 · notebook step B2.7)

```text
Data: rnaseq_log2norm.tsv.gz (expressed genes × 12 samples), rnaseq_metadata.tsv (condition, batch), de_results.tsv.
Question: do the 40 top DE genes separate samples by condition or by batch?
Z-score each gene; cluster rows and columns with 1 − Pearson r and average linkage; show dendrograms.
Diverging blue–white–orange centred at 0 (−2.5 to 2.5), labelled colour bar; annotation bars for condition and batch.
No jet/rainbow. Assert row means 0, SD 1; print column order.
```

**Principles enforced:** scaling · distance stated · diverging at 0 · annotation · no rainbow

**Check against:** Every row has mean 0 and SD 1 after z-scoring; the first split of the column dendrogram separates control from treated (CHECK B2.heatmap_columns_split_by_condition).

**Reference code:** `python solutions/figures.py 3`

*Q. Why ask for annotation bars for both condition and batch?* So the reader can see whether clusters follow biology (condition) or a technical factor (batch).

## Card 4 · PCA and UMAP  (slide 70 · notebook step B2.3 (part A) · B4 (part B))

```text
A: rnaseq_log2norm.tsv.gz: PCA on the 1,000 most variable genes, centred; PC1 vs PC2 with % variance on the axes;
colour = condition, marker shape = batch; add a scree plot. Print variance for PC1–5.
B: join scrna_cells.tsv (umap1, umap2) to results/B4_clusters.tsv on cell_id; label each cluster on the map.
Remove UMAP axis numbers and add a note: distances between clusters are not interpretable.
Report n_neighbors, min_dist, input PCs and seed.
```

**Principles enforced:** % variance · batch visible · direct labels · no distance claims · settings reported

**Check against:** A: PC1 38.9%, PC2 23.0%. B: 2,940 cells in 8 named clusters; 98.5% of cells named correctly in B4.

**Reference code:** `python solutions/figures.py 4`

*Q. Why remove the axis numbers from a UMAP?* UMAP coordinates have no units or meaning; numbers invite readers to measure distances.

## Card 5 · Manhattan + QQ  (slide 71 · notebook step B3)

```text
Data: gwas_sumstats.tsv.gz (SNP, CHR 1–22, BP hg38 1-based, A1, A2, MAF, BETA, SE, P, P_no_ancestry_PCs); 299,185 SNPs.
Manhattan: cumulative position with gaps, ticks at chromosome centres, two alternating blues, dashed line at 5e-8, rasterized points.
QQ: expected vs observed −log10 p for both P columns, identity line, λ = median χ² / 0.4549 in the legend.
Print λ for both and the number of SNPs with P < 5e-8; assert 0.95 < λ(P) < 1.1.
```

**Principles enforced:** build stated · threshold drawn · λ reported · rasterized · assertion

**Check against:** λ(P) = 1.03, λ(P_no_ancestry_PCs) = 1.33; 53 SNPs with P < 5e-8.

**Reference code:** `python solutions/figures.py 5`

*Q. Which line in this prompt would catch using the wrong p-value column?* The assertion 0.95 < λ(P) < 1.1: the uncorrected column gives λ = 1.33 and fails.

## Card 6 · Genome tracks  (slide 72 · notebook step B5)

```text
Data: cnv_chr22_10kb.tsv (chrom, start, end, gc, mappable, depth_patient, depth_reference; BED-style 0-based starts).
Question: where is the deletion, how large, how strong is the evidence?
log2 ratio = log2[(patient ÷ its median) ÷ (reference ÷ its median)] on mappable bins; segment it with dviz.segment.
Tracks on one x-axis in Mb: depth, log2 + segments, genes; lines at 0 and −1.
Print the call (1-based start), its length, median log2 inside/outside; shade it on every track.
```

**Principles enforced:** shared axis · convention stated · expected value · segment shown · numbers printed

**Check against:** Call chr22:18,940,001–21,110,000 (2.17 Mb); median log2 -0.96 inside, 0.05 outside.

**Reference code:** `python solutions/figures.py 6`

*Q. Why include the coordinate convention in a plotting prompt?* Off-by-one errors between BED (0-based) and 1-based labels are silent; stating it prevents wrong breakpoints.

## Card 7 · UpSet + lollipop  (slide 73 · notebook step B6 · B8)

```text
Data: variant_callsets.tsv (GATK, DeepVariant, FreeBayes, Strelka2, in_truth_set as 0/1: convert with .astype(bool)).
UpSet: top 10 intersections sorted by size, dot matrix, set sizes. Do NOT draw a four-set Venn.
Lollipop: tp53_variants.tsv (one row per tumour) + tp53_domains.tsv; x = position 1–393, stem = tumours, colour = class (Okabe–Ito), domains as boxes; label ≥ 30 tumours.
Print each caller’s precision against in_truth_set.
```

**Principles enforced:** sets as bars · sorted · no Venn · domains shown · truth check

**Check against:** Precision: GATK 0.969, DeepVariant 0.990, FreeBayes 0.949, Strelka2 0.972. Lollipop labels (≥ 30 tumours): R175H (71), R248Q/W (68), R273H/C (60), R282W (40), G245S (38).

**Reference code:** `python solutions/figures.py 7`

*Q. Why ask for precision against the truth set as well as the UpSet?* Callers can agree on the same errors; only a truth set measures correctness.

## Card 8 · ROC, PR, calibration  (slide 74 · notebook step B7)

```text
Data: classifier_scores.tsv (split, label_pathogenic, score_modelA, score_modelB). Fit calibration on 'calibration'; evaluate on 'test' only.
Question: which model to use when ~3.5% of variants are pathogenic?
Panels: ROC with AUROC; PR with AP and a baseline at prevalence; calibration, 10 bins (skip < 20): B raw vs isotonic.
Print prevalence, AUROC and AP per model; assert AP > prevalence.
```

**Principles enforced:** test split only · PR baseline · calibration · prevalence stated · assertion

**Check against:** Test split: prevalence 3.5%; AUROC A 0.937, B 0.967; AP A 0.551, B 0.833.

**Reference code:** `python solutions/figures.py 8`

*Q. Why draw the prevalence as a baseline on the precision–recall plot?* A random classifier’s precision equals the prevalence (3.5%), so it shows how much better the model is.
