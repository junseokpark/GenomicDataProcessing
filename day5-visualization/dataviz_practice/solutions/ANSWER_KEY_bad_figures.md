# Answer key · Exercise 1 · Spot the misleader (M1–M8)

Open this **after** you have written your own answer for each figure.
All eight figures are made from the practice data by `build/make_bad_figures.py`. They deliberately use naive defaults.
For each figure: the main principle it breaks, the evidence you can see, the fix, and where the practice shows the fix.

| Figure | Main principle broken | Where the fix is shown |
|---|---|---|
| M1 | Bar of means hides n = 4 and the spread; star without a named test | `practice.ipynb` Step B1 |
| M2 | Bars do not start at zero | Step A3 |
| M3 | Rainbow colour map; rows not scaled; no colour bar or labels | Step B2.7 (heatmap) |
| M4 | Draw order hides the expressing cells | `electives/B4_single_cell.ipynb` |
| M5 | Pie with eight slices; angle and area for a precise comparison | Step A4 |
| M6 | Red–green; raw p-values; no thresholds; untested genes drawn as p = 1 | Steps B2.5–B2.6 (volcano) |
| M7 | Distance claim from a UMAP | `electives/B4_single_cell.ipynb` |
| M8 | Dual y-axes make two series look related; causal title | Step A4 (time course) |

---

### M1 · "GENE_X is higher in cases" (bar of means ± error, with a star)
* **Evidence.** Two bars of means with error bars and a star. You cannot see that each group has only **n = 4**
  samples, that the error bars are standard errors (SE, not stated), or which test gave the star. One Case sample
  (150 TPM) pulls the Case mean up to 100.2, while the Case median is 89.8.
* **Why it misleads.** With n = 4, the SE (17.3) is half the SD (34.7), so the bars look tight. The star is not
  robust: the Mann–Whitney U test gives p = 0.029 (with 4 vs 4 this is the smallest p the test can give), but Welch's
  t-test gives p = 0.107. Without the points you cannot judge which to believe.
* **Fix.** Show every sample as a point with the median; write n under each group; name the test and give the exact p;
  if you show error bars, say whether they are SD, SE or 95% CI. (The 95% CI of the Case mean is ±55.1 TPM.)

### M2 · "Pipeline C maps far more reads!" (truncated bar axis)
* **Evidence.** The y-axis starts at 95.8%. Bar C (98.0%) looks **7.3 times** as tall as bar A (96.1%); the true
  ratio is 1.02 — a difference of 1.9 percentage points. Green highlights the "winner".
* **Why it misleads.** A bar encodes its value by length from zero; cutting the axis breaks that.
* **Fix.** Start bars at zero, or use dots (one per library) with the mean on a zoomed axis that is clearly labelled.
  Decide whether 1.9 points matters for the analysis; the title should state the difference, not "far more".

### M3 · "Top genes" (rainbow heatmap)
* **Evidence.** Jet/rainbow colours, no colour bar, row numbers instead of gene names, sample numbers instead of names,
  no annotation for condition or batch, rows not clustered.
* **Why it misleads.** Rainbow maps are not perceptually ordered (yellow and cyan look like sharp boundaries) and fail
  for colour-vision deficiency. Raw log2 counts are dominated by how highly each gene is expressed, so the pattern
  shows expression level more than the treatment effect.
* **Fix.** Z-score each gene (row mean 0, SD 1); use a diverging map centred at 0 with a labelled colour bar; add
  condition and batch annotation bars; state the distance and linkage used for clustering.

### M4 · "FCGR3A is rarely expressed" (UMAP feature plot)
* **Evidence.** Most cells look pale. The cells with high FCGR3A were drawn **first**, and the many cells with zero
  counts were drawn on top of them. No colour bar, so you cannot tell what the colours mean.
* **Numbers.** FCGR3A has at least one UMI in 23.2% of all cells, but in 93.7% of FCGR3A monocytes and 90.2% of NK cells (FCGR3A encodes CD16, carried by both). It is not rare; it is specific, and the draw order hid it (computed from data/scrna_counts.npz with the truth labels).
* **Fix.** Draw cells in order of expression (zeros first, highest last), add a colour bar with the unit
  (log-normalised counts), or use a dot plot (fraction of cells expressing × mean expression per cell type).

### M5 · "Cell types" (exploded 3-D-style pie)
* **Evidence.** Eight slices, shadows and "exploded" slices, a legend instead of labels, no numbers.
* **Why it misleads.** People compare angles and areas poorly. You cannot tell whether CD14 monocytes (19.0%) are
  more common than memory CD4 T cells (16.3%). Shadows and gaps change the apparent areas.
* **Fix.** Sorted horizontal bars with the percentage written at the end of each bar and n = 2,940 cells in the
  title; for per-donor comparisons, one small bar chart per cell type rather than stacked bars.

### M6 · Volcano plot in red and green
* **Evidence.** Red (up) and green (down) only; y-axis labelled "p" (it is −log10 of the **raw** p-value); x-axis
  labelled "FC" (it is log2 fold change, not shrunken); no threshold lines; no counts.
* **Hidden problems.** 758 genes have raw p < 0.05; if no gene changed you would still expect about 370 by chance among
  the 7,392 tested genes. After Benjamini–Hochberg correction 396 genes have padj < 0.05. The 608 genes that DESeq2 did
  not test (padj = NA) were silently drawn as p = 1 (y = 0) instead of being removed and counted.
* **Fix.** y = −log10(padj); x = shrunken log2 fold change; dashed lines at padj 0.05 and |log2FC| 1; up in vermilion,
  down in blue (Okabe–Ito), the rest grey; write "177 up / 106 down"; report the 608 untested genes in the legend.
  Check it with the deuteranopia simulation: red and green become nearly the same colour.

### M7 · "B cells are the most different (farthest from everything)" (UMAP)
* **Evidence.** A line drawn between cluster centres on a UMAP, and a claim about distance.
* **Why it misleads.** UMAP keeps neighbourhoods but not distances between clusters or cluster sizes; the same data
  with another seed or `min_dist` gives a different layout.
  In the practice data the most isolated cluster in the 30 PCs is NK, not B; B is the most isolated island only on this UMAP (15 neighbours). On t-SNE it is CD8 T, on the UMAP with 5 or 200 neighbours NK (results/B4.json, slide_numbers).
* **Fix.** Use the UMAP only to show which cells group together; label clusters directly, remove the axis numbers,
  state the settings (n_neighbors, min_dist, input PCs, seed). To compare cell types, use marker genes or distances in
  the original (PCA) space.

### M8 · "TSC22D3 drives FKBP5" (dual y-axes)
* **Evidence.** Two lines with **two different y-axes**; the right axis range was picked so the curves end together.
  No axis labels or units, no legend, no error bars, and the title makes a causal claim.
* **Why it misleads.** Both series are already in the same unit (log2 fold change vs 0 h), so a second axis only
  serves to make them look aligned. Over the six time points the two means correlate (r = 0.86), but correlation in a
  time course is not causation: both genes respond to dexamethasone.
* **Fix.** One shared axis with direct labels, mean ± SD with the three wells as points (n = 3 per time point), and a
  descriptive title ("TSC22D3 rises within 1 h; FKBP5 keeps rising to 24 h").
