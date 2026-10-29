# Figure review checklist

Use it on your own figure before you present it, and on your partner's figure when you review it.
Sixteen items in four groups. Tick an item only when you can point to the evidence (a printed number, a label, a line of code).

## 1. Message

| ✓ | Item | How to check |
|---|------|--------------|
| ☐ | **One question; the finding is stated** | You can say the question in one sentence and the title or legend states the answer with a number ("177 up, 106 down at FDR 5%"). |
| ☐ | **The chart type fits the question** | Comparing groups → points + box or dot plot. Change over time → lines. Parts of a whole → sorted bars. Thousands of genes → volcano / MA. Positions on a genome → tracks. |

## 2. Data honesty

| ✓ | Item | How to check |
|---|------|--------------|
| ☐ | **Points shown when n is small; n stated** | With fewer than ~20 values per group, every value is visible. n is written on the figure or in the legend, with its unit (donors, wells, cells, libraries). |
| ☐ | **Bars start at zero; no dual y-axes** | Bar length must be proportional to the value. If you need to zoom, use dots on a labelled axis. Two series with the same unit share one axis. |
| ☐ | **Test, correction and thresholds named** | The legend names the test ("Mann–Whitney U"), the correction ("Benjamini–Hochberg"), and every cut-off line has its value ("padj 0.05", "5 × 10⁻⁸"). |
| ☐ | **Nothing hidden by overplotting or draw order** | Dense scatter → transparency, hexbin or density. Feature plots → expressing cells drawn last. Missing values are reported, not drawn as 0 or p = 1. |
| ☐ | **No distance or size claims from UMAP / t-SNE** | Claims about how far apart or how big clusters are come from a test or from the original data, not from the embedding. |

## 3. Encoding and colour

| ✓ | Item | How to check |
|---|------|--------------|
| ☐ | **The key comparison is on position** | The numbers your question compares are read along a common axis (not from angles, areas or colour shades). |
| ☐ | **Perceptual or diverging colour maps** | Ordered values → viridis-type map with a colour bar. Values above and below a centre (z-scores, log2 ratios) → diverging map centred at 0. No rainbow / jet. |
| ☐ | **Works in greyscale and for colour-vision deficiency** | Run `dviz.show_cvd("results/figures/<your figure>.png")`. Groups stay distinct because of shape, labels or lightness, not hue alone. |

## 4. Labels and integrity

| ✓ | Item | How to check |
|---|------|--------------|
| ☐ | **Units, genome build and coordinate convention** | Axis labels have units ("TPM", "Mb, hg38"). Interval files say BED (0-based starts); positions say 1-based. |
| ☐ | **Text ≥ 7 pt at the final size** | Save at the final width (89 mm single column, 183 mm double) and check the smallest label. |
| ☐ | **Script, seed and vector output** | The figure is made by a script or notebook cell you can rerun; random steps have a seed; you saved SVG/PDF as well as PNG. |
| ☐ | **Every number on the figure matches the printed output** | Each count, p-value or percentage in the title, legend or labels appears in a `CHECK ·` cell output. |
| ☐ | **AI use logged; no generative-AI images of data** | If an assistant wrote code or text, `ai_use_log.md` has an entry. Images of data or results are never made by an image generator. |
| ☐ | **Limitations said out loud** | You can name one thing the figure cannot show (n, batch, simulation, model assumptions). |

---
*The four groups follow the lecture (Parts 2–5). The CVD simulation is in `dviz.py` (Machado et al., 2009).*
