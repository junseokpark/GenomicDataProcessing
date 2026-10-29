# Data visualization practice · Genomic Data Processing, 10/27 (D5)

Hands-on files for Part 6 of the lecture *Data visualization in genetics and genomics*. Every dataset is
**simulated with real anchors** (real gene names, genome coordinates, published loci and protein domains) and has a
**planted answer**, so you can check whether your analysis and your figure find it. There is no patient data.

* **Practice A · simple data → simple figures.** Small tables you can read by eye (scatter, distributions, bars vs dots,
  composition, time course, 60,000-point scatter).
* **Practice B · realistic data → test, check, then draw.** Every exercise has four beats:
  **Question** (one sentence) → **Test** (run the method a real analysis would use) → **Check** (compare with the
  truth in a `CHECK` cell) → **Visualize** (draw the figure that shows the result honestly).

## Quick start

**Google Colab.** Upload `dataviz_practice_dataset.zip`, then in a code cell:

```python
!unzip -q dataviz_practice_dataset.zip
%cd dataviz_practice
!pip install -q umap-learn adjustText
```

Open `practice.ipynb` (File → Open notebook → Upload, or the file browser) and run the cells from the top.

**Your own computer.** With conda: `conda env create -f environment.yml`, `conda activate dataviz-practice`, then
`jupyter lab` in this folder. With pip: `pip install -r requirements.txt`. Open `practice.ipynb`.

Run the notebook **from top to bottom**. If something looks wrong, use *Restart kernel and run all*. The whole core
notebook runs in about 10 seconds; each elective in under a minute.

## Plan for the session (about 90 minutes, in pairs)

| Part | Minutes | Where |
|---|---|---|
| Setup and file check | 5 | `practice.ipynb` · Step 0 |
| Exercise 1 · diagnose eight misleading figures (M1–M8) | 8 | `bad_figures/`, `checklist.md` |
| Practice A · simple data (A1–A5) | 20 | `practice.ipynb` |
| Practice B · B1 two groups, B2 bulk RNA-seq | 27 | `practice.ipynb` |
| One elective per pair (B3–B8) | 20 | `electives/` |
| Defend one figure: verify, accessibility, optional AI remake | 10 | `practice.ipynb` · D1–D3 |
| Present one figure | in class | `checklist.md` |

**Electives** (choose one): `B3_gwas` (GWAS with population structure: λ and ancestry PCs), `B4_single_cell`
(normalise, PCA, k-means, UMAP and t-SNE: what the maps keep), `B5_cnv_tracks` (find the 22q11.2 deletion with
genome tracks), `B6_variant_callers` (four callers against a truth set: UpSet), `B7_classifier` (ROC, precision–recall
and calibration at 3.5% prevalence), `B8_tp53_lollipop` (where TP53 mutations fall; about 12 minutes).

## Folder layout

```
practice.ipynb        everyone: Setup, Exercise 1, Practice A (A1–A5), B1, B2, Defend (D1–D3)
electives/            B3_gwas · B4_single_cell · B5_cnv_tracks · B6_variant_callers · B7_classifier · B8_tp53_lollipop
data/                 all inputs (table below); data/truth/ = the planted answers, open them only in CHECK steps
dviz.py               shared helpers: style and colours, BH adjustment, λ, segmentation, calibration, colour-vision check, UpSet, check
bad_figures/          M1–M8, the misleading figures of Exercise 1
checklist.md          16-item figure checklist (Exercise 1, Defend, Present)
prompts/              prompt_cards.md (cards 1–8 of the lecture) and critique_prompt.md (optional AI)
ai_use_log.md         record of any AI use (tool, prompt, edits, checks)
solutions/            figures.py (reference code for the demos and prompt cards), ANSWER_KEY_bad_figures.md,
                      executed/ (every notebook with its outputs)
results/              written by the notebooks: <exercise>.json and figures/ (empty at the start)
fonts/                SUIT font (SIL Open Font License 1.1), so figures match the slides
docs/                 DATA_DICTIONARY.md: every column of every file, how it was simulated
INSTRUCTOR_NOTES.md   for instructors: timing, expected values, common mistakes, answers
```

Instructors: `INSTRUCTOR_NOTES.md` and `solutions/` contain the answers. Remove them before sharing the folder with students if you want to release the answers after the session.

## Conventions

* **Genome build hg38** everywhere. **Intervals** are BED-style (`chrom`, `start`, `end`; 0-based, half-open);
  **single positions** (SNP `BP`, CpG `pos`, variant `pos`, protein `position`) are 1-based.
* **Missing values** are empty fields, never 0. In `de_results.tsv` an empty `padj` means *not tested*: exclude
  and count those genes; never draw them as significant or at p = 1.
* **CHECK cells** print `CHECK ✓ name = value (reference …)` when your result matches the planted answer, `✗` when
  it does not, and `CHECK · name = value (recorded)` for values that are only recorded. All values are saved to
  `results/<exercise>.json`. A ✗ is information: find out why before you change anything.
* **Truth files** (`data/truth/`) are opened only inside CHECK steps. A real analysis has no truth file; that is
  why every method here is first checked on data where the answer is known.
* **Seeds**: every generator uses `20261027 + k`; rerunning a notebook gives the same numbers (run times excepted).

## Data dictionary (short)

| File (`data/`) | Rows × columns | Size | What it is | Used in |
|---|---|---|---|---|
| `anscombe.tsv` | 44 × 3 | 468 B | Anscombe’s quartet (published values, Anscombe 1973): 4 sets × 11 points | A1 |
| `mapping_rates.tsv` | 12 × 3 | 148 B | % reads mapped, 3 pipeline versions × 4 libraries | A3, M2 |
| `expression_by_group.tsv` | 88 × 4 | 2 kB | one gene (GENE_X, TPM), Control vs Case: `small_n` (4 + 4) and `large_n` (40 + 40) | A2, B1, M1 |
| `timecourse.tsv` | 54 × 4 | 873 B | FKBP5, TSC22D3, PER1 after dexamethasone: 6 times × 3 wells (log2 FC vs 0 h) | A4, M8 |
| `celltype_composition.tsv` | 32 × 4 | 629 B | 8 blood cell types × 4 donors (counts and %), the same 2,940 cells as B4 | A4, M5 |
| `methylation_cpg.tsv.gz` | 60,000 × 6 | 527 kB | 60,000 CpGs (chr1–22, 1-based `pos`), β by bisulfite (WGBS) and nanopore, with depths | A5 |
| `rnaseq_counts.tsv.gz` | 8,000 × 15 | 245 kB | bulk RNA-seq counts, 8,000 genes × 12 samples (6 control, 6 dexamethasone; batches A and B) | B2 |
| `rnaseq_metadata.tsv` | 12 × 4 | 285 B | sample, condition, batch, RIN | B2 |
| `rnaseq_log2norm.tsv.gz` | 6,940 × 13 | 273 kB | log2(median-of-ratios normalised count + 1), 6,940 expressed genes (mean ≥ 10) | B2, cards 3–4 |
| `de_results.tsv` | 8,000 × 8 | 562 kB | PyDESeq2, design ~ batch + condition: baseMean, log2FC, apeGLM-shrunken log2FC, lfcSE, p, BH padj (empty = not tested) | B2, card 2 |
| `de_results_no_batch.tsv` | 8,000 × 7 | 498 kB | the same without batch in the design (no shrunken column) | B2 |
| `gwas_sumstats.tsv.gz` | 299,185 × 10 | 8.3 MB | GWAS summary statistics, AMD-like disease, 1,500 cases + 1,500 controls; P with 4 ancestry PCs, P_no_ancestry_PCs without | B3, card 5 |
| `gwas_genotypes_subset.npz` | matrix | 7.3 MB | genotypes (0/1/2) of 3,000 people × 10,097 SNPs for the PCA and the trend test | B3 |
| `gwas_snps_subset.tsv` | 10,097 × 8 | 440 kB | the 10,097 SNPs of the genotype subset (hg38, 1-based BP) | B3 |
| `gwas_samples.tsv` | 3,000 × 5 | 51 kB | IID, case (0/1), site, sex, age | B3 |
| `gwas_genes_hg38.tsv` | 18,352 × 5 | 593 kB | gene positions for labelling (MANE Select, BED-style `start0`) | B3 |
| `scrna_counts.npz` | matrix | 863 kB | UMI counts, 2,940 cells × 2,000 genes (scipy.sparse CSR; load with `scipy.sparse.load_npz`) | B4 |
| `scrna_cells.tsv` | 2,940 × 6 | 135 kB | cell_id (barcode), donor, n_counts, n_genes, umap1, umap2 (rows = matrix rows) | B4, card 4 |
| `scrna_genes.tsv` | 2,000 × 2 | 50 kB | gene symbol and Ensembl ID (columns = matrix columns) | B4 |
| `scrna_markers.tsv` | 8 × 2 | 283 B | marker genes per cell type from the literature (prior knowledge, not the truth) | B4 |
| `cnv_chr22_10kb.tsv` | 5,082 × 7 | 233 kB | read depth of a patient and a reference panel in 10-kb bins, chr22 (BED-style 0-based starts), GC, mappable | B5, card 6 |
| `cnv_genes_chr22.tsv` | 82 × 8 | 7 kB | genes on chr22 (MANE Select, BED-style) | B5 |
| `cnv_segdups_chr22.tsv` | 590 × 9 | 37 kB | real segmental duplications, chr22 17.98–22.41 Mb (UCSC genomicSuperDups, hg38), LCR22A–D marked | B5 |
| `cnv_genome_1mb.tsv` | 2,887 × 5 | 109 kB | genome-wide log2 ratios in 1-Mb bins for the patient and a tumour | B5 |
| `variant_callsets.tsv` | 5,210 × 12 | 320 kB | 5,210 candidate variants on chr20 with calls by GATK, DeepVariant, FreeBayes, Strelka2 and the truth (0/1) | B6, card 7 |
| `classifier_scores.tsv` | 16,000 × 6 | 661 kB | 16,000 variants: split (train/calibration/test), label_pathogenic, scores of models A and B | B7, card 8 |
| `tp53_variants.tsv` | 1,200 × 4 | 31 kB | 1,200 TP53 variants, one per tumour: protein change, residue (1–393), class | B8, card 7 |
| `tp53_domains.tsv` | 4 × 6 | 450 B | TP53 domains (UniProt P04637): TAD, PRD, DBD, TET | B8, card 7 |

Truth files in `data/truth/`: `expression_by_group_truth.tsv` (responders among the 40 large-study cases),
`rnaseq_truth.tsv` (is_de, true log2 FC, batch shift, known glucocorticoid targets), `gwas_truth.tsv/.json` (8 planted
loci with odds ratios), `gwas_ancestry.tsv`, `scrna_truth.tsv/.json` (cell types, simulation parameters),
`cnv_truth.tsv` (planted copy-number events), `callsets_truth.json`, `classifier_truth.tsv/.json` (true
probabilities), `tp53_truth.tsv` (expected counts). Full column-by-column descriptions, units and how each file was
simulated: `docs/DATA_DICTIONARY.md`.

## Checking your work

1. The **CHECK cells** in each notebook: your numbers against the planted answer.
2. `results/<exercise>.json`: every CHECK value, for your legend and for the D1 table on the slides.
3. `python solutions/figures.py` (or `python solutions/figures.py 2 5`): the reference code for demos and prompt
   cards 1–8; it prints the numbers behind each figure (card 4B needs `results/B4_clusters.tsv`, written by B4).
4. `solutions/executed/`: every notebook with its outputs, to compare after you have tried.

## Using AI (optional)

AI assistants are allowed but not needed. `prompts/prompt_cards.md` has the eight prompt cards from the lecture;
`prompts/critique_prompt.md` asks a model to critique a figure against `checklist.md`. Whatever writes the code,
the CHECK cells decide whether the figure is right. Record what you did in `ai_use_log.md`. Paste only code, column
names and printed numbers; never paste patient data or unpublished results into a public tool.

## Tested with

* Python 3.12: numpy 2.5, pandas 3.1, scipy 1.18, matplotlib 3.11, scikit-learn 1.9, umap-learn 0.5.12,
  adjustText 1.4, pydeseq2 0.5.4.
* Python 3.11: numpy 1.26, pandas 2.2, scipy 1.11, matplotlib 3.8, scikit-learn 1.4, umap-learn 0.5.5,
  adjustText 1.2 (without pydeseq2).

Every CHECK passed with both sets. Numbers that depend on t-SNE and UMAP (B4) move slightly between versions (for example, the share of 15 nearest neighbours a UMAP keeps was 15.8% and 15.2%); those CHECKs use ranges.

pydeseq2 is needed only for the optional re-run in Step B2.5b; every other step reads `data/de_results.tsv`.
Without umap-learn, B4 uses the stored `umap1`/`umap2` columns (computed with the same settings).

## Sources and licences

* Data: simulated for teaching (seeds above); no patient data. Real anchors: gene symbols and coordinates from
  MANE Select (NCBI/EMBL-EBI, hg38); glucocorticoid-response genes as in Himes et al. (2014) *PLoS ONE* 9:e99625;
  AMD loci as in Fritsche et al. (2016) *Nat Genet* 48:134–143; the 22q11.2 region and LCR22A–D (McDonald-McGinn
  et al., 2015, *Nat Rev Dis Primers* 1:15071); TP53 domains from UniProt P04637; benchmark design after Genome in
  a Bottle (Zook et al., 2019, *Nat Biotechnol* 37:561–566); Anscombe (1973) *Am Stat* 27:17–21.
* Font: SUIT, SIL Open Font License 1.1 (`fonts/OFL.txt`).
* The generator scripts that made every file (byte-identical on rerun) are in `build/` of the course repository
  (github.com/junseokpark/GenomicDataProcessing, folder `day5-visualization/dataviz_practice/build/`).
