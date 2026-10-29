# build/ · how the Day 5 practice was made

Every file in `data/`, the notebooks and the lecture-slide figures come from the scripts in this folder. Students do not
need it. Run every command from `dataviz_practice/` (the folder that contains `dviz.py`), after:

```bash
export MPLBACKEND=Agg PYTHONPATH="$PWD:$PWD/build" NUMBA_CACHE_DIR=/tmp/numba_d5
python -m pip install pydeseq2      # only for make_rnaseq.py (the provided DESeq2 results)
```

## 1. Data

`bash build/make_data.sh` runs the ten generators in dependency order (about 75 s; `make_gwas.py` needs about 3 GB RAM).
All of them use fixed seeds, so a rerun rewrites the 40 files in `data/` byte for byte.

| Generator | Writes (in `data/`) | Used in |
|---|---|---|
| `make_simple.py` | anscombe.tsv, expression_by_group.tsv, mapping_rates.tsv, methylation_cpg.tsv.gz, timecourse.tsv, truth/expression_by_group_truth.tsv | Practice A, B1 |
| `make_cell_labels.py` | celltype_composition.tsv, truth/scrna_truth.tsv | Practice A, B4 |
| `make_rnaseq.py` | de_results.tsv, de_results_no_batch.tsv, rnaseq_counts.tsv.gz, rnaseq_metadata.tsv, truth/rnaseq_truth.tsv (and build/rnaseq_build_summary.json) | B2 |
| `make_log2norm.py` | rnaseq_log2norm.tsv.gz | B2, prompt cards 3–4 |
| `make_scrna.py` | scrna_cells.tsv, scrna_counts.npz, scrna_genes.tsv, scrna_markers.tsv, truth/scrna_truth.json | B4 |
| `make_gwas.py` | gwas_genes_hg38.tsv, gwas_genotypes_subset.npz, gwas_samples.tsv, gwas_snps_subset.tsv, gwas_sumstats.tsv.gz, truth/gwas_ancestry.tsv, truth/gwas_truth.json, truth/gwas_truth.tsv | B3 |
| `make_cnv.py` | cnv_chr22_10kb.tsv, cnv_genes_chr22.tsv, cnv_genome_1mb.tsv, cnv_segdups_chr22.tsv, truth/cnv_truth.tsv | B5 |
| `make_tp53.py` | tp53_domains.tsv, tp53_variants.tsv, truth/tp53_truth.tsv | B8 |
| `make_callsets.py` | truth/callsets_truth.json, variant_callsets.tsv | B6 |
| `make_classifier.py` | classifier_scores.tsv, truth/classifier_truth.json, truth/classifier_truth.tsv | B7 |

`make_cell_labels.py` must run before `make_scrna.py`, and `make_rnaseq.py` before `make_log2norm.py`.

## 2. Notebooks

The builders write `practice.ipynb` and `electives/*.ipynb`. Numbers quoted in their text are read from `results/*.json`
of the previous run, so build and execute twice; the first pass shows "…" where a number is not known yet.
`answers_B4.py` computes the numbers behind the folded B4 answers and must run first.

```bash
python build/answers_B4.py                     # -> results/B4_yourturn.json
for pass in 1 2; do
  python build/build_practice_notebook.py
  python build/run_notebook.py practice.ipynb practice_executed.ipynb
  for X in B3_gwas:B3 B4_single_cell:B4 B5_cnv_tracks:B5 B6_variant_callers:B6 B7_classifier:B7 B8_tp53_lollipop:B8; do
    NB=${X%%:*}; EX=${X##*:}
    python build/build_${EX}_notebook.py
    python build/run_notebook.py electives/$NB.ipynb electives/${NB}_executed.ipynb
  done
done
```

`run_notebook.py` executes the cells in-process (no Jupyter server needed). Every `# CHECK ·` line must print ✓. The
executed copies are written next to the notebooks and are ignored by git; copy them to `solutions/executed/` when they
are final.

## 3. Figures

```bash
python build/make_bad_figures.py      # bad_figures/M1–M8 (Exercise 1) and the slide panel image-84-1
python build/slides_core.py           # lecture-slide figures from the Practice A, B1 and B2 data
python build/slides_practice.py       # practice_B1.png and practice_B2.png (Part 6 slides)
python build/slides_gwas.py
python build/slides_single_cell.py    # reads results/B4.json: run the B4 notebook first
python build/slides_tracks.py         # reads results/B5.json and results/B8.json
python build/slides_evaluation.py
```

Slide figures go to `build/slide_figs/`, each at the exact pixel size of its frame on the slide
(`build/picture_frames.tsv`). That folder is ignored by git.

## Helpers and reference files

| File | What it is |
|---|---|
| `nbbuild.py`, `nb_part_a.py`, `nb_part_b.py` | cells and text of the notebooks |
| `scrna_helpers.py` | B4 analysis helpers, copied into the B4 notebook so notebook, generator and slides run the same code |
| `slidefig.py` | saves a figure at the exact pixel size of its slide frame |
| `rnaseq_build_summary.json` | DESeq2 summary written by `make_rnaseq.py`, quoted in the practice notebook |
| `ref/hg38_chrom_sizes.tsv` | hg38 chromosome lengths |
| `ref/genes_hg38_mane_select.tsv.gz` | MANE Select genes on hg38 (0-based start, BED style) |
| `ref/d4_tracks_hg38/` | public UCSC annotation tracks downloaded for Day 4 (segmental duplications, ClinGen, ClinVar, DGV, gnomAD-SV, ENCODE cCREs, phastCons, pLI, MANE); retrieval date and regions in `provenance_tracks_hg38.json` |
