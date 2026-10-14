# Practice 1 · From reads to variants
### Calling SNVs and indels in a synthetic trio (GRCh38)
*Genomic Data Processing · Day 4 (10/15) · Understanding genetic variants*

A child has sickle cell disease and both parents are healthy carriers. You will take the family's sequencing reads all the way to a named variant: read QC → alignment → duplicate marking → coverage → pileup and genotype likelihoods → joint variant calling → normalisation → filtering → benchmarking against the known truth → trio checks → the clinical (HGVS) name **NM_000518.5(HBB):c.20A>T p.(Glu7Val)**.

> **All data are synthetic.** The reads were simulated from the public GRCh38 reference sequence with a fixed random seed. No real person's DNA or genotypes are included.

---

## 1. Set up (once)

You need **conda** (or the faster **mamba**). We recommend **Miniforge**, which includes both: <https://github.com/conda-forge/miniforge>.

**macOS and Linux** – open a terminal in this folder and run:
```bash
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate d4-practice
```

> **One environment for both Day 4 practices.** The same `d4-practice` environment also runs Practice 2. If it already exists, skip `conda env create` and just run `conda activate d4-practice`.
On Apple-silicon Macs, if conda cannot find a bioconda package for your machine, create the environment for Intel Macs instead (it runs through Rosetta):
```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
```

**Windows** – the bioinformatics tools do not run natively on Windows. Use **WSL2** (Windows Subsystem for Linux):
1. In PowerShell (as administrator): `wsl --install`, restart, and open the *Ubuntu* app.
2. Inside Ubuntu, install Miniforge for Linux (`curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh` then `bash Miniforge3-Linux-x86_64.sh`).
3. Copy this folder into your Linux home directory (for example `cp -r /mnt/c/Users/<you>/Downloads/practice1_variant_calling ~/`). Working under `/mnt/c/...` also works but is much slower.
4. Run the two commands of the macOS/Linux section.

Check that the tools are found:
```bash
samtools --version | head -n 1     # samtools 1.24
bcftools --version | head -n 1     # bcftools 1.24
```

## 2. Run the practice

```bash
cd practice1_variant_calling
conda activate d4-practice
jupyter lab
```
Open **`practice1_variant_calling_student.ipynb`** and run the cells from top to bottom (Shift + Enter). Every step ends with a **Check** cell: if it prints numbers and no red error, you are on track. Questions marked **Your turn** have their answer hidden under *Answer* – try first, then click.

Running all cells takes well under a minute (the optional internet query in Step 11 usually adds 10–20 seconds, longer if the Ensembl server is busy). If you get stuck, `practice1_variant_calling_executed.ipynb` shows every output as it should look.

To run everything without opening Jupyter:
```bash
jupyter nbconvert --to notebook --execute practice1_variant_calling_student.ipynb --output my_run.ipynb
```

## 3. What is in this folder

```
practice1_variant_calling/
├── README.md                                    this file
├── INSTRUCTOR_NOTES.md                          expected values, timings, common errors, answers
├── environment.yml                              conda environment (tool versions)
├── practice1_variant_calling_student.ipynb      the exercise (no outputs)
├── practice1_variant_calling_executed.ipynb     the same notebook with all outputs
├── data/
│   ├── ref/chr11_GRCh38_betaglobin.fa.gz (+ .fai, .gzi)   reference: chr11, only the beta-globin window kept
│   ├── fastq/S01_R1.fastq.gz … S03_R2.fastq.gz             paired-end reads (2 x 150 bp, ~30x)
│   ├── trio.ped                                            family file: S01 child, S02 father, S03 mother
│   ├── truth/truth_trio.vcf.gz (+ .tbi)                    the planted variants (normalised)
│   ├── truth/truth_raw_unnormalised.vcf.gz                 the same variants, written the "raw" way (Step 7)
│   ├── truth/planted_events.tsv                            answer key: every planted variant and artefact
│   ├── annotation/HBB_MANE_NM_000518.5_hg38.tsv            HBB exons and coding sequence (MANE Select)
│   ├── annotation/genes_window_hg38.bed                    MANE genes in the window
│   ├── provenance/                                         downloaded source files and their versions
│   └── MD5SUMS.txt                                         checksums of the data files
├── scripts/
│   ├── generate_data.py                         rebuilds every file in data/ (seed 20261015)
│   └── make_slide_figures.py                    redraws the lecture figures from your results
└── results/
    ├── p1_results.json                          key numbers of the reference run
    └── figures/p1_fig_*.png                     the six lecture figures
```
The notebook writes its own outputs (BAM, VCF, QC reports, logs) into `results/`.

## 4. About the data

**Reference genome.** `data/ref/chr11_GRCh38_betaglobin.fa.gz` contains chromosome 11 of GRCh38 at its full length (135,086,622 bp), but only `chr11:5,200,001–5,300,000` – the β-globin locus with *HBB*, *HBD*, *HBG1*, *HBG2* and *HBE1* – holds real sequence. All other bases are `N`. This keeps the files small (the reference is about 0.5 MB) while every position is a real GRCh38 position, so you can compare results with genome browsers, ClinVar or Ensembl. The sequence was downloaded from the UCSC Genome Browser API (hg38) on 2026-10-14 (see `data/provenance/provenance.json`).

**Annotation.** Exon and coding-sequence coordinates of *HBB* come from the UCSC hg38 `mane` track (MANE v1.5: NCBI RefSeq release GCF_000001405.40-RS_2025_08, Ensembl 116).

**Reads.** About 10,400 read pairs per person (≈ 30× depth over the window), fragment length ≈ 380 bp, NovaSeq-like binned base qualities (`F` = 37, `:` = 25, `,` = 11, `#` = 2), sequencing errors drawn from those qualities, and about 4 % PCR duplicates. A few problems were planted on purpose; you will discover them during the practice (the answer key is `data/truth/planted_events.tsv`).

**Rebuilding the data.** Everything in `data/` (except the downloaded files in `data/provenance/`) can be rebuilt:
```bash
python scripts/generate_data.py              # uses the saved UCSC downloads, no internet needed
python scripts/generate_data.py --download   # downloads the sequence and annotation again first
```
With the same package versions the files are identical; compare with `data/MD5SUMS.txt` (`md5sum -c data/MD5SUMS.txt` on Linux/WSL).

## 5. Common problems

| message | what to do |
|---|---|
| `command not found: samtools` (or bcftools, minimap2 …) | the environment is not active: `conda activate d4-practice`, then restart Jupyter |
| `AssertionError: Start Jupyter inside the practice1_variant_calling folder` | start `jupyter lab` from inside this folder |
| `NameError: name '...' is not defined` | a cell was skipped; run the notebook from the top (*Run → Run All Cells*) |
| `[W::bcf_hdr_check_sanity] MQ should be declared as Type=Float` | harmless warning, see the note in Step 6 |
| Step 11 prints "Could not reach the Ensembl REST API" | you are offline or the public server is busy; the step is optional |
