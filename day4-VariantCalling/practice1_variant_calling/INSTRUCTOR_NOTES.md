# Practice 1 · Instructor notes
*From reads to variants: calling SNVs and indels in a synthetic trio (GRCh38) · Day 4 (10/15)*

All numbers below come from the reference run (`practice1_variant_calling_executed.ipynb`, `results/p1_results.json`). The data are deterministic (seed 20261015), so students using the pinned environment should see **identical** values. With other tool versions small differences are possible; the `assert` statements in the Check cells use ranges or structural conditions, not exact numbers, wherever that is reasonable.

## 1. At a glance
* **Audience:** beginners; no command-line experience assumed. Every step: *Why* → code → *Check* → *What changed? / Your turn*.
* **Compute:** < 1 min for all 49 code cells (≈ 7 s without the optional VEP query; ≈ 17 s with it on the build machine). Memory < 1 GB. The run writes ≈ 11 MB into `results/` (BAM, VCF, QC reports, depth table).
* **Data:** 3 × ~10,400 read pairs (2 × 150 bp) over GRCh38 chr11:5,200,001–5,300,000; full-length chr11 reference with all other bases set to `N`.
* **Story:** child S01 with sickle cell disease (HBB c.20A>T homozygous), carrier parents S02 (father) and S03 (mother).

## 2. Suggested timing (75 min)
| block | notebook steps | minutes | key idea to land |
|---|---|---|---|
| Inputs and QC | 0, 1a–1c | 10 | coordinates are real GRCh38; Phred scale; QC measures, it does not change reads |
| Alignment, duplicates, coverage | 2–4 | 15 | read groups; duplicates are one molecule seen many times; 8× vs 30× |
| Pileup, likelihoods, calling | 5–6 | 15 | genotype = which hypothesis explains the reads best; QUAL vs GQ |
| Normalise, filter, benchmark | 7–9 | 15 | same variant, different spellings; filters trade recall for precision |
| Trio and naming | 10–11 | 15 | Mendelian logic; minus strand; predictors ≠ classification |
| Wrap-up | 12 | 5 | exit questions |

## 3. Expected Check values per step
| step | what students should see |
|---|---|
| 0 | samtools 1.24, bcftools 1.24, htslib 1.24, minimap2 2.31-r1302, fastp 1.4.0, pysam 0.24.1, python 3.12 |
| 1a | chr11 length 135,086,622; base at chr11:5,227,002 = **T**; chr11:5,226,991–5,227,010 = `CAGACTTCTCCTCAGGAGTC`; 0 N inside the window |
| 1b | first record `@SYNTH_S01_000001 1:N:0:1`, 150 bp; read pairs S01 10,400, S02 10,400, S03 10,415 (S03 includes the 15-copy jackpot) |
| 1c | Q30 88.7 / 88.8 / 88.7 %; Q20 96.8 %; GC 38.8–38.9 %; 0 reads fail fastp filters; mean quality of the last 10 cycles R1 33.4 > R2 31.9; fastp duplication ≈ 0.1 % (see answer in notebook) |
| 2 | 100.0 % mapped, 100.0 % properly paired in all samples; mean fragment length 377.7 / 378.9 / 378.9 bp; example read `SYNTH_S01_004055` FLAG 99, MAPQ 60, 150M |
| 3 | duplicate reads 3.86 % (S01), 3.86 % (S02), 3.99 % (S03); markdup stats S03: `DUPLICATE PAIR: 832`; jackpot chr11:5,214,982 (ref G): before {C 15, G 37}, after {C 1, G 37}; 15 pairs share read start 5,214,922 / mate start 5,215,173 |
| 4 | mean depth 30.0× in each sample; ≥ 20×: 96.8 / 97.3 / 96.2 % of the window; low-coverage interval chr11:5,291,001–5,292,000: 6.9 / 8.5 / 9.0× |
| 5 | pileup at 5,227,002: S01 T 0 / A 25 (VAF 1.00); S02 T 17 / A 12 (0.41); S03 T 20 / A 17 (0.46). Worked example (S02, ε = 0.001): log10 L = −36.0 / −8.7 / −51.0; PL by hand 273 / **0** / 423; bcftools PL 227 / **0** / 255 |
| 6 | 33 variant records; HBB record QUAL 700.3, INFO/DP 94, MQ 60; S01 1/1 DP 25 AD 0,25 GQ 74 PL 255,75,0; S02 0/1 DP 29 AD 17,12 GQ 127; S03 0/1 DP 37 AD 20,17 GQ 127 |
| 7 | norm: `33/1/0/2/0/0/0` (total/split/joined/realigned/…) → 34 records. Deletion `CTTTTTTTT>CTTTTTTT` → `CT>C` at 5,239,407; insertion `TCACACACACACACAC>TCACACACACACACACAC` → `T>TCA` at 5,230,685; multi-allelic `G>A,C` at 5,245,200 (S01 1/2, S02 0/1, S03 0/2) → two records. Truth demo: raw vs raw 31 shared / 2 / 2; normalised calls vs raw truth 31 / 3 / 2; both normalised 34 / 0 / 0. Truth norm: `33/1/0/1/0/0/0` |
| 8 | FILTER: 33 PASS, 1 LowQual = chr11:5,291,568 A>T (QUAL 95.1, INFO/DP 23; S03 0/1, DP 8, AD 3,5, GQ 78) → S03 genotype `./.`; one genotype per person set to missing (all at that site) |
| 9 | pooled TP/FP/FN: main 71/0/0 (precision 1.000, recall 1.000); no markdup 71/1/0 (0.986, 1.000; FP S03 5,214,982 G>C); `-Q 0` 71/1/0 (0.986, 1.000; FP S02 5,257,081 G>T, QUAL 19.8, AD 14,18); filters 70/0/1 (1.000, 0.986; FN S03 5,291,568). Per sample TP: S01 28, S02 21, S03 22. Genotype concordance 1.0 everywhere. Low-quality site: bases ≥ Q13 {G 14}; all bases {G 14, T 18} |
| 10 | mendelian2: sites_good 32, sites_merr 1, sites_missing 1; error site chr11:5,272,758 G>A (S01 0/1 AD 16,9; S02 0/0 AD 27,0; S03 0/0 AD 19,0; every read in the parents: G only); recessive filter → only chr11:5,227,002 |
| 11 | CDS 444 bp, starts ATG, ends TAA, protein MVHLTPEEKS…; `NM_000518.5(HBB):c.20A>T p.(Glu7Val)` (GAG>GTG); `NM_000518.5(HBB):c.9T>C p.(His3=)` (CAT>CAC); `NC_000011.10:g.5227002T>A`. VEP (Ensembl 116, REST 15.12, retrieved 2026-10-14): missense_variant, gAg/gTg, CADD 14.4, AlphaMissense 0.2229 likely_benign, SpliceAI Δ 0, colocated rs334 |

## 4. Answer key: planted events (`data/truth/planted_events.tsv`)
| event | position (GRCh38) | change | S01 | S02 | S03 | what it teaches | observed |
|---|---|---|---|---|---|---|---|
| HBB c.20A>T (rs334) | chr11:5,227,002 | T>A | 1/1 | 0/1 | 0/1 | recessive disease variant; minus strand naming | called, recessive hit |
| HBB c.9T>C (rs713040) | chr11:5,227,013 | A>G | 0/1 | 0/1 | 0/0 | synonymous variant; in *cis* with c.20 on the father's transmitted copy | called |
| 1-bp deletion in T×8 | chr11:5,239,407 | CT>C | 0/1 | 0/1 | 0/0 | left-alignment (raw truth stores it as 5,239,414 TT>T) | called |
| 2-bp insertion in (CA)×7 | chr11:5,230,685 | T>TCA | 0/0 | 0/0 | 0/1 | indel in a repeat; untransmitted maternal allele | called |
| 4-bp deletion, unique | chr11:5,206,900 | TGACC>T | 0/1 | 0/0 | 0/1 | simple indel | called |
| multi-allelic SNV | chr11:5,245,200 | G>A,C | 1/2 | 0/1 | 0/2 | splitting with `norm -m -any` | called, split |
| de novo SNV | chr11:5,272,758 | G>A | 0/1 | 0/0 | 0/0 | Mendelian error | the only Mendelian error |
| low-coverage SNV | chr11:5,291,568 | A>T | 0/0 | 0/0 | 0/1 | depth filters remove true variants | called, then filtered |
| 25 background SNVs | 5,202,015–5,296,652 | various | – | – | – | ordinary inherited variation; none fits the recessive pattern | all called |
| **artefact:** PCR jackpot | chr11:5,214,982 | G>C (reads only) | – | – | 15 copies | duplicate marking | false 0/1 in S03 without markdup |
| **artefact:** systematic low-quality error | chr11:5,257,081 | G>T (reads only) | – | 20/34 pairs, Q12 | – | minimum base quality | false 0/1 in S02 with `-Q 0` (or bcftools' default `-Q 1`) |
| low-coverage interval | chr11:5,291,001–5,292,000 | 27 % of fragments kept | | | | depth plot | 6.9–9.0× |

Truth totals: 33 variant sites (34 records after splitting); non-reference genotypes S01 28, S02 21, S03 22 (71 in total).

## 5. Design notes (read before teaching)
1. **bcftools' default minimum base quality is 1, not 13.** In bcftools 1.24 `mpileup -Q` defaults to 1 (see `bcftools mpileup` help); `samtools mpileup` uses 13, and some online manuals still list 13 for bcftools. The course pipeline therefore sets `-Q 13` explicitly. Running the main command **without** `-Q` gives the same false call as `-Q 0`.
2. **The low-quality artefact uses base quality 12, not 8.** In tests bcftools' likelihood model never called quality-8 errors with `-Q 0`, even at 70 % of reads, because it already gives such bases little weight. To make the `-Q 0` lesson work, the systematic error has quality 12 (still below 13) in 20 of 34 S02 read pairs (≈ 59 %, AD 14,18 with `-Q 0`). Its QUAL is only 19.8, so the QUAL < 30 filter would also catch it – a good discussion point ("two safety nets").
3. **The low-coverage true variant is in the mother only** (S03 0/1, untransmitted), so exactly one genotype is lost by the filters. The generator chose its exact position from a dry run of the fragment sampling so that S03 has ~8 reads there.
4. **Main pipeline is perfect by design**, so each scenario shows exactly one planted problem. Real data are never this clean – say so.
5. **Overlapping read pairs.** samtools/bcftools count an overlapping pair once; in the Step 5 pileup an overlapping pair appears as a single base with summed quality (`k` = Q74).
6. **Executed notebook.** The build machine's sandbox did not allow Jupyter kernels to open sockets, so the executed notebook was produced by running every cell in order inside one IPython process (same outputs as `jupyter nbconvert --execute`). Students run it normally in JupyterLab.

## 6. Common errors and fixes
| symptom | cause | fix |
|---|---|---|
| `command not found: minimap2` (or samtools, bcftools, fastp) in a `%%bash` cell | Jupyter started outside the conda environment | `conda activate d4-practice`, restart `jupyter lab`; the first cell also adds the environment's `bin` to PATH |
| `AssertionError: Start Jupyter inside the practice1_variant_calling folder` | notebook opened from another working folder | start Jupyter in the package folder |
| `NameError` (e.g. `pile`, `NORM`, `base_counts`) | cells skipped or kernel restarted | *Run → Run All Cells* |
| minimap2 output has no `@RG` line / `SM:` missing | `\t` lost when the command is retyped with other quotes | keep the double quotes: `-R "@RG\tID:${S}\tSM:${S}..."` |
| `[W::bcf_hdr_check_sanity] MQ should be declared as Type=Float` | htslib 1.24 header check | harmless; explained in Step 6; such messages are sent to `results/logs/warnings.log` |
| `[mpileup] 3 samples in 3 input files` in red | samtools information message | harmless |
| `samtools markdup: error, no MC tag` | `fixmate -m` step skipped | run the four-command recipe in order |
| `Could not reach the Ensembl REST API` / HTTP 500 | offline or public server busy | optional step; re-run the cell later |
| very slow on Windows | files under `/mnt/c/...` in WSL | copy the folder to the WSL home directory |
| bioconda package not found on Apple silicon | missing osx-arm64 build | `CONDA_SUBDIR=osx-64 conda env create -f environment.yml` |
| different numbers from the table above | other tool versions | use `environment.yml`; compare `data/MD5SUMS.txt` |

## 7. Discussion notes beyond the in-notebook answers
* **QUAL vs GQ.** QUAL is about the site (anyone variant?), GQ about one person's genotype. A site can have very high QUAL and still contain a person with a low GQ (e.g. low depth in that person).
* **Why bcftools' PL are smaller than the hand calculation.** Per-read qualities, quality caps, an error-dependency model (errors at the same site tend to co-occur), and PL values stored up to 255. Ranking matters more than the exact numbers.
* **The child's GQ (74).** Limited by depth: PL(0/1) ≈ −10·log10(0.5²⁵) ≈ 75. Ask students what depth would be needed for GQ 99 (≈ 33 reads).
* **Multi-allelic splitting.** After `-m -any`, the child's `1/2` becomes `1/0` and `0/1`; the AD of each split record only counts REF and that ALT (S01 AD 0,12 and 0,7).
* **VCF left-alignment vs HGVS 3′ rule.** Opposite directions for plus-strand genes; the same direction for minus-strand genes such as *HBB*.
* **Real benchmarking.** GIAB HG002/HG003/HG004 (a trio), stratified by region; hap.py / RTG vcfeval compare haplotypes, so representation differences do not count as errors.
* **De novo rate.** Roughly 1–1.5 × 10⁻⁸ per base per generation, i.e. ~50–100 new SNVs per genome; in a 100-kb window a true de novo is very unlikely, which is why such calls need confirmation.
* **HbS and predictors.** AlphaMissense (0.22, likely benign) and CADD (14.4) assess the typical impact of an amino-acid change; HbS acts through polymerisation of deoxygenated haemoglobin. ClinVar: pathogenic for sickle cell anaemia; the same allele protects carriers against severe malaria (VEP lists "protective" among rs334's clinical significances).
* **p.Glu7Val vs "β6 Glu→Val".** HGVS counts the initiator methionine; classical haemoglobin numbering does not.
* **Phase.** In the father, c.9T>C and c.20A>T are 11 bp apart and on the same reads (in *cis*); visible in IGV.

## 8. Rebuilding data and figures
* Data: `python scripts/generate_data.py` (offline, from `data/provenance/`) or `--download` (re-fetch from the UCSC API). Seed 20261015; output must match `data/MD5SUMS.txt` with the pinned environment.
* Slide figures: run the notebook, then `python scripts/make_slide_figures.py --font-dir <folder with SUIT-Regular.ttf and SUIT-Bold.ttf>` (falls back to the default sans-serif font).
* Provenance: UCSC hg38 sequence and `mane` track retrieved 2026-10-14 (`downloadTime` fields in `data/provenance/`); MANE v1.5 (NCBI GCF_000001405.40-RS_2025_08, Ensembl 116), track dataTime 2025-12-08.
