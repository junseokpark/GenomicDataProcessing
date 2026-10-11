# Practice session package — Genomic Data Processing, Day 1 (10/12)

Everything the 30-minute hands-on block needs. Unpack anywhere and open
`genomic_data_processing_practice.ipynb` from the directory **above** `practice_session/`
(or edit the `BASE` variable in the first code cell).

## Requirements

| Needed | Why |
|---|---|
| Python 3.9+, `pandas` | the whole notebook |
| `pysam` *(optional)* | only for reading BAM; every BAM ships with a `.sam` text twin, so the notebook runs without it |
| IGV *(optional)* | section 5; session files are generated for you |

No network access, no reference download, no aligner index build. Nothing in the notebook
takes more than a second to run.

## What is real and what is synthetic

**Real** — the `chr11` sequence in both mini-references is UCSC sequence
(hg19 chr11:5,240,001–5,300,000 and the byte-identical hg38 chr11:5,218,771–5,278,770);
the gene and exon coordinates in `genes.hg38.bed` are NCBI RefSeq Select transcripts for the
β-globin cluster; the `HBB_checkpoint` correspondence hg19 `chr11:5,248,232` ↔ hg38
`chr11:5,227,002` is independently documented (ClinVar RCV005357130,
`NM_000518.5(HBB):c.20A>T`) and is used **for its coordinates only**.

**Synthetic** — samples S01–S03, all reads, all genotypes, all depths, the `toy1` contig, and
the chain file. No result produced here is a clinical statement about anyone.

Both mini-references are real chr11 sequence **N-padded up to the real coordinate**, so
positions in IGV, `samtools` and the notebook are genuine hg19/hg38 addresses. They are
teaching mini-chromosomes, not genome builds: their contig lengths (5,300,000 and 5,280,000)
are deliberately *not* the real chr11 lengths, which is what makes the §1 length check work.

## Layout

```
practice_session/
├── manifest.tsv                        three samples, their builds and their files
├── collaborator_share/                 what the hg19 collaborator sent us
│   ├── reference.hg19.fa (+ .fai)      hg19 teaching mini-chr11
│   ├── regions.S0{1,2,3}.hg19.tsv      legacy 1-based closed interval tables
│   ├── S0{1,2,3}.hg19.vcf              variant calls (S03 has the planted build conflict)
│   ├── S01.hg19.sorted.bam (+ .bai, .sam)
│   ├── S01_R{1,2}.fastq.gz             the original reads, 80 bp, 91 pairs
│   ├── toy1.fa / toy1.vcf              the ATC>A and A>G examples from slides 51–54
│   └── genes.toy1.bed                  the [200,250) gene for the boundary test
├── target/                             the hg38 analysis environment
│   ├── reference.hg38.fa (+ .fai)
│   ├── genes.hg38.bed                  HBB, HBD, HBG1, HBG2, HBE1 + exons + promoters
│   └── hg19ToHg38.teaching.over.chain  single-block teaching chain (see below)
├── checkpoints/                        verified outputs, for when a cell will not run
│   └── S01.hg38.sorted.bam (+ .bai, .sam)
└── outputs/                            everything the notebook writes
```

### The chain file
`hg19ToHg38.teaching.over.chain` is **not** the UCSC chain. It is one ungapped block over
hg19 `chr11:5,246,500–5,280,000`, where the two assemblies are byte-identical, so
hg38 = hg19 − 21,230 holds exactly. The notebook proves this from the two FASTAs rather than
asking you to believe it. Records outside the block fail to map **on purpose** — the
`unmapped` output is part of the lesson. For real work:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver input.hg19.bed hg19ToHg38.over.chain.gz mapped.hg38.bed unmapped.hg19.bed
```

## How the notebook maps onto the slides

| Notebook section | Slides | Practice step |
|---|---|---|
| 1 Reference / history | 59–62, 84 | 01–02 |
| 2 hg19 coordinates | 47–50, 63–64, 85 | 03 |
| 3 VCF semantics | 51–55, 65–68, 86 | 04 |
| 4 hg19 → hg38 mapping | 69–72, 87 | 05 |
| 5 Assembly-matched IGV | 46, 73–74, 88 | 06 |
| 6 Read realignment | 41–45, 75–78, 89–90 | 07–08 |
| 7 Batch and QC | 56–57, 79–80, 91 | 09 |
| 8 Evidence card | 81–82, 92–93 | 10 |
| Appendix | 20–25 | — |

## Deliberate failures (do not "fix" them)

| Where | What | Correct response |
|---|---|---|
| S03 | manifest and header both claim hg38; the contig length says otherwise | `fail`, record the conflict, convert nothing |
| `r05` | interval with `start > end` | `fail`; never silently swap the two numbers |
| `r04`, `S01_offchain` | outside the chain block | `unmapped`, keep source coordinates |
| `S01_symdel`, `S01_bnd1` | symbolic `<DEL>` and a breakend | `review`; do not force them into an ordinary interval |
| `S01_ins1` | insertion — no reference base is lost | `review`; an affected-base interval is the wrong question |

A justified stop is a valid result, and these five rows are where most of the marks are.

## Instructor notes

* Total runtime of every cell is a few seconds; the 30 minutes are for reading and deciding.
* Section 7's last cell is the one to put on the board: the un-converted coordinate returns a
  clean empty result for a single base, and the **wrong gene (HBG1)** for a gene span.
* The five presentation questions on slide 82 map to checkpoints 2, 3, 4, 5 and 7 in order.
* Rebuild the whole package from scratch with `build_lab.py` and `make_notebook.py`; the only
  network calls are the two UCSC sequence fetches, cached in `lab_src/`.
