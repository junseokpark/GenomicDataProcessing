#!/usr/bin/env bash
# Rebuild every public input and every synthetic file of Practice 2 (run from the package root).
# Seed: 20261015 (scripts/config.py). Needs internet access to api.genome.ucsc.edu and hgdownload.soe.ucsc.edu.
set -euo pipefail
python scripts/fetch_public_data.py --stage all          # reference windows, chain, gene tables, hg19/hg38 tracks
python scripts/simulate_dna.py --samples P1 P2 P3 C1 --outdir data/fastq   # shipped FASTQs
python scripts/build_pon.py                              # C1 (from data/fastq) + C2-C5 (temporary) -> panel of normals
python scripts/make_truth.py                             # truth BED/TSV + ClinVar exact matches
python scripts/simulate_rna.py                           # synthetic RNA counts (hg38 MANE genes)
echo "done: data/ rebuilt"
