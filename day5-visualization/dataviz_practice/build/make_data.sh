#!/usr/bin/env bash
# Rebuild every file in data/ with the ten generators, in dependency order (fixed seeds: the files come out
# byte-identical to the ones in the repository). Needs pydeseq2 for make_rnaseq.py; make_gwas.py needs about 3 GB RAM.
# Run from dataviz_practice/ :   bash build/make_data.sh
set -euo pipefail
[ -f dviz.py ] || { echo "Run this from dataviz_practice/ (the folder that contains dviz.py)." >&2; exit 1; }
export MPLBACKEND=Agg PYTHONPATH="$PWD:$PWD/build" NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$PWD/.numba_cache}"
for g in make_simple make_cell_labels make_rnaseq make_log2norm make_scrna make_gwas make_cnv make_tp53 make_callsets make_classifier; do
  t0=$(date +%s)
  python "build/$g.py" > /dev/null
  echo "build/$g.py  $(( $(date +%s) - t0 )) s"
done
