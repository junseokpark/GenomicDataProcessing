"""Simulate PBMC-like 10x UMI counts for Practice B4 (single cell).

Cells: exactly the 2,940 cells of data/truth/scrna_truth.tsv (fixed labels, not changed here).
Genes: 2,000 real protein-coding genes (MANE Select, hg38): canonical PBMC markers, ribosomal and
housekeeping genes, and randomly drawn background genes (most of them barely expressed in blood).

Model (per cell c of type t, gene g):
  expected UMIs   mu[c, g] = L[c] * p[t(c), g] * donor[d(c), g] * cell_state[c, g]
  observed UMIs   counts ~ Poisson(mu * Gamma(shape = 1/phi_g, scale = phi_g))   (gamma-Poisson = NB)
  L[c]            library size, log-normal (median ~2,200 UMIs, clipped to 550-7,600), x cell-type factor
  p[t, g]         expression profile of type t: marker levels set by hand to resemble 10x PBMC 3k
                  (approximate UMIs per cell), normalised to sum to 1
Naive <-> memory CD4 T and CD14 <-> FCGR3A monocytes are continua (a per-cell Beta score mixes the
two profiles in log space), so each pair is close but separable. CD8 T cells have an effector score.
NKG7 is exactly 0 in B cells and CD14 monocytes (used by the violin-plot slide).

Writes data/scrna_counts.npz, data/scrna_genes.tsv, data/scrna_cells.tsv (incl. a fallback UMAP),
data/scrna_markers.tsv. Deterministic: seed 20261027 + 5; run twice -> byte-identical files.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import scrna_helpers as H  # noqa: E402

SEED = 20261027 + 5
N_GENES = 2000
TYPES = H.CELL_TYPES          # Naive CD4 T, Memory CD4 T, CD8 T, NK, B, CD14 Mono, FCGR3A Mono, DC
MEDIAN_UMIS = 2200
COLS = ["NV", "MEM", "CD8", "NK", "B", "C14", "F16", "DC"]   # TABLE columns, same order as TYPES

# Approximate UMIs per cell (at ~2,200 UMIs) for each gene and cell type, set by hand to resemble
# 10x PBMC data (Seurat PBMC 3k tutorial markers). Columns follow TYPES. Small values = ambient RNA.
TABLE = """
gene      NV    MEM   CD8   NK    B     C14   F16   DC    group
ACTB      25    28    28    30    25    40    40    35    hk
B2M       60    65    70    70    45    45    55    35    hk
GAPDH     4     5     5     6     4     10    9     9     hk
TMSB4X    50    55    60    65    40    65    70    55    hk
TMSB10    12    14    14    16    12    15    18    14    hk
FTL       6     8     8     10    8     70    65    45    hk
FTH1      10    12    12    14    10    45    45    30    hk
EEF1A1    45    40    40    32    40    25    22    28    hk
TPT1      35    32    30    24    30    18    16    20    hk
MYL6      6     8     8     10    6     12    12    10    hk
PTMA      12    10    10    8     10    6     6     8     hk
PFN1      6     7     8     10    6     12    14    10    hk
CFL1      4     5     6     7     5     10    10    8     hk
HLA-A     6     7     8     8     6     6     7     6     hk
HLA-B     12    14    16    16    10    12    14    10    hk
HLA-C     6     7     8     8     5     6     7     6     hk
HLA-E     3     4     4     5     3     4     5     4     hk
UBA52     5     5     5     5     5     4     4     4     hk
UBB       3     3     3     3     3     4     4     4     hk
UBC       3     3     3     3     3     4     4     4     hk
SRGN      3     4     6     8     2     12    10    6     hk
ARHGDIB   4     5     6     7     4     6     8     6     hk
CORO1A    3     4     4     4     3     4     4     4     hk
LAPTM5    2     2     2     3     3     5     6     4     hk
PTPRC     2     2     2     2     2     2     2     2     hk
EIF1      6     6     6     6     6     8     8     7     hk
NACA      3     3     3     3     3     2     2     2     hk
EEF2      3     3     3     3     3     3     3     3     hk
H3-3B     3     3     3     3     3     4     4     4     hk
FAU       8     7     7     6     7     5     5     5     hk
SH3BGRL3  2     3     3     4     2     6     6     5     hk
LCP1      2     2     2     2     2     4     4     3     hk
CD52      5     6     6     5     4     3     3     2     hk
CD53      2     2     2     2     2     3     3     3     hk
S100A6    1.5   4     4     6     1.5   12    10    6     hk
CALM2     2     2     2     2     2     3     3     3     hk
ARPC3     2     2     2     3     2     4     4     3     hk
ARPC2     2     2     2     3     2     4     5     3     hk
GMFG      1     1     1.5   2     1     2     2     2     hk
RAC2      2     2     2.5   3     2     3     3     3     hk
MYL12A    2     2.5   3     4     2     4     5     4     hk
CD3D      5     5     5     0.4   0.05  0.05  0.05  0.05  t
CD3E      4     4     4     0.3   0.05  0.05  0.05  0.05  t
CD3G      2     2     2     0.2   0.02  0.02  0.02  0.02  t
CD2       1.5   2.5   2.5   1.5   0.05  0.05  0.05  0.05  t
CD247     0.8   1     1.5   2     0.02  0.02  0.02  0.02  t
LCK       1.2   1.2   1.5   0.8   0.05  0.02  0.02  0.02  t
BCL11B    0.8   0.8   0.8   0.05  0.01  0.01  0.01  0.01  t
TRAT1     0.8   0.6   0.5   0.02  0.01  0.01  0.01  0.01  t
CD5       0.6   0.7   0.4   0.02  0.05  0.01  0.01  0.01  t
CD6       0.6   0.7   0.5   0.05  0.01  0.01  0.01  0.01  t
CD27      1.4   1.0   1.0   0.1   0.4   0.02  0.02  0.02  t
CD28      0.5   0.6   0.3   0.02  0.01  0.01  0.01  0.01  t
IL32      2.5   9     7     2     0.3   0.3   0.3   0.3   t
LTB       12    12    3     0.8   6     0.6   0.6   1     t
LDHB      6     4     4     1.5   3     0.4   0.4   0.6   t
CD4       0.3   0.4   0.02  0.02  0.01  0.6   0.6   0.5   cd4
MAL       1.5   1.2   0.6   0.05  0.05  0.02  0.02  0.02  cd4
CD40LG    0.4   0.6   0.05  0.02  0.02  0.01  0.01  0.01  cd4
IL7R      2.5   4     1     0.3   0.1   0.1   0.05  0.1   cd4
CCR7      3     0.5   0.8   0.05  0.8   0.05  0.05  0.2   naive
LEF1      2.5   0.4   1     0.05  0.2   0.02  0.02  0.02  naive
TCF7      2     0.6   0.9   0.2   0.05  0.02  0.02  0.02  naive
SELL      2.2   0.35  0.7   0.4   1     2     0.3   0.3   naive
NOSIP     1.5   0.7   0.8   0.5   0.6   0.4   0.4   0.4   naive
PIK3IP1   1.6   0.6   0.9   0.4   0.6   0.2   0.2   0.2   naive
LRRN3     0.6   0.1   0.2   0.01  0.01  0.01  0.01  0.01  naive
SATB1     1.2   0.6   0.7   0.3   0.3   0.05  0.05  0.05  naive
TXK       0.6   0.3   0.3   0.4   0.01  0.01  0.01  0.01  naive
PRKCA     0.5   0.2   0.2   0.05  0.05  0.02  0.02  0.02  naive
FHIT      0.8   0.3   0.4   0.05  0.1   0.05  0.05  0.05  naive
NELL2     1.2   0.3   0.8   0.05  0.05  0.02  0.02  0.02  naive
TSHZ2     0.8   0.2   0.1   0.01  0.01  0.01  0.01  0.01  naive
BACH2     0.9   0.3   0.5   0.1   0.4   0.02  0.02  0.02  naive
SERINC5   0.8   0.25  0.4   0.1   0.2   0.05  0.05  0.05  naive
ABLIM1    0.9   0.3   0.4   0.05  0.3   0.02  0.02  0.05  naive
PECAM1    0.7   0.2   0.3   0.3   0.05  0.3   0.3   0.1   naive
SOX4      0.6   0.1   0.2   0.02  0.1   0.02  0.02  0.02  naive
CHI3L2    0.7   0.1   0.2   0.01  0.01  0.01  0.01  0.01  naive
AK5       0.8   0.2   0.3   0.02  0.02  0.01  0.01  0.01  naive
ACTN1     0.9   0.3   0.4   0.1   0.05  0.05  0.05  0.05  naive
IL6ST     0.9   0.4   0.4   0.05  0.1   0.05  0.05  0.05  naive
NUCB2     1.0   0.4   0.5   0.3   0.2   0.05  0.05  0.05  naive
S100A4    0.8   7     5     6     1     20    15    5     memory
GPR183    0.6   2     0.6   0.3   2     1     0.5   0.8   memory
ANXA1     0.4   3.5   2     2     0.3   2.5   1     0.6   memory
KLRB1     0.1   2.5   1.5   1.5   0.05  0.05  0.05  0.05  memory
AQP3      0.4   2.5   0.5   0.1   0.1   0.05  0.05  0.05  memory
ITGB1     0.4   1.5   0.6   0.5   0.4   0.3   0.3   0.3   memory
CRIP1     0.5   2     1.5   1     0.3   0.3   0.5   0.3   memory
LMNA      0.2   1.2   0.6   0.4   0.1   0.6   0.5   0.4   memory
VIM       2     4     3     3     1.5   6     6     4     memory
S100A11   0.3   1.5   1     1     0.3   8     8     3     memory
S100A10   0.4   2     1.5   0.8   0.3   1     1     0.6   memory
CD58      0.2   0.8   0.6   0.6   0.2   0.5   0.5   0.4   memory
FAS       0.1   0.6   0.4   0.2   0.1   0.05  0.05  0.05  memory
ICOS      0.1   0.7   0.1   0.05  0.01  0.01  0.01  0.01  memory
LGALS1    0.3   2     1.5   3     0.3   4     5     2     memory
EMP3      0.5   2     1.5   1.5   0.5   2     3     1     memory
TIMP1     0.3   1.5   0.8   0.5   0.2   3     3     1     memory
MYADM     0.3   1.2   0.8   0.6   0.2   1     1     0.6   memory
ANXA2     0.2   1.5   1     0.8   0.1   2     2     1     memory
KLF6      0.5   1.5   1.2   1.2   0.3   2.5   1.5   1     memory
TNFRSF4   0.05  0.8   0.05  0.02  0.01  0.01  0.01  0.01  memory
MAF       0.1   0.8   0.3   0.1   0.05  0.05  0.1   0.05  memory
RORA      0.3   1.2   0.6   0.5   0.05  0.05  0.05  0.05  memory
PRDM1     0.05  0.8   0.5   0.3   0.02  0.02  0.02  0.02  memory
CXCR3     0.1   0.9   0.6   0.3   0.02  0.02  0.02  0.02  memory
CCR6      0.05  0.8   0.1   0.02  0.3   0.01  0.01  0.01  memory
CCR4      0.1   0.9   0.1   0.02  0.01  0.01  0.01  0.01  memory
CTLA4     0.05  0.6   0.05  0.02  0.01  0.01  0.01  0.01  memory
IL2RA     0.05  0.6   0.05  0.02  0.01  0.01  0.01  0.01  memory
CD8A      0.02  0.02  2     0.15  0.01  0.01  0.01  0.01  cd8
CD8B      0.02  0.02  1.6   0.05  0.01  0.01  0.01  0.01  cd8
CCL5      0.1   1     15    10    0.1   0.1   0.2   0.1   cytotoxic
GZMK      0.02  0.3   2.5   0.6   0.01  0.01  0.01  0.01  cytotoxic
GZMA      0.05  0.4   2.5   6     0.02  0.02  0.05  0.02  cytotoxic
GZMH      0.01  0.2   2.5   4     0.01  0.01  0.01  0.01  cytotoxic
CST7      0.05  0.3   2.2   4     0.02  0.05  0.05  0.02  cytotoxic
KLRG1     0.02  0.4   1.2   1.2   0.01  0.01  0.01  0.01  cytotoxic
NKG7      0.02  0.4   4     20    0     0     0.3   0.05  cytotoxic
GNLY      0.05  0.3   3     40    0.05  0.05  0.1   0.05  cytotoxic
GZMB      0.01  0.02  0.4   6     0.02  0.01  0.02  0.1   cytotoxic
PRF1      0.02  0.1   0.8   4     0.01  0.01  0.02  0.01  cytotoxic
KLRD1     0.01  0.05  0.5   2.5   0.01  0.01  0.02  0.01  cytotoxic
CTSW      0.1   0.4   2     5     0.02  0.02  0.05  0.02  cytotoxic
FGFBP2    0.01  0.05  0.6   4     0.01  0.01  0.01  0.01  cytotoxic
CCL4      0.05  0.2   1.5   3     0.05  0.1   0.1   0.05  cytotoxic
KLRF1     0.01  0.01  0.05  1.5   0.01  0.01  0.01  0.01  nk
SPON2     0.01  0.02  0.3   2     0.01  0.01  0.02  0.01  nk
XCL2      0.01  0.02  0.2   2     0.01  0.01  0.01  0.01  nk
XCL1      0.01  0.02  0.1   1     0.01  0.01  0.01  0.01  nk
KLRC1     0.01  0.01  0.1   1     0.01  0.01  0.01  0.01  nk
HOPX      0.05  0.3   0.6   2     0.01  0.01  0.01  0.01  nk
IL2RB     0.05  0.2   0.6   2     0.01  0.01  0.01  0.01  nk
CD7       0.8   0.8   1.5   3     0.02  0.02  0.02  0.02  nk
CLIC3     0.01  0.05  0.4   2     0.01  0.01  0.01  0.01  nk
AKR1C3    0.01  0.02  0.05  1.5   0.01  0.01  0.01  0.01  nk
SH2D1B    0.01  0.01  0.05  0.8   0.01  0.01  0.01  0.01  nk
TYROBP    0.2   0.3   0.6   10    0.2   12    15    8     nk_myeloid
FCER1G    0.1   0.1   0.3   6     0.2   10    15    6     nk_myeloid
FCGR3A    0.01  0.02  0.1   3     0.01  0.15  6     0.1   nk_myeloid
MS4A1     0.02  0.02  0.02  0.02  4     0.01  0.01  0.05  b
CD79A     0.02  0.02  0.02  0.02  8     0.02  0.02  0.05  b
CD79B     0.05  0.05  0.1   0.1   4     0.05  0.05  0.1   b
BANK1     0.01  0.01  0.01  0.01  1.2   0.01  0.01  0.01  b
CD19      0.01  0.01  0.01  0.01  0.6   0.01  0.01  0.01  b
CD22      0.01  0.01  0.01  0.01  0.8   0.01  0.01  0.02  b
TCL1A     0.01  0.01  0.01  0.01  5     0.01  0.01  0.05  b
VPREB3    0.01  0.01  0.01  0.01  1.2   0.01  0.01  0.01  b
FCRLA     0.01  0.01  0.01  0.01  1     0.01  0.01  0.01  b
CD37      1     1.2   1     0.6   6     0.4   0.6   0.5   b
FCER2     0.01  0.01  0.01  0.01  0.8   0.01  0.01  0.01  b
PAX5      0.01  0.01  0.01  0.01  0.5   0.01  0.01  0.01  b
BLK       0.01  0.01  0.01  0.01  0.8   0.01  0.01  0.02  b
SPIB      0.01  0.01  0.01  0.01  0.6   0.01  0.01  0.2   b
CD72      0.01  0.01  0.01  0.01  0.6   0.02  0.02  0.05  b
TNFRSF13C 0.01  0.01  0.01  0.01  0.5   0.01  0.01  0.01  b
CXCR5     0.01  0.05  0.01  0.01  0.6   0.01  0.01  0.01  b
HVCN1     0.01  0.01  0.01  0.01  0.8   0.01  0.01  0.05  b
BLNK      0.01  0.01  0.01  0.01  0.5   0.01  0.01  0.05  b
CD24      0.01  0.01  0.01  0.01  0.4   0.01  0.01  0.01  b
FCRL1     0.01  0.01  0.01  0.01  0.5   0.01  0.01  0.01  b
CD74      3     4     3     2     40    25    25    60    mhc2
HLA-DRA   0.3   0.5   0.4   0.3   25    12    15    50    mhc2
HLA-DRB1  0.2   0.4   0.3   0.2   12    6     8     25    mhc2
HLA-DQA1  0.05  0.05  0.05  0.05  3     1     1.5   8     mhc2
HLA-DQB1  0.05  0.05  0.05  0.05  3     1     1.5   6     mhc2
HLA-DPA1  0.2   0.3   0.2   0.2   6     4     5     15    mhc2
HLA-DPB1  0.3   0.4   0.3   0.3   8     6     7     20    mhc2
HLA-DMA   0.05  0.05  0.05  0.05  1.5   0.6   0.6   2.5   mhc2
HLA-DMB   0.05  0.05  0.05  0.05  1.5   0.4   0.4   2     mhc2
LYZ       0.4   0.4   0.4   0.4   1     60    18    35    myeloid
CST3      0.1   0.1   0.1   0.1   0.2   12    10    25    myeloid
AIF1      0.1   0.1   0.1   0.2   0.1   6     12    5     myeloid
LST1      0.1   0.2   0.2   0.3   0.1   4     10    2     myeloid
CTSS      0.3   0.3   0.3   0.5   0.4   8     8     5     myeloid
PSAP      0.5   0.6   0.6   0.6   0.5   8     8     6     myeloid
SAT1      0.6   0.8   0.8   1     0.6   8     10    4     myeloid
LGALS3    0.05  0.1   0.1   0.1   0.05  4     3     2     myeloid
CD68      0.05  0.05  0.05  0.05  0.05  3     3     1     myeloid
SPI1      0.05  0.05  0.05  0.05  0.4   2     2.5   1     myeloid
CFD       0.02  0.02  0.02  0.02  0.02  4     3     1     myeloid
SERPINA1  0.01  0.01  0.01  0.01  0.01  2.5   4     0.5   myeloid
CD33      0.01  0.01  0.01  0.01  0.01  0.8   0.3   0.5   myeloid
CYBB      0.05  0.05  0.05  0.05  0.2   2     2     1.5   myeloid
NCF2      0.02  0.02  0.02  0.05  0.05  2     1.5   1     myeloid
TYMP      0.1   0.1   0.1   0.1   0.1   3     2     1.5   myeloid
IFI30     0.05  0.05  0.05  0.05  0.5   4     4     4     myeloid
GRN       0.05  0.05  0.05  0.05  0.05  2     2     1     myeloid
CPVL      0.02  0.02  0.02  0.02  0.02  1.5   0.6   1.5   myeloid
CTSZ      0.05  0.05  0.05  0.05  0.05  1.5   1.5   1     myeloid
CD86      0.01  0.01  0.01  0.01  0.05  0.8   0.8   1.2   myeloid
CLEC7A    0.01  0.01  0.01  0.01  0.01  1     0.4   0.6   myeloid
CD14      0.02  0.02  0.02  0.02  0.02  2.5   0.3   0.2   cd14
S100A8    0.3   0.3   0.3   0.3   0.3   25    3     3     cd14
S100A9    0.5   0.5   0.5   0.5   0.5   40    8     5     cd14
S100A12   0.02  0.02  0.02  0.02  0.02  4     0.3   0.3   cd14
VCAN      0.01  0.01  0.01  0.01  0.01  3     0.3   0.4   cd14
FCN1      0.02  0.02  0.02  0.02  0.02  8     4     1     cd14
LGALS2    0.02  0.02  0.02  0.02  0.02  3     0.3   2     cd14
CSF3R     0.01  0.01  0.01  0.01  0.01  1.5   0.3   0.1   cd14
MNDA      0.02  0.02  0.02  0.02  0.05  3     1.5   0.5   cd14
CD36      0.01  0.01  0.01  0.01  0.01  2     0.5   0.3   cd14
FPR1      0.01  0.01  0.01  0.01  0.01  1     0.2   0.05  cd14
CEBPD     0.05  0.05  0.05  0.05  0.05  2     0.6   0.4   cd14
CEBPB     0.05  0.05  0.05  0.05  0.05  2.5   1.5   0.6   cd14
MS4A6A    0.01  0.01  0.01  0.01  0.01  2     1     1.2   cd14
CD163     0.01  0.01  0.01  0.01  0.01  0.6   0.1   0.05  cd14
RBP7      0.01  0.01  0.01  0.01  0.01  0.8   0.05  0.05  cd14
CYP1B1    0.01  0.01  0.01  0.01  0.01  1     0.1   0.1   cd14
CD93      0.01  0.01  0.01  0.01  0.01  0.8   0.1   0.05  cd14
PLBD1     0.01  0.01  0.01  0.01  0.01  1     0.1   0.3   cd14
CLEC4E    0.01  0.01  0.01  0.01  0.01  0.8   0.1   0.05  cd14
NCF1      0.02  0.02  0.02  0.02  0.05  2     1.2   0.8   cd14
ITGAM     0.01  0.01  0.01  0.05  0.01  0.6   0.2   0.1   cd14
MS4A7     0.01  0.01  0.01  0.01  0.01  0.7   5     0.1   fcgr3a
CDKN1C    0.01  0.01  0.01  0.01  0.01  0.05  3     0.05  fcgr3a
IFITM2    0.5   0.6   1     2     0.4   2     6     1     fcgr3a
IFITM3    0.3   0.3   0.3   0.5   0.2   2     6     1     fcgr3a
RHOC      0.05  0.05  0.1   0.3   0.05  0.4   3     0.3   fcgr3a
HES4      0.01  0.01  0.01  0.05  0.01  0.05  1.5   0.05  fcgr3a
CKB       0.01  0.01  0.01  0.01  0.01  0.05  1     0.05  fcgr3a
CX3CR1    0.01  0.01  0.1   0.5   0.01  0.1   1.5   0.05  fcgr3a
SIGLEC10  0.01  0.01  0.01  0.01  0.01  0.05  1     0.05  fcgr3a
TCF7L2    0.02  0.02  0.02  0.02  0.02  0.2   1.5   0.1   fcgr3a
VMO1      0.01  0.01  0.01  0.01  0.01  0.05  1.5   0.05  fcgr3a
CSF1R     0.01  0.01  0.01  0.01  0.01  0.3   2     0.4   fcgr3a
MTSS1     0.05  0.05  0.05  0.05  0.05  0.3   1.5   0.2   fcgr3a
LYPD2     0.01  0.01  0.01  0.01  0.01  0.01  0.8   0.01  fcgr3a
C1QA      0.01  0.01  0.01  0.01  0.01  0.05  1     0.05  fcgr3a
FCER1A    0.01  0.01  0.01  0.01  0.01  0.05  0.02  6     dc
CLEC10A   0.01  0.01  0.01  0.01  0.01  0.1   0.05  2.5   dc
CD1C      0.01  0.01  0.01  0.01  0.05  0.02  0.02  1.5   dc
ENHO      0.01  0.01  0.01  0.01  0.01  0.02  0.02  1.2   dc
PKIB      0.01  0.01  0.01  0.01  0.01  0.02  0.02  1     dc
FLT3      0.01  0.01  0.01  0.01  0.01  0.02  0.02  0.8   dc
NDRG2     0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.8   dc
CD1E      0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.8   dc
IRF8      0.01  0.01  0.01  0.05  0.3   0.3   0.6   1.5   dc
CLEC4C    0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.2   dc
LILRA4    0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.3   dc
PLD4      0.01  0.01  0.01  0.01  0.05  0.02  0.05  0.5   dc
"""

# Prior-knowledge marker lists (literature / Seurat PBMC 3k tutorial) - written to data/scrna_markers.tsv
MARKERS = {
    "Naive CD4 T": ["CCR7", "LEF1", "TCF7", "SELL", "IL7R"],
    "Memory CD4 T": ["S100A4", "GPR183", "ANXA1", "KLRB1", "IL7R"],
    "CD8 T": ["CD8A", "CD8B", "CCL5", "GZMK"],
    "NK": ["GNLY", "NKG7", "GZMB", "PRF1", "KLRD1"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7", "CDKN1C", "LST1"],
    "DC": ["FCER1A", "CST3", "CLEC10A"],
}

N_RIBO = 40                                   # ribosomal protein genes (always detected)
RIBO_FACTOR = dict(zip(TYPES, [1.25, 1.10, 1.10, 0.85, 1.0, 0.60, 0.55, 0.75]))
LIB_FACTOR = dict(zip(TYPES, [0.90, 1.00, 1.00, 1.00, 1.05, 1.35, 1.30, 1.40]))
LIB_SDLOG, LIB_CLIP = 0.36, (550, 7600)
BG_SHARE = 0.015                              # share of UMIs that go to the background genes
BG_SDLOG, TYPE_SDLOG, DONOR_SDLOG = 1.4, 0.25, 0.08
PHI = {"hk": 0.04, "ribo": 0.04, "marker": 0.30, "background": 0.20}
STATE = {"Naive CD4 T": (1.3, 10.0), "Memory CD4 T": (10.0, 1.3),       # Beta(a, b) of the memory score
         "CD14 Mono": (1.3, 7.0), "FCGR3A Mono": (7.0, 1.3)}         # Beta(a, b) of the FCGR3A score
EFFECTOR = (2.4, 1.8)        # CD8 effector score s: cytotoxic genes x exp(2.4 (s - 0.5)), naive genes x exp(-1.8 (s - 0.5))
NAIVE_LIKE_EXTRA = ["LTB", "IL7R", "CD27", "CD28"]          # also lost by effector CD8 T cells (TEMRA-like)


def gene_panel(rng):
    mane = pd.read_csv(HERE / "ref" / "genes_hg38_mane_select.tsv.gz", sep="\t")
    table = pd.read_csv(io.StringIO(TABLE), sep=r"\s+")
    missing = sorted(set(table.gene) - set(mane.gene))
    assert not missing, f"not in MANE Select: {missing}"
    ribo = mane.loc[mane.gene.str.fullmatch(r"RP[LS]\d+[AX]?|RPLP[012]") & (mane.chrom != "chrY"), "gene"]
    ribo = sorted(ribo)
    ribo = list(rng.choice(ribo, N_RIBO, replace=False))
    used = set(table.gene) | set(ribo)
    pool = mane.loc[~mane.gene.isin(used) & (mane.chrom != "chrY") & ~mane.gene.str.startswith(("RPL", "RPS")),
                    "gene"].sort_values().tolist()
    n_bg = N_GENES - len(table) - len(ribo)
    background = list(rng.choice(pool, n_bg, replace=False))
    return mane, table, sorted(ribo), sorted(background)


def profiles(rng, table, ribo, background):
    """Expected share of UMIs per gene (columns) for each cell type (rows) before cell-level effects."""
    E = pd.DataFrame(0.0, index=TYPES, columns=list(table.gene) + ribo + background)
    E.loc[:, table.gene] = table[COLS].to_numpy().T
    rib = np.exp(rng.normal(np.log(17.0), 0.45, len(ribo)))
    for t in TYPES:
        E.loc[t, ribo] = rib * RIBO_FACTOR[t]
    bg = np.exp(rng.normal(0.0, BG_SDLOG, len(background)))
    bg_type = np.exp(rng.normal(0.0, TYPE_SDLOG, (len(TYPES), len(background))))  # many small differences
    for i, t in enumerate(TYPES):
        named = E.loc[t].sum()
        w = bg * bg_type[i]
        E.loc[t, background] = w / w.sum() * named * BG_SHARE / (1 - BG_SHARE)
    return E.div(E.sum(axis=1), axis=0)          # each row sums to 1


def main():
    rng = np.random.default_rng(SEED)
    mane, table, ribo, background = gene_panel(rng)
    P = profiles(rng, table, ribo, background)
    group = dict(zip(table.gene, table.group))
    kind = np.array(["hk" if group.get(g) == "hk" else "marker" if g in group else
                     "ribo" if g in ribo else "background" for g in P.columns])
    phi = np.array([PHI[k] for k in kind]) * np.exp(rng.normal(0, 0.3, len(kind)))

    cells = pd.read_csv(ROOT / "data" / "truth" / "scrna_truth.tsv", sep="\t")
    assert len(cells) == 2940 and set(cells.cell_type) == set(TYPES)
    cells = cells.sort_values("cell_id").reset_index(drop=True)      # output order: by barcode
    n, G = len(cells), P.shape[1]
    logP = np.log(np.maximum(P.to_numpy(), 1e-300))
    zero = P.to_numpy() == 0
    tix = cells.cell_type.map({t: i for i, t in enumerate(TYPES)}).to_numpy()

    # per-cell expression profile: own type, or a mix of two types along a continuum
    logp = logP[tix].copy()
    pairs = {"Naive CD4 T": "Memory CD4 T", "Memory CD4 T": "Naive CD4 T",
             "CD14 Mono": "FCGR3A Mono", "FCGR3A Mono": "CD14 Mono"}
    state = np.zeros(n)
    for t, (a, b) in STATE.items():
        m = (cells.cell_type == t).to_numpy()
        s = rng.beta(a, b, m.sum())                           # 0 = naive / CD14 end, 1 = memory / FCGR3A end
        state[m] = s
        lo, hi = (TYPES.index(t), TYPES.index(pairs[t])) if t in ("Naive CD4 T", "CD14 Mono") else \
                 (TYPES.index(pairs[t]), TYPES.index(t))
        both = ~(zero[lo] | zero[hi])
        mix = np.outer(1 - s, logP[lo]) + np.outer(s, logP[hi])
        logp[np.ix_(m, both)] = mix[:, both]
    cd8 = (cells.cell_type == "CD8 T").to_numpy()
    eff = rng.beta(2.0, 2.0, cd8.sum()) - 0.5                # CD8 effector score (centred)
    state[cd8] = eff + 0.5
    cyto = np.isin(P.columns, table.gene[table.group == "cytotoxic"])
    naive = np.isin(P.columns, list(table.gene[table.group == "naive"]) + NAIVE_LIKE_EXTRA)
    logp[np.ix_(cd8, cyto)] += EFFECTOR[0] * eff[:, None]
    logp[np.ix_(cd8, naive)] -= EFFECTOR[1] * eff[:, None]

    donor_eff = rng.normal(0, DONOR_SDLOG, (4, G))
    dix = cells.donor.str[1].astype(int).to_numpy() - 1
    logp += donor_eff[dix]
    p = np.exp(logp)
    p[zero[tix]] = 0.0                                        # exact zeros stay zero (NKG7 in B, CD14 Mono)
    p /= p.sum(axis=1, keepdims=True)

    lib = MEDIAN_UMIS * np.exp(rng.normal(0, LIB_SDLOG, n)) * cells.cell_type.map(LIB_FACTOR).to_numpy()
    lib = np.clip(lib / np.median(lib) * MEDIAN_UMIS, *LIB_CLIP)
    mu = lib[:, None] * p
    lam = rng.gamma(shape=1.0 / phi, scale=phi * mu)       # gamma-Poisson (negative binomial)
    counts = rng.poisson(lam).astype(np.int32)
    return dict(cells=cells, P=P, table=table, ribo=ribo, background=background, mane=mane,
                counts=counts, state=state, phi=phi, kind=kind)


def fallback_umap(counts):
    """The notebook's pipeline with the same settings: genes in >= 3 cells, 10,000 + log1p,
    1,000 HVGs, 30 PCs, UMAP (15 neighbours, min_dist 0.3, seed 0)."""
    keep = (counts > 0).sum(axis=0) >= 3
    lognorm, scaled = H.normalise_log1p(counts[:, keep])
    hvg, _ = H.highly_variable_genes(scaled, n_top=1000)
    pcs, _ = H.run_pca(lognorm[:, hvg], n_pcs=30)
    emb, _ = H.run_umap(pcs, n_neighbors=15, min_dist=0.3, seed=0)
    return emb


def write(sim):
    data = ROOT / "data"
    mane, cells, counts = sim["mane"], sim["cells"], sim["counts"]
    genes = np.array(sim["P"].columns)
    pos = mane.set_index("gene").loc[genes]
    chrom_rank = pos["chrom"].str.replace("chr", "").map(lambda c: int(c) if c.isdigit() else 23)
    order = np.lexsort((pos["start0"].to_numpy(), chrom_rank.to_numpy()))   # genomic order, like 10x
    counts, genes = counts[:, order], genes[order]
    pd.DataFrame({"gene": genes, "gene_id": pos["gene_id"].to_numpy()[order]}).to_csv(
        data / "scrna_genes.tsv", sep="\t", index=False)
    sp.save_npz(data / "scrna_counts.npz", sp.csr_matrix(counts, dtype=np.int32), compressed=True)
    emb = np.round(fallback_umap(counts), 4) + 0.0                         # + 0.0 turns -0.0 into 0.0
    out = cells[["cell_id", "donor"]].assign(n_counts=counts.sum(axis=1), n_genes=(counts > 0).sum(axis=1),
                                             umap1=emb[:, 0], umap2=emb[:, 1])
    out.to_csv(data / "scrna_cells.tsv", sep="\t", index=False)
    pd.DataFrame({"cell_type": list(MARKERS), "markers": [",".join(v) for v in MARKERS.values()]}).to_csv(
        data / "scrna_markers.tsv", sep="\t", index=False)
    table = sim["table"]
    params = {
        "description": "Simulation parameters of data/scrna_counts.npz (Practice B4). SIMULATED data.",
        "seed": SEED, "n_cells": int(counts.shape[0]), "n_genes": int(counts.shape[1]),
        "model": "counts ~ Poisson(mu * Gamma(1/phi, phi)); mu = library size x cell-type profile x donor x state",
        "library_size": {"median_umis": MEDIAN_UMIS, "sdlog": LIB_SDLOG, "clip": list(LIB_CLIP),
                         "cell_type_factor": LIB_FACTOR},
        "background_share_of_umis": BG_SHARE, "background_sdlog": BG_SDLOG,
        "type_sdlog_background": TYPE_SDLOG, "donor_sdlog": DONOR_SDLOG, "phi_by_class": PHI,
        "n_ribosomal_genes": N_RIBO, "ribosomal_factor": RIBO_FACTOR,
        "continua_beta": {k: list(v) for k, v in STATE.items()},
        "cd8_effector": f"score s ~ Beta(2, 2); cytotoxic genes x exp({EFFECTOR[0]} (s - 0.5)); naive genes and "
                        f"{', '.join(NAIVE_LIKE_EXTRA)} x exp(-{EFFECTOR[1]} (s - 0.5))",
        "gene_class_counts": {k: int((sim["kind"] == k).sum()) for k in ["hk", "ribo", "marker", "background"]},
        "marker_table_umis_per_cell": {r["gene"]: {"group": r["group"], **{t: float(r[c]) for t, c in zip(TYPES, COLS)}}
                                       for r in table.to_dict(orient="records")},
        "ribosomal_genes": sim["ribo"],
    }
    (data / "truth" / "scrna_truth.json").write_text(json.dumps(params, indent=1) + "\n")
    return out


if __name__ == "__main__":
    sim = main()
    out = write(sim)
    c = sim["counts"]
    print(f"cells {c.shape[0]:,}  genes {c.shape[1]:,}  median UMIs {np.median(c.sum(1)):,.0f}  "
          f"range {c.sum(1).min():,}-{c.sum(1).max():,}  zeros {100 * (c == 0).mean():.2f}%  "
          f"median genes/cell {np.median((c > 0).sum(1)):.0f}")
