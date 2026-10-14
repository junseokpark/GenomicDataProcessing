#!/usr/bin/env python
"""Write the truth set of planted CNVs (hg19) and look up matching public ClinVar records.

Outputs
  data/truth/planted_cnvs_hg19.bed   BED (0-based start), name = cnv_id, score = copy number
  data/truth/planted_cnvs_hg19.tsv   one row per planted event with source record, expected CN and notes
ClinVar matches are exact-coordinate matches in the UCSC hg19 clinvarCnv track
(data/tracks_hg19/clinvarCnv_hg19_windows.tsv.gz, downloaded by fetch_public_data.py).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402


def main():
    cv = pd.read_csv("data/tracks_hg19/clinvarCnv_hg19_windows.tsv.gz", sep="\t")
    rows = []
    for ev in C.PLANTED:
        want = "copy number loss" if ev["type"] == "DEL" else "copy number gain"
        m = cv[(cv.chrom == ev["chrom"]) & (cv.chromStart + 1 == ev["start_1based"]) & (cv.chromEnd == ev["end_1based"])]
        m_same = m[m["type"].str.lower().isin([want, "deletion" if ev["type"] == "DEL" else "duplication"])]
        matches = "; ".join(f"{r['name']} ({r['type']}, {r.clinSign})" for _, r in m_same.iterrows()) or "none (exact coordinates)"
        rows.append(dict(cnv_id=ev["cnv_id"], sample=ev["sample"], chrom=ev["chrom"],
                         start0=ev["start_1based"] - 1, end=ev["end_1based"],
                         start_1based=ev["start_1based"], end_1based=ev["end_1based"],
                         size_bp=ev["end_1based"] - ev["start_1based"] + 1, type=ev["type"],
                         expected_cn=ev["cn"], expected_log2=(round(float(np.log2(ev["cn"] / 2)), 3) if ev["cn"] > 0 else "-inf (CN 0)"),
                         zygosity=ev["zygosity"], source_record=ev["source"],
                         clinvar_hg19_exact_matches=matches, public_grch38=ev["public_grch38"],
                         notes=ev["note"] + " SYNTHETIC sample; coordinates from public records."))
    df = pd.DataFrame(rows)
    os.makedirs("data/truth", exist_ok=True)
    df.to_csv("data/truth/planted_cnvs_hg19.tsv", sep="\t", index=False)
    df[["chrom", "start0", "end", "cnv_id", "expected_cn"]].assign(strand=".", sample=df["sample"], type=df["type"]) \
        .to_csv("data/truth/planted_cnvs_hg19.bed", sep="\t", index=False, header=False)
    print(df[["cnv_id", "chrom", "start_1based", "end_1based", "expected_cn", "clinvar_hg19_exact_matches"]].to_string(index=False))


if __name__ == "__main__":
    main()
