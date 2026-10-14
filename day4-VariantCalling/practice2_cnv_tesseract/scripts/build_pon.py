#!/usr/bin/env python
"""Build the panel of normals (PON): 5-kb bin fragment counts for the five SYNTHETIC controls.

C1 is read from data/fastq/ (the same FASTQ the students process). C2-C5 are simulated into a
temporary folder with simulate_dna.py (same seeds/model), mapped with exactly the commands used in
the notebook (minimap2 -ax sr | samtools fixmate -m | sort | markdup), counted with
cnvlib.count_fragments_in_bins (MAPQ >= 20, proper pairs, primary, non-duplicate, fragment
midpoint), and deleted again. Output: data/pon/pon_bincounts_hg19_5kb.tsv.gz (+ pon_build_info.json)
"""
import json
import os
import shutil
import subprocess
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cnvlib  # noqa: E402
import config as C  # noqa: E402
import simulate_dna as S  # noqa: E402

REF = "data/ref/hg19_cnv_windows.fa.gz"


def main(threads=4):
    tmp = "data/_pon_tmp"
    os.makedirs(tmp, exist_ok=True)
    os.makedirs("data/pon", exist_ok=True)
    W = S.load_windows(REF)
    bins = cnvlib.make_bins("data/ref/windows_hg19.bed", C.BIN_SIZE, REF)
    info = {"controls": {}, "bin_size": C.BIN_SIZE, "min_mapq": cnvlib.MIN_MAPQ}
    for s in C.CONTROLS:
        fq = "data/fastq" if os.path.exists(f"data/fastq/{s}_R1.fastq.gz") else tmp
        if fq == tmp:
            info["controls"][s] = S.simulate_sample(s, W, tmp, C.MEAN_DEPTH * C.DEPTH_FACTOR[s])
        bam = os.path.join(tmp, f"{s}.bam")
        cnvlib.map_sample(s, f"{fq}/{s}_R1.fastq.gz", f"{fq}/{s}_R2.fastq.gz", REF, bam, threads=threads, tmpdir=tmp)
        bins[s] = cnvlib.count_fragments_in_bins(bam, bins)
        info["controls"].setdefault(s, {})["fragments_counted"] = int(bins[s].sum())
        print(s, int(bins[s].sum()), flush=True)
    bins.to_csv("data/pon/pon_bincounts_hg19_5kb.tsv.gz", sep="\t", index=False)
    info["tool_versions"] = {
        "minimap2": subprocess.run(["minimap2", "--version"], capture_output=True, text=True).stdout.strip(),
        "samtools": subprocess.run(["samtools", "--version"], capture_output=True, text=True).stdout.split("\n")[0]}
    info["note"] = "SYNTHETIC controls. C2 carries a heterozygous UGT2B17 deletion (population variation)."
    json.dump(info, open("data/pon/pon_build_info.json", "w"), indent=1)
    shutil.rmtree(tmp)


if __name__ == "__main__":
    main()
