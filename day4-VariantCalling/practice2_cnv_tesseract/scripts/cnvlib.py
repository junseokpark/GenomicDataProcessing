"""Functions shared by the data generator (build_pon.py) and mirrored in the notebook.

The notebook re-implements count_fragments_in_bins() in a visible cell; the CHECK in Step 2
verifies that the student's C1 counts are identical to the C1 column of the panel of normals
built here, so the two implementations must stay identical.
"""
import subprocess

import numpy as np
import pandas as pd
import pysam

MIN_MAPQ = 20


def make_bins(windows_bed, bin_size, fasta_path):
    """5-kb bins tiling each window, with GC fraction and N fraction from the reference."""
    fa = pysam.FastaFile(fasta_path)
    rows = []
    for line in open(windows_bed):
        chrom, s0, e = line.split("\t")[:3]
        s0, e = int(s0), int(e)
        seq = fa.fetch(chrom, s0, e).upper()
        for b in range(s0, e, bin_size):
            sub = seq[b - s0:b - s0 + bin_size]
            acgt = sum(sub.count(x) for x in "ACGT")
            gc = (sub.count("G") + sub.count("C")) / acgt if acgt else np.nan
            rows.append((chrom, b, min(b + bin_size, e), round(gc, 4), round(1 - acgt / len(sub), 4)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gc", "n_frac"])


def count_fragments_in_bins(bam_path, bins, min_mapq=MIN_MAPQ):
    """Count DNA fragments (read pairs) per bin.

    One fragment is counted once, at its midpoint, using the LEFT read of each properly
    paired fragment (template length > 0). Reads are skipped if unmapped, secondary,
    supplementary, QC-fail or marked duplicate, or if MAPQ < min_mapq.
    """
    counts = np.zeros(len(bins), dtype=np.int64)
    with pysam.AlignmentFile(bam_path) as bam:
        for chrom, grp in bins.groupby("chrom", sort=False):
            first = grp.index[0]
            w0, w1 = int(grp.start.min()), int(grp.end.max())
            size = int(grp.end.iloc[0] - grp.start.iloc[0])
            for r in bam.fetch(chrom, w0, w1):
                if (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail
                        or r.is_duplicate or not r.is_proper_pair):
                    continue
                if r.template_length <= 0 or r.mapping_quality < min_mapq:
                    continue
                mid = r.reference_start + r.template_length // 2
                if w0 <= mid < w1:
                    counts[first + (mid - w0) // size] += 1
    return counts


def map_sample(sample, r1, r2, ref, out_bam, threads=4, tmpdir="."):
    """minimap2 -ax sr (+@RG) | samtools fixmate -m | samtools sort | samtools markdup ; index."""
    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}_lib1"
    sorted_bam = out_bam.replace(".bam", ".sorted.bam")
    cmd = (f"minimap2 -ax sr -t {threads} -R '{rg}' {ref} {r1} {r2} 2> {out_bam}.minimap2.log "
           f"| samtools fixmate -m -@ 2 - - "
           f"| samtools sort -@ 2 -T {tmpdir}/{sample}.tmp -o {sorted_bam} - "
           f"&& samtools markdup -@ 2 -s {sorted_bam} {out_bam} 2> {out_bam}.markdup.log "
           f"&& samtools index {out_bam} && rm -f {sorted_bam}")
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
    return out_bam
