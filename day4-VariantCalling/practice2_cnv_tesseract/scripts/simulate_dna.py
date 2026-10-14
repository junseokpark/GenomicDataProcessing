#!/usr/bin/env python
"""Simulate paired-end short-read DNA sequencing for the SYNTHETIC Practice 2 cohort.

Everything produced here is synthetic (no real patient data). Reads are sampled from REAL hg19
sequence (data/ref/hg19_cnv_windows.fa.gz, downloaded by fetch_public_data.py) after planting the
CNVs listed in config.PLANTED into per-sample diploid haplotypes:
  * heterozygous deletion -> the interval is removed from haplotype 1
  * homozygous deletion   -> removed from both haplotypes (CN 0)
  * tandem duplication    -> a second copy is inserted right after the interval on haplotype 1

Read model (identical for every sample):
  * fragment length ~ Normal(350, 50), truncated to 200-600 bp; 2 x 100 bp reads (FR orientation)
  * fragments containing N are never produced (gaps get no reads)
  * mild GC bias: a fragment with GC fraction g is kept with probability exp(-0.5*((g-0.45)/0.15)^2)
    (same curve for every sample, so the ratio to the controls removes it)
  * binned base qualities (2/12/23/37) that drop towards the 3' end; sequencing errors are drawn
    with probability 10^(-Q/10) (Q2 bases are written as N, as Illumina does)
  * 3% of fragments get one duplicate copy (PCR/optical duplicate, independent sequencing errors)
  * read order is shuffled

Seed: numpy default_rng([20261015, sample_index]) with sample_index from config.SAMPLES order.

Depth: config.MEAN_DEPTH (3.0x) x config.DEPTH_FACTOR[sample] (0.88-1.15).

Usage: python scripts/simulate_dna.py --samples P1 P2 P3 C1 --outdir data/fastq [--depth 3.0]
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402

LUT = np.frombuffer(b"ACGTN", dtype=np.uint8)
COMP = np.array([3, 2, 1, 0, 4], dtype=np.uint8)
ENC = np.full(256, 4, dtype=np.uint8)
for i, b in enumerate(b"ACGT"):
    ENC[b] = i
    ENC[ord(chr(b).lower())] = i


def gc_keep_prob(gc):
    return np.exp(-0.5 * ((gc - 0.45) / 0.15) ** 2)


def load_windows(fasta):
    fa = pysam.FastaFile(fasta)
    out = []
    for chrom, s0, e, name in C.windows_bed():
        seq = np.frombuffer(fa.fetch(chrom, s0, e).encode(), dtype=np.uint8)
        out.append((chrom, s0, e, ENC[seq]))
    return out


def sample_haplotypes(sample, windows):
    """Return list of haplotype code arrays (2 per window) with the sample's CNVs planted."""
    events = [ev for ev in C.PLANTED if ev["sample"] == sample]
    haps = []
    for chrom, w0, w1, ref in windows:
        h = [ref.copy(), ref.copy()]
        for ev in events:
            if ev["chrom"] != chrom:
                continue
            s, e = ev["start_1based"] - 1 - w0, ev["end_1based"] - w0   # window-local, half-open
            assert 0 <= s < e <= len(ref), ev
            targets = [0, 1] if ev["zygosity"] == "homozygous" else [0]
            for t in targets:
                if ev["type"] == "DEL":
                    h[t] = np.concatenate([h[t][:s], h[t][e:]])
                elif ev["type"] == "DUP":
                    h[t] = np.concatenate([h[t][:e], h[t][s:e], h[t][e:]])
        haps += h
    return haps


def sample_fragments(hap, n_target, rng):
    """GC-biased rejection sampling of fragment (start, length) on one haplotype."""
    L = len(hap)
    is_gc = ((hap == 1) | (hap == 2)).astype(np.int32)
    is_n = (hap == 4).astype(np.int32)
    cgc = np.concatenate([[0], np.cumsum(is_gc)])
    cn = np.concatenate([[0], np.cumsum(is_n)])
    starts, lens = [], []
    have = 0
    while have < n_target:
        m = int((n_target - have) * 1.8) + 100
        ln = np.clip(np.round(rng.normal(C.FRAG_MEAN, C.FRAG_SD, m)), C.FRAG_MIN, C.FRAG_MAX).astype(np.int64)
        st = (rng.random(m) * (L - ln)).astype(np.int64)
        ok = (cn[st + ln] - cn[st]) == 0
        gc = (cgc[st + ln] - cgc[st]) / ln
        ok &= rng.random(m) < gc_keep_prob(gc)
        starts.append(st[ok]); lens.append(ln[ok]); have += int(ok.sum())
    st = np.concatenate(starts)[:n_target]
    ln = np.concatenate(lens)[:n_target]
    return st, ln


def fragment_reads(hap, st, ln, rng):
    ar = np.arange(C.READ_LEN)
    A = hap[st[:, None] + ar]                                        # left end, + strand
    B = COMP[hap[(st + ln - C.READ_LEN)[:, None] + ar]][:, ::-1]     # right end, reverse complement
    plus = rng.random(len(st)) < 0.5
    R1 = np.where(plus[:, None], A, B)
    R2 = np.where(plus[:, None], B, A)
    return R1, R2


def add_quals_and_errors(R, rng, read_number):
    n, L = R.shape
    pos = np.arange(L) / (L - 1)
    p37 = (0.93 - 0.18 * pos ** 2) if read_number == 1 else (0.90 - 0.25 * pos ** 2)
    rest = 1 - p37
    th = np.stack([p37, p37 + rest * 0.75, p37 + rest * 0.98])        # cumulative: Q37, Q23, Q12, Q2
    u = rng.random((n, L))
    qidx = (u >= th[0]).astype(np.int8) + (u >= th[1]) + (u >= th[2])  # 0:Q37 1:Q23 2:Q12 3:Q2
    qvals = np.array([37, 23, 12, 2], dtype=np.int64)[qidx]
    perr = 10.0 ** (-qvals / 10.0)
    err = (rng.random((n, L)) < perr) & (qidx < 3)
    R = R.copy()
    R[err] = (R[err] + rng.integers(1, 4, int(err.sum()))) % 4
    R[qidx == 3] = 4                                                   # Q2 -> N
    R[R > 4] = 4
    return R, (qvals + 33).astype(np.uint8)


def write_fastq(path, names, S, Q):
    seqs = LUT[S]
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0) as fh:
        buf = []
        for i, nm in enumerate(names):
            buf.append(b"@" + nm + b"\n" + seqs[i].tobytes() + b"\n+\n" + Q[i].tobytes() + b"\n")
            if len(buf) == 20000:
                fh.write(b"".join(buf)); buf = []
        fh.write(b"".join(buf))


def simulate_sample(sample, windows, outdir, depth):
    idx = [s["sample"] for s in C.SAMPLES].index(sample)
    rng = np.random.default_rng([C.SEED, idx])
    haps = sample_haplotypes(sample, windows)
    R1s, R2s = [], []
    for hap in haps:
        n_target = int(round((depth / 2.0) * len(hap) / (2 * C.READ_LEN)))
        st, ln = sample_fragments(hap, n_target, rng)
        r1, r2 = fragment_reads(hap, st, ln, rng)
        R1s.append(r1); R2s.append(r2)
    R1 = np.concatenate(R1s); R2 = np.concatenate(R2s)
    n = len(R1)
    dup = rng.random(n) < C.DUP_RATE
    R1 = np.concatenate([R1, R1[dup]]); R2 = np.concatenate([R2, R2[dup]])
    order = rng.permutation(len(R1))
    R1, R2 = R1[order], R2[order]
    S1, Q1 = add_quals_and_errors(R1, rng, 1)
    S2, Q2 = add_quals_and_errors(R2, rng, 2)
    names = [f"{sample}_{i + 1:07d}".encode() for i in range(len(R1))]
    os.makedirs(outdir, exist_ok=True)
    write_fastq(os.path.join(outdir, f"{sample}_R1.fastq.gz"), names, S1, Q1)
    write_fastq(os.path.join(outdir, f"{sample}_R2.fastq.gz"), names, S2, Q2)
    return {"sample": sample, "fragments": int(n), "duplicate_copies": int(dup.sum()), "read_pairs": int(len(R1))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", default=C.SHIPPED_FASTQ)
    ap.add_argument("--outdir", default="data/fastq")
    ap.add_argument("--depth", type=float, default=C.MEAN_DEPTH, help="base depth; multiplied by config.DEPTH_FACTOR[sample]")
    ap.add_argument("--fasta", default="data/ref/hg19_cnv_windows.fa.gz")
    a = ap.parse_args()
    W = load_windows(a.fasta)
    for s in a.samples:
        print(simulate_sample(s, W, a.outdir, a.depth * C.DEPTH_FACTOR[s]))
