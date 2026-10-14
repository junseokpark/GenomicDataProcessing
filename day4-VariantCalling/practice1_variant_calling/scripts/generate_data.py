#!/usr/bin/env python3
"""Rebuild every synthetic input file of Practice 1 from the public GRCh38 download.

SYNTHETIC DATA ONLY. No real person's reads or genotypes are contained in this package.
The family (S01 child, S02 father, S03 mother) is invented for teaching.

What the script does
  1. Reads GRCh38 chr11:5,200,001-5,300,000 (UCSC API response cached in data/provenance/,
     or re-downloaded with --download) and the UCSC hg38 'mane' track for the same window.
  2. Writes a full-length chr11 FASTA (135,086,622 bp) in which every base outside the
     window is N. Files stay small, but every coordinate is a real GRCh38 coordinate.
     The FASTA is bgzip-compressed and indexed with samtools faidx (.fai + .gzi).
  3. Plants truth variants in four parental haplotypes (S02: P1/P2, S03: M1/M2). The child
     S01 inherits P1 and M1 unchanged (no recombination in 100 kb) plus one de novo SNV.
  4. Simulates 2 x 150 bp paired-end reads: fragment length ~N(380, 60) clipped to
     [220, 700]; 30x of unique fragments per sample; NovaSeq-like binned base qualities
     (37/25/11/2 -> F : , #) that fall towards the read end and are slightly worse in R2;
     substitution errors drawn from the stated qualities; ~4% random PCR duplicates.
     Planted problems: a duplicate 'jackpot' in S03, a systematic low-quality error in S02
     (base quality 12: skipped when the pipeline sets bcftools mpileup -Q 13, as the course
     pipeline does; NOTE that bcftools 1.24's own default is -Q 1, which keeps these bases),
     and a 1-kb interval with reduced fragment sampling (~8x) in all samples.
     Two positions are chosen after a dry run of the fragment sampling (same random streams):
     the low-coverage SNV sits where the carrier has ~8x, and the low-quality error sits where
     S02 has ~30-36x. Moving an SNV does not change haplotype lengths, so the dry run and the
     real run sample identical fragments.
  5. Writes truth VCFs, a table of planted events, the PED file, annotation tables,
     provenance and MD5 checksums.

Usage (run from the package folder, inside the conda environment):
    python scripts/generate_data.py              # use the cached UCSC responses (offline)
    python scripts/generate_data.py --download   # re-download from api.genome.ucsc.edu first

Random seed: 20261015 (numpy SeedSequence; one child stream per purpose, see SEED below).
Same seed + same package versions (environment.yml) -> identical data (see data/MD5SUMS.txt).
"""
import argparse
import copy
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import urllib.request
from collections import Counter
from datetime import datetime, timezone

import numpy as np

# ----------------------------------------------------------------------------- constants
SEED = 20261015
CHROM, CHROM_LEN = "chr11", 135_086_622
WIN_START, WIN_END = 5_200_001, 5_300_000            # 1-based, inclusive (GRCh38)
WLEN = WIN_END - WIN_START + 1
UCSC_SEQ_URL = ("https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr11;"
                "start=5200000;end=5300000")
UCSC_MANE_URL = ("https://api.genome.ucsc.edu/getData/track?genome=hg38;track=mane;"
                 "chrom=chr11;start=5200000;end=5300000")
UCSC_MANE_SCHEMA_URL = "https://api.genome.ucsc.edu/list/schema?genome=hg38;track=mane"
MANE_VERSION_URL = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/README_versions.txt"
CHECK_20MER = ("chr11:5226991-5227010", "CAGACTTCTCCTCAGGAGTC")   # verified on UCSC hg38

READ_LEN = 150
FRAG_MEAN, FRAG_SD, FRAG_MIN, FRAG_MAX = 380, 60, 220, 700
TARGET_DEPTH = 30            # mean depth of UNIQUE fragments (duplicates come on top)
DUP_FRACTION = 0.04          # extra copies, as a fraction of unique fragments
QBINS = np.array([37, 25, 11, 2])
LOWCOV = (5_291_001, 5_292_000)   # 1-based inclusive, reduced fragment sampling
LOWCOV_KEEP = 0.27                # keep 27% of fragments whose reads touch it (~8x)
JACKPOT_COPIES = 15               # 1 original + 14 PCR duplicates with identical ends
JACKPOT_R1_OFFSET = 60            # error base sits at R1 position 61
LOWBQ_FRACTION = 0.60             # fraction of S02 read pairs (at that position) showing the error
LOWBQ_Q = 12                      # its base quality (Phred 12 -> '-'); below the pipeline's explicit -Q 13
                                  # (bcftools mpileup's own default is -Q 1, which would keep these bases)
N_BACKGROUND_SNV = 25
SAMPLES = ["S01", "S02", "S03"]
BASES = "ACGT"
B2I = {b: i for i, b in enumerate(BASES)}

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
D = lambda *p: os.path.join(PKG, "data", *p)


def log(*a):
    print("[generate_data]", *a, flush=True)


# ----------------------------------------------------------------------------- download / load
def fetch(url, path):
    with urllib.request.urlopen(url, timeout=60) as r:
        data = r.read()
    with open(path, "wb") as fh:
        fh.write(data)


def load_inputs(download):
    os.makedirs(D("provenance"), exist_ok=True)
    seq_json = D("provenance", "ucsc_hg38_chr11_5200001-5300000_sequence.json")
    mane_json = D("provenance", "ucsc_hg38_mane_chr11_5200001-5300000.json")
    schema_json = D("provenance", "ucsc_hg38_mane_schema.json")
    mane_ver = D("provenance", "ucsc_mane_README_versions.txt")
    if download or not os.path.exists(seq_json):
        log("downloading from the UCSC API ...")
        fetch(UCSC_SEQ_URL, seq_json)
        fetch(UCSC_MANE_URL, mane_json)
        fetch(UCSC_MANE_SCHEMA_URL, schema_json)
        fetch(MANE_VERSION_URL, mane_ver)
    seq = json.load(open(seq_json))
    mane = json.load(open(mane_json))
    schema = json.load(open(schema_json))
    version_txt = open(mane_ver).read()
    dna = seq["dna"]
    assert len(dna) == WLEN, len(dna)
    return dna, mane, schema, version_txt, seq


# ----------------------------------------------------------------------------- sequence helpers
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


class Window:
    """The 100-kb reference window; index i (0-based) <-> GRCh38 position WIN_START + i."""

    def __init__(self, dna, exons):
        self.soft = dna                      # UCSC soft-masked (lower case = repeat)
        self.seq = dna.upper()
        self.lower = np.array([c.islower() for c in dna])
        assert set(self.seq) <= set("ACGT")
        k = 40
        cnt = Counter()
        for s in (self.seq, revcomp(self.seq)):
            for i in range(len(s) - k + 1):
                cnt[s[i:i + k]] += 1
        nonuniq = np.array([cnt[self.seq[i:i + k]] > 1 for i in range(WLEN - k + 1)], int)
        self._cs = np.concatenate([[0], np.cumsum(nonuniq)])
        self._k = k
        self.inexon = np.zeros(WLEN, bool)
        for a, b in exons:                   # 1-based inclusive
            lo, hi = max(a, WIN_START), min(b, WIN_END)
            if lo <= hi:
                self.inexon[lo - WIN_START:hi - WIN_START + 1] = True

    def base(self, pos):
        return self.seq[pos - WIN_START]

    def sub(self, a, b):                     # 1-based inclusive
        return self.seq[a - WIN_START:b - WIN_START + 1]

    def unique_context(self, pos, flank=100):
        """True if every 40-mer within +/- flank bp of pos occurs once in the window (both strands)."""
        i = pos - WIN_START
        a, b = max(0, i - flank), min(WLEN - self._k, i + flank - self._k)
        return self._cs[b + 1] - self._cs[a] == 0

    def near_homopolymer(self, pos, run=4, dist=3):
        s = self.sub(pos - dist - run, pos + dist + run)
        return re.search(r"A{%d,}|C{%d,}|G{%d,}|T{%d,}" % ((run,) * 4), s) is not None

    def eligible_snv(self, pos):
        i = pos - WIN_START
        return (2000 <= i < WLEN - 2000 and not self.lower[i - 20:i + 21].any()
                and self.unique_context(pos) and not self.inexon[i - 20:i + 21].any()
                and not self.near_homopolymer(pos))


def pick_near(win, anchor, used, min_gap=500, extra=lambda p: True):
    """First position >= anchor that passes the SNV rules and is >= min_gap from used sites."""
    p = anchor
    while p < WIN_END - 2000:
        if win.eligible_snv(p) and extra(p) and all(abs(p - u) >= min_gap for u in used):
            return p
        p += 1
    raise RuntimeError(f"no eligible site after {anchor}")


def other_base(ref, rng=None, prefer=None):
    if prefer is not None:
        assert prefer != ref
        return prefer
    choices = [b for b in BASES if b != ref]
    return choices[int(rng.integers(0, 3))]


# ----------------------------------------------------------------------------- annotation
def parse_mane(mane):
    rows = []
    for r in mane["mane"]:
        starts = [int(x) for x in r["chromStarts"].strip(",").split(",")]
        sizes = [int(x) for x in r["blockSizes"].strip(",").split(",")]
        ex = [(r["chromStart"] + s + 1, r["chromStart"] + s + z) for s, z in zip(starts, sizes)]
        rows.append(dict(gene=r["geneName2"], enst=r["name"], refseq=r["ncbiId"],
                         strand=r["strand"], tx_start=r["chromStart"] + 1, tx_end=r["chromEnd"],
                         cds_start=r["thickStart"] + 1, cds_end=r["thickEnd"], exons=ex,
                         status=r["maneStat"]))
    return rows


def write_annotation(genes, retrieved, version_txt, schema):
    os.makedirs(D("annotation"), exist_ok=True)
    mv = dict(l.split("\t") for l in version_txt.strip().splitlines())
    with open(D("annotation", "genes_window_hg38.bed"), "w") as fh:
        fh.write("# MANE transcripts overlapping chr11:5200001-5300000 (GRCh38); UCSC hg38 track 'mane', "
                 f"MANE v{mv.get('MANE Version')}; retrieved {retrieved}. BED is 0-based, half-open.\n")
        fh.write("# columns: chrom start end gene score strand MANE_Ensembl MANE_RefSeq "
                 "(OR51V1 starts 266 bp before the window)\n")
        for g in genes:
            fh.write(f"{CHROM}\t{g['tx_start'] - 1}\t{g['tx_end']}\t{g['gene']}\t0\t{g['strand']}"
                     f"\t{g['enst']}\t{g['refseq']}\n")
    hbb = [g for g in genes if g["gene"] == "HBB"][0]
    assert hbb["refseq"] == "NM_000518.5" and hbb["strand"] == "-"
    exons_tx = sorted(hbb["exons"], key=lambda e: -e[0])   # transcript order (minus strand)
    with open(D("annotation", "HBB_MANE_NM_000518.5_hg38.tsv"), "w") as fh:
        fh.write(f"# HBB MANE Select {hbb['refseq']} / {hbb['enst']} on GRCh38 (UCSC hg38 track 'mane', "
                 f"MANE v{mv.get('MANE Version')}, NCBI {mv.get('NCBI RefSeq Annotation Release')}, "
                 f"Ensembl {mv.get('Ensembl Release')}; track dataTime {schema.get('dataTime')}; "
                 f"retrieved {retrieved})\n")
        fh.write("# coordinates are 1-based inclusive; exons are listed in transcript order "
                 "(5' -> 3'); the gene is on the minus strand\n")
        fh.write("transcript\tensembl\tgene\tchrom\tstrand\texon\texon_start\texon_end\t"
                 "cds_start\tcds_end\tcds_length\n")
        for n, (a, b) in enumerate(exons_tx, 1):
            ca, cb = max(a, hbb["cds_start"]), min(b, hbb["cds_end"])
            cl = cb - ca + 1 if ca <= cb else 0
            fh.write(f"{hbb['refseq']}\t{hbb['enst']}\tHBB\t{CHROM}\t-\t{n}\t{a}\t{b}\t"
                     f"{ca if cl else 'NA'}\t{cb if cl else 'NA'}\t{cl}\n")
    return hbb


# ----------------------------------------------------------------------------- reference FASTA
def write_reference(win):
    os.makedirs(D("ref"), exist_ok=True)
    out = D("ref", "chr11_GRCh38_betaglobin.fa.gz")
    full = bytearray(b"N" * CHROM_LEN)
    full[WIN_START - 1:WIN_END] = win.seq.encode()
    header = (f">{CHROM} GRCh38 chr11 (UCSC hg38); only {WIN_START}-{WIN_END} kept, "
              "all other bases set to N (teaching reference)\n").encode()
    mv = memoryview(full)
    with open(out, "wb") as fh:
        p = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=fh)
        p.stdin.write(header)
        step = 60 * 100_000
        for a in range(0, CHROM_LEN, step):
            chunk = mv[a:min(a + step, CHROM_LEN)]
            p.stdin.write(b"\n".join(chunk[i:i + 60] for i in range(0, len(chunk), 60)) + b"\n")
        p.stdin.close()
        assert p.wait() == 0
    for ext in (".fai", ".gzi"):
        if os.path.exists(out + ext):
            os.remove(out + ext)
    subprocess.run(["samtools", "faidx", out], check=True)
    chk = subprocess.run(["samtools", "faidx", out, CHECK_20MER[0]], capture_output=True,
                         text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "")
    assert chk == CHECK_20MER[1], chk
    return out


# ----------------------------------------------------------------------------- variants
class Variant:
    def __init__(self, eid, pos, ref, alt, haps, cls, context, purpose, expected, vid="."):
        self.eid, self.pos, self.ref, self.alt = eid, pos, ref, alt
        self.haps = set(haps)        # which haplotypes carry ALT (P1, P2, M1, M2, CP)
        self.cls, self.context, self.purpose, self.expected, self.vid = cls, context, purpose, expected, vid

    def gt(self, sample):
        h = {"S01": ("CP", "M1"), "S02": ("P1", "P2"), "S03": ("M1", "M2")}[sample]
        return int(h[0] in self.haps) + int(h[1] in self.haps)


def plant_variants(win, rng):
    used = []
    V = []
    # (a) HBB c.20A>T (HbS) and (b) HBB c.9T>C, fixed checkpoint positions
    assert win.base(5_227_002) == "T" and win.base(5_227_013) == "A"
    V.append(Variant("HBB_c.20A>T", 5_227_002, "T", "A", {"P1", "M1"}, "SNV",
                     "HBB exon 1; NM_000518.5:c.20A>T p.(Glu7Val), HbS (minus-strand gene)",
                     "recessive disease variant: child 1/1, both parents 0/1",
                     "called in all three; recessive filter finds it", vid="rs334"))
    V.append(Variant("HBB_c.9T>C", 5_227_013, "A", "G", {"P1"}, "SNV",
                     "HBB exon 1; NM_000518.5:c.9T>C p.(His3=), synonymous; in cis with c.20A>T on P1",
                     "second coding variant for the naming exercise", "called in S01 and S02",
                     vid="rs713040"))
    used += [5_227_002, 5_227_013]
    # (d1) 1-bp deletion inside a homopolymer run >= 6
    hp_start = 5_239_408
    m = re.match(r"(A+|C+|G+|T+)", win.seq[hp_start - WIN_START:])
    run_len, run_base = len(m.group(1)), m.group(1)[0]
    assert run_len >= 6 and win.base(hp_start - 1) != run_base
    anchor = win.base(hp_start - 1)
    V.append(Variant("DEL1_homopolymer", hp_start - 1, anchor + run_base, anchor, {"P1"}, "DEL",
                     f"{run_base}x{run_len} homopolymer at {CHROM}:{hp_start}-{hp_start + run_len - 1}",
                     "indel in a repeat: shows left-alignment by bcftools norm",
                     "called in S01 and S02; truth_raw stores it right-shifted"))
    V[-1].run = (hp_start, run_len, run_base)
    used.append(hp_start)
    # (d2) 2-bp insertion in a dinucleotide repeat
    di_start = 5_230_686
    m = re.match(r"((?:CA)+)", win.seq[di_start - WIN_START:])
    units = len(m.group(1)) // 2
    assert units >= 5
    anchor = win.base(di_start - 1)
    V.append(Variant("INS2_dinucleotide", di_start - 1, anchor, anchor + "CA", {"M2"}, "INS",
                     f"(CA)x{units} repeat at {CHROM}:{di_start}-{di_start + 2 * units - 1}",
                     "insertion in a repeat; untransmitted maternal allele",
                     "called in S03 only (child 0/0)"))
    used.append(di_start)
    # (d3) 4-bp deletion in unique sequence (no equivalent shifted representation)

    def del4_ok(p):
        s = win.sub(p - 1, p + 4)
        return s[1] != s[5] and s[0] != s[4] and len(set(s[1:5])) >= 3
    p = pick_near(win, 5_206_900, used, extra=del4_ok)
    V.append(Variant("DEL4_unique", p - 1, win.sub(p - 1, p + 3), win.base(p - 1), {"M1"}, "DEL",
                     "unique sequence; only one possible representation",
                     "simple indel, maternal, transmitted", "called in S01 and S03"))
    used.append(p)
    # (e) multi-allelic SNV: ALT1 on P1 (father), ALT2 on M1 (mother) -> child 1/2
    p = pick_near(win, 5_245_200, used)
    r = win.base(p)
    a1, a2 = [b for b in BASES if b != r][:2]
    V.append(Variant("MULTI_ALT1", p, r, a1, {"P1"}, "SNV", "multi-allelic site, allele 1 (paternal)",
                     "multi-allelic record: split by bcftools norm -m -any", "S01 1/2, S02 0/1, S03 0/2"))
    V.append(Variant("MULTI_ALT2", p, r, a2, {"M1"}, "SNV", "multi-allelic site, allele 2 (maternal)",
                     "multi-allelic record: split by bcftools norm -m -any", "S01 1/2, S02 0/1, S03 0/2"))
    used.append(p)
    # (f) de novo SNV in the child only
    p = pick_near(win, 5_272_400, used)
    V.append(Variant("DENOVO", p, win.base(p), other_base(win.base(p), rng), {"CP"}, "SNV",
                     "non-coding; present only in the child",
                     "de novo candidate: the one Mendelian inconsistency", "S01 0/1, parents 0/0"))
    used.append(p)
    # (g) true heterozygous SNV inside the low-coverage interval (S03 only, untransmitted)
    p = pick_near(win, (LOWCOV[0] + LOWCOV[1]) // 2 - 20, used)
    assert LOWCOV[0] + 150 < p < LOWCOV[1] - 150
    V.append(Variant("LOWCOV_HET", p, win.base(p), other_base(win.base(p), rng), {"M2"}, "SNV",
                     f"inside the low-coverage interval {CHROM}:{LOWCOV[0]}-{LOWCOV[1]} (~8x)",
                     "true variant that a depth filter removes",
                     "called (S03 0/1) but removed by FMT/DP>=10 and INFO/DP>=30 filters"))
    used.append(p)
    # artifacts (not true variants) - positions only, applied during read simulation
    pj = pick_near(win, 5_213_000, used)
    used.append(pj)
    pq = pick_near(win, 5_257_100, used)              # provisional; final choice after the dry run
    used.append(pq)
    # (c) background non-coding SNVs with Mendelian-consistent genotypes
    cand = [q for q in range(WIN_START + 2000, WIN_END - 2000) if win.eligible_snv(q)]
    rng.shuffle(cand)
    bg = []
    for q in cand:
        if len(bg) == N_BACKGROUND_SNV:
            break
        if all(abs(q - u) >= 1000 for u in used + bg):
            bg.append(q)
    bg.sort()
    for k, q in enumerate(bg, 1):
        while True:
            carry = {h for h in ("P1", "P2", "M1", "M2") if rng.random() < 0.35}
            if not carry:
                continue
            if carry == {"P1", "M1"}:      # would mimic the recessive pattern (child 1/1, parents 0/1)
                continue
            break
        if "P1" in carry:
            carry.add("CP")
        V.append(Variant(f"SNV{k:02d}", q, win.base(q), other_base(win.base(q), rng), carry, "SNV",
                         "non-coding background SNV", "ordinary inherited variation",
                         "called; Mendelian-consistent"))
    # haplotype CP (child paternal) = P1 + de novo
    for v in V:
        if "P1" in v.haps:
            v.haps.add("CP")
    V.sort(key=lambda v: (v.pos, v.alt))
    arts = dict(jackpot=dict(pos=pj, ref=win.base(pj), alt=other_base(win.base(pj), rng)),
                lowbq=dict(pos=pq, ref=win.base(pq), alt=other_base(win.base(pq), rng)))
    return V, arts


class Haplotype:
    """A haplotype sequence plus a map from haplotype index to reference window index."""

    def __init__(self, name, win, variants):
        self.name = name
        seq, h2r = [], []
        i = 0
        vs = sorted([v for v in variants if name in v.haps], key=lambda v: v.pos)
        for v in vs:
            j = v.pos - WIN_START
            assert j >= i, f"overlapping variants on {name}"
            seq.append(win.seq[i:j]); h2r.extend(range(i, j))
            assert win.seq[j:j + len(v.ref)] == v.ref, (v.eid, v.ref)
            # anchor/first base(s) of ALT map to the REF positions they replace
            seq.append(v.alt)
            for t in range(len(v.alt)):
                h2r.append(j + min(t, len(v.ref) - 1))
            i = j + len(v.ref)
        seq.append(win.seq[i:]); h2r.extend(range(i, WLEN))
        s = "".join(seq)
        self.seq = np.frombuffer(s.encode(), np.uint8).copy()
        self.codes = np.array([B2I[c] for c in s], np.uint8)
        self.h2r = np.array(h2r, np.int64)
        assert len(self.h2r) == len(self.codes)
        r2h = np.full(WLEN, -1, np.int64)
        for hi, ri in enumerate(self.h2r):
            if r2h[ri] == -1:
                r2h[ri] = hi
        self.r2h = r2h


# ----------------------------------------------------------------------------- read simulation
def quality_profile(read):
    """Per-cycle probabilities of the four quality bins (37, 25, 11, 2)."""
    t = np.linspace(0, 1, READ_LEN)
    worse = 1.0 if read == 1 else 1.5
    p2 = 0.002 + 0.020 * worse * t ** 4
    p11 = 0.008 + 0.040 * worse * t ** 2
    p25 = 0.030 + 0.120 * worse * t ** 2
    # the first two cycles are slightly worse as well
    p25[:2] += 0.03
    p11[:2] += 0.01
    p37 = 1 - p2 - p11 - p25
    return np.vstack([p37, p25, p11, p2]).T          # (150, 4)


def draw_quals(rng, n, read):
    prof = quality_profile(read)
    cum = np.cumsum(prof, axis=1)                     # (150, 4)
    u = rng.random((n, READ_LEN))
    idx = (u[:, :, None] > cum[None, :, :]).sum(axis=2)
    idx = np.minimum(idx, 3)
    return QBINS[idx].astype(np.int16)


def add_errors(rng, codes, quals):
    pe = 10.0 ** (-quals / 10.0)
    err = rng.random(codes.shape) < pe
    shift = rng.integers(1, 4, size=codes.shape).astype(np.uint8)
    out = codes.copy()
    out[err] = (out[err] + shift[err]) % 4
    return out, int(err.sum())


def rc_codes(a):
    return (3 - a[::-1]).astype(np.uint8)


def sample_fragments(haps, rng):
    """Draw unique fragments (hap, start, length, strand); reads touching LOWCOV are thinned."""
    lc0, lc1 = LOWCOV[0] - WIN_START, LOWCOV[1] - WIN_START
    n_unique = int(round(TARGET_DEPTH * WLEN / (2 * READ_LEN)))
    frags = []
    n_rejected = 0
    while len(frags) < n_unique:
        m = 4000
        h = rng.integers(0, 2, m)
        L = np.clip(np.rint(rng.normal(FRAG_MEAN, FRAG_SD, m)), FRAG_MIN, FRAG_MAX).astype(int)
        u = rng.random(m)
        strand = rng.integers(0, 2, m)
        keep_u = rng.random(m)
        for k in range(m):
            hap = haps[h[k]]
            s = int(u[k] * (len(hap.codes) - L[k] + 1))
            e = s + L[k]
            a0, a1 = hap.h2r[s], hap.h2r[s + READ_LEN - 1]
            b0, b1 = hap.h2r[e - READ_LEN], hap.h2r[e - 1]
            touch = (a0 <= lc1 and a1 >= lc0) or (b0 <= lc1 and b1 >= lc0)
            if touch and keep_u[k] > LOWCOV_KEEP:
                n_rejected += 1
                continue
            frags.append((int(h[k]), s, int(L[k]), int(strand[k])))
            if len(frags) == n_unique:
                break
    return frags, n_unique, n_rejected


def fragment_depth(haps, frags):
    """Unique-fragment depth per window position (a read pair counts once), total and per haplotype."""
    diff = np.zeros((2, WLEN + 1), np.int64)
    for h, s, L, _ in frags:
        hap = haps[h]
        a = (int(hap.h2r[s]), int(hap.h2r[s + READ_LEN - 1]))
        b = (int(hap.h2r[s + L - READ_LEN]), int(hap.h2r[s + L - 1]))
        spans = [a, b] if b[0] > a[1] else [(a[0], max(a[1], b[1]))]
        for x, y in spans:
            diff[h, x] += 1
            diff[h, y + 1] -= 1
    per_hap = np.cumsum(diff, axis=1)[:, :WLEN]
    return per_hap.sum(axis=0), per_hap


def simulate_sample(sample, haps, rng, arts):
    frags, n_unique, n_rejected = sample_fragments(haps, rng)
    n_dup = int(round(DUP_FRACTION * n_unique))
    dup_idx = rng.integers(0, n_unique, n_dup)
    templates = [(f, False) for f in frags] + [(frags[i], False) for i in dup_idx]
    info = dict(sample=sample, unique_fragments=n_unique, random_duplicates=n_dup,
                fragments_rejected_lowcov=n_rejected)
    # the S03 jackpot: one fragment with a PCR error, copied JACKPOT_COPIES times
    jp = None
    if sample == "S03":
        hap_i = 1                                       # M2 carries no planted variant here
        hap = haps[hap_i]
        hi = hap.r2h[arts["jackpot"]["pos"] - WIN_START]
        L = int(np.clip(np.rint(rng.normal(FRAG_MEAN, FRAG_SD)), FRAG_MIN, FRAG_MAX))
        s = int(hi - JACKPOT_R1_OFFSET)
        jp = (hap_i, s, L, 0)
        templates += [(jp, True)] * JACKPOT_COPIES
        info.update(jackpot_copies=JACKPOT_COPIES, jackpot_fragment_start=int(hap.h2r[s] + WIN_START),
                    jackpot_fragment_length=L)
    n = len(templates)
    R1 = np.empty((n, READ_LEN), np.uint8)
    R2 = np.empty((n, READ_LEN), np.uint8)
    cover = []                                          # per template: (r1_hap_idx0, r1_rev, r2_idx0, r2_rev)
    for j, ((hi_, s, L, strand), is_jp) in enumerate(templates):
        hap = haps[hi_]
        frag = hap.codes[s:s + L].copy()
        if is_jp:
            frag[JACKPOT_R1_OFFSET] = B2I[arts["jackpot"]["alt"]]
        left, right = frag[:READ_LEN], rc_codes(frag[L - READ_LEN:])
        if strand == 0:
            R1[j], R2[j] = left, right
        else:
            R1[j], R2[j] = right, left
        cover.append((hi_, s, L, strand))
    Q1, Q2 = draw_quals(rng, n, 1), draw_quals(rng, n, 2)
    R1e, e1 = add_errors(rng, R1, Q1)
    R2e, e2 = add_errors(rng, R2, Q2)
    info["substitution_errors"] = e1 + e2
    if jp is not None:                                  # PCR error base: high quality, no extra error
        jrows = [j for j, (_, is_jp) in enumerate(templates) if is_jp]
        for j in jrows:
            Q1[j, JACKPOT_R1_OFFSET] = 37
            R1e[j, JACKPOT_R1_OFFSET] = B2I[arts["jackpot"]["alt"]]
    if sample == "S02":                                 # systematic low-quality error
        pq = arts["lowbq"]["pos"] - WIN_START
        alt = B2I[arts["lowbq"]["alt"]]
        hits = {}                                       # template -> [(read 1/2, offset, reverse?)]
        for j, (hi_, s, L, strand) in enumerate(cover):
            hp = haps[hi_].r2h[pq]
            if s <= hp < s + READ_LEN:                  # covered by the forward-strand read
                hits.setdefault(j, []).append((1 if strand == 0 else 2, hp - s, False))
            if s + L - READ_LEN <= hp < s + L:          # covered by the reverse-strand read
                hits.setdefault(j, []).append((2 if strand == 0 else 1, (s + L - 1) - hp, True))
        tlist = sorted(hits)
        n_pick = int(round(LOWBQ_FRACTION * len(tlist)))
        pick = rng.choice(len(tlist), n_pick, replace=False)
        for t in sorted(pick):
            for rd, off, rev in hits[tlist[t]]:
                j = tlist[t]
                code = (3 - alt) if rev else alt
                if rd == 1:
                    R1e[j, off], Q1[j, off] = code, LOWBQ_Q
                else:
                    R2e[j, off], Q2[j, off] = code, LOWBQ_Q
        info.update(lowbq_pairs_covering=len(tlist), lowbq_pairs_with_error=n_pick)
    order = rng.permutation(n)
    return R1e[order], Q1[order], R2e[order], Q2[order], info


def write_fastq(path, sample, mate, codes, quals):
    seq_lut = np.frombuffer(BASES.encode(), np.uint8)
    seqs = seq_lut[codes]
    qchr = (quals + 33).astype(np.uint8)
    lines = []
    for i in range(codes.shape[0]):
        lines.append(b"@SYNTH_%s_%06d %d:N:0:1\n%s\n+\n%s\n"
                     % (sample.encode(), i + 1, mate, seqs[i].tobytes(), qchr[i].tobytes()))
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=6) as gz:
            gz.write(b"".join(lines))


# ----------------------------------------------------------------------------- truth files
VCF_HEADER = """##fileformat=VCFv4.2
##source=practice1 generate_data.py (SYNTHETIC trio; seed {seed})
##reference=data/ref/chr11_GRCh38_betaglobin.fa.gz (GRCh38 chr11; bases outside {a}-{b} are N)
##contig=<ID=chr11,length={L},assembly=GRCh38>
##INFO=<ID=EVENT,Number=.,Type=String,Description="Planted event ID (see data/truth/planted_events.tsv)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS01\tS02\tS03
"""

GTSTR = {0: "0/0", 1: "0/1", 2: "1/1"}


def write_truth(V, win, ref_fa):
    os.makedirs(D("truth"), exist_ok=True)
    hdr = VCF_HEADER.format(seed=SEED, a=WIN_START, b=WIN_END, L=CHROM_LEN)
    recs = []
    multi = [v for v in V if v.eid.startswith("MULTI")]
    for v in V:
        if v.eid.startswith("MULTI"):
            continue
        g = [GTSTR[v.gt(s)] for s in SAMPLES]
        pos, ref, alt = v.pos, v.ref, v.alt
        if v.eid == "DEL1_homopolymer":            # raw file: right-shifted representation
            hp, rl, rb = v.run
            pos, ref, alt = hp + rl - 2, rb + rb, rb
        recs.append((pos, v.vid, ref, alt, v.eid, g))
    m1, m2 = multi
    gts = []
    for s in SAMPLES:
        h = {"S01": ("CP", "M1"), "S02": ("P1", "P2"), "S03": ("M1", "M2")}[s]
        al = sorted(1 if x in m1.haps else 2 if x in m2.haps else 0 for x in h)
        gts.append(f"{al[0]}/{al[1]}")
    recs.append((m1.pos, ".", m1.ref, f"{m1.alt},{m2.alt}", "MULTI_ALT1,MULTI_ALT2", gts))
    recs.sort()
    raw_txt = D("truth", "truth_raw_unnormalised.vcf")
    with open(raw_txt, "w") as fh:
        fh.write(hdr)
        for pos, vid, ref, alt, eid, g in recs:
            fh.write(f"{CHROM}\t{pos}\t{vid}\t{ref}\t{alt}\t.\tPASS\tEVENT={eid}\tGT\t" + "\t".join(g) + "\n")
    raw_gz = D("truth", "truth_raw_unnormalised.vcf.gz")
    subprocess.run(["bcftools", "view", "--no-version", "-Oz", "-o", raw_gz, raw_txt], check=True)
    os.remove(raw_txt)
    norm_gz = D("truth", "truth_trio.vcf.gz")
    p = subprocess.run(["bcftools", "norm", "--no-version", "-f", ref_fa, "-m", "-any", "-Oz",
                        "-o", norm_gz, raw_gz], capture_output=True, text=True, check=True)
    subprocess.run(["bcftools", "index", "-f", "-t", norm_gz], check=True)
    # the normalised file must contain exactly the left-aligned, split records we planted
    out = subprocess.run(["bcftools", "query", "-f", "%POS\t%REF\t%ALT\t[%GT ]\n", norm_gz],
                         capture_output=True, text=True, check=True).stdout.split("\n")
    got = sorted(tuple(l.split("\t")[:3]) for l in out if l)
    want = sorted((str(v.pos), v.ref, v.alt) for v in V)
    assert got == want, (set(got) ^ set(want))
    return raw_gz, norm_gz, p.stderr.strip()


def write_ped():
    with open(D("trio.ped"), "w") as fh:
        fh.write("#family\tindividual\tfather\tmother\tsex(1=male,2=female)\tphenotype(1=unaffected,2=affected)\n")
        fh.write("FAM01\tS01\tS02\tS03\t2\t2\n")
        fh.write("FAM01\tS02\t0\t0\t1\t1\n")
        fh.write("FAM01\tS03\t0\t0\t2\t1\n")


def write_events(V, arts, sim_info):
    rows = []
    for v in V:
        rows.append([v.eid, "true_variant", v.cls, CHROM, v.pos, v.ref, v.alt,
                     GTSTR[v.gt("S01")], GTSTR[v.gt("S02")], GTSTR[v.gt("S03")],
                     v.context, v.purpose, v.expected])
    s2 = [i for i in sim_info if i["sample"] == "S02"][0]
    s3 = [i for i in sim_info if i["sample"] == "S03"][0]
    j, q = arts["jackpot"], arts["lowbq"]
    rows.append(["ARTIFACT_JACKPOT", "artifact", "PCR_duplicate_jackpot", CHROM, j["pos"], j["ref"], j["alt"],
                 "0/0", "0/0", "0/0 (artifact in reads)",
                 f"S03: one fragment with a PCR error, {s3['jackpot_copies']} copies with identical ends, "
                 f"error at R1 cycle {JACKPOT_R1_OFFSET + 1} with base quality 37",
                 "shows why duplicate marking matters",
                 "false 0/1 call in S03 without markdup; disappears after samtools markdup"])
    rows.append(["ARTIFACT_LOWBQ", "artifact", "systematic_low_quality_error", CHROM, q["pos"], q["ref"], q["alt"],
                 "0/0", "0/0 (artifact in reads)", "0/0",
                 f"S02: {s2['lowbq_pairs_with_error']} of {s2['lowbq_pairs_covering']} read pairs show the same "
                 f"wrong base, always with base quality {LOWBQ_Q}",
                 "shows what the minimum base-quality filter protects against",
                 "skipped by the course pipeline (bcftools mpileup -Q 13, set explicitly); becomes a weak "
                 "false 0/1 call with -Q 0 and also with bcftools' own default (-Q 1)"])
    rows.append(["LOWCOV_INTERVAL", "coverage_feature", "reduced_sampling", CHROM, LOWCOV[0], ".", ".",
                 ".", ".", ".", f"{CHROM}:{LOWCOV[0]}-{LOWCOV[1]}; {int(LOWCOV_KEEP * 100)}% of fragments kept (~8x)",
                 "depth plot feature; hosts LOWCOV_HET", "visible dip in the depth plot"])
    cols = ["event_id", "category", "type", "chrom", "pos", "ref", "alt", "S01", "S02", "S03",
            "context", "purpose", "expected_behaviour"]
    with open(D("truth", "planted_events.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--download", action="store_true", help="re-download the UCSC inputs first")
    args = ap.parse_args()
    dna, mane, schema, version_txt, seqjson = load_inputs(args.download)
    retrieved = seqjson.get("downloadTime", "?").replace(":", "-", 2)
    genes = parse_mane(mane)
    exons = [e for g in genes for e in g["exons"]]
    win = Window(dna, exons)
    assert win.sub(5_226_991, 5_227_010) == CHECK_20MER[1]
    log("window loaded:", WLEN, "bp; repeat-masked fraction %.3f" % win.lower.mean())
    ref_fa = write_reference(win)
    log("reference written:", os.path.relpath(ref_fa, PKG))
    hbb = write_annotation(genes, retrieved, version_txt, schema)
    ss = np.random.SeedSequence(SEED)
    rng_var, rng_s1, rng_s2, rng_s3 = [np.random.default_rng(s) for s in ss.spawn(4)]
    V, arts = plant_variants(win, rng_var)
    rngs = dict(zip(SAMPLES, (rng_s1, rng_s2, rng_s3)))

    def build(V):
        haps = {h: Haplotype(h, win, V) for h in ("P1", "P2", "M1", "M2", "CP")}
        return {"S01": (haps["CP"], haps["M1"]), "S02": (haps["P1"], haps["P2"]),
                "S03": (haps["M1"], haps["M2"])}
    sample_haps = build(V)
    # dry run of the fragment sampling on copies of the random streams -> depth profiles
    depth = {}
    for s in SAMPLES:
        frags, _, _ = sample_fragments(sample_haps[s], copy.deepcopy(rngs[s]))
        depth[s] = fragment_depth(sample_haps[s], frags)
    swap = {"A": "C", "C": "A", "G": "T", "T": "G"}
    # (g) place the low-coverage SNV where S03 has 7-9 fragments (3-5 from M2) and S01/S02 have 6-10
    lv = [v for v in V if v.eid == "LOWCOV_HET"][0]
    centre = (LOWCOV[0] + LOWCOV[1]) // 2
    cands = []
    for p in range(LOWCOV[0] + 150, LOWCOV[1] - 150 + 1):
        i = p - WIN_START
        if not win.eligible_snv(p):
            continue
        d1, d2, d3 = (depth[s][0][i] for s in SAMPLES)
        alt3 = depth["S03"][1][1][i]                     # S03 haplotype index 1 = M2
        if 7 <= d3 <= 9 and 3 <= alt3 <= 5 and 6 <= d1 <= 10 and 6 <= d2 <= 10:
            cands.append((abs(p - centre), p))
    p = min(cands)[1]
    lv.pos, lv.ref = p, win.base(p)
    lv.alt = lv.alt if lv.alt != lv.ref else swap[lv.ref]
    # (i) place the low-quality error where S02 has 32-36 fragments, >= 1 kb from every other event
    taken = [v.pos for v in V] + [arts["jackpot"]["pos"]]
    cands = []
    for p in range(WIN_START + 2000, WIN_END - 2000):
        if any(abs(p - u) < 1000 for u in taken) or not win.eligible_snv(p):
            continue
        if 32 <= depth["S02"][0][p - WIN_START] <= 36:
            cands.append((abs(p - 5_257_100), p))
    p = min(cands)[1]
    old_alt = arts["lowbq"]["alt"]
    arts["lowbq"] = dict(pos=p, ref=win.base(p), alt=old_alt if old_alt != win.base(p) else swap[win.base(p)])
    V.sort(key=lambda v: (v.pos, v.alt))
    sample_haps = build(V)
    os.makedirs(D("fastq"), exist_ok=True)
    sim_info = []
    for s in SAMPLES:
        R1, Q1, R2, Q2, info = simulate_sample(s, sample_haps[s], rngs[s], arts)
        write_fastq(D("fastq", f"{s}_R1.fastq.gz"), s, 1, R1, Q1)
        write_fastq(D("fastq", f"{s}_R2.fastq.gz"), s, 2, R2, Q2)
        info["read_pairs"] = int(R1.shape[0])
        info["dry_run_fragment_depth_at_lowcov_het"] = int(depth[s][0][lv.pos - WIN_START])
        info["dry_run_fragment_depth_at_lowbq_site"] = int(depth[s][0][arts["lowbq"]["pos"] - WIN_START])
        sim_info.append(info)
        log(s, info)
    raw_gz, norm_gz, norm_msg = write_truth(V, win, ref_fa)
    log("truth VCFs written; bcftools norm said:", norm_msg)
    write_ped()
    write_events(V, arts, sim_info)
    prov = dict(
        description="Provenance of the public inputs used to build the SYNTHETIC Practice 1 data",
        seed=SEED, numpy_version=np.__version__, python=sys.version.split()[0],
        generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        sources=[
            dict(what="GRCh38 chr11:5200001-5300000 sequence", url=UCSC_SEQ_URL, genome="hg38",
                 downloadTime=seqjson.get("downloadTime"),
                 md5_uppercase_window=hashlib.md5(win.seq.encode()).hexdigest()),
            dict(what="MANE transcripts in the window (UCSC hg38 track 'mane')", url=UCSC_MANE_URL,
                 downloadTime=mane.get("downloadTime"), track_dataTime=schema.get("dataTime"),
                 dataVersion_file=MANE_VERSION_URL, dataVersion=version_txt.strip().replace("\t", ": ").split("\n")),
            dict(what="chr11 length", value=CHROM_LEN, url="https://api.genome.ucsc.edu/list/chromosomes?genome=hg38"),
        ],
        simulation=dict(read_length=READ_LEN, fragment_length=f"N({FRAG_MEAN},{FRAG_SD}) clipped to [{FRAG_MIN},{FRAG_MAX}]",
                        target_depth_unique=TARGET_DEPTH, duplicate_fraction=DUP_FRACTION,
                        quality_bins={"37": "F", "25": ":", "11": ",", "2": "#"},
                        lowcov_interval=f"{CHROM}:{LOWCOV[0]}-{LOWCOV[1]}", lowcov_keep_fraction=LOWCOV_KEEP,
                        artifacts=arts, per_sample=sim_info),
    )
    with open(D("provenance", "provenance.json"), "w") as fh:
        json.dump(prov, fh, indent=2)
    files = [os.path.join(dp, f) for dp, _, fs in os.walk(D()) for f in fs
             if not f.startswith(".") and f != "MD5SUMS.txt" and "provenance" not in dp]
    with open(D("MD5SUMS.txt"), "w") as fh:
        for f in sorted(files):
            fh.write(f"{md5(f)}  {os.path.relpath(f, PKG)}\n")
    log("done.")


if __name__ == "__main__":
    main()
