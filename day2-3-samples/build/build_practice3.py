#!/usr/bin/env python3
"""Build data/practice3/ and checkpoints/practice3_*.tsv for Mini practice 3.

Mini practice 3 . One person, six datasets

One synthetic person (P3-01), one synthetic 120-kb locus (chrS), six libraries: DNA and RNA,
each read by Illumina, PacBio and Oxford Nanopore. Everything is synthetic teaching data: the
read models approximate public platform specifications but are not calibrated to any real run.
Deterministic (fixed seed), no network. Writes about 6 MB.

Coordinates in this file are 1-based inclusive unless a name ends in 0 (0-based index).

Usage
  python build/build_practice3.py
"""
import gzip
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D3 = ROOT / "data" / "practice3"
CHECK = ROOT / "checkpoints"
SEED = 3033
CHROM, LOCUS_LEN = "chrS", 120_000
PERSON = "P3-01"

# ---------------------------------------------------------------- the truth (spec section 4)
CGI_ALTX = (18_801, 20_000)
ALTX_EXONS = {"E1": (20_001, 20_200), "E2a": (22_001, 22_120), "E2b": (23_001, 23_120),
              "E3": (25_001, 25_900), "E4a": (27_001, 27_150), "E4b": (28_001, 28_150),
              "E5": (30_001, 30_800)}
TRANSCRIPTS = {
    "ALTX-201": ("ALTX", [ALTX_EXONS[e] for e in ("E1", "E2a", "E3", "E4a", "E5")]),
    "ALTX-202": ("ALTX", [ALTX_EXONS[e] for e in ("E1", "E2b", "E3", "E4b", "E5")]),
    "ASE1-201": ("ASE1", [(38_001, 38_300), (40_001, 40_200), (44_001, 44_250), (47_001, 48_500)]),
    "OFF1-201": ("OFF1", [(66_001, 66_250), (69_001, 69_200), (72_001, 72_180), (75_001, 75_900)]),
    "LOW1-201": ("LOW1", [(94_001, 94_300), (96_001, 96_400), (98_001, 99_200)]),
}
GENES = {"ALTX": (20_001, 30_800), "ASE1": (38_001, 48_500), "OFF1": (66_001, 75_900),
         "LOW1": (94_001, 99_200)}
ABUNDANCE = {"ASE1": 10, "ALTX": 4, "LOW1": 1, "OFF1": 0}          # molecules in blood
ISOFORM_SHARE = {"ALTX-201": 0.60, "ALTX-202": 0.40}
REP1 = (34_001, 34_060)                    # (CAG)20 in the reference
REP_UNITS = {1: 25, 2: 180}                # hap1, hap2
SNV1 = dict(name="SNV1", pos=44_101, ref="G", alt="A", hap=1, note="exon 3 of ASE1; heterozygous, on hap1")
SNV2 = dict(name="SNV2", pos=41_101, ref="C", alt="T", hap=2, note="intron 2 of ASE1; heterozygous, on hap2")
ASE_HAP1_SHARE = 0.85                      # hap1 makes 85% of ASE1 RNA
HP1 = (52_001, 52_010)                     # A x 10, both haplotypes
DUP_A, DUP_B = (56_001, 62_000), (84_001, 90_000)
DUP_OFFSETS = (1_500, 3_000, 4_500)        # DUP_B differs from DUP_A only here
OFF1_PROMOTER = (65_201, 66_000)
METHYLATION = [("CGI_ALTX", CGI_ALTX, 0.05), ("OFF1_PROMOTER", OFF1_PROMOTER, 0.90),
               ("ALTX_body", GENES["ALTX"], 0.80)]
METH_DEFAULT = 0.70
# background GC in 5-kb blocks (spec: between ~30% and ~65%)
GC_BLOCKS = [0.42, 0.62, 0.36, 0.47, 0.52, 0.40, 0.33, 0.48, 0.45, 0.57, 0.44, 0.46,
             0.43, 0.50, 0.64, 0.38, 0.31, 0.46, 0.61, 0.42, 0.65, 0.35, 0.49, 0.60]

LIBRARIES = ["DNA-ILMN", "DNA-HIFI", "DNA-ONT", "RNA-ILMN", "RNA-KINNEX", "RNA-ONT"]

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
LUT = np.zeros(256, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
COMP = np.zeros(256, dtype=np.uint8)
for _a, _b in zip(b"ACGTN", b"TGCAN"):
    COMP[_a] = _b
M, I, D, N, S = 0, 1, 2, 3, 4
OPS = "MIDNS"


def b(s):
    return np.frombuffer(s.encode(), dtype=np.uint8).copy()


# ---------------------------------------------------------------- the reference and haplotypes
def random_seq(rng, n, gc):
    return BASES[rng.choice(4, size=n, p=[(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])]


def gc_before_depletion(target, keep):
    """Depletion removes a fraction (1 - keep) of CG pairs, one G/C base each: GC drops from g to
    about g - c*g^2 with c = (1 - keep)/4. Solve for the g that lands on the target."""
    c = (1 - keep) / 4
    return (1 - math.sqrt(1 - 4 * c * target)) / (2 * c)


def deplete_cpg(rng, seq, start0, end0, keep):
    """Turn most CG dinucleotides into TG or CA, as in human DNA outside CpG islands."""
    s = seq[start0:end0]
    cg = np.flatnonzero((s[:-1] == ord("C")) & (s[1:] == ord("G")))
    hit = cg[rng.random(len(cg)) > keep]
    flip = rng.random(len(hit)) < 0.5
    s[hit[flip]] = ord("T")
    s[hit[~flip] + 1] = ord("A")


def build_reference(rng):
    seq = np.concatenate([random_seq(rng, 5_000, gc_before_depletion(gc, 0.25)) for gc in GC_BLOCKS])
    deplete_cpg(rng, seq, 0, LOCUS_LEN, keep=0.25)
    for _, (_, exons) in TRANSCRIPTS.items():
        for s, e in exons:
            seq[s - 1:e] = random_seq(rng, e - s + 1, gc_before_depletion(0.50, 0.45))
            deplete_cpg(rng, seq, s - 1, e, keep=0.45)
    s, e = CGI_ALTX
    seq[s - 1:e] = random_seq(rng, e - s + 1, 0.72)                      # CpG-dense island
    s, e = OFF1_PROMOTER
    seq[s - 1:e] = random_seq(rng, e - s + 1, gc_before_depletion(0.65, 0.80))
    deplete_cpg(rng, seq, s - 1, e, keep=0.80)
    # REP1: (CAG)20 between flanks that contain no CAG
    s, e = REP1
    seq[s - 1 - 20:s - 1] = b("GTTTGACTTGTTTATGTTGT")
    seq[s - 1:e] = b("CAG" * 20)
    seq[e:e + 20] = b("TTGTAATGTTTGTGATTGTA")
    # HP1: A x 10 with non-A flanks
    s, e = HP1
    seq[s - 1 - 5:s - 1] = b("GTCTC")
    seq[s - 1:e] = b("A" * 10)
    seq[e:e + 5] = b("GCTGT")
    # SNV reference bases, with neighbours that keep CpGs out of the picture
    seq[SNV1["pos"] - 2:SNV1["pos"]] = b("TG")
    seq[SNV2["pos"] - 1:SNV2["pos"] + 1] = b("CA")
    # DUP_B is DUP_A with three transitions
    a0, b0 = DUP_A[0] - 1, DUP_B[0] - 1
    size = DUP_A[1] - DUP_A[0] + 1
    seq[b0:b0 + size] = seq[a0:a0 + size]
    transition = {ord("A"): ord("G"), ord("G"): ord("A"), ord("C"): ord("T"), ord("T"): ord("C")}
    for off in DUP_OFFSETS:
        seq[b0 + off] = transition[int(seq[b0 + off])]
    return seq


def build_haplotype(ref, hap):
    """Haplotype sequence and, for every haplotype base, its 0-based reference position
    (-1 for bases that are not in the reference: the extra CAG units)."""
    seq = ref.copy()
    for snv in (SNV1, SNV2):
        if snv["hap"] == hap:
            seq[snv["pos"] - 1] = ord(snv["alt"])
    s, e = REP1
    extra = 3 * (REP_UNITS[hap] - 20)
    seq = np.concatenate([seq[:e], b("CAG" * (REP_UNITS[hap] - 20)), seq[e:]])
    refpos = np.concatenate([np.arange(e), np.full(extra, -1), np.arange(e, LOCUS_LEN)])
    rep0 = (s - 1, s - 1 + 3 * REP_UNITS[hap])       # the repeat on this haplotype, 0-based [ )
    return seq, refpos, rep0


# ---------------------------------------------------------------- read simulation
def homopolymer_runs(tmpl, kmin):
    change = np.flatnonzero(np.diff(tmpl) != 0) + 1
    starts = np.r_[0, change]
    lens = np.diff(np.r_[starts, len(tmpl)])
    keep = lens >= kmin
    return starts[keep], lens[keep]


def draw_errors(rng, tmpl, p_sub, p_ins, p_del, hp=None):
    """Per template base: kind (0 kept, 1 substituted, 3 deleted), the substituted base, how
    many bases are inserted after it and which base. hp = (kmin, rate per base of run,
    probabilities of -2, -1, +1, +2) adds homopolymer length errors."""
    n = len(tmpl)
    u = rng.random(n)
    p_sub = np.broadcast_to(np.asarray(p_sub, dtype=float), (n,))
    kind = np.zeros(n, dtype=np.int8)
    kind[u < p_sub] = 1
    kind[(u >= p_sub) & (u < p_sub + p_del)] = 3
    ins_n = (rng.random(n) < p_ins).astype(np.int64)
    ins_base = BASES[rng.integers(0, 4, n)]
    sub_base = BASES[(LUT[tmpl] + rng.integers(1, 4, n)) % 4]
    if hp is not None:
        kmin, rate, probs = hp
        starts, lens = homopolymer_runs(tmpl, kmin)
        if len(starts):
            p_change = np.minimum(0.6, rate * lens)
            hit = rng.random(len(starts)) < p_change
            what = rng.choice(4, size=len(starts), p=probs)          # -2, -1, +1, +2
            for st, ln, h, w in zip(starts, lens, hit, what):
                if not h:
                    continue
                last = st + ln - 1
                if w <= 1:                                            # shorter run
                    kind[last] = 3
                    if w == 0:
                        kind[last - 1] = 3
                else:                                                 # longer run
                    ins_n[last] += 1 if w == 2 else 2
                    ins_base[last] = tmpl[last]
    return kind, sub_base, ins_n, ins_base


def project(tmpl, refpos, ref, kind, sub_base, ins_n, ins_base, spliced=False):
    """Turn a simulated read into a reference alignment.
    Returns (pos1, cigar, seq, nm, n_aligned) or None when nothing aligns."""
    present = kind != 3
    base = np.where(kind == 1, sub_base, tmpl)
    counts = present.astype(np.int64) + ins_n
    first = np.cumsum(counts) - counts
    seq = np.empty(int(counts.sum()), dtype=np.uint8)
    seq[first[present]] = base[present]
    for t in np.flatnonzero(ins_n):
        off = first[t] + int(present[t])
        seq[off:off + ins_n[t]] = ins_base[t]
    on_ref = refpos >= 0
    # fast path: no indels, no insertion relative to the reference, no intron
    if (ins_n.sum() == 0 and present.all() and on_ref.all()
            and refpos[-1] - refpos[0] == len(refpos) - 1):
        mm = int((base != ref[refpos]).sum())
        return int(refpos[0]) + 1, f"{len(tmpl)}M", seq, mm, len(tmpl)
    prim = np.where(present, np.where(on_ref, M, I), np.where(on_ref, D, -1))
    prim_len = np.where(prim >= 0, 1, 0)
    gap = np.zeros(len(tmpl), dtype=np.int64)
    if spliced:
        gap[:-1] = np.where(on_ref[:-1] & on_ref[1:], refpos[1:] - refpos[:-1] - 1, 0)
    codes = np.stack([prim, np.full(len(tmpl), I), np.full(len(tmpl), N)], axis=1).ravel()
    lens = np.stack([prim_len, ins_n, gap], axis=1).ravel()
    keep = lens > 0
    codes, lens = codes[keep], lens[keep]
    if len(codes) == 0:
        return None
    cut = np.flatnonzero(np.diff(codes) != 0) + 1
    starts = np.r_[0, cut]
    ops = [[int(c), int(n)] for c, n in zip(codes[starts], np.add.reduceat(lens, starts))]
    # read ends: anything before the first M or after the last M is an insertion (soft clip)
    # or a deletion/skip (dropped)
    if not any(c == M for c, _ in ops):
        return None
    first_m = next(k for k, (c, _) in enumerate(ops) if c == M)
    last_m = max(k for k, (c, _) in enumerate(ops) if c == M)
    head = sum(n for c, n in ops[:first_m] if c in (I, S))
    tail = sum(n for c, n in ops[last_m + 1:] if c in (I, S))
    ops = ([[S, head]] if head else []) + ops[first_m:last_m + 1] + ([[S, tail]] if tail else [])
    if not any(c == M for c, _ in ops):
        return None
    m_mask = present & on_ref
    if not m_mask.any():
        return None
    pos1 = int(refpos[np.flatnonzero(m_mask)[0]]) + 1
    mm = int((base[m_mask] != ref[refpos[m_mask]]).sum())
    nm = mm + sum(n for c, n in ops if c in (I, D))
    n_aligned = sum(n for c, n in ops if c == M)
    return pos1, "".join(f"{n}{OPS[c]}" for c, n in ops), seq, nm, n_aligned


def cigar_ref_len(cigar):
    total, num = 0, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            if ch in "MDN":
                total += int(num)
            num = ""
    return total


def illumina_error_profile(n_cycles):
    """Per-cycle error 0.1% at base 1 rising to 0.8% at base 150; QUAL = Phred of that."""
    cyc = np.arange(n_cycles)
    p = 0.001 + 0.007 * cyc / 149
    q = np.round(-10 * np.log10(p)).astype(int)
    return p, "".join(chr(33 + x) for x in q)


def dup_twin(pos1, end1):
    """If [pos1, end1] lies inside one copy of the duplicated block and avoids the three
    differing bases, the read fits both copies equally: return the shift to the other copy."""
    for (s, e), shift in ((DUP_A, DUP_B[0] - DUP_A[0]), (DUP_B, DUP_A[0] - DUP_B[0])):
        if s <= pos1 and end1 <= e:
            diffs = [s + off for off in DUP_OFFSETS]
            if not any(pos1 <= d <= end1 for d in diffs):
                return shift
    return None


class Library:
    def __init__(self, lib):
        self.lib = lib
        self.records = []          # dicts: qname flag pos mapq cigar rnext pnext tlen seq qual nm
        self.meth = None


def revcomp_codes(x):
    return COMP[x][::-1]


def add_pair(L, rng, qname, r1, r2, qual1, qual2, mapq_rule, max_span=1_000):
    """r1/r2: (alignment or None, strand_reverse, seq_if_unmapped). Writes two SAM records.
    A proper pair faces inward; for DNA it must also span <= max_span bases. RNA pairs may span
    an intron (max_span=None), as a spliced aligner would accept."""
    recs = []
    for (aln, rev, raw), mate_bit, qual in ((r1, 0x40, qual1), (r2, 0x80, qual2)):
        rec = dict(qname=qname, flag=0x1 | mate_bit | (0x10 if rev else 0), qual=qual)
        if aln is None:
            rec.update(flag=rec["flag"] | 0x4, pos=0, mapq=0, cigar="*", seq=raw, nm=None)
        else:
            pos1, cigar, seq, nm, _ = aln
            end1 = pos1 + cigar_ref_len(cigar) - 1
            mapq = 60
            shift = dup_twin(pos1, end1) if mapq_rule else None
            if shift is not None:
                mapq = 0
                if rng.random() < 0.5:
                    pos1 += shift
            rec.update(pos=pos1, mapq=mapq, cigar=cigar, seq=seq, nm=nm, end=pos1 + cigar_ref_len(cigar) - 1)
        recs.append(rec)
    a, c = recs
    for me, mate in ((a, c), (c, a)):
        if mate["flag"] & 0x10:
            me["flag"] |= 0x20
        if mate["cigar"] == "*":
            me["flag"] |= 0x8
    both_mapped = a["cigar"] != "*" and c["cigar"] != "*"
    if both_mapped:
        left, right = sorted((a, c), key=lambda r: r["pos"])
        span = max(a["end"], c["end"]) - left["pos"] + 1
        proper = ((max_span is None or span <= max_span)
                  and (left["flag"] & 0x10) == 0 and (right["flag"] & 0x10) != 0)
        for me, mate in ((a, c), (c, a)):
            me.update(rnext="=", pnext=mate["pos"], tlen=span if me is left else -span)
            if proper:
                me["flag"] |= 0x2
    else:
        mapped = [r for r in (a, c) if r["cigar"] != "*"]
        if mapped:
            m = mapped[0]
            for me in (a, c):
                me.update(rnext="=", pnext=m["pos"], tlen=0)
                if me["cigar"] == "*":
                    me["pos"] = m["pos"]
        else:
            for me in (a, c):
                me.update(rnext="*", pnext=0, tlen=0)
    L.records += [a, c]


def read_qual(err, n):
    """One Phred value per long read: the read's simulated error rate (real files vary per base)."""
    q = int(np.clip(round(-10 * math.log10(err)), 1, 60))
    return chr(33 + q) * n


def add_single(L, rng, qname, aln, rev, mapq_rule, qual):
    pos1, cigar, seq, nm, _ = aln
    end1 = pos1 + cigar_ref_len(cigar) - 1
    mapq = 60
    shift = dup_twin(pos1, end1) if mapq_rule else None
    if shift is not None:
        mapq = 0
        if rng.random() < 0.5:
            pos1 += shift
    L.records.append(dict(qname=qname, flag=0x10 if rev else 0, pos=pos1, mapq=mapq, cigar=cigar,
                          rnext="*", pnext=0, tlen=0, seq=seq, qual=qual, nm=nm,
                          end=pos1 + cigar_ref_len(cigar) - 1))


def sim_dna_illumina(rng, ref, haps, n_pairs=12_000, read_len=150):
    L = Library("DNA-ILMN")
    p_cyc, qual = illumina_error_profile(read_len)
    gc_cum = {h: np.r_[0, np.cumsum((s == ord("G")) | (s == ord("C")))] for h, (s, _, _) in haps.items()}
    made = 0
    while made < n_pairs:
        h = 1 + int(rng.random() < 0.5)
        hseq, hpos, rep0 = haps[h]
        ins = int(np.clip(np.rint(rng.normal(350, 50)), 200, 600))
        s = int(rng.integers(0, len(hseq) - ins + 1))
        gc = (gc_cum[h][s + ins] - gc_cum[h][s]) / ins
        if rng.random() > math.exp(-0.5 * ((gc - 0.45) / 0.18) ** 2):   # PCR under-amplifies extreme GC
            continue
        made += 1
        frag_rev = rng.random() < 0.5
        mates = []
        for mate in (0, 1):
            fwd = (mate == 0) != frag_rev                  # R1 reads the fragment's own strand
            a0 = s if fwd else s + ins - read_len
            tmpl, tpos = hseq[a0:a0 + read_len], hpos[a0:a0 + read_len]
            p = p_cyc if fwd else p_cyc[::-1]
            kind, sub_base, ins_n, ins_base = draw_errors(rng, tmpl, p, 0.0, 0.0)
            inside_repeat = rep0[0] <= a0 and a0 + read_len <= rep0[1]
            aln = None if inside_repeat else project(tmpl, tpos, ref, kind, sub_base, ins_n, ins_base)
            if aln is not None and aln[4] < 30:            # too little unique sequence to place
                aln = None
            raw = np.where(kind == 1, sub_base, tmpl)
            mates.append((aln, not fwd, raw, qual if fwd else qual[::-1]))
        add_pair(L, rng, f"DNA-ILMN_{made:06d}", mates[0][:3], mates[1][:3], mates[0][3], mates[1][3], True)
    return L


def long_read_lengths(rng, kind):
    if kind == "hifi":
        while True:
            x = rng.normal(16_000, 3_000)
            if 8_000 <= x <= 30_000:
                return int(x)
    sigma = 0.9                                            # ONT: log-normal, N50 about 22 kb
    mu = math.log(22_000) - sigma ** 2
    while True:
        x = rng.lognormal(mu, sigma)
        if 1_000 <= x <= 100_000:
            return int(x)


HIFI_HP = (4, 0.004, [0.05, 0.55, 0.37, 0.03])
ONT_HP = (5, 0.035, [0.12, 0.45, 0.35, 0.08])


MEAN_LONG_LEN = {"hifi": 16_000, "ont": 14_500}     # means of the two length models above


def sim_dna_long(rng, ref, haps, lib, target_x=20, pad=100_000):
    """PacBio HiFi or ONT reads from either haplotype. The haplotype is padded with sequence from
    'the rest of the chromosome' so reads keep their full length at the locus edges; the part
    outside chrS becomes a soft clip, as it would in a whole-genome alignment.

    Read starts sit on a jittered grid and alternate between the two haplotypes. A real run
    scatters reads at random (Poisson-like coverage); with only ~20 reads over any base, that
    noise alone could hide a lesson (for example a heterozygous site reading 75:25), so the
    teaching data spread reads evenly instead."""
    L = Library(lib)
    kind = "hifi" if lib == "DNA-HIFI" else "ont"
    padded = {}
    for h, (hseq, hpos, rep0) in haps.items():
        left, right = random_seq(rng, pad, 0.42), random_seq(rng, pad, 0.42)
        padded[h] = (np.concatenate([left, hseq, right]),
                     np.concatenate([np.full(pad, -1), hpos, np.full(pad, -1)]))
    spacing = MEAN_LONG_LEN[kind] / target_x
    n_grid = int((pad + LOCUS_LEN) / spacing)
    i = 0
    meth_n = np.zeros(LOCUS_LEN, dtype=np.int64)
    meth_m = np.zeros(LOCUS_LEN, dtype=np.int64)
    for k in range(n_grid):
        h = 1 + k % 2
        tseq, tpos = padded[h]
        ln = long_read_lengths(rng, kind)
        s = int((k + rng.random()) * spacing)
        s = min(s, len(tseq) - ln)
        tmpl, rpos = tseq[s:s + ln], tpos[s:s + ln]
        overlap = int((rpos >= 0).sum())
        if overlap < 1_000:
            continue
        if kind == "hifi":
            err = 10 ** np.clip(rng.normal(math.log10(5e-4), 0.25), -4, math.log10(3e-3))
            ps, pi, pdl = 0.2 * err, 0.4 * err, 0.4 * err
            hp = HIFI_HP
        else:
            err = 10 ** np.clip(rng.normal(-2.0, 0.18), math.log10(0.003), math.log10(0.03))
            ps, pi, pdl = 0.35 * err, 0.20 * err, 0.45 * err
            hp = ONT_HP
        rev = rng.random() < 0.5
        kd, sb, inn, ib = draw_errors(rng, tmpl, ps, pi, pdl, hp)
        aln = project(tmpl, rpos, ref, kd, sb, inn, ib)
        if aln is None:
            continue
        i += 1
        add_single(L, rng, f"{lib}_{i:05d}", aln, rev, True, read_qual(err, len(aln[2])))
        # methylation calls from the same molecule
        on = rpos[rpos >= 0]
        lo, hi = int(on.min()), int(on.max())
        meth_n[lo:hi + 1] += 1
        call_error = 0.03 if kind == "hifi" else 0.05
        cs = CPGS[(CPGS >= lo) & (CPGS <= hi)]
        state = rng.random(len(cs)) < CPG_P[cs]                      # methylated on this molecule?
        state ^= rng.random(len(cs)) < call_error                    # platform call errors
        np.add.at(meth_m, cs, state.astype(np.int64))
    L.meth = (meth_n, meth_m)
    return L


def transcript_template(ref, tid, hap):
    gene, exons = TRANSCRIPTS[tid]
    pos0 = np.concatenate([np.arange(s - 1, e) for s, e in exons])
    seq = ref[pos0].copy()
    if gene == "ASE1" and hap == SNV1["hap"]:
        seq[pos0 == SNV1["pos"] - 1] = ord(SNV1["alt"])
    return seq, pos0


def pre_mrna_template(ref, gene, hap):
    s, e = GENES[gene]
    pos0 = np.arange(s - 1, e)
    seq = ref[pos0].copy()
    for snv in (SNV1, SNV2):
        if snv["hap"] == hap and s <= snv["pos"] <= e:
            seq[snv["pos"] - s] = ord(snv["alt"])
    return seq, pos0


def pick_molecule(rng, by_length):
    """A transcript and a haplotype. by_length: short-read fragments are drawn in proportion to
    molecules x length; long reads count molecules."""
    tids, w = [], []
    for tid, (gene, exons) in TRANSCRIPTS.items():
        share = ISOFORM_SHARE.get(tid, 1.0)
        length = sum(e - s + 1 for s, e in exons)
        tids.append(tid)
        w.append(ABUNDANCE[gene] * share * (length if by_length else 1))
    w = np.array(w) / np.sum(w)
    tid = tids[int(rng.choice(len(tids), p=w))]
    gene = TRANSCRIPTS[tid][0]
    if gene == "ASE1":
        hap = 1 if rng.random() < ASE_HAP1_SHARE else 2
    else:
        hap = 1 + int(rng.random() < 0.5)
    return tid, gene, hap


def sim_rna_illumina(rng, ref, n_frags=10_000, read_len=100):
    L = Library("RNA-ILMN")
    p_cyc, qual = illumina_error_profile(read_len)
    for i in range(1, n_frags + 1):
        tid, gene, hap = pick_molecule(rng, by_length=True)
        if rng.random() < 0.03:                                  # unspliced pre-mRNA
            tseq, tpos = pre_mrna_template(ref, gene, hap)
        else:
            tseq, tpos = transcript_template(ref, tid, hap)
        T = len(tseq)
        fl = int(np.clip(np.rint(rng.normal(250, 50)), 120, min(500, T)))
        a = 0.2                                                  # mild 3' bias: 1.2x at the 3' end
        x = (-1 + math.sqrt(1 + 2 * a * rng.random() * (1 + a / 2))) / a
        s = int(round(x * (T - fl)))
        mates = []
        for mate in (0, 1):                                      # dUTP: read 1 is antisense
            fwd = mate == 1
            a0 = s if fwd else s + fl - read_len
            tmpl, rpos = tseq[a0:a0 + read_len], tpos[a0:a0 + read_len]
            p = p_cyc if fwd else p_cyc[::-1]
            k, sb, inn, ib = draw_errors(rng, tmpl, p, 0.0, 0.0)
            aln = project(tmpl, rpos, ref, k, sb, inn, ib, spliced=True)
            mates.append((aln, not fwd, np.where(k == 1, sb, tmpl), qual if fwd else qual[::-1]))
        add_pair(L, rng, f"RNA-ILMN_{i:06d}", mates[0][:3], mates[1][:3], mates[0][3], mates[1][3], False,
                 max_span=None)
    return L


def sim_rna_long(rng, ref, lib, n_reads):
    L = Library(lib)
    i = 0
    while i < n_reads:
        tid, gene, hap = pick_molecule(rng, by_length=False)
        tseq, tpos = transcript_template(ref, tid, hap)
        T = len(tseq)
        if lib == "RNA-KINNEX":                                  # full-length cDNA, 15% 5'-truncated
            start = 0 if rng.random() < 0.85 else int(rng.uniform(0.05, 0.8) * T)
            end = T
            err = 10 ** np.clip(rng.normal(math.log10(5e-4), 0.25), -4, math.log10(3e-3))
            ps, pi, pdl, hp = 0.2 * err, 0.4 * err, 0.4 * err, HIFI_HP
        else:                                                    # direct RNA: starts at the 3' end
            ln = int(min(T, 100 + rng.exponential(1_750)))
            start, end = T - ln, T
            err = 10 ** np.clip(rng.normal(math.log10(0.03), 0.10), math.log10(0.015), math.log10(0.06))
            ps, pi, pdl, hp = 0.30 * err, 0.20 * err, 0.50 * err, ONT_HP
        tmpl, rpos = tseq[start:end], tpos[start:end]
        k, sb, inn, ib = draw_errors(rng, tmpl, ps, pi, pdl, hp)
        aln = project(tmpl, rpos, ref, k, sb, inn, ib, spliced=True)
        if aln is None:
            continue
        i += 1
        add_single(L, rng, f"{lib}_{i:05d}", aln, False, False, read_qual(err, len(aln[2])))
    return L


# ---------------------------------------------------------------- writing
SAMPLES = [
    ("DNA-ILMN", "DNA", "EDTA blood", "Illumina NovaSeq X", "2 x 150 bp paired-end",
     "paired-end 2 x 150 bp; insert ~350 bp", "PCR-amplified library (shows GC bias)", "~30x"),
    ("DNA-HIFI", "DNA", "EDTA blood", "PacBio Revio", "HiFi reads, 15-20 kb library",
     "single long reads", "PCR-free; 5mC from polymerase kinetics", "~20x"),
    ("DNA-ONT", "DNA", "EDTA blood", "Oxford Nanopore PromethION",
     "R10.4.1 flow cell, Ligation Sequencing Kit V14", "single long reads",
     "PCR-free; 5mC from the current signal", "~20x"),
    ("RNA-ILMN", "RNA", "PAXgene blood RNA", "Illumina NovaSeq X", "2 x 100 bp, stranded poly(A) mRNA",
     "paired-end 2 x 100 bp; read 1 antisense", "cDNA fragments ~250 bp", "~10,000 fragments"),
    ("RNA-KINNEX", "RNA", "PAXgene blood RNA", "PacBio Revio",
     "Kinnex full-length RNA (MAS-Seq concatenation of Iso-Seq cDNA; reads segmented back to single cDNAs)",
     "single full-length cDNA reads", "full-length cDNA, poly(A) trimmed", "~1,500 reads"),
    ("RNA-ONT", "RNA", "PAXgene blood RNA", "Oxford Nanopore PromethION",
     "Direct RNA, SQK-RNA004, RNA flow cell", "single native RNA reads, read 3' to 5'",
     "native RNA, no cDNA; written in DNA letters (U as T)", "~2,000 reads"),
]
PLATFORM_TAG = {"DNA-ILMN": "ILLUMINA", "RNA-ILMN": "ILLUMINA", "DNA-HIFI": "PACBIO",
                "RNA-KINNEX": "PACBIO", "DNA-ONT": "ONT", "RNA-ONT": "ONT"}


def write_gz_text(path, text):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(text.encode())


def write_sam(L, path, sample_row):
    lib, molecule, tube, platform, chem, layout, prep, amount = sample_row
    head = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:{CHROM}\tLN:{LOCUS_LEN}",
            f"@RG\tID:{lib}\tSM:{PERSON}\tPL:{PLATFORM_TAG[lib]}\tLB:{lib}\tDS:{molecule}; {chem}",
            f"@PG\tID:build_practice3\tPN:build_practice3.py\tVN:1\tDS:synthetic alignments projected "
            f"from simulated haplotypes (teaching data, not a real aligner)",
            "@CO\tSynthetic teaching data for Mini practice 3. Not calibrated to any real run."]
    recs = sorted(L.records, key=lambda r: (r["pos"] == 0, r["pos"], r["qname"], r["flag"] & 0xC0))
    lines = []
    for r in recs:
        rname = CHROM if r["pos"] else "*"
        seq = r["seq"].tobytes().decode() if isinstance(r["seq"], np.ndarray) else r["seq"]
        fields = [r["qname"], str(r["flag"]), rname, str(r["pos"]), str(r["mapq"]), r["cigar"],
                  r["rnext"], str(r["pnext"]), str(r["tlen"]), seq, r["qual"]]
        if r["nm"] is not None:
            fields.append(f"NM:i:{r['nm']}")
        lines.append("\t".join(fields))
    write_gz_text(path, "\n".join(head + lines) + "\n")


def write_annotation(ref):
    fa = D3 / "chrS.fa"
    s = ref.tobytes().decode()
    fa.write_text(f">{CHROM} synthetic teaching locus, 120 kb\n" + "\n".join(s[i:i + 60] for i in range(0, len(s), 60)) + "\n")
    offset = len(f">{CHROM} synthetic teaching locus, 120 kb\n")
    (D3 / "chrS.fa.fai").write_text(f"{CHROM}\t{LOCUS_LEN}\t{offset}\t60\t61\n")
    gtf = []
    for tid, (gene, exons) in TRANSCRIPTS.items():
        gs, ge = GENES[gene]
        attrs = f'gene_id "{gene}"; transcript_id "{tid}";'
        gtf.append(f"{CHROM}\tsynthetic\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t+\t.\t{attrs}")
        names = ({v: k for k, v in ALTX_EXONS.items()} if gene == "ALTX"
                 else {x: f"E{k}" for k, x in enumerate(exons, 1)})
        for k, (es, ee) in enumerate(exons, 1):
            gtf.append(f"{CHROM}\tsynthetic\texon\t{es}\t{ee}\t.\t+\t.\t{attrs} exon_number \"{k}\"; "
                       f"exon_id \"{gene}-{names[(es, ee)]}\";")
    (D3 / "genes.gtf").write_text("\n".join(gtf) + "\n")
    regions = [
        ("CGI_ALTX", *CGI_ALTX, "CpG island", "GC-rich, CpG-dense promoter of ALTX"),
        ("REP1", *REP1, "short tandem repeat", "(CAG)n repeat; 20 units in the reference"),
        ("HP1", *HP1, "homopolymer", "ten A bases in a row"),
        ("DUP_A", *DUP_A, "segmental duplication", "6-kb block with a near-identical copy (DUP_B)"),
        ("DUP_B", *DUP_B, "segmental duplication", "copy of DUP_A that differs at 3 bases"),
        ("OFF1_PROMOTER", *OFF1_PROMOTER, "promoter", "GC-rich promoter of OFF1"),
    ]
    pd.DataFrame(regions, columns=["name", "start", "end", "kind", "description"]).to_csv(
        D3 / "regions.tsv", sep="\t", index=False)
    pd.DataFrame([(v["name"], v["pos"], v["ref"], v["alt"], v["note"]) for v in (SNV1, SNV2)],
                 columns=["name", "pos", "ref", "alt", "note"]).to_csv(D3 / "sites.tsv", sep="\t", index=False)


CPGS = None
CPG_P = None


def cpg_table(ref):
    c0 = np.flatnonzero((ref[:-1] == ord("C")) & (ref[1:] == ord("G")))
    p = np.full(LOCUS_LEN, METH_DEFAULT)
    for _, (s, e), val in reversed(METHYLATION):        # promoters win over the gene body
        p[s - 1:e] = val
    return c0, p


# ---------------------------------------------------------------- reference CHECK numbers
# An independent re-implementation of the notebook's measurements, computed from the SAM files
# that were just written. checkpoints/practice3_checks.tsv must equal the notebook's output.
import re  # noqa: E402

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
TX_LEN = {g: sum(e - s + 1 for s, e in TRANSCRIPTS[t][1]) for g, t in
          (("ALTX", "ALTX-201"), ("ASE1", "ASE1-201"), ("LOW1", "LOW1-201"), ("OFF1", "OFF1-201"))}
GENE_ORDER = ["ASE1", "ALTX", "LOW1", "OFF1"]
DNA_LIBS, RNA_LIBS = LIBRARIES[:3], LIBRARIES[3:]


class Rec:
    __slots__ = ("qname", "flag", "pos", "mapq", "cigar", "seq", "nm", "ops", "blocks", "end",
                 "n_m", "n_i", "n_d")

    def __init__(self, f):
        self.qname, self.flag, self.pos, self.mapq = f[0], int(f[1]), int(f[3]), int(f[4])
        self.cigar, self.seq = f[5], f[9]
        self.nm = next((int(x[5:]) for x in f[11:] if x.startswith("NM:i:")), None)
        self.ops = [(int(n), o) for n, o in CIG.findall(self.cigar)] if self.cigar != "*" else []
        self.blocks, r, q = [], self.pos, 0
        self.n_m = self.n_i = self.n_d = 0
        for n, o in self.ops:
            if o == "M":
                self.blocks.append((r, r + n - 1, q))
                r += n
                q += n
                self.n_m += n
            elif o in "DN":
                r += n
                if o == "D":
                    self.n_d += n
            elif o in "IS":
                q += n
                if o == "I":
                    self.n_i += n
        self.end = r - 1

    @property
    def mapped(self):
        return not self.flag & 0x4

    def qidx(self, refp):
        for rs, re_, qs in self.blocks:
            if rs <= refp <= re_:
                return qs + refp - rs
        return None


def load(lib):
    out = []
    with gzip.open(D3 / f"{lib}.sam.gz", "rt") as fh:
        for line in fh:
            if not line.startswith("@"):
                out.append(Rec(line.rstrip("\n").split("\t")))
    return out


def depth(recs, min_mapq=0):
    diff = np.zeros(LOCUS_LEN + 2, dtype=np.int64)
    for r in recs:
        if r.mapped and r.mapq >= min_mapq:
            for rs, re_, _ in r.blocks:
                diff[rs] += 1
                diff[re_ + 1] -= 1
    return np.cumsum(diff)[1:LOCUS_LEN + 1]          # index i = position i + 1


def n50(lengths):
    s = sorted(lengths, reverse=True)
    half, run = sum(s) / 2, 0
    for x in s:
        run += x
        if run >= half:
            return x


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def cag_identity(s):
    best = 0.0
    for motif in ("CAG", "CTG"):
        if s.count(motif) < 0.25 * len(s):           # quick reject
            continue
        for f in range(3):
            pattern = (motif * (len(s) // 3 + 3))[f:f + len(s)]
            best = max(best, sum(a == c for a, c in zip(s, pattern)) / len(s))
    return best


def exon_positions(gene):
    pos = set()
    for tid, (g, exons) in TRANSCRIPTS.items():
        if g == gene:
            for s, e in exons:
                pos.update(range(s, e + 1))
    return np.array(sorted(pos))


def reference_checks():
    data = {lib: load(lib) for lib in LIBRARIES}
    rows = []

    def put(fig, lib, metric, value):
        rows.append((fig, lib, metric, value))

    dep = {lib: depth(data[lib]) for lib in LIBRARIES}
    # F01 read length
    for lib in LIBRARIES:
        lens = [len(r.seq) for r in data[lib]]
        put("F01", lib, "n_reads", len(lens))
        put("F01", lib, "median_length", round(float(np.median(lens)), 1))
        put("F01", lib, "read_n50", n50(lens))
    # F02 errors outside true variants, and HP1
    skip = lambda r: r.pos <= REP1[1] and r.end >= REP1[0] or any(r.pos <= v["pos"] <= r.end for v in (SNV1, SNV2))
    for lib in LIBRARIES:
        use = [r for r in data[lib] if r.mapped and not skip(r)]
        al = sum(r.n_m for r in use)
        ins = sum(r.n_i for r in use)
        dele = sum(r.n_d for r in use)
        mm = sum(r.nm - r.n_i - r.n_d for r in use)
        put("F02", lib, "reads_used", len(use))
        put("F02", lib, "mismatch_per_kb", round(1000 * mm / al, 3))
        put("F02", lib, "insertion_per_kb", round(1000 * ins / al, 3))
        put("F02", lib, "deletion_per_kb", round(1000 * dele / al, 3))
        put("F02", lib, "accuracy_pct", round(100 * (1 - (mm + ins + dele) / al), 3))
    for lib in DNA_LIBS:
        hl = []
        for r in data[lib]:
            if r.mapped:
                a, z = r.qidx(HP1[0] - 1), r.qidx(HP1[1] + 1)
                if a is not None and z is not None:
                    hl.append(z - a - 1)
        hl = np.array(hl)
        put("F02", lib, "hp1_reads", len(hl))
        for k in (9, 10, 11):
            put("F02", lib, f"hp1_pct_{k}A", round(100 * float((hl == k).mean()), 1))
    # F03 DNA vs RNA coverage
    for lib in LIBRARIES:
        for g in GENE_ORDER:
            ex = exon_positions(g)
            put("F03", lib, f"{g}_exon_depth", round(float(dep[lib][ex - 1].mean()), 2))
        allex = np.concatenate([exon_positions(g) for g in GENE_ORDER])
        put("F03", lib, "pct_bases_in_exons", round(100 * float(dep[lib][allex - 1].sum() / dep[lib].sum()), 1))
    # F04 allele fractions at SNV1
    for lib in LIBRARIES:
        ref_n = alt_n = 0
        for r in data[lib]:
            if r.mapped:
                q = r.qidx(SNV1["pos"])
                if q is not None:
                    ref_n += r.seq[q] == SNV1["ref"]
                    alt_n += r.seq[q] == SNV1["alt"]
        lo, hi = wilson(alt_n, ref_n + alt_n)
        put("F04", lib, "n_ref", ref_n)
        put("F04", lib, "n_alt", alt_n)
        put("F04", lib, "alt_fraction", round(alt_n / (ref_n + alt_n), 3))
        put("F04", lib, "ci_low", round(lo, 3))
        put("F04", lib, "ci_high", round(hi, 3))
    # F05 junctions and exon chains
    E = ALTX_EXONS
    junctions = [("E1", "E2a"), ("E1", "E2b"), ("E2a", "E3"), ("E2b", "E3"),
                 ("E3", "E4a"), ("E3", "E4b"), ("E4a", "E5"), ("E4b", "E5")]
    jcount = {j: 0 for j in junctions}
    lookup = {(E[a][1], E[c][0]): (a, c) for a, c in junctions}
    for r in data["RNA-ILMN"]:
        if not r.mapped:
            continue
        for k, (n, o) in enumerate(r.ops):
            if o == "N" and 0 < k < len(r.ops) - 1:
                left = sum(m for m, oo in r.ops[:k] if oo in "MDN")
                donor = r.pos + left - 1
                acceptor = donor + n + 1
                if (donor, acceptor) in lookup and r.ops[k - 1][1] == "M" and r.ops[k + 1][1] == "M" \
                        and r.ops[k - 1][0] >= 8 and r.ops[k + 1][0] >= 8:
                    jcount[lookup[(donor, acceptor)]] += 1
    for a, c in junctions:
        put("F05", "RNA-ILMN", f"junction_{a}_{c}", jcount[(a, c)])
    for lib in RNA_LIBS:
        cats = {"E2a+E4a": 0, "E2b+E4b": 0, "E2a+E4b": 0, "E2b+E4a": 0, "incomplete": 0}
        for r in data[lib]:
            if not r.mapped or r.end < GENES["ALTX"][0] or r.pos > GENES["ALTX"][1]:
                continue
            has = {}
            for ex in ("E2a", "E2b", "E4a", "E4b"):
                s, e = E[ex]
                has[ex] = sum(max(0, min(re_, e) - max(rs, s) + 1) for rs, re_, _ in r.blocks) >= 20
            two = [x for x in ("E2a", "E2b") if has[x]]
            four = [x for x in ("E4a", "E4b") if has[x]]
            if len(two) == 1 and len(four) == 1:
                cats[f"{two[0]}+{four[0]}"] += 1
            else:
                cats["incomplete"] += 1
        for k, v in cats.items():
            put("F05", lib, "chain_" + k.replace("+", "_"), v)
    # F06 coverage along the ASE1 transcript
    exons = TRANSCRIPTS["ASE1-201"][1]
    tx_index = {}
    for s, e in exons:
        for x in range(s, e + 1):
            tx_index[x] = len(tx_index)
    T = len(tx_index)
    for lib in RNA_LIBS:
        cov = np.zeros(T)
        n_gene = reach = 0
        for r in data[lib]:
            if not r.mapped or r.end < GENES["ASE1"][0] or r.pos > GENES["ASE1"][1]:
                continue
            idx = [tx_index[x] for rs, re_, _ in r.blocks for x in range(rs, re_ + 1) if x in tx_index]
            if not idx:
                continue
            np.add.at(cov, idx, 1)
            n_gene += 1
            reach += min(idx) < 50
        k = int(T * 0.2)
        put("F06", lib, "ratio_3to5", round(float(cov[-k:].mean() / cov[:k].mean()), 2))
        put("F06", lib, "pct_reads_reach_5prime", "NA" if lib == "RNA-ILMN" else round(100 * reach / n_gene, 1))
    # F07 the CAG repeat
    for lib in DNA_LIBS:
        units = []
        for r in data[lib]:
            if r.mapped and all(r.qidx(x) is not None for x in (REP1[0] - 10, REP1[0] - 1, REP1[1] + 1, REP1[1] + 10)):
                units.append((r.qidx(REP1[1] + 1) - r.qidx(REP1[0] - 1) - 1) / 3)
        units = np.array(units)
        put("F07", lib, "spanning_reads", len(units))
        for name, sel in (("short", units < 100), ("long", units >= 100)):
            u = units[sel]
            put("F07", lib, f"{name}_n", len(u))
            if len(u):
                vals, cnt = np.unique(np.floor(u + 0.5).astype(int), return_counts=True)
                put("F07", lib, f"{name}_mode", int(vals[np.argmax(cnt)]))
                put("F07", lib, f"{name}_sd", round(float(u.std()), 2))
            else:
                put("F07", lib, f"{name}_mode", "NA")
                put("F07", lib, f"{name}_sd", "NA")
        put("F07", lib, "all_cag_reads", sum(len(r.seq) >= 30 and cag_identity(r.seq) >= 0.95 for r in data[lib]))
    # F08 the duplicated block
    flank = np.r_[np.arange(53_001, 56_001), np.arange(62_001, 65_001)]
    inside = np.arange(DUP_A[0], DUP_A[1] + 1)
    for lib in DNA_LIBS:
        d20 = depth(data[lib], 20)
        a, f = float(d20[inside - 1].mean()), float(d20[flank - 1].mean())
        put("F08", lib, "dupA_depth", round(a, 2))
        put("F08", lib, "flank_depth", round(f, 2))
        put("F08", lib, "ratio_inside_flank", round(a / f, 3))
    # F09 phasing
    for lib in DNA_LIBS:
        combo = {"alt_ref": 0, "ref_alt": 0, "ref_ref": 0, "alt_alt": 0, "other": 0}
        for r in data[lib]:
            if r.mapped:
                q1, q2 = r.qidx(SNV1["pos"]), r.qidx(SNV2["pos"])
                if q1 is not None and q2 is not None:
                    a1 = "alt" if r.seq[q1] == SNV1["alt"] else "ref" if r.seq[q1] == SNV1["ref"] else None
                    a2 = "alt" if r.seq[q2] == SNV2["alt"] else "ref" if r.seq[q2] == SNV2["ref"] else None
                    combo[f"{a1}_{a2}" if a1 and a2 else "other"] += 1
        put("F09", lib, "reads_covering_both", sum(combo.values()))
        for k, v in combo.items():
            put("F09", lib, k, v)
    # F10 GC and coverage, 1-kb windows
    fasta = "".join(x.strip() for x in (D3 / "chrS.fa").read_text().splitlines()[1:])
    gc = np.array([(fasta[i:i + 1000].count("G") + fasta[i:i + 1000].count("C")) / 1000
                   for i in range(0, LOCUS_LEN, 1000)])
    for lib in DNA_LIBS:
        w = dep[lib].reshape(-1, 1000).mean(axis=1)
        norm = w / np.median(w)
        hi, mid = gc >= 0.60, (gc >= 0.40) & (gc < 0.50)
        put("F10", lib, "n_windows_gc60", int(hi.sum()))
        put("F10", lib, "n_windows_gc40_50", int(mid.sum()))
        put("F10", lib, "norm_depth_gc60", round(float(norm[hi].mean()), 3))
        put("F10", lib, "norm_depth_gc40_50", round(float(norm[mid].mean()), 3))
        put("F10", lib, "ratio", round(float(norm[hi].mean() / norm[mid].mean()), 3))
        put("F10", lib, "sd_gc40_50", round(float(norm[mid].std()), 3))
    # F11 methylation
    for lib in ("DNA-HIFI", "DNA-ONT"):
        tab = pd.read_csv(D3 / f"methylation_{lib}.tsv", sep="\t")
        tab = tab[tab.n_reads >= 5]
        region = np.full(len(tab), "other", dtype=object)
        for name, (s, e), _ in reversed(METHYLATION):
            region[(tab.pos >= s).to_numpy() & (tab.pos <= e).to_numpy()] = name
        for name in ("CGI_ALTX", "OFF1_PROMOTER", "ALTX_body", "other"):
            sel = region == name
            put("F11", lib, f"n_cpg_{name}", int(sel.sum()))
            put("F11", lib, f"mean_{name}", round(float(tab.fraction.to_numpy()[sel].mean()), 3))
    # F12 expression
    gene_exons = {g: sorted({ex for tid, (gg, exs) in TRANSCRIPTS.items() if gg == g for ex in exs})
                  for g in GENE_ORDER}
    for lib in RNA_LIBS:
        counts = {g: 0 for g in GENE_ORDER}
        for r in data[lib]:
            if not r.mapped or (lib == "RNA-ILMN" and not r.flag & 0x40):
                continue
            for g in GENE_ORDER:
                if any(min(re_, e) >= max(rs, s) for rs, re_, _ in r.blocks for s, e in gene_exons[g]):
                    counts[g] += 1
                    break
        per = {g: counts[g] / (TX_LEN[g] if lib == "RNA-ILMN" else 1) for g in GENE_ORDER}
        tot = sum(per.values())
        for g in GENE_ORDER:
            share = per[g] / tot
            c = counts[g]
            if c:
                lo, hi = share * max(0.0, 1 - 1.96 / math.sqrt(c)), share * (1 + 1.96 / math.sqrt(c))
            else:
                lo, hi = 0.0, (3 / (TX_LEN[g] if lib == "RNA-ILMN" else 1)) / tot
            put("F12", lib, f"{g}_reads", c)
            put("F12", lib, f"{g}_share", round(share, 4))
            put("F12", lib, f"{g}_share_lo", round(lo, 4))
            put("F12", lib, f"{g}_share_hi", round(hi, 4))
    checks = pd.DataFrame(rows, columns=["figure", "library", "metric", "value"])
    return checks, data


def reference_scorecard(checks):
    v = {(f, l, m): x for f, l, m, x in checks.itertuples(index=False)}
    rows = []
    for lib in LIBRARIES:
        dna = lib.startswith("DNA")
        row = dict(library=lib, molecule="DNA" if dna else "RNA",
                   platform={"ILMN": "Illumina", "HIFI": "PacBio", "KINNEX": "PacBio", "ONT": "ONT"}[lib.split("-")[1]],
                   read_n50=v[("F01", lib, "read_n50")],
                   accuracy_pct=round(float(v[("F02", lib, "accuracy_pct")]), 2))
        if dna:
            row["sizes_expanded_repeat"] = "yes" if v[("F07", lib, "long_n")] >= 3 else "no"
            row["resolves_duplication"] = "yes" if v[("F08", lib, "ratio_inside_flank")] >= 0.5 else "no"
            row["phases_snv1_snv2"] = "yes" if v[("F09", lib, "alt_ref")] + v[("F09", lib, "ref_alt")] >= 3 else "no"
            row["methylation_from_reads"] = "no" if lib == "DNA-ILMN" else "yes"
            row["allele_specific_expression"] = row["full_length_isoforms"] = row["gene_expression"] = "n/a"
        else:
            lo, hi = v[("F04", lib, "ci_low")], v[("F04", lib, "ci_high")]
            row["allele_specific_expression"] = "yes" if lo > 0.5 or hi < 0.5 else "no"
            full = sum(v[("F05", lib, f"chain_{k}")] for k in ("E2a_E4a", "E2b_E4b", "E2a_E4b", "E2b_E4a"))
            row["full_length_isoforms"] = "yes" if full >= 10 else "no"
            share = [v[("F12", lib, f"{g}_share")] for g in ("ASE1", "ALTX", "LOW1")]
            row["gene_expression"] = "yes" if share[0] > share[1] > share[2] and v[("F12", lib, "OFF1_reads")] == 0 else "no"
            for k in ("sizes_expanded_repeat", "resolves_duplication", "phases_snv1_snv2", "methylation_from_reads"):
                row[k] = "n/a"
        rows.append(row)
    cols = ["library", "molecule", "platform", "read_n50", "accuracy_pct", "sizes_expanded_repeat",
            "resolves_duplication", "phases_snv1_snv2", "methylation_from_reads",
            "allele_specific_expression", "full_length_isoforms", "gene_expression"]
    return pd.DataFrame(rows)[cols]


# ---------------------------------------------------------------- the teaching facts must stay visible
def teaching_facts(checks):
    """Each fact the figures teach, as a test on the reference CHECK numbers. Returns
    [(description, passed)]; the build stops if any fails."""
    v = {(f, l, m): x for f, l, m, x in checks.itertuples(index=False)}
    g = lambda f, l, m: float(v[(f, l, m)])
    facts = [
        ("F01 short reads have fixed lengths (150 / 100 bp)",
         g("F01", "DNA-ILMN", "read_n50") == 150 and g("F01", "RNA-ILMN", "read_n50") == 100),
        ("F01 long DNA reads: N50 14-30 kb; ONT N50 well above its median",
         14_000 <= g("F01", "DNA-HIFI", "read_n50") <= 30_000 and 14_000 <= g("F01", "DNA-ONT", "read_n50") <= 30_000
         and g("F01", "DNA-ONT", "read_n50") > 1.5 * g("F01", "DNA-ONT", "median_length")),
        ("F01 long RNA reads are no longer than the transcripts (N50 <= 2,250 nt)",
         g("F01", "RNA-KINNEX", "read_n50") <= 2_250 and g("F01", "RNA-ONT", "read_n50") <= 2_250),
        ("F02 Illumina and HiFi >= 99.5% / 99.85% accurate; ONT DNA 98.5-99.5%; ONT RNA least accurate",
         g("F02", "DNA-ILMN", "accuracy_pct") >= 99.5 and g("F02", "DNA-HIFI", "accuracy_pct") >= 99.85
         and 98.5 <= g("F02", "DNA-ONT", "accuracy_pct") <= 99.5
         and g("F02", "RNA-ONT", "accuracy_pct") < min(g("F02", l, "accuracy_pct") for l in LIBRARIES[:5])),
        ("F02 ONT errors are mostly deletions; Illumina makes no indels",
         all(g("F02", l, "deletion_per_kb") > max(g("F02", l, "mismatch_per_kb"), g("F02", l, "insertion_per_kb"))
             for l in ("DNA-ONT", "RNA-ONT"))
         and g("F02", "DNA-ILMN", "insertion_per_kb") == 0 and g("F02", "DNA-ILMN", "deletion_per_kb") == 0),
        ("F02 at HP1, ONT reads 9 or 11 A far more often than HiFi or Illumina",
         g("F02", "DNA-ONT", "hp1_pct_9A") + g("F02", "DNA-ONT", "hp1_pct_11A") >= 25
         and g("F02", "DNA-ILMN", "hp1_pct_10A") >= 95 and g("F02", "DNA-HIFI", "hp1_pct_10A") >= 80),
        ("F03 OFF1 has normal DNA depth and zero RNA depth; RNA bases sit in exons",
         all(g("F03", l, "OFF1_exon_depth") >= 10 for l in DNA_LIBS)
         and all(g("F03", l, "OFF1_exon_depth") == 0 and g("F03", l, "pct_bases_in_exons") >= 95 for l in RNA_LIBS)
         and all(g("F03", l, "pct_bases_in_exons") <= 15 for l in DNA_LIBS)),
        ("F04 SNV1: every DNA interval contains 0.5; every RNA interval excludes 0.5 and contains 0.85",
         all(g("F04", l, "ci_low") <= 0.5 <= g("F04", l, "ci_high") for l in DNA_LIBS)
         and all(g("F04", l, "ci_low") > 0.5 and g("F04", l, "ci_low") <= 0.85 <= g("F04", l, "ci_high")
                 for l in RNA_LIBS)),
        ("F05 junctions ~60/40 at both alternative exons; no complete chain in short reads; no mixed chain",
         all(0.5 <= g("F05", "RNA-ILMN", f"junction_{a}") / (g("F05", "RNA-ILMN", f"junction_{a}")
                                                             + g("F05", "RNA-ILMN", f"junction_{b_}")) <= 0.7
             for a, b_ in (("E1_E2a", "E1_E2b"), ("E2a_E3", "E2b_E3"), ("E3_E4a", "E3_E4b"), ("E4a_E5", "E4b_E5")))
         and sum(g("F05", "RNA-ILMN", f"chain_{k}") for k in ("E2a_E4a", "E2b_E4b", "E2a_E4b", "E2b_E4a")) == 0
         and all(g("F05", l, "chain_E2a_E4b") == 0 and g("F05", l, "chain_E2b_E4a") == 0
                 and g("F05", l, "chain_E2a_E4a") > g("F05", l, "chain_E2b_E4b") > 30 for l in ("RNA-KINNEX", "RNA-ONT"))),
        ("F06 direct RNA is 3'-biased (ratio >= 2, 20-40% reach the 5' end); Kinnex and Illumina near flat",
         g("F06", "RNA-ONT", "ratio_3to5") >= 2 and 20 <= g("F06", "RNA-ONT", "pct_reads_reach_5prime") <= 40
         and g("F06", "RNA-KINNEX", "ratio_3to5") <= 1.3 and g("F06", "RNA-KINNEX", "pct_reads_reach_5prime") >= 75
         and 1.05 <= g("F06", "RNA-ILMN", "ratio_3to5") <= 1.4),
        ("F07 HiFi and ONT size both alleles (25 and ~180 units), ONT wider; Illumina sees only 25 plus all-CAG reads",
         all(g("F07", l, "short_mode") == 25 and 177 <= g("F07", l, "long_mode") <= 182
             and g("F07", l, "short_n") >= 3 and g("F07", l, "long_n") >= 3 for l in ("DNA-HIFI", "DNA-ONT"))
         and g("F07", "DNA-ONT", "long_sd") > g("F07", "DNA-HIFI", "long_sd")
         and g("F07", "DNA-ILMN", "long_n") == 0 and g("F07", "DNA-ILMN", "short_n") >= 3
         and g("F07", "DNA-ILMN", "all_cag_reads") >= 10),
        ("F08 at MAPQ >= 20 Illumina loses DUP_A (ratio < 0.3); long reads keep it (>= 0.8)",
         g("F08", "DNA-ILMN", "ratio_inside_flank") < 0.3
         and all(g("F08", l, "ratio_inside_flank") >= 0.8 for l in ("DNA-HIFI", "DNA-ONT"))),
        ("F09 no Illumina read covers both SNVs; long reads show only alt-ref and ref-alt (trans)",
         g("F09", "DNA-ILMN", "reads_covering_both") == 0
         and all(g("F09", l, "alt_ref") >= 3 and g("F09", l, "ref_alt") >= 3
                 and g("F09", l, "ref_ref") + g("F09", l, "alt_alt") == 0 for l in ("DNA-HIFI", "DNA-ONT"))),
        ("F10 Illumina loses depth in GC >= 60% windows (ratio < 0.8, drop > 3 x scatter); long reads within 0.85-1.15",
         g("F10", "DNA-ILMN", "ratio") < 0.8
         and 1 - g("F10", "DNA-ILMN", "ratio") > 3 * g("F10", "DNA-ILMN", "sd_gc40_50")
         and all(0.85 <= g("F10", l, "ratio") <= 1.15 for l in ("DNA-HIFI", "DNA-ONT"))
         and g("F10", "DNA-ILMN", "n_windows_gc60") >= 10),
        ("F11 ALTX promoter low (< 0.15) and OFF1 promoter high (> 0.75) in both long-read tables",
         all(g("F11", l, "mean_CGI_ALTX") < 0.15 and g("F11", l, "mean_OFF1_PROMOTER") > 0.75
             for l in ("DNA-HIFI", "DNA-ONT"))),
        ("F12 all RNA datasets rank ASE1 > ALTX > LOW1 > OFF1 = 0, shares within 0.05 of 10:4:1",
         all(g("F12", l, "ASE1_share") > g("F12", l, "ALTX_share") > g("F12", l, "LOW1_share") > 0
             and g("F12", l, "OFF1_reads") == 0
             and abs(g("F12", l, "ASE1_share") - 10 / 15) < 0.05 and abs(g("F12", l, "ALTX_share") - 4 / 15) < 0.05
             for l in RNA_LIBS)),
    ]
    return facts


def write_truth(ref, data, checks, facts):
    fasta = ref.tobytes().decode()
    gc5 = [round((fasta[i:i + 5000].count("G") + fasta[i:i + 5000].count("C")) / 5000, 3)
           for i in range(0, LOCUS_LEN, 5000)]
    region_gc = {name: round((fasta[s - 1:e].count("G") + fasta[s - 1:e].count("C")) / (e - s + 1), 3)
                 for name, (s, e) in (("CGI_ALTX", CGI_ALTX), ("OFF1_PROMOTER", OFF1_PROMOTER))}
    diffs = [dict(offset=off, pos_DUP_A=DUP_A[0] + off, base_DUP_A=fasta[DUP_A[0] + off - 1],
                  pos_DUP_B=DUP_B[0] + off, base_DUP_B=fasta[DUP_B[0] + off - 1]) for off in DUP_OFFSETS]
    depth_mean = {lib: round(float(depth(data[lib]).mean()), 2) for lib in DNA_LIBS}
    truth = {
        "about": "Answer key for Mini practice 3 (instructor use). Synthetic teaching data: one person, "
                 "one 120-kb locus, six libraries. Read models approximate public platform specifications "
                 "but are not calibrated to any real run.",
        "coordinates": "1-based inclusive",
        "seed": SEED,
        "person": PERSON,
        "locus": {"chrom": CHROM, "length": LOCUS_LEN,
                  "gc_5kb_blocks_target": GC_BLOCKS, "gc_5kb_blocks_built": gc5,
                  "gc_regions_built": region_gc, "exon_gc_target": 0.50,
                  "cpg_depletion_keep_fraction": {"background": 0.25, "exons": 0.45, "OFF1_PROMOTER": 0.80,
                                                  "CGI_ALTX": "not depleted"}},
        "genes": {g: {"start": GENES[g][0], "end": GENES[g][1], "strand": "+",
                      "relative_abundance_in_blood": ABUNDANCE[g],
                      "transcripts": {tid: {"exons": [list(x) for x in ex], "length_nt": sum(e - s + 1 for s, e in ex),
                                            "share_of_gene": ISOFORM_SHARE.get(tid, 1.0)}
                                      for tid, (gg, ex) in TRANSCRIPTS.items() if gg == g}}
                  for g in ("ALTX", "ASE1", "OFF1", "LOW1")},
        "ALTX_exons": {k: list(x) for k, x in ALTX_EXONS.items()},
        "ALTX_mixed_chains": "E2a+E4b and E2b+E4a do not exist (0%)",
        "variants": {"SNV1": SNV1, "SNV2": SNV2,
                     "phase": "trans: SNV1 alt on hap1, SNV2 alt on hap2 (3,000 bp apart)"},
        "allele_specific_expression": {"gene": "ASE1", "hap1_share_of_RNA": ASE_HAP1_SHARE,
                                       "expected_RNA_alt_fraction_at_SNV1": ASE_HAP1_SHARE,
                                       "expected_DNA_alt_fraction_at_SNV1": 0.5},
        "REP1": {"start": REP1[0], "end": REP1[1], "motif": "CAG", "reference_units": 20,
                 "hap1_units": REP_UNITS[1], "hap2_units": REP_UNITS[2],
                 "hap2_extra_bases": 3 * (REP_UNITS[2] - 20)},
        "HP1": {"start": HP1[0], "end": HP1[1], "base": "A", "length": 10, "genotype": "10/10",
                "flanks": "GTCTC ... GCTGT"},
        "duplication": {"DUP_A": list(DUP_A), "DUP_B": list(DUP_B), "differences": diffs,
                        "illumina_rule": "a read inside one copy that avoids the 3 differing bases gets MAPQ 0 "
                                         "and is placed in either copy at random"},
        "methylation": {"CGI_ALTX": 0.05, "OFF1_PROMOTER": 0.90, "ALTX_body": 0.80, "other_CpGs": METH_DEFAULT,
                        "per_CpG_call_error": {"DNA-HIFI": 0.03, "DNA-ONT": 0.05}},
        "libraries": {
            "DNA-ILMN": {"pairs": 12_000, "read_length": 150, "insert": "N(350, 50) clipped to 200-600",
                         "error": "substitutions only; 0.1% at base 1 rising to 0.8% at base 150",
                         "qual": "per cycle, Phred of the error model",
                         "gc_bias": "fragment kept with probability exp(-0.5*((gc-0.45)/0.18)^2)",
                         "mean_depth": depth_mean["DNA-ILMN"]},
            "DNA-HIFI": {"reads": int(checks.query("figure=='F01' and library=='DNA-HIFI' and metric=='n_reads'").value.iloc[0]),
                         "length": "N(16 kb, 3 kb) truncated to 8-30 kb",
                         "error": "per-read rate log-normal around 5e-4 (Q33); 20% substitutions, 40% insertions, "
                                  "40% deletions; extra homopolymer length errors",
                         "qual": "one value per read (Phred of its error rate)", "mean_depth": depth_mean["DNA-HIFI"]},
            "DNA-ONT": {"reads": int(checks.query("figure=='F01' and library=='DNA-ONT' and metric=='n_reads'").value.iloc[0]),
                        "length": "log-normal, sigma 0.9, N50 about 22 kb, 1-100 kb",
                        "error": "per-read rate log-normal around 1% (Q20), 0.3-3%; 35% substitutions, "
                                 "20% insertions, 45% deletions; homopolymer length errors (runs >= 5)",
                        "qual": "one value per read (Phred of its error rate)", "mean_depth": depth_mean["DNA-ONT"]},
            "RNA-ILMN": {"fragments": 10_000, "read_length": 100, "fragment": "N(250, 50)",
                         "pre_mRNA_fraction": 0.03, "three_prime_bias": "linear, 1.2x at the 3' end",
                         "strand": "dUTP: read 1 antisense", "error": "as DNA-ILMN (first 100 cycles)"},
            "RNA-KINNEX": {"reads": 1_500, "full_length_fraction": 0.85,
                           "truncated_start": "uniform 5-80% of the transcript", "error": "as DNA-HIFI"},
            "RNA-ONT": {"reads": 2_000, "length": "min(transcript, 100 + Exponential(mean 1,750)) from the 3' end",
                        "error": "per-read rate log-normal around 3% (Q15), 1.5-6%; 30% substitutions, "
                                 "20% insertions, 50% deletions; homopolymer length errors"},
        },
        "sampling": "Long DNA reads start on a jittered grid and alternate between haplotypes, so that "
                    "~20x over one locus is not dominated by sampling noise. Real runs scatter reads at random.",
        "simplifications": [
            "Alignments are projections of the simulated haplotypes, not the output of an aligner.",
            "Illumina reads carry substitutions only; real Illumina data have rare indels.",
            "Long-read QUAL holds one value per read; real HiFi and nanopore files vary base by base.",
            "Methylation tables are per-CpG summaries computed from simulated per-molecule calls.",
            "All four genes are on the + strand; there is no antisense transcription.",
        ],
        "teaching_facts_checked_by_build": [{"fact": d, "passed": bool(ok)} for d, ok in facts],
    }
    with open(CHECK / "practice3_truth.json", "w") as fh:
        json.dump(truth, fh, indent=1)
        fh.write("\n")


def main():
    global CPGS, CPG_P
    D3.mkdir(parents=True, exist_ok=True)
    CHECK.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    ref = build_reference(rng)
    assert len(ref) == LOCUS_LEN
    assert ref[SNV1["pos"] - 1] == ord("G") and ref[SNV2["pos"] - 1] == ord("C")
    CPGS, CPG_P = cpg_table(ref)
    haps = {h: build_haplotype(ref, h) for h in (1, 2)}
    write_annotation(ref)
    libs = {
        "DNA-ILMN": sim_dna_illumina(np.random.default_rng(SEED + 1), ref, haps),
        "DNA-HIFI": sim_dna_long(np.random.default_rng(SEED + 2), ref, haps, "DNA-HIFI"),
        "DNA-ONT": sim_dna_long(np.random.default_rng(SEED + 3), ref, haps, "DNA-ONT"),
        "RNA-ILMN": sim_rna_illumina(np.random.default_rng(SEED + 4), ref),
        "RNA-KINNEX": sim_rna_long(np.random.default_rng(SEED + 5), ref, "RNA-KINNEX", 1_500),
        "RNA-ONT": sim_rna_long(np.random.default_rng(SEED + 6), ref, "RNA-ONT", 2_000),
    }
    rows = []
    for row in SAMPLES:
        lib = row[0]
        write_sam(libs[lib], D3 / f"{lib}.sam.gz", row)
        n_reads = len(libs[lib].records)
        rows.append((lib, PERSON, row[1], row[2], row[3], row[4], row[5], row[6], n_reads,
                     f"{row[7]}; synthetic teaching data, not calibrated to a real run"))
    pd.DataFrame(rows, columns=["library_id", "person", "molecule", "specimen_tube", "platform", "chemistry",
                                "read_layout", "library_prep", "reads_in_locus", "notes"]).to_csv(
        D3 / "samples.tsv", sep="\t", index=False)
    for lib in ("DNA-HIFI", "DNA-ONT"):
        n, m = libs[lib].meth
        c0 = CPGS[n[CPGS] > 0]
        tab = pd.DataFrame({"chrom": CHROM, "pos": c0 + 1, "n_reads": n[c0], "n_methylated": m[c0],
                            "fraction": np.round(m[c0] / n[c0], 3)})
        tab.to_csv(D3 / f"methylation_{lib}.tsv", sep="\t", index=False)
    checks, data = reference_checks()
    facts = teaching_facts(checks)
    failed = [d for d, ok in facts if not ok]
    assert not failed, "teaching facts not visible: " + "; ".join(failed)
    checks.to_csv(CHECK / "practice3_checks.tsv", sep="\t", index=False)
    reference_scorecard(checks).to_csv(CHECK / "practice3_scorecard.tsv", sep="\t", index=False)
    write_truth(ref, data, checks, facts)
    print(f"data/practice3: {sum(f.stat().st_size for f in D3.iterdir()) / 1e6:.1f} MB; "
          f"{len(checks)} CHECK numbers; {len(facts)} teaching facts hold")
    return ref, haps, libs, checks


if __name__ == "__main__":
    main()
