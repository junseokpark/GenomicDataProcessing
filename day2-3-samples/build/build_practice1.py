#!/usr/bin/env python3
"""Build data/practice1/ and checkpoints/practice1_*.tsv for Mini practice 1.

Mini practice 1 . Is this sample fit for this platform?

Real inputs (shipped in data/practice1/, verified by SHA-256 before anything is simulated):
  hbb_window.hg38.fa   UCSC hg38 chr11:5,218,771-5,278,770 (60,000 bp), the Day 1 window
  genes.hg38.bed       NCBI RefSeq Select beta-globin gene models, copied from Day 1

Everything else is synthetic and deterministic (fixed seeds). No network access.

Usage
  python build/build_practice1.py
  python build/build_practice1.py --day1-seq PATH --day1-genes PATH   # bootstrap from Day 1 files
"""
import argparse
import gzip
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D1 = ROOT / "data" / "practice1"
READS = D1 / "reads"
CHECK = ROOT / "checkpoints"

CHROM, CHROM_LEN = "chr11", 135_086_622            # hg38 chr11 length (used in the PAF)
WIN_START0, WIN_END = 5_218_770, 5_278_770          # BED [start, end) of the 60-kb window
WIN_NAME = "chr11:5218771-5278770"                  # 1-based closed, genome-browser style
TARGET = (CHROM, 5_226_000, 5_234_000, "HBB_HBD_window")
SEQ_SHA256 = "74b7c7effcc64b4d7083d89d4cb774e659547fdfda1e27f757af9e734464205f"
GENES_SHA256 = "4192b42e27eebf0ebffa125858c9f0bf53ffeccf2fa1e9d3cf98a3eb7ed8cac9"
DAY1_18MER = (5_226_995, 5_227_012, "CTTCTCCTCAGGAGTCAG")   # Day 1 slide-46 bases, 1-based

NA = np.nan

# --------------------------------------------------------------------------------------
# Part A tables
# --------------------------------------------------------------------------------------
SHEET_COLUMNS = ["sample_id", "subject_id", "specimen", "source_detail", "molecule",
                 "conc_ng_per_ul", "volume_ul", "a260_280", "a260_230", "din", "rin", "dv200",
                 "modal_fragment_kb", "cell_count", "requested_assay", "study_question",
                 "matched_normal_id", "consent_genomic_research"]

SAMPLES = [
    ("S01", "P001", "blood (EDTA)", "EDTA whole blood, high-molecular-weight extraction", "DNA",
     45, 50, 1.86, 2.10, 9.2, NA, NA, 48, NA, "PACBIO_HIFI_WGS",
     "Resolve structural variants and phase haplotypes", NA, "yes"),
    ("S02", "P002", "saliva (Oragene)", "Oragene saliva kit, column extraction", "DNA",
     12, 40, 1.80, 1.45, 7.6, NA, NA, 25, NA, "ILMN_PCRFREE_WGS",
     "Genome-wide SNVs and indels", NA, "yes"),
    ("S03", "P003", "FFPE tumour", "archival formalin-fixed paraffin-embedded tumour block", "DNA",
     8, 30, 1.90, 1.95, 2.6, NA, NA, 0.35, NA, "PACBIO_HIFI_WGS",
     "Somatic structural variants", NA, "yes"),
    ("S04", "P004", "fresh-frozen tumour", "tumour biopsy, snap-frozen", "DNA",
     60, 40, 1.88, 2.15, 8.4, NA, NA, 35, NA, "ILMN_WES",
     "Somatic coding mutations", "S05", "yes"),
    ("S05", "P004", "blood", "EDTA whole blood; matched normal for S04 (same subject)", "DNA",
     40, 30, 1.84, 2.08, 9.0, NA, NA, 40, NA, "ILMN_WES",
     "Germline background for S04", NA, "yes"),
    ("S06", "P006", "plasma cfDNA", "cell-free DNA extracted from plasma", "DNA",
     0.4, 50, 1.75, 1.60, NA, NA, NA, 0.167, NA, "ONT_LIGATION_LONG",
     "Phase variants across 50 kb", NA, "yes"),
    ("S07", "P007", "lymphoblastoid cell line", "cultured cells, frozen pellet; DNA not extracted",
     "cells", NA, NA, NA, NA, NA, NA, NA, NA, 6_500_000, "ONT_ULTRALONG",
     "Assemble a complex region (e.g. MHC)", NA, "yes"),
    ("S08", "P008", "blood (EDTA)", "EDTA whole blood", "DNA",
     30, 20, 1.83, 2.05, 8.9, NA, NA, 22, NA, "PACBIO_HIFI_WGS",
     "Size a tandem-repeat expansion", NA, "yes"),
    ("S09", "P009", "fresh-frozen brain cortex", "post-mortem cortex, snap-frozen", "RNA",
     120, 20, 2.05, 2.10, NA, 8.7, NA, NA, NA, "RNA_POLYA",
     "Gene expression in cortex", NA, "yes"),
    ("S10", "P010", "FFPE tumour", "archival FFPE tumour block", "RNA",
     15, 15, 1.95, 1.90, NA, 2.3, 62, NA, NA, "RNA_POLYA",
     "Tumour gene expression", NA, "yes"),
    ("S11", "P011", "buccal swab", "cheek swab", "DNA",
     6, 50, 1.82, 2.02, 7.0, NA, NA, 18, NA, "PACBIO_HIFI_WGS",
     "Structural variants", NA, "yes"),
    ("S12", "P012", "blood", "EDTA whole blood", "DNA",
     50, 40, 1.85, 2.12, 9.1, NA, NA, 45, NA, "ONT_LIGATION_LONG",
     "Methylation and structural variants", NA, "no"),
]

REQ_COLUMNS = ["assay", "platform", "kit_or_chemistry", "molecule", "min_input", "input_unit",
               "min_fragment_kb", "integrity_rule", "purity_rule", "notes", "source"]
DNA_PURITY = "A260/280 1.8-2.0; A260/230 2.0-2.2"
RNA_PURITY = "A260/280 1.9-2.1; A260/230 >= 1.8 (teaching range)"

REQUIREMENTS = [
    ("ILMN_PCRFREE_WGS", "Illumina short-read WGS (2x150)", "Illumina DNA PCR-Free Prep", "DNA",
     25, "ng", 1, "modal_fragment_kb >= 1", DNA_PURITY,
     "Product page: input 25-300 ng; reference guide: protocol accepts 25-2000 ng "
     "(>= 300 ng gives normalized libraries). Not an FFPE workflow. "
     "min_fragment_kb is a teaching threshold.",
     "Illumina DNA PCR-Free Prep product page and reference guide (1000000086922); "
     "min_fragment_kb: teaching threshold"),
    ("ILMN_WES", "Illumina short-read exome (2x150)", "Illumina DNA Prep with Enrichment (exome panel)",
     "DNA", 100, "ng", 0.15, "modal_fragment_kb >= 0.15", DNA_PURITY,
     "Recommended minimum input for human DNA 100-500 ng. Accepts extracted FFPE DNA with the "
     "FFPE workflow, hence the short-fragment teaching threshold.",
     "Illumina DNA Prep reference guide (1000000025416); min_fragment_kb: teaching threshold"),
    ("PACBIO_HIFI_WGS", "PacBio HiFi long-read WGS", "PacBio Revio, SPRQ-Nx chemistry", "DNA",
     500, "ng", 15, "modal_fragment_kb >= 15", DNA_PURITY,
     "500 ng DNA per sample. Human WGS libraries are 15-20 kb, so the DNA must be at least as "
     "long as the library (teaching rule).",
     "PacBio Revio product page / specification sheet; length rule: teaching threshold"),
    ("ONT_LIGATION_LONG", "Oxford Nanopore long-read", "Ligation Sequencing Kit V14 (SQK-LSK114), "
     "R10.4.1 flow cell", "DNA", 1000, "ng", 10, "modal_fragment_kb >= 10", DNA_PURITY,
     "1 ug (1,000 ng) gDNA for long-fragment (> 10 kb) libraries. Read length = fragment length.",
     "ONT store page / chemistry technical document (SQK-LSK114)"),
    ("ONT_ULTRALONG", "Oxford Nanopore ultra-long", "Ultra-Long DNA Sequencing Kit V14 (SQK-ULK114), "
     "R10.4.1 flow cell", "cells", 6_000_000, "cells", NA, NA, NA,
     "6 million cells (PBMCs, cultured cells or cells from frozen tissue). The protocol extracts "
     "the DNA itself because ordinary extraction breaks long DNA. Read N50 50-100+ kb, "
     "reads up to 4+ Mb.",
     "ONT store page (SQK-ULK114)"),
    ("RNA_POLYA", "Illumina short-read RNA-seq", "poly(A) selection mRNA-seq", "RNA",
     100, "ng", NA, "rin >= 7", RNA_PURITY,
     "Poly(A) selection pulls transcripts out by their 3' tail; on fragmented RNA only pieces "
     "near the tail are kept (3' bias). RIN >= 7 and 100 ng: teaching thresholds "
     "(common core-facility practice).",
     "teaching threshold"),
    ("RNA_CAPTURE_FFPE", "Illumina short-read RNA-seq", "capture-based RNA-seq for degraded/FFPE RNA "
     "(e.g. Illumina RNA exome, TruSeq RNA Access)", "RNA", 20, "ng", NA, "dv200 >= 30", RNA_PURITY,
     "DV200 >= 30% (< 30% not recommended). Input by DV200: > 70% 20 ng; 50-70% 20-40 ng; "
     "30-50% 40-100 ng. min_input is the lowest band; the S10 re-check uses the band.",
     "Illumina technical note 'Evaluating RNA Quality from FFPE Samples' (470-2014-001)"),
]

ALTERNATIVE = {
    "S02": "Bead clean-up, re-measure A260/230, then run ILMN_PCRFREE_WGS as requested",
    "S03": "ILMN_WES or a panel with an FFPE-compatible prep, plus a matched normal "
           "(e.g. blood) from P003",
    "S06": "Mutation detection: Illumina cfDNA targeted panel. Phasing: collect blood and "
           "extract long DNA",
    "S10": "RNA_CAPTURE_FFPE: DV200 62% (medium quality) needs 20-40 ng; see the re-check",
    "S11": "Re-quantify (fluorometer) and concentrate; or recollect as blood; or switch to "
           "ILMN_PCRFREE_WGS (resolves fewer SVs)",
    "S12": "Do not process; refer back to the study team (consent)",
}

EXPECTED = {"S01": "PASS", "S02": "REVIEW", "S03": "FAIL", "S04": "PASS", "S05": "PASS",
            "S06": "FAIL", "S07": "PASS", "S08": "PASS", "S09": "PASS", "S10": "FAIL",
            "S11": "REVIEW", "S12": "BLOCKED"}

# purity windows (low, high); high None = no upper limit. Same numbers as the purity_rule column.
PURITY = {"DNA": [("a260_280", 1.8, 2.0), ("a260_230", 2.0, 2.2)],
          "RNA": [("a260_280", 1.9, 2.1), ("a260_230", 1.8, None)]}


def reference_decision(s, r, by_id):
    """Independent re-implementation of the notebook's decide(); returns (decision, reasons)."""
    if s["consent_genomic_research"] != "yes":
        return "BLOCKED", ["BLOCKED: no consent for genomic research "
                           "(do not process, do not assess quality)"]
    if s["molecule"] != r["molecule"]:
        return "FAIL", [f"FAIL: assay needs {r['molecule']} but the tube holds {s['molecule']}"]
    found = []
    unit = r["input_unit"]
    have = s["cell_count"] if unit == "cells" else s["conc_ng_per_ul"] * s["volume_ul"]
    need = r["min_input"]
    if have < need:
        lvl = "REVIEW" if have >= 0.5 * need else "FAIL"
        found.append(f"{lvl}: input {have:,.0f} {unit} is {have / need:.0%} of the "
                     f"{need:,.0f} {unit} minimum")
    rule = r["integrity_rule"]
    if isinstance(rule, str):
        metric, thr = rule.split(" >= ")
        thr = float(thr)
        v = s[metric]
        if pd.isna(v):
            found.append(f"REVIEW: {metric} not measured")
        elif v < thr:
            found.append(f"FAIL: {metric} {v:g} < {thr:g} required (cannot be repaired)")
    for ratio, low, high in PURITY.get(s["molecule"], []):
        v = s[ratio]
        if pd.isna(v):
            found.append(f"REVIEW: {ratio} not measured")
        elif v < low:
            found.append(f"REVIEW: {ratio} {v:.2f} below {low:.2f} (clean up, re-measure)")
        elif high is not None and v > high:
            found.append(f"REVIEW: {ratio} {v:.2f} above {high:.2f} (clean up, re-measure)")
    if "tumour" in s["specimen"] and "somatic" in s["study_question"].lower():
        mn = s["matched_normal_id"]
        if pd.isna(mn):
            found.append("REVIEW: tumour sample for a somatic question but no matched normal")
        elif mn not in by_id or by_id[mn]["subject_id"] != s["subject_id"]:
            found.append(f"REVIEW: matched normal {mn} is not from the same subject")
    levels = [f.split(":")[0] for f in found]
    decision = "FAIL" if "FAIL" in levels else ("REVIEW" if "REVIEW" in levels else "PASS")
    return decision, found


# --------------------------------------------------------------------------------------
# Part B simulation
# --------------------------------------------------------------------------------------
LUT = np.full(256, -1, dtype=np.int64)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
RC_TABLE = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(RC_TABLE)[::-1]


def run_lengths(codes):
    change = np.flatnonzero(np.diff(codes) != 0) + 1
    starts, ends = np.r_[0, change], np.r_[change, len(codes)]
    return np.repeat(ends - starts, ends - starts)


def q_to_str(q):
    return (np.asarray(q, dtype=np.uint8) + 33).tobytes().decode()


def simulate_long_read(rng, ref, hp, start, length, q_fn, hp_probs, non_probs):
    """One long read from ref[start:start+length]. Quality is drawn first, then each reference
    base is wrong with probability 10^(-Q/10); error type depends on homopolymer context.
    Returns (bases_on_forward_strand, quals, n_match, aln_len)."""
    seg = ref[start:start + length]
    in_hp = hp[start:start + length] >= 4
    n = len(seg)
    q = q_fn(rng, n, in_hp)
    err = rng.random(n) < 10.0 ** (-q / 10.0)
    u = rng.random(n)
    t_hp = 1 + np.searchsorted(np.cumsum(hp_probs), u, side="right")
    t_non = 1 + np.searchsorted(np.cumsum(non_probs), u, side="right")
    kind = np.where(err, np.where(in_hp, t_hp, t_non), 0)        # 0 ok, 1 sub, 2 ins, 3 del
    kind = np.minimum(kind, 3)
    counts = np.where(kind == 3, 0, np.where(kind == 2, 2, 1))
    out = np.repeat(seg, counts)
    outq = np.repeat(q, counts)
    first = np.cumsum(counts) - counts
    sub = np.flatnonzero(kind == 1)
    out[first[sub]] = BASES[(LUT[seg[sub]] + rng.integers(1, 4, len(sub))) % 4]
    ins = np.flatnonzero(kind == 2)
    ins_non = ins[~in_hp[ins]]
    out[first[ins_non] + 1] = BASES[rng.integers(0, 4, len(ins_non))]  # HP insertions copy the base
    outq[first[ins] + 1] = rng.integers(4, 15, len(ins))
    n_ins = len(ins)
    n_match = int((kind == 0).sum() + n_ins)
    aln_len = int(n + n_ins)
    return out.tobytes().decode(), outq, n_match, aln_len


def hifi_q(rng, n, in_hp):
    q = np.clip(np.rint(rng.normal(45, 6, n)), 30, 60)
    q_hp = np.clip(np.rint(rng.normal(24, 4, n)), 12, 35)
    return np.where(in_hp, q_hp, q).astype(np.int64)


def make_ont_q(mu):
    def ont_q(rng, n, in_hp):
        q = np.clip(np.rint(rng.normal(mu, 5, n)), 3, 40)
        q_hp = np.clip(np.rint(rng.normal(mu - 7, 4, n)), 3, 40)
        return np.where(in_hp, q_hp, q).astype(np.int64)
    return ont_q


def write_fastq_gz(path, records):
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            for name, bases, quals in records:
                gz.write(f"@{name}\n{bases}\n+\n{q_to_str(quals)}\n".encode())


def paf_line(qname, qlen, strand, w_start, w_end, n_match, aln_len):
    # minimal 12-column PAF; target coordinates are real chr11 (0-based, half-open like BED);
    # MAPQ 255 = "not available": these alignments come from the simulation, not an aligner
    return (f"{qname}\t{qlen}\t0\t{qlen}\t{strand}\t{CHROM}\t{CHROM_LEN}\t"
            f"{WIN_START0 + w_start}\t{WIN_START0 + w_end}\t{n_match}\t{aln_len}\t255")


def simulate_illumina(rng, ref, n_pairs=2000, read_len=150):
    win = len(ref)
    insert = np.clip(np.rint(rng.normal(400, 25, n_pairs)), 350, 450).astype(int)
    start = (rng.random(n_pairs) * (win - insert + 1)).astype(int)
    frag_plus = rng.random(n_pairs) < 0.5
    cyc = np.arange(read_len) / (read_len - 1)
    r1, r2, paf = [], [], []
    for mate, mean_q in ((1, 38 - 8 * cyc ** 2), (2, 36 - 9 * cyc ** 2)):
        q = np.clip(np.rint(mean_q + rng.normal(0, 2.5, (n_pairs, read_len))), 2, 41)
        dip = rng.random((n_pairs, read_len)) < 0.004
        q = np.where(dip, rng.integers(8, 20, (n_pairs, read_len)), q).astype(np.int64)
        err = rng.random((n_pairs, read_len)) < 10.0 ** (-q / 10.0)
        shift = rng.integers(1, 4, (n_pairs, read_len))
        for i in range(n_pairs):
            s, e = start[i], start[i] + insert[i]
            # mate 1 reads the fragment's own strand from its start; mate 2 reads the other strand
            forward = frag_plus[i] if mate == 1 else not frag_plus[i]
            if forward:
                w0, w1, strand = s, s + read_len, "+"
                seq = ref[w0:w1].copy()
            else:
                w0, w1, strand = e - read_len, e, "-"
                seq = np.frombuffer(revcomp(ref[w0:w1].tobytes().decode()).encode(),
                                    dtype=np.uint8).copy()
            bad = np.flatnonzero(err[i])
            seq[bad] = BASES[(LUT[seq[bad]] + shift[i, bad]) % 4]
            name = f"S08_ilmn_{i + 1:05d}/{mate}"
            (r1 if mate == 1 else r2).append((name, seq.tobytes().decode(), q[i]))
            paf.append((i, mate, paf_line(name, read_len, strand, w0, w1,
                                          read_len - len(bad), read_len)))
    paf.sort(key=lambda x: (x[0], x[1]))
    return r1, r2, [p[2] for p in paf]


def simulate_long(rng, ref, hp, lengths, prefix, q_fns, hp_probs, non_probs):
    win = len(ref)
    recs, paf = [], []
    for i, L in enumerate(lengths):
        start = int(rng.integers(0, win - L + 1))
        seq, q, n_match, aln_len = simulate_long_read(rng, ref, hp, start, int(L), q_fns[i],
                                                      hp_probs, non_probs)
        strand = "+" if rng.random() < 0.5 else "-"
        if strand == "-":
            seq, q = revcomp(seq), q[::-1]
        name = f"S08_{prefix}_{i + 1:04d}"
        recs.append((name, seq, q))
        paf.append(paf_line(name, len(seq), strand, start, start + int(L), n_match, aln_len))
    return recs, paf


def hifi_lengths(rng, n=40):
    return np.clip(np.rint(rng.normal(16_000, 3_000, n)), 8_000, 25_000).astype(int)


def ont_lengths(rng, n=50, n50=20_000, sigma=0.8, max_len=58_000, min_len=1_000):
    mu = np.log(n50) - sigma ** 2            # length-weighted median of a log-normal = N50
    out = []
    while len(out) < n:
        x = int(round(rng.lognormal(mu, sigma)))
        if min_len <= x <= max_len:
            out.append(x)
    return np.array(out)


# --------------------------------------------------------------------------------------
# reference statistics, computed by re-reading the files that were written
# --------------------------------------------------------------------------------------
def n50(lengths):
    s = np.sort(np.asarray(lengths))[::-1]
    return int(s[np.searchsorted(np.cumsum(s), s.sum() / 2)])


def fastq_stats(paths):
    lengths, err_sum = [], 0.0
    for p in paths:
        with gzip.open(p, "rt") as fh:
            lines = fh.read().splitlines()
        assert len(lines) % 4 == 0
        for j in range(0, len(lines), 4):
            assert lines[j].startswith("@") and lines[j + 2] == "+"
            q = np.frombuffer(lines[j + 3].encode(), dtype=np.uint8).astype(np.int64) - 33
            assert len(q) == len(lines[j + 1])
            lengths.append(len(q))
            err_sum += float((10.0 ** (-q / 10.0)).sum())
    return np.array(lengths), err_sum


def read_stats_table(win_len):
    paf = pd.read_csv(READS / "alignments.paf", sep="\t", header=None,
                      names=["qname", "qlen", "qstart", "qend", "strand", "tname", "tlen",
                             "tstart", "tend", "nmatch", "alnlen", "mapq"])
    plat = paf.qname.str.extract(r"^S08_(ilmn|hifi|ont)_")[0]
    spans = ((paf.tstart <= TARGET[1]) & (paf.tend >= TARGET[2])).groupby(plat).sum()
    files = {"illumina_2x150": (["illumina_R1.fastq.gz", "illumina_R2.fastq.gz"], "ilmn"),
             "pacbio_hifi": (["pacbio_hifi.fastq.gz"], "hifi"),
             "ont_ligation": (["ont_ligation.fastq.gz"], "ont")}
    rows = []
    for platform, (names, key) in files.items():
        L, e = fastq_stats([READS / n for n in names])
        rows.append(dict(platform=platform, reads=len(L), total_bases=int(L.sum()),
                         mean_length=round(L.mean(), 1), max_length=int(L.max()),
                         read_n50=n50(L), mean_error_rate=round(e / L.sum(), 6),
                         mean_accuracy_pct=round(100 * (1 - e / L.sum()), 3),
                         approx_depth_x=round(L.sum() / win_len, 2),
                         reads_spanning_target=int(spans.get(key, 0))))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------
def load_or_bootstrap(args):
    fa, fai, genes = D1 / "hbb_window.hg38.fa", D1 / "hbb_window.hg38.fa.fai", D1 / "genes.hg38.bed"
    if args.day1_seq:
        seq = Path(args.day1_seq).read_text().strip().upper()
        lines = [seq[i:i + 60] for i in range(0, len(seq), 60)]
        fa.write_text(f">{WIN_NAME}\n" + "\n".join(lines) + "\n")
    if args.day1_genes:
        genes.write_bytes(Path(args.day1_genes).read_bytes())
    if not fa.exists() or not genes.exists():
        sys.exit("missing data/practice1/hbb_window.hg38.fa or genes.hg38.bed; "
                 "bootstrap with --day1-seq/--day1-genes (the Day 1 package files)")
    text = fa.read_text().splitlines()
    assert text[0] == f">{WIN_NAME}", text[0]
    seq = "".join(text[1:])
    assert hashlib.sha256(seq.encode()).hexdigest() == SEQ_SHA256, "window sequence changed"
    assert hashlib.sha256(genes.read_bytes()).hexdigest() == GENES_SHA256, "gene models changed"
    assert len(seq) == WIN_END - WIN_START0 == 60_000
    s1, e1, kmer = DAY1_18MER
    assert seq[s1 - 1 - WIN_START0:e1 - WIN_START0] == kmer       # same bases as Day 1
    offset = len(text[0]) + 1
    fai.write_text(f"{WIN_NAME}\t{len(seq)}\t{offset}\t60\t61\n")
    return seq


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--day1-seq")
    ap.add_argument("--day1-genes")
    args = ap.parse_args()
    READS.mkdir(parents=True, exist_ok=True)
    CHECK.mkdir(exist_ok=True)
    seq = load_or_bootstrap(args)

    # ---- Part A tables ------------------------------------------------------------------
    sheet = pd.DataFrame(SAMPLES, columns=SHEET_COLUMNS)
    sheet["cell_count"] = sheet["cell_count"].astype("Int64")
    sheet.to_csv(D1 / "sample_sheet.tsv", sep="\t", index=False, na_rep="NA")
    req = pd.DataFrame(REQUIREMENTS, columns=REQ_COLUMNS)
    req.to_csv(D1 / "platform_requirements.tsv", sep="\t", index=False, na_rep="NA")

    # re-read exactly as a student would, then decide with the reference implementation
    sheet = pd.read_csv(D1 / "sample_sheet.tsv", sep="\t")
    req = pd.read_csv(D1 / "platform_requirements.tsv", sep="\t").set_index("assay")
    by_id = {r["sample_id"]: r for _, r in sheet.iterrows()}
    rows = []
    for _, s in sheet.iterrows():
        r = req.loc[s["requested_assay"]]
        for ratio, low, high in PURITY.get(r["molecule"], []):    # table and code agree
            assert f"{low:.1f}" in r["purity_rule"]
        decision, reasons = reference_decision(s, r, by_id)
        total = s["conc_ng_per_ul"] * s["volume_ul"]
        rows.append(dict(sample_id=s["sample_id"], requested_assay=s["requested_assay"],
                         total_ng=round(total, 1) if pd.notna(total) else np.nan,
                         decision=decision, reasons="; ".join(reasons) if reasons else "-",
                         suggested_alternative=ALTERNATIVE.get(s["sample_id"], "-")))
    dec = pd.DataFrame(rows)
    got = dict(zip(dec.sample_id, dec.decision))
    assert got == EXPECTED, {k: (v, EXPECTED[k]) for k, v in got.items() if v != EXPECTED[k]}
    assert dec.decision.value_counts().to_dict() == {"PASS": 6, "FAIL": 3, "REVIEW": 2, "BLOCKED": 1}
    assert set(ALTERNATIVE) == set(dec.sample_id[dec.decision != "PASS"])
    s10 = by_id["S10"]
    assert reference_decision(s10, req.loc["RNA_CAPTURE_FFPE"], by_id)[0] == "PASS"
    assert reference_decision(by_id["S03"], req.loc["ILMN_WES"], by_id)[0] == "REVIEW"
    dec.to_csv(CHECK / "practice1_decisions.tsv", sep="\t", index=False, na_rep="NA")

    # ---- Part B reads -------------------------------------------------------------------
    target = "\t".join(map(str, TARGET)) + "\n"
    (D1 / "target_window.bed").write_text(target)
    ref = np.frombuffer(seq.encode(), dtype=np.uint8)
    hp = run_lengths(ref)

    r1, r2, paf_il = simulate_illumina(np.random.default_rng(1013), ref)
    write_fastq_gz(READS / "illumina_R1.fastq.gz", r1)
    write_fastq_gz(READS / "illumina_R2.fastq.gz", r2)

    rng = np.random.default_rng(2026)
    hifi, paf_hifi = simulate_long(rng, ref, hp, hifi_lengths(rng), "hifi", [hifi_q] * 40,
                                   hp_probs=(0.05, 0.45, 0.50), non_probs=(0.60, 0.20, 0.20))
    write_fastq_gz(READS / "pacbio_hifi.fastq.gz", hifi)

    rng = np.random.default_rng(ONT_SEED)
    lengths = ont_lengths(rng)
    q_fns = [make_ont_q(mu) for mu in rng.uniform(22, 27, len(lengths))]
    ont, paf_ont = simulate_long(rng, ref, hp, lengths, "ont", q_fns,
                                 hp_probs=(0.10, 0.30, 0.60), non_probs=(0.45, 0.20, 0.35))
    write_fastq_gz(READS / "ont_ligation.fastq.gz", ont)

    (READS / "alignments.paf").write_text("\n".join(paf_il + paf_hifi + paf_ont) + "\n")

    stats = read_stats_table(len(seq))
    st = stats.set_index("platform")
    assert st.loc["illumina_2x150", "reads_spanning_target"] == 0
    assert st.loc["pacbio_hifi", "reads_spanning_target"] >= 3, st
    assert st.loc["ont_ligation", "reads_spanning_target"] >= 3, st
    assert st.loc["ont_ligation", "max_length"] <= 58_000
    assert 17_000 <= st.loc["ont_ligation", "read_n50"] <= 23_000, st
    assert st.loc["pacbio_hifi", "mean_accuracy_pct"] >= 99.9, st
    assert 98.0 <= st.loc["ont_ligation", "mean_accuracy_pct"] <= 99.5, st
    stats.to_csv(CHECK / "practice1_read_stats.tsv", sep="\t", index=False)
    print(dec[["sample_id", "requested_assay", "decision"]].to_string(index=False))
    print(stats.to_string(index=False))


ONT_SEED = 326

if __name__ == "__main__":
    main()
