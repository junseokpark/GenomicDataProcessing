"""Build the Day-1 practice-session data package.

Real sequence: UCSC chr11 windows (hg19 5,240,001-5,300,000 and the identical
hg38 5,218,771-5,278,770).  Everything else -- samples, reads, genotypes -- is
synthetic and labelled as such.
"""
import gzip, os, random, shutil, textwrap
import pysam

random.seed(20261012)
OUT = "practice_session"
D = lambda *p: os.path.join(OUT, *p)
OFFSET = 21230                      # hg19 pos - OFFSET = hg38 pos, verified over the window
H19_START, H19_END = 5240001, 5300000
H38_START, H38_END = H19_START - OFFSET, H19_END - OFFSET
H19_LEN, H38_LEN = 5300000, 5280000   # mini-chromosome lengths (N-padded)

for sub in ("", "collaborator_share", "target", "checkpoints", "igv", "templates"):
    os.makedirs(D(sub), exist_ok=True)

h19_win = open("lab_src/hg19_5240001_5300000.txt").read().strip().upper()
h38_win = open("lab_src/hg38_5218771_5278770.txt").read().strip().upper()
assert h19_win == h38_win and len(h19_win) == 60000

RC = str.maketrans("ACGTN", "TGCAN")
rc = lambda s: s[::-1].translate(RC)


# ---------------------------------------------------------------- references
def write_fasta(path, name, seq, width=60):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")
    # .fai: name, length, offset of first base, bases per line, bytes per line
    with open(path + ".fai", "w") as fh:
        fh.write(f"{name}\t{len(seq)}\t{len(name) + 2}\t{width}\t{width + 1}\n")


def padded(seq, start, total):
    """Place `seq` at 1-based `start` inside an N-padded contig of length `total`."""
    assert start - 1 + len(seq) <= total
    return "N" * (start - 1) + seq + "N" * (total - (start - 1) - len(seq))


hg19_mini = padded(h19_win, H19_START, H19_LEN)
hg38_mini = padded(h38_win, H38_START, H38_LEN)
write_fasta(D("collaborator_share", "reference.hg19.fa"), "chr11", hg19_mini)
write_fasta(D("target", "reference.hg38.fa"), "chr11", hg38_mini)

base19 = lambda p: hg19_mini[p - 1]
ref19 = lambda s, e: hg19_mini[s - 1:e]           # 1-based inclusive

# toy1: the format-teaching contig used by slides 49-54 (ATC at 100-102, A at 300)
toy = list("".join(random.choice("ACGT") for _ in range(1000)))
toy[99:102] = list("ATC")
toy[299] = "A"
toy1 = "".join(toy)
write_fasta(D("collaborator_share", "toy1.fa"), "toy1", toy1)
assert toy1[99:102] == "ATC" and toy1[299] == "A"

# ------------------------------------------------------------ hg38 features
# Real NCBI RefSeq Select transcripts for the beta-globin cluster (UCSC, hg38).
GENES = [  # name, txStart(1-based), txEnd, strand, exons [(start1,end1)]
    ("HBB",  5225464, 5227071, "-", [(5225464, 5225726), (5225891, 5226114), (5226929, 5227071)]),
    ("HBD",  5232838, 5234483, "-", [(5232838, 5233100), (5233265, 5233487), (5234341, 5234483)]),
    ("HBG1", 5248269, 5249857, "-", [(5248269, 5248531), (5248696, 5248918), (5249715, 5249857)]),
    ("HBG2", 5253188, 5254781, "-", [(5253188, 5253450), (5253615, 5253837), (5254639, 5254781)]),
    ("HBE1", 5268345, 5269945, "-", [(5268345, 5268607), (5268772, 5268994), (5269803, 5269945)]),
]
with open(D("target", "genes.hg38.bed"), "w") as fh:
    fh.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
    for nm, s, e, st, exons in GENES:
        fh.write(f"chr11\t{s-1}\t{e}\t{nm}\t0\t{st}\n")
        ordered = exons[::-1] if st == "-" else exons
        for i, (es, ee) in enumerate(ordered, 1):
            fh.write(f"chr11\t{es-1}\t{ee}\t{nm}_exon{i}\t0\t{st}\n")
        # 500 bp promoter window upstream of the TSS (minus strand -> after txEnd)
        ps, pe = (e + 1, e + 500) if st == "-" else (s - 500, s - 1)
        fh.write(f"chr11\t{ps-1}\t{pe}\t{nm}_promoter\t0\t{st}\n")

# toy1 gene file for the boundary exercise of slide 85
with open(D("collaborator_share", "genes.toy1.bed"), "w") as fh:
    fh.write("#chrom\tchromStart\tchromEnd\tname\tscore\tstrand\n")
    fh.write("toy1\t200\t250\ttoyGeneA\t0\t+\n")
    fh.write("toy1\t400\t600\ttoyGeneB\t0\t-\n")

# ----------------------------------------------------------------- the chain
# Single ungapped block; verified byte-identical over the declared span.
CHAIN_H19_S, CHAIN_H19_E = 5246500, 5280000        # 0-based half-open, hg19
with open(D("target", "hg19ToHg38.teaching.over.chain"), "w") as fh:
    fh.write("# TEACHING CHAIN -- NOT the UCSC hg19ToHg38.over.chain.gz file.\n")
    fh.write("# One ungapped block over chr11:%d-%d (hg19). Within this span the two\n"
             % (CHAIN_H19_S, CHAIN_H19_E))
    fh.write("# assemblies are byte-identical, so hg38 = hg19 - %d exactly.\n" % OFFSET)
    fh.write("# Anything outside the block is intentionally unmappable.\n")
    size = CHAIN_H19_E - CHAIN_H19_S
    fh.write(f"chain {size*2} chr11 {H19_LEN} + {CHAIN_H19_S} {CHAIN_H19_E} "
             f"chr11 {H38_LEN} + {CHAIN_H19_S-OFFSET} {CHAIN_H19_E-OFFSET} 1\n")
    fh.write(f"{size}\n\n")

# ------------------------------------------------------------- legacy tables
CHECKPOINT = 5248232                                # hg19 1-based, HBB_checkpoint
assert base19(CHECKPOINT) == "T"
DEL_POS, DEL_LEN = 5248100, 3                       # padded deletion anchor
INS_POS = 5248400                                   # padded insertion anchor

regions = {  # sample -> rows of (region_id, chrom, start1, end1, note)
    "S01": [
        ("r01", "chr11", CHECKPOINT, CHECKPOINT, "HBB_checkpoint"),
        ("r02", "chr11", 5248160, 5248251, "HBB_exon1_cds"),
        ("r03", "chr11", 5248252, 5248252, "single_base_next_to_exon"),
        ("r04", "chr11", 5246010, 5246100, "outside_teaching_chain"),
        ("r05", "chr11", 5248500, 5248401, "invalid_start_gt_end"),
    ],
    "S02": [
        ("r06", "chr11", 5247000, 5247099, "intronic_block"),
        ("r07", "chr11", 5248232, 5248232, "HBB_checkpoint"),
    ],
    "S03": [
        ("r08", "chr11", 5248200, 5248299, "declared_hg38_actually_hg19"),
    ],
}
for smp, rows in regions.items():
    with open(D("collaborator_share", f"regions.{smp}.hg19.tsv"), "w") as fh:
        fh.write("region_id\tchrom\tstart_1based\tend_1based\tnote\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
shutil.copy(D("collaborator_share", "regions.S01.hg19.tsv"),
            D("collaborator_share", "regions.hg19.tsv"))

# ----------------------------------------------------------- reads and BAMs
READLEN, FRAGLEN, NPAIRS = 80, 300, 90
REG_S = 5247900                                     # first fragment start (hg19)
STEP = 6


def apply_variants(seq, seq_start, carry):
    """Apply the designed variants to a reference slice; return (sequence, full CIGAR ops)."""
    end = seq_start + len(seq) - 1
    events = []
    if carry.get("snv") and seq_start <= CHECKPOINT <= end:
        events.append(("snv", CHECKPOINT))
    if carry.get("del") and seq_start <= DEL_POS and DEL_POS + DEL_LEN <= end:
        events.append(("del", DEL_POS))
    if carry.get("ins") and seq_start <= INS_POS < end:
        events.append(("ins", INS_POS))
    events.sort(key=lambda e: e[1])
    out, ops, cur = [], [], 0
    for kind, pos in events:
        off = pos - seq_start
        if kind == "snv":
            out.append(seq[cur:off] + "A")
            cur = off + 1
        elif kind == "del":
            out.append(seq[cur:off + 1])
            ops.append((off + 1 - cur, "M")); ops.append((DEL_LEN, "D"))
            cur = off + 1 + DEL_LEN
        elif kind == "ins":
            out.append(seq[cur:off + 1] + "TA")
            ops.append((off + 1 - cur, "M")); ops.append((2, "I"))
            cur = off + 1
    out.append(seq[cur:])
    full = "".join(out)
    used = sum(n for n, o in ops if o in "MIS")
    ops = ops + [(len(full) - used, "M")]
    return full, ops


def trim(seq, ops, n):
    """Cut a read (and its CIGAR) to exactly n read bases."""
    kept, used = [], 0
    for cnt, op in ops:
        if op in "DN":
            if used < n:
                kept.append((cnt, op))
            continue
        take = min(cnt, n - used)
        if take:
            kept.append((take, op)); used += take
        if used == n:
            break
    while kept and kept[-1][1] in "DN":
        kept.pop()
    return seq[:n], [(c, o) for c, o in kept if c > 0]


reads = []   # (qname, r1seq, r1pos, r1cigar, r2seq, r2pos, r2cigar)
PAD = 8
for k in range(NPAIRS):
    fstart = REG_S + k * STEP
    carry = {"snv": k % 2 == 0, "del": k % 2 == 0, "ins": k % 2 == 0}
    r2_start = fstart + FRAGLEN - READLEN
    for (name_pos, pos) in (("r1", fstart), ("r2", r2_start)):
        pass
    s1, o1 = trim(*apply_variants(ref19(fstart, fstart + READLEN + PAD - 1),
                                  fstart, carry), READLEN)
    s2, o2 = trim(*apply_variants(ref19(r2_start, r2_start + READLEN + PAD - 1),
                                  r2_start, carry), READLEN)
    reads.append((f"S01_read{k+1:03d}", s1, fstart, o1, s2, r2_start, o2))

# one deliberately soft-clipped read for the CIGAR discussion
sc_start = 5248450
sc_seq = "GGCTTAGTCA" + ref19(sc_start, sc_start + 69)
reads.append(("S01_read_softclip", sc_seq, sc_start, [(10, "S"), (70, "M")],
              ref19(sc_start + 220, sc_start + 299), sc_start + 220, [(80, "M")]))

CIG = {"M": 0, "I": 1, "D": 2, "S": 4}


def write_bam(path, contig, clen, shift, sample):
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"},
           "SQ": [{"SN": contig, "LN": clen}],
           "RG": [{"ID": sample, "SM": sample, "PL": "SYNTHETIC",
                   "DS": "simulated teaching reads"}]}
    recs = []
    for qname, s1, p1, c1, s2, p2, c2 in reads:
        for mate, (sq, pos, cig, other) in enumerate(
                [(s1, p1, c1, p2), (s2, p2, c2, p1)]):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = sq          # stored on the forward reference strand
            a.flag = (99 if mate == 0 else 147)
            a.reference_id = 0
            a.reference_start = pos - 1 - shift        # 0-based
            a.mapping_quality = 60
            a.cigar = [(CIG[o], n) for n, o in cig]
            a.next_reference_id = 0
            a.next_reference_start = other - 1 - shift
            a.template_length = (FRAGLEN if mate == 0 else -FRAGLEN)
            a.query_qualities = pysam.qualitystring_to_array("I" * len(sq))
            a.set_tag("RG", sample)
            recs.append(a)
    recs.sort(key=lambda a: a.reference_start)
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for a in recs:
            out.write(a)
    pysam.index(path)
    with pysam.AlignmentFile(path) as bam, open(path[:-4] + ".sam", "w") as sam:
        sam.write(str(bam.header))
        for a in bam:
            sam.write(a.to_string() + "\n")


write_bam(D("collaborator_share", "S01.hg19.sorted.bam"), "chr11", H19_LEN, 0, "S01")
write_bam(D("checkpoints", "S01.hg38.sorted.bam"), "chr11", H38_LEN, OFFSET, "S01")
shutil.copy(D("collaborator_share", "S01.hg19.sorted.bam"),
            D("collaborator_share", "S01.sorted.bam"))
shutil.copy(D("collaborator_share", "S01.hg19.sorted.bam.bai"),
            D("collaborator_share", "S01.sorted.bam.bai"))
shutil.copy(D("collaborator_share", "S01.hg19.sorted.sam"),
            D("collaborator_share", "S01.sorted.sam"))   # every BAM keeps its .sam twin

# FASTQ -- the original observations, identical for both alignment runs
with gzip.open(D("collaborator_share", "S01_R1.fastq.gz"), "wt") as f1, \
     gzip.open(D("collaborator_share", "S01_R2.fastq.gz"), "wt") as f2:
    for qname, s1, p1, c1, s2, p2, c2 in reads:
        f1.write(f"@{qname}/1\n{s1}\n+\n{'I'*len(s1)}\n")
        f2.write(f"@{qname}/2\n{rc(s2)}\n+\n{'I'*len(s2)}\n")

# observed depth at each variant site, taken from the simulated BAM
_bamdepth = {}
with pysam.AlignmentFile(D("collaborator_share", "S01.hg19.sorted.bam")) as _b:
    for _col in _b.pileup("chr11", 5245000, 5250000, min_base_quality=0):
        _bamdepth[_col.reference_pos + 1] = _col.nsegments
depth = lambda p: _bamdepth.get(p, 0)

# ------------------------------------------------------------------- the VCFs
def vcf_header(assembly, contig, length, sample):
    return textwrap.dedent(f"""\
        ##fileformat=VCFv4.2
        ##reference={assembly}
        ##contig=<ID={contig},length={length}>
        ##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
        ##INFO=<ID=NOTE,Number=1,Type=String,Description="Teaching annotation (synthetic data)">
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
        ##SYNTHETIC=These genotypes are teaching fiction; only the coordinates are real.
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
        """)


del_ref = ref19(DEL_POS, DEL_POS + DEL_LEN)          # anchor + 3 deleted bases
ins_ref = base19(INS_POS)
S01_RECS = [
    (CHECKPOINT, "HBB_checkpoint", "T", "A", 60, "PASS",
     "NOTE=public_coordinate_anchor", "0/1", depth(CHECKPOINT)),
    (DEL_POS, "S01_del1", del_ref, del_ref[0], 55, "PASS",
     "NOTE=padded_deletion", "0/1", depth(DEL_POS)),
    (INS_POS, "S01_ins1", ins_ref, ins_ref + "TA", 52, "PASS",
     "NOTE=padded_insertion", "0/1", depth(INS_POS)),
    (5246050, "S01_offchain", base19(5246050), "G", 44, "PASS",
     "NOTE=outside_teaching_chain", "0/1", depth(5246050)),
    (5248600, "S01_symdel", base19(5248600), "<DEL>", 30, "PASS",
     "SVTYPE=DEL;END=5248900;NOTE=symbolic_allele", "0/1", depth(5248600)),
    (5249000, "S01_bnd1", base19(5249000), f"{base19(5249000)}[chr11:5252000[", 25,
     "PASS", "SVTYPE=BND;NOTE=breakend_unsupported_by_exercise", "0/1", depth(5249000)),
]
assert S01_RECS[0][2] == "T" and S01_RECS[3][3] != S01_RECS[3][2]

with open(D("collaborator_share", "S01.hg19.vcf"), "w") as fh:
    fh.write(vcf_header("hg19/GRCh37 (teaching mini-chr11)", "chr11", H19_LEN, "S01"))
    for pos, vid, r, a, q, flt, info, gt, dp in S01_RECS:
        if dp == 0:   # no read support in this small prepared subset
            info, dp = info + ";NOTE2=no_read_coverage_in_prepared_subset", "."
        fh.write(f"chr11\t{pos}\t{vid}\t{r}\t{a}\t{q}\t{flt}\t{info}\tGT:DP\t{gt}:{dp}\n")

with open(D("collaborator_share", "S02.hg19.vcf"), "w") as fh:
    fh.write(vcf_header("hg19/GRCh37 (teaching mini-chr11)", "chr11", H19_LEN, "S02"))
    fh.write(f"chr11\t{CHECKPOINT}\tHBB_checkpoint\tT\tA\t58\tPASS\tNOTE=public_coordinate_anchor\tGT:DP\t1/1:{depth(CHECKPOINT)}\n")
    fh.write(f"chr11\t5247500\tS02_snv1\t{base19(5247500)}\t{'A' if base19(5247500)!='A' else 'G'}\t47\tPASS\tNOTE=synthetic\tGT:DP\t0/1:33\n")

# S03: header declares GRCh38 but the coordinates are hg19 -- the controlled failure
with open(D("collaborator_share", "S03.hg19.vcf"), "w") as fh:
    fh.write(vcf_header("hg38/GRCh38  <-- DECLARED, but positions are hg19", "chr11",
                        H38_LEN, "S03"))
    fh.write(f"chr11\t{CHECKPOINT}\tS03_site1\tT\tA\t51\tPASS\tNOTE=synthetic\tGT:DP\t0/1:29\n")

# toy1 VCF for the pure-semantics exercises (slides 51-54)
with open(D("collaborator_share", "toy1.vcf"), "w") as fh:
    fh.write(vcf_header("toy1 synthetic contig", "toy1", 1000, "TOY"))
    fh.write("toy1\t100\td1\tATC\tA\t60\tPASS\tNOTE=padded_deletion\tGT:DP\t0/1:30\n")
    fh.write("toy1\t300\tsnv1\tA\tG\t60\tPASS\tNOTE=substitution\tGT:DP\t0/1:30\n")

# ------------------------------------------------------------ the manifest
with open(D("manifest.tsv"), "w") as fh:
    cols = ["sample_id", "source_assembly", "target_assembly", "source_fasta",
            "target_fasta", "chain", "regions", "vcf", "bam", "annotation_release",
            "expected_status", "note"]
    fh.write("\t".join(cols) + "\n")
    rows = [
        ["S01", "hg19", "hg38", "collaborator_share/reference.hg19.fa",
         "target/reference.hg38.fa", "target/hg19ToHg38.teaching.over.chain",
         "collaborator_share/regions.S01.hg19.tsv", "collaborator_share/S01.hg19.vcf",
         "collaborator_share/S01.hg19.sorted.bam", "RefSeq Select (UCSC hg38, 2026-10)",
         "pass", "full teaching sample"],
        ["S02", "hg19", "hg38", "collaborator_share/reference.hg19.fa",
         "target/reference.hg38.fa", "target/hg19ToHg38.teaching.over.chain",
         "collaborator_share/regions.S02.hg19.tsv", "collaborator_share/S02.hg19.vcf",
         "", "RefSeq Select (UCSC hg38, 2026-10)", "pass", "variants only, no BAM"],
        ["S03", "hg38", "hg38", "target/reference.hg38.fa",
         "target/reference.hg38.fa", "", "collaborator_share/regions.S03.hg19.tsv",
         "collaborator_share/S03.hg19.vcf", "", "RefSeq Select (UCSC hg38, 2026-10)",
         "fail", "manifest says the source is already hg38 - check before believing it"],
    ]
    for r in rows:
        fh.write("\t".join(r) + "\n")

print("built", sum(len(f) for _, _, f in os.walk(OUT)), "files")
