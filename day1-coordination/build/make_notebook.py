import nbformat as nbf

nb = nbf.v4.new_notebook()
C, M = [], []
cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
co = lambda s: cells.append(nbf.v4.new_code_cell(s.strip("\n")))

md(r"""
# Genomic Data Processing — Day 1 practice session

**Data Processing for Big Science · visiting lecture · 10/12**

You have 30 minutes. Most of the code is written for you. Your job is to look at the
data, decide what each step means, and keep one piece of evidence from each checkpoint
for your 5-minute presentation.

### The story
A collaborator sends you results from a study processed on **hg19 (GRCh37)**.
Your annotation workflow runs on **hg38 (GRCh38)**. Both files say `chr11` and both
contain numbers that look like positions. The numbers do not mean the same thing.

### What is real and what is synthetic — read this before you trust anything
| Part of the package | Status |
|---|---|
| `chr11` sequence in both mini-references | **Real** UCSC sequence (hg19 chr11:5,240,001–5,300,000 and the identical hg38 chr11:5,218,771–5,278,770) |
| Gene and exon coordinates in `genes.hg38.bed` | **Real** NCBI RefSeq Select transcripts for the β-globin cluster |
| `HBB_checkpoint` correspondence hg19 5,248,232 ↔ hg38 5,227,002 | **Real**, independently documented (ClinVar RCV005357130, NM_000518.5(HBB):c.20A>T) |
| Samples S01–S03, their reads, genotypes and depths | **Synthetic teaching fiction** |
| `hg19ToHg38.teaching.over.chain` | **Teaching excerpt**, not the UCSC chain file (see §4) |
| `toy1` contig | **Synthetic**, exists only to make the VCF examples readable |

Nothing here is a clinical result. We use a public variant site **only as a coordinate
anchor**, exactly as slide 101 says.

### Notebook sections (these match slide 83)
1. Reference / history 2. hg19 coordinates 3. VCF semantics 4. hg19 → hg38 mapping
5. Assembly-matched IGV 6. Read realignment 7. Batch / QC 8. Evidence card
""")

co(r"""
import os, gzip, re, textwrap
from pathlib import Path
import pandas as pd

BASE = Path("practice_session")           # change this if you unpacked elsewhere
SHARE, TARGET = BASE / "collaborator_share", BASE / "target"
CHECK, OUT = BASE / "checkpoints", BASE / "outputs"
OUT.mkdir(exist_ok=True)

try:
    import pysam
    HAVE_PYSAM = True
except ImportError:
    HAVE_PYSAM = False      # everything still works: we ship .sam twins of every .bam

print("pysam available:", HAVE_PYSAM)
print(sorted(p.name for p in SHARE.iterdir()))
""")

# ------------------------------------------------------------------ section 1
md(r"""
---
## 1 · Reference / history — *practice steps 01–02*

> **Slide 59 asked:** you receive a folder of genomic files. Which file do you open first?

Not the biggest one. Open the **manifest** and the **headers**, because those are the only
places where a file says what reference it belongs to. A position without an assembly is
not an address; it is just a number.
""")

co(r"""
manifest = pd.read_csv(BASE / "manifest.tsv", sep="\t").fillna("")
manifest[["sample_id", "source_assembly", "target_assembly", "expected_status", "note"]]
""")

md(r"""
### 1.1 Read-only header checks

On a real folder you would run

```bash
samtools view -H S01.sorted.bam     # @SQ lines = the sequence dictionary
bcftools view -h  S01.vcf.gz        # ##reference and ##contig lines
```

Here we do the same thing in Python so the notebook runs anywhere. **Read only** — we are
establishing facts, not fixing anything yet.
""")

co(r"""
def fasta_dict(fai_path):
    '''Sequence dictionary from a .fai index: {name: length}.'''
    return {l.split("\t")[0]: int(l.split("\t")[1])
            for l in open(fai_path) if l.strip()}

def vcf_header(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return [l.rstrip("\n") for l in fh if l.startswith("#")]

def sam_header(path):
    '''Works on the .sam twin, so no pysam needed.'''
    with open(path) as fh:
        return [l.rstrip("\n") for l in fh if l.startswith("@")]

print("hg19 mini-reference :", fasta_dict(SHARE / "reference.hg19.fa.fai"))
print("hg38 mini-reference :", fasta_dict(TARGET / "reference.hg38.fa.fai"))
print()
for line in vcf_header(SHARE / "S01.hg19.vcf"):
    if line.startswith(("##reference", "##contig", "##SYNTHETIC")):
        print("S01 VCF  |", line)
print()
for line in sam_header(SHARE / "S01.hg19.sorted.sam"):
    print("S01 BAM  |", line[:110])
""")

md(r"""
**Stop and read the output above.** Two things should bother you:

* The contig is called `chr11` in *every* file. Identical names, two different coordinate
  systems. Nothing in the name protects you.
* The two mini-references have **different lengths** (5,300,000 vs 5,280,000). A length
  mismatch is one of the cheapest reference checks you can run, and it is the check that
  catches the problem below.

### 1.2 The compatibility table
We now produce the small table slide 60 asks for: one row per input, with `pass`,
`review` or `fail` **and a reason**. A tool that silently drops the reason is a tool that
loses the only information you needed.
""")

co(r"""
def compat_row(sample, declared, vcf_path, fai_path):
    hdr = "\n".join(vcf_header(vcf_path))
    declared_in_file = re.search(r"##reference=(\S+)", hdr)
    declared_in_file = declared_in_file.group(1) if declared_in_file else "ABSENT"
    contig_len = int(re.search(r"##contig=<ID=chr11,length=(\d+)>", hdr).group(1))
    fasta_len = fasta_dict(fai_path)["chr11"]
    if contig_len != fasta_len:
        status, reason = "fail", (f"VCF declares chr11 length {contig_len:,} but the "
                                  f"{declared} FASTA is {fasta_len:,}")
    elif declared_in_file.split("/")[0] != declared:
        status, reason = "review", f"manifest says {declared}, header says {declared_in_file}"
    else:
        status, reason = "pass", "-"
    return dict(sample_id=sample, manifest_says=declared,
                header_says=declared_in_file.split()[0], vcf_contig_len=contig_len,
                fasta_len=fasta_len, status=status, reason=reason)

rows = []
for _, m in manifest.iterrows():
    fai = BASE / (m.source_fasta + ".fai")
    rows.append(compat_row(m.sample_id, m.source_assembly, BASE / m.vcf, fai))
compat = pd.DataFrame(rows)
compat.to_csv(OUT / "compatibility_table.tsv", sep="\t", index=False)
compat
""")

md(r"""
### Checkpoint 1 — keep this for your presentation

**S03 fails, and it fails for a reason you can state in one sentence.** Its manifest row and
its VCF header both claim the data are already on hg38, but the declared contig length is
the hg19 mini-chromosome's. The positions inside are hg19 positions wearing an hg38 label.

This is the controlled failure from slide 80. Notice what we did **not** do: we did not
rename anything, and we did not convert S03 "just to see". A file whose source build is
unknown or contradicted has no defensible starting point for conversion (slide 62).

> **Q.** The manifest and the header disagree. Which one do you believe?
> **A.** Neither, until a third piece of evidence breaks the tie. Record the conflict and stop.
""")

# ------------------------------------------------------------------ section 2
md(r"""
---
## 2 · hg19 coordinates — *practice step 03*

Two conventions, five bases, one lifetime of off-by-one errors:

```
1-based position   1   2   3   4   5
base               A   C   G   T   A
BED boundary     0   1   2   3   4   5
```

* **1-based closed** `[2, 4]` → counts the first base as 1 and includes both ends → `C G T`
* **0-based half-open (BED)** `[1, 4)` → counts the *gaps between* bases and excludes the end → `C G T`

Same three bases, different arithmetic. The conversion is `start0 = start1 - 1`, `end0 = end1`
— subtract one from the start **only**. The length check is what makes it true:
`end1 - start1 + 1 == end0 - start0`.
""")

co(r"""
def one_based_to_bed(start1, end1):
    '''1-based closed [start1, end1] -> 0-based half-open [start0, end0).'''
    if start1 > end1:
        raise ValueError(f"invalid interval: start {start1} > end {end1}")
    return start1 - 1, end1

# the two cases you worked by hand on slide 63
assert one_based_to_bed(101, 200) == (100, 200)
assert one_based_to_bed(101, 101) == (100, 101)
assert (200 - 101 + 1) == (200 - 100) == 100        # length preserved
assert (101 - 101 + 1) == (101 - 100) == 1          # single base survives
print("conversion assertions passed")
""")

md(r"""
Those assertions are deliberately boring, and that is the point. A single-base interval is
the cheapest test that exposes an off-by-one: if you had subtracted one from *both* ends you
would have produced `[100, 100)`, an interval containing **zero** bases, and nothing would
have crashed.

### 2.1 Convert the collaborator's legacy interval table
Note the two rows we treat differently from the rest.
""")

co(r"""
def convert_region_table(path):
    df = pd.read_csv(path, sep="\t")
    out = []
    for _, r in df.iterrows():
        rec = dict(region_id=r.region_id, chrom=r.chrom,
                   start_1based=r.start_1based, end_1based=r.end_1based, note=r.note)
        try:
            s0, e0 = one_based_to_bed(int(r.start_1based), int(r.end_1based))
            rec.update(bed_start=s0, bed_end=e0, length=e0 - s0, qc_status="pass", reason="-")
        except ValueError as err:
            # slide 85: do NOT silently swap start and end
            rec.update(bed_start=pd.NA, bed_end=pd.NA, length=pd.NA,
                       qc_status="fail", reason=str(err))
        out.append(rec)
    return pd.DataFrame(out)

s01_bed = convert_region_table(SHARE / "regions.S01.hg19.tsv")
s01_bed
""")

md(r"""
`r05` has `start > end`. Swapping the two numbers would make the table look clean and would
invent a region nobody asked for. We keep it, mark it `fail`, and carry the reason forward.

### 2.2 The boundary test
Half-open intervals **touch without overlapping**. `[199, 200)` ends at boundary 200;
`[200, 250)` starts there. They share no base.
""")

co(r"""
def overlap_bp(a, b):
    '''Shared bases between two non-empty half-open intervals.'''
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

print("[199,200) vs [200,250) ->", overlap_bp((199, 200), (200, 250)), "bp shared")
print("[199,201) vs [200,250) ->", overlap_bp((199, 201), (200, 250)), "bp shared")
print("[100,200) vs [200,250) ->", overlap_bp((100, 200), (200, 250)), "bp shared")
""")

md(r"""
### Checkpoint 2 — keep this
A **one-base** shift flips a gene-overlap answer from *no* to *yes* at a boundary, and no
error message appears either way. This is the evidence Student 1 presents.
""")

# ------------------------------------------------------------------ section 3
md(r"""
---
## 3 · VCF semantics — *practice step 04*

A BED interval says *where*. A VCF record says *where, what was there, what replaced it, and
in whom*. Converting the second into the first always throws something away — the exercise is
to throw away the right things, on purpose, and write down what you did.

### 3.1 The padded deletion
```
POS=100  REF=ATC  ALT=A
reference  100:A  101:T  102:C
alternate  100:A                 <- T and C are gone, A is still there
```
`POS` is **not** the first deleted base. It is the anchor: a base kept in both alleles so
that neither REF nor ALT is an empty string. Two different intervals can be derived, and they
answer different questions:

| interval | BED | answers |
|---|---|---|
| reference span incl. anchor | `[99, 102)` | "which reference bases does this record describe?" |
| deleted bases only | `[100, 102)` | "which bases actually disappeared?" |

Our overlap exercise uses **deleted bases**. We say so in the output schema rather than hoping
the next person guesses.
""")

co(r"""
def read_vcf(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    recs = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            recs.append(dict(chrom=f[0], pos=int(f[1]), vid=f[2], ref=f[3], alt=f[4],
                             qual=f[5], filt=f[6], info=info,
                             fmt=f[8] if len(f) > 8 else "",
                             sample=f[9] if len(f) > 9 else ""))
    return recs

def vcf_to_interval(rec):
    '''Classify a record and return the affected-base BED interval, or a review status.'''
    ref, alt, pos = rec["ref"], rec["alt"], rec["pos"]
    if alt.startswith("<") or "[" in alt or "]" in alt:
        kind = "symbolic_del" if alt == "<DEL>" else "unsupported"
        if kind == "symbolic_del" and "END" in rec["info"]:
            return dict(interval_kind="deleted_bases_symbolic", start0=pos,
                        end0=int(rec["info"]["END"]), status="review",
                        reason="symbolic allele: END is taken from INFO, not from sequence")
        return dict(interval_kind=None, start0=pd.NA, end0=pd.NA, status="review",
                    reason=f"unsupported allele representation: {alt[:24]}")
    if len(ref) == 1 and len(alt) == 1:
        return dict(interval_kind="substituted_base", start0=pos - 1, end0=pos,
                    status="mapped", reason="-")
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:
        return dict(interval_kind="deleted_bases", start0=pos, end0=pos + len(ref) - 1,
                    status="mapped", reason="-")
    if len(alt) > 1 and len(ref) == 1 and ref[0] == alt[0]:
        # an insertion adds bases between two reference positions: zero reference bases
        # are affected, so an "affected-base" interval is the wrong question to ask
        return dict(interval_kind="insertion_point", start0=pos - 1, end0=pos,
                    status="review",
                    reason="insertion: the anchor base is reported, no reference base is lost")
    return dict(interval_kind=None, start0=pd.NA, end0=pd.NA, status="review",
                reason="complex allele: outside this exercise's policy")

toy = read_vcf(SHARE / "toy1.vcf")
pd.DataFrame([{**{k: r[k] for k in ("vid", "pos", "ref", "alt")}, **vcf_to_interval(r)}
              for r in toy])
""")

md(r"""
Check the two toy rows against what you worked out on slides 52–54:
`d1` → `[100, 102)` (the anchor at 100 is *excluded*), `snv1` → `[299, 300)` (one base, not
the empty `[299, 299)`).

### 3.2 Now the human records — and verify REF against the actual sequence
A VCF can claim anything in its REF column. The reference FASTA is the arbiter. This check
costs three lines and catches wrong-build files, wrong-strand files and off-by-ones at once.
""")

co(r"""
class MiniFasta:
    '''Minimal random-access FASTA reader driven by the .fai index.'''
    def __init__(self, path):
        self.path = Path(path)
        name, length, offset, bpl, bypl = open(str(path) + ".fai").readline().split("\t")
        self.name, self.length = name, int(length)
        self.offset, self.bpl, self.bypl = int(offset), int(bpl), int(bypl)
        self.fh = open(path, "rb")

    def fetch(self, start1, end1):
        '''1-based inclusive.'''
        s, e = start1 - 1, end1
        self.fh.seek(self.offset + s // self.bpl * self.bypl + s % self.bpl)
        need = e - s
        raw = self.fh.read(need + need // self.bpl + 2).decode()
        return raw.replace("\n", "")[:need].upper()

hg19 = MiniFasta(SHARE / "reference.hg19.fa")
hg38 = MiniFasta(TARGET / "reference.hg38.fa")

# sanity: the 18 reference bases drawn on slide 46
assert hg38.fetch(5226995, 5227012) == "CTTCTCCTCAGGAGTCAG"
print("hg38 chr11:5,226,995-5,227,012 =", hg38.fetch(5226995, 5227012))

s01_vcf = read_vcf(SHARE / "S01.hg19.vcf")
rows = []
for r in s01_vcf:
    obs = hg19.fetch(r["pos"], r["pos"] + len(r["ref"]) - 1)
    rows.append(dict(vid=r["vid"], pos=r["pos"], ref=r["ref"], alt=r["alt"][:20],
                     ref_in_fasta=obs, ref_ok=(obs == r["ref"]), **vcf_to_interval(r)))
s01_variants = pd.DataFrame(rows)
s01_variants[["vid", "pos", "ref", "alt", "ref_in_fasta", "ref_ok",
              "interval_kind", "start0", "end0", "status", "reason"]]
""")

md(r"""
### Checkpoint 3 — keep this
Every S01 REF allele matches the hg19 sequence, so the file really is hg19. Three of six
records convert cleanly; the insertion, the symbolic `<DEL>` and the breakend are kept with a
`review` status and their original fields. They are **not** forced into an ordinary interval
just so the output table looks complete.

> **Q.** Can you rebuild the original VCF from the exported BED?
> **A.** No. `chr11 5248231 5248232 HBB_checkpoint` has lost REF, ALT, the genotype, the
> sample and the build. That is why we keep `source_POS`, `source_REF`, `source_ALT` and a
> stable ID next to every interval. This is Student 2's evidence.
""")

# ------------------------------------------------------------------ section 4
md(r"""
---
## 4 · hg19 → hg38 mapping — *practice step 05*

Now the real conversion. Two separate transformations, and students mix them up every year:

1. **convention**: 1-based position → BED interval, *on hg19* (no biology involved)
2. **assembly**: BED interval on hg19 → BED interval on hg38 (biology and alignment involved)

Step 2 cannot be done with arithmetic you invent. It needs a **chain file**: a stored
alignment between the two assemblies, block by block.

### About the chain we ship
`hg19ToHg38.teaching.over.chain` is **not** the UCSC file. It is a single-block excerpt over
`chr11:5,246,500–5,280,000` (hg19). Within that span the two assemblies are byte-identical,
which the next cell proves from the two FASTAs, so hg38 = hg19 − 21,230 holds **exactly
here**. Outside the block, mapping fails on purpose — that is the `unmapped` output you need
to see. For anything real, use:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
liftOver input.hg19.bed hg19ToHg38.over.chain.gz mapped.hg38.bed unmapped.hg19.bed
```
""")

co(r"""
print(open(TARGET / "hg19ToHg38.teaching.over.chain").read())

# prove the claim rather than asserting it: compare 60 kb of both mini-references
a = hg19.fetch(5240001, 5300000)
b = hg38.fetch(5240001 - 21230, 5300000 - 21230)
print("60,000 bp identical between the two assemblies here:", a == b)
""")

co(r"""
def read_chain(path):
    '''Parse a (simple, ungapped-block) chain file into mapping blocks.'''
    blocks = []
    with open(path) as fh:
        hdr = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("chain"):
                f = line.split()
                hdr = dict(t_name=f[2], t_start=int(f[5]), t_end=int(f[6]),
                           q_name=f[7], q_start=int(f[10]), q_end=int(f[11]))
            elif hdr:
                size = int(line.split()[0])
                blocks.append(dict(src_chrom=hdr["t_name"], src_start=hdr["t_start"],
                                   src_end=hdr["t_start"] + size, dst_chrom=hdr["q_name"],
                                   dst_start=hdr["q_start"]))
                hdr = None
    return blocks

CHAIN = read_chain(TARGET / "hg19ToHg38.teaching.over.chain")

def lift(chrom, start0, end0, chain=CHAIN):
    '''Map a half-open interval. Returns (chrom, start, end) or (None, reason).'''
    for b in chain:
        if b["src_chrom"] != chrom:
            continue
        if start0 >= b["src_start"] and end0 <= b["src_end"]:
            shift = b["dst_start"] - b["src_start"]
            return (b["dst_chrom"], start0 + shift, end0 + shift), None
        if start0 < b["src_end"] and end0 > b["src_start"]:
            return None, "partially outside the chain block (split mapping)"
    return None, "no chain block covers this interval"

# the documented checkpoint: hg19 chr11:5,248,232 -> hg38 chr11:5,227,002
got, err = lift("chr11", 5248231, 5248232)
print("lifted:", got, "| error:", err)
assert got == ("chr11", 5227001, 5227002), got
print("matches the ClinVar-documented hg38 position 5,227,002 ->", got[2])
""")

md(r"""
That `assert` is the heart of the whole session. We did not check that the command *ran*; we
checked the answer against a coordinate that was documented **independently of our pipeline**
(ClinVar RCV005357130: `NC_000011.9:g.5248232T>A` ↔ `NC_000011.10:g.5227002T>A`). A pipeline
that agrees with an outside source at a known point has earned a little trust. It has not
earned trust everywhere — see the unmapped rows below.

### 4.1 Lift the interval table and the variants, keeping both outputs
""")

co(r"""
def lift_table(df, source_sample):
    out = []
    for _, r in df.iterrows():
        rec = dict(sample_id=source_sample, source_id=r.region_id, source_assembly="hg19",
                   target_assembly="hg38", source_start=r.bed_start, source_end=r.bed_end,
                   note=r.note)
        if r.qc_status == "fail":
            rec.update(target_start=pd.NA, target_end=pd.NA,
                       mapping_status="not_attempted", reason=r.reason)
        else:
            got, err = lift(r.chrom, int(r.bed_start), int(r.bed_end))
            if got:
                rec.update(target_start=got[1], target_end=got[2],
                           mapping_status="mapped", reason="-")
            else:
                rec.update(target_start=pd.NA, target_end=pd.NA,
                           mapping_status="unmapped", reason=err)
        out.append(rec)
    return pd.DataFrame(out)

lifted = lift_table(s01_bed, "S01")

with open(OUT / "mapped.hg38.bed", "w") as fh:
    for _, r in lifted[lifted.mapping_status == "mapped"].iterrows():
        fh.write(f"chr11\t{r.target_start}\t{r.target_end}\t{r.source_id}\n")
with open(OUT / "unmapped.hg19.bed", "w") as fh:
    for _, r in lifted[lifted.mapping_status != "mapped"].iterrows():
        fh.write(f"#{r.mapping_status}\t{r.reason}\n")
        fh.write(f"chr11\t{r.source_start}\t{r.source_end}\t{r.source_id}\n")

print("mapped.hg38.bed\n" + open(OUT / "mapped.hg38.bed").read())
print("unmapped.hg19.bed\n" + open(OUT / "unmapped.hg19.bed").read())
lifted
""")

md(r"""
### 4.2 Variants need allele-aware conversion, not just a new number

Moving `POS` is the easy part. After the move you must re-check **REF against the target
reference**, because the target genome does not have to agree with the source. In production
this is `picard LiftoverVcf` with the chain *and* the target FASTA:

```bash
java -jar picard.jar LiftoverVcf I=S01.hg19.vcf O=small.hg38.vcf \
     CHAIN=hg19ToHg38.over.chain.gz REJECT=rejected.hg19.vcf R=hg38.fa
```

Below we do the same bookkeeping by hand so you can watch the rejects accumulate.
""")

co(r"""
accepted, rejected = [], []
for r in s01_vcf:
    got, err = lift("chr11", r["pos"] - 1, r["pos"] - 1 + len(r["ref"]))
    if not got:
        rejected.append((r, err)); continue
    new_pos = got[1] + 1
    ref_on_target = hg38.fetch(new_pos, new_pos + len(r["ref"]) - 1)
    if ref_on_target != r["ref"]:
        rejected.append((r, f"REF {r['ref']} != hg38 sequence {ref_on_target}")); continue
    if r["alt"].startswith("<") or "[" in r["alt"]:
        rejected.append((r, "allele type outside this exercise's supported set")); continue
    accepted.append((r, new_pos))

with open(OUT / "small.hg38.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n##reference=hg38 teaching mini-chr11\n")
    fh.write("##contig=<ID=chr11,length=%d>\n" % hg38.length)
    fh.write("##INFO=<ID=SRC,Number=1,Type=String,Description=\"source hg19 position\">\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS01\n")
    for r, np_ in accepted:
        fh.write(f"chr11\t{np_}\t{r['vid']}\t{r['ref']}\t{r['alt']}\t{r['qual']}\t"
                 f"{r['filt']}\tSRC=chr11:{r['pos']}\t{r['fmt']}\t{r['sample']}\n")
with open(OUT / "rejected.hg19.vcf", "w") as fh:
    fh.write("##fileformat=VCFv4.2\n##reference=hg19 (SOURCE coordinates retained)\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for r, why in rejected:
        fh.write(f"chr11\t{r['pos']}\t{r['vid']}\t{r['ref']}\t{r['alt']}\t{r['qual']}\t"
                 f"{r['filt']}\tREJECT_REASON={why.replace(' ', '_')}\n")

print(open(OUT / "small.hg38.vcf").read())
print(open(OUT / "rejected.hg19.vcf").read())
""")

md(r"""
### Checkpoint 4 — keep this
Three records converted, three were rejected **with their source coordinates intact**. The
rejected file is not a failure of the workflow; it is part of the result. A conversion that
reports 100% success on this input would be hiding three records from you.

This is Student 3's evidence: *two* transformations, one independent verification, and a
rejects file that still points back at the original data.
""")

# ------------------------------------------------------------------ section 5
md(r"""
---
## 5 · Assembly-matched IGV — *practice step 06*

Rule, no exceptions: **one session per assembly.** An hg19 BAM under an hg38 ruler will load
happily, draw neatly, and be wrong.

The next cell writes two IGV session files. Open them with *File → Open Session* after
loading each mini-reference as a custom genome (*Genomes → Load Genome from File*).
""")

co(r"""
def igv_session(path, genome_fa, locus, tracks):
    t = "\n".join(f'        <Resource path="{p}"/>' for p in tracks)
    xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<Session genome="{genome_fa}" locus="{locus}" version="8">
    <Resources>
{t}
    </Resources>
</Session>
'''
    Path(path).write_text(xml)
    return xml

rel = lambda p: os.path.relpath(p, OUT)
igv_session(OUT / "session_hg19.xml", rel(SHARE / "reference.hg19.fa"),
            "chr11:5,248,200-5,248,270",
            [rel(SHARE / "S01.hg19.vcf"), rel(SHARE / "S01.hg19.sorted.bam")])
igv_session(OUT / "session_hg38.xml", rel(TARGET / "reference.hg38.fa"),
            "chr11:5,226,970-5,227,040",
            [rel(OUT / "small.hg38.vcf"), rel(OUT / "mapped.hg38.bed"),
             rel(CHECK / "S01.hg38.sorted.bam"), rel(TARGET / "genes.hg38.bed")])
print("wrote", OUT / "session_hg19.xml", "and", OUT / "session_hg38.xml")
""")

md(r"""
### 5.1 The empty-track exercise (slide 73)
You open the hg38 session, type `chr11:5,248,232` into the locus box, and the VCF track is
empty. Before concluding anything biological, reproduce the situation here.
""")

co(r"""
def features_in_view(vcf_path, start1, end1):
    return [r["vid"] for r in read_vcf(vcf_path) if start1 <= r["pos"] <= end1]

print("hg38 VCF at the hg19 address chr11:5,248,223-5,248,240 ->",
      features_in_view(OUT / "small.hg38.vcf", 5248223, 5248240))
print("hg38 VCF at the hg38 address chr11:5,226,995-5,227,012 ->",
      features_in_view(OUT / "small.hg38.vcf", 5226995, 5227012))
""")

md(r"""
The track was never broken. You were standing 21,230 bases away. The troubleshooting order
that gets you there fastest: **which genome is loaded → which locus am I at → is the track
loaded and visible → is the file the one I think it is.** Only after all four do you start
asking biological questions.

### Checkpoint 5 — keep this
A labelled before/after pair: same record, same stable ID `HBB_checkpoint`, two addresses.
Put the assembly **and** the locus on every screenshot you paste into a slide. This is
Student 4's evidence.
""")

# ------------------------------------------------------------------ section 6
md(r"""
---
## 6 · Read realignment — *practice steps 07–08*

**Lifting a VCF does not move a single read.** The BAM still points at hg19. To get read
placements on hg38 you must align the original reads again, against the new reference.

Production command (do **not** run this in class — nobody builds a human index in four
minutes):

```bash
bwa mem -t 8 hg38.fa S01_R1.fastq.gz S01_R2.fastq.gz | samtools sort -o S01.hg38.sorted.bam
samtools index S01.hg38.sorted.bam
```

What we run instead is a 20-line seed-and-extend aligner restricted to our 60 kb window.
Slide 78 is explicit about what that means: **a tiny search space is an approximation, not
full-reference mapping.** A read that is unique in 60 kb may have twenty equally good homes in
3.1 billion. We use it to show the mechanics; we do not report its MAPQ as if it were real.
""")

co(r"""
SEED = 25
WIN_S, WIN_E = 5247000, 5249500                      # hg19 window we align into
win19 = hg19.fetch(WIN_S, WIN_E)
win38 = hg38.fetch(WIN_S - 21230, WIN_E - 21230)

def build_index(seq, offset1, k=SEED):
    # k-mer -> 1-based position; k-mers that occur more than once are dropped
    idx = {}
    for i in range(len(seq) - k + 1):
        idx.setdefault(seq[i:i + k], []).append(offset1 + i)
    return {kmer: v[0] for kmer, v in idx.items() if len(v) == 1}

IDX = {"hg19": build_index(win19, WIN_S), "hg38": build_index(win38, WIN_S - 21230)}
REF = {"hg19": (win19, WIN_S), "hg38": (win38, WIN_S - 21230)}
print({k: f"{len(v):,} unique {SEED}-mers" for k, v in IDX.items()})

def mini_align(read, asm):
    # anchor the first and last seed that occur exactly once, then describe the
    # middle with a single indel operation. Teaching aligner: 2.5 kb search space.
    idx, (ref, roff) = IDX[asm], REF[asm]
    head = next(((j, idx[read[j:j + SEED]]) for j in range(len(read) - SEED + 1)
                 if read[j:j + SEED] in idx), None)
    tail = next(((j, idx[read[j:j + SEED]]) for j in range(len(read) - SEED, -1, -1)
                 if read[j:j + SEED] in idx), None)
    if head is None or tail is None:
        return None
    jh, ph = head
    jt, pt = tail
    start = ph - jh                                   # if the prefix aligns too
    prefix_mm = (sum(1 for i in range(jh) if read[i] != ref[start - roff + i])
                 if start >= roff else 99)
    clip = 0 if prefix_mm <= 2 else jh                # otherwise soft-clip it
    if clip:
        start = ph
    body = len(read) - clip
    delta = (pt + SEED - start) - body                # reference span minus read span
    # first mismatch after the anchor = where the indel sits
    split = 0
    while split < body - abs(delta) and read[clip + split] == ref[start - roff + split]:
        split += 1
    pre = [(clip, "S")] if clip else []
    if delta == 0:
        cig = pre + [(body, "M")]
    elif delta > 0:
        cig = pre + [(split, "M"), (delta, "D"), (body - split, "M")]
    else:
        cig = pre + [(split, "M"), (-delta, "I"), (body + delta - split, "M")]
    cig = [(n, o) for n, o in cig if n > 0]
    return dict(pos=start, cigar="".join(f"{n}{o}" for n, o in cig))

def read_fastq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            yield h[1:].strip(), fh.readline().strip(), fh.readline(), fh.readline().strip()

rows = []
for name, seq, _, _ in read_fastq(SHARE / "S01_R1.fastq.gz"):
    a19, a38 = mini_align(seq, "hg19"), mini_align(seq, "hg38")
    if a19 and a38:
        rows.append(dict(read_id=name, hg19_pos=a19["pos"], hg19_cigar=a19["cigar"],
                         hg38_pos=a38["pos"], hg38_cigar=a38["cigar"],
                         pos_delta=a19["pos"] - a38["pos"],
                         cigar_changed=a19["cigar"] != a38["cigar"]))
cmp_df = pd.DataFrame(rows)
cmp_df.to_csv(OUT / "read_comparison.tsv", sep="\t", index=False)
print(f"{len(cmp_df)} R1 reads placed on both assemblies")
print("position shift values observed:", sorted(cmp_df.pos_delta.unique()))
print("reads whose CIGAR changed:", int(cmp_df.cigar_changed.sum()))
print("CIGAR shapes seen:", sorted(cmp_df.hg19_cigar.unique())[:8])
cmp_df.head(8)
""")

md(r"""
Read that output carefully, because it is the single most-misunderstood point of the day:

* **Every read moved by exactly 21,230 bases.** The sequences never changed; the reference did.
* **No CIGAR changed.** In this window the two assemblies are identical, so the alignment
  *shape* is preserved and only the address differs. A changed coordinate does not imply a
  changed alignment, and an unchanged CIGAR does not imply nothing happened.

One honest caveat about the toy aligner: a few reads that carry the inserted bases come back
as soft-clipped (`21S59M`) where `bwa mem` would open a 2 bp insertion. The CIGARs are
internally consistent and identical on both assemblies, so the comparison above still holds,
but this is exactly the kind of artefact a 2.5 kb search space produces. Never quote a MAPQ or
a clipping pattern from a toy index as if it came from a genome-wide alignment.

Compare that with the records we lifted in §4: those were *variant coordinates*, converted
with a chain file. These are *reads*, placed again by an aligner. Different operation,
different failure modes, different evidence.

### 6.1 Walk one CIGAR by hand (slides 45, 75–76)
""")

co(r"""
def walk_cigar(pos1, cigar):
    ops = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    read_used = ref_used = 0
    ref_at = pos1
    print(f"{'op':>6} {'read +':>7} {'ref +':>7}   reference positions")
    for n, o in ops:
        n = int(n)
        dr = n if o in "MIS=X" else 0
        dref = n if o in "MDN=X" else 0
        span = f"{ref_at}-{ref_at + dref - 1}" if dref else "(no reference)"
        print(f"{n}{o:>5} {dr:>7} {dref:>7}   {span}")
        read_used += dr; ref_used += dref; ref_at += dref
    print(f"\ntotal read bases consumed      : {read_used}")
    print(f"total reference bases consumed : {ref_used}")
    print(f"inclusive alignment end        : {pos1} + {ref_used} - 1 = {pos1 + ref_used - 1}")

walk_cigar(101, "10M2I5M3D8M")
""")

md(r"""
25 read bases, 26 reference bases, end at 126 — and the deleted reference positions 116–118
are inside the alignment span while containing no read base at all. That is why an interval
derived from a BAM record is not a statement that every base in it was observed.

### Checkpoint 6 — keep this
`read_comparison.tsv`: one read ID, two assemblies, position changed by 21,230, CIGAR
unchanged, and one sentence naming the limitation (60 kb search space, simulated reads).
""")

# ------------------------------------------------------------------ section 7
md(r"""
---
## 7 · Batch and QC — *practice step 09*

One sample teaches you the file format. A cohort is where the science is — and where silent
failures become invisible. The rule from slide 91: **a failed conversion must never become a
zero gene count.**

* `0` = the overlap check ran against a declared annotation and found nothing. A result.
* `NA` = the check did not run, or the annotation was unavailable. Not a result.
""")

co(r"""
def load_bed(path):
    feats = []
    for line in open(path):
        if line.startswith(("#", "track")):
            continue
        f = line.rstrip("\n").split("\t")
        feats.append((f[0], int(f[1]), int(f[2]), f[3]))
    return feats

GENES38 = load_bed(TARGET / "genes.hg38.bed")

def genes_overlapping(chrom, s0, e0, feats=GENES38):
    return sorted({n for c, fs, fe, n in feats
                   if c == chrom and overlap_bp((s0, e0), (fs, fe)) > 0})

features = []
for _, m in manifest.iterrows():
    sample = m.sample_id
    row_status = compat.loc[compat.sample_id == sample, "status"].iloc[0]
    tbl = convert_region_table(BASE / m.regions)
    if row_status == "fail":
        for _, r in tbl.iterrows():
            features.append(dict(sample_id=sample, source_id=r.region_id,
                                 mapping_status="not_attempted", target_start=pd.NA,
                                 target_end=pd.NA, n_distinct_genes=pd.NA,
                                 genes="NA", qc_status="fail",
                                 reason="source assembly unresolved - see compatibility table"))
        continue
    for _, r in lift_table(tbl, sample).iterrows():
        if r.mapping_status != "mapped":
            features.append(dict(sample_id=sample, source_id=r.source_id,
                                 mapping_status=r.mapping_status, target_start=pd.NA,
                                 target_end=pd.NA, n_distinct_genes=pd.NA, genes="NA",
                                 qc_status="review", reason=r.reason))
        else:
            g = genes_overlapping("chr11", int(r.target_start), int(r.target_end))
            features.append(dict(sample_id=sample, source_id=r.source_id,
                                 mapping_status="mapped", target_start=int(r.target_start),
                                 target_end=int(r.target_end), n_distinct_genes=len(g),
                                 genes=",".join(g) if g else "-", qc_status="pass",
                                 reason="-"))

features = pd.DataFrame(features)
features.insert(1, "source_assembly", "hg19")
features.insert(2, "target_assembly", "hg38")
features["annotation_release"] = "RefSeq Select (UCSC hg38, retrieved 2026-10)"
features["interval_kind"] = "region_table_interval"
features.to_csv(OUT / "features.hg38.tsv", sep="\t", index=False)
features[["sample_id", "source_id", "mapping_status", "target_start", "target_end",
          "n_distinct_genes", "genes", "qc_status", "reason"]]
""")

md(r"""
Three samples in, and the status column already carries three different meanings: `pass`,
`review` (no chain block) and `fail` (unresolved source build). Count them separately and
report all three. A summary that shows only the `pass` rows is not a summary, it is a
selection.

### 7.1 What skipping the conversion would have cost you

This is the punchline of the whole lecture, and it is worth running twice. We take the
checkpoint coordinate and ask the hg38 annotation about it **with** and **without** the
conversion.
""")

co(r"""
# (a) the single checkpoint base
converted = genes_overlapping("chr11", 5227001, 5227002)     # correct hg38 address
naive     = genes_overlapping("chr11", 5248231, 5248232)     # hg19 number used as-is

# (b) the whole legacy HBB transcript span, hg19 chr11:5,246,694-5,248,301
converted_span = genes_overlapping("chr11", 5246693 - 21230, 5248301 - 21230)
naive_span     = genes_overlapping("chr11", 5246693, 5248301)

print("single base, converted  chr11:5,227,002        ->", converted)
print("single base, un-converted                      ->", naive)
print()
print("HBB span,   converted   chr11:5,225,464-5,227,071 ->", converted_span)
print("HBB span,   un-converted chr11:5,246,694-5,248,301 ->", naive_span)
""")

md(r"""
Read both comparisons, because they fail in two different ways.

**The single base.** The un-converted coordinate returns an empty list. No error, no warning —
a perfectly well-formed `0`. Downstream it is indistinguishable from the legitimate result
"this site overlaps no annotated gene", and it would be written into a results table as a
fact. Slide 91 asked when zero is a justified count: *only when the check ran on the right
input.* This zero ran cleanly on the wrong input, which is worse than a missing value,
because `NA` at least admits it does not know.

**The gene span.** Here the shortcut is not merely empty, it is confidently wrong. The legacy
hg19 *HBB* transcript span, read as if it were hg38, overlaps **HBG1** — fetal haemoglobin
instead of adult haemoglobin, a different gene about 21 kb away, in the same cluster, on the
same strand. Every one of those properties makes the wrong answer *look* right. Nothing
downstream can rescue you.

That is the entire argument of today's lecture in two lines of output: convert deliberately,
then verify against a coordinate documented outside your own pipeline.

### Checkpoint 7 — keep this
`features.hg38.tsv`, the three-way status count, and the HBB-vs-HBG1 comparison above. This
is Student 5's evidence, and it is the strongest single slide anyone can show today.
""")

# ------------------------------------------------------------------ section 8
md(r"""
---
## 8 · Evidence card — *practice step 10*

Four minutes, one minute of questions, at most three slides. Fill the template below with
*your* assigned question (slide 82) and the checkpoint output you saved.

> *"I used **[operation]** because **[problem]**. I checked **[evidence]**. The result supports
> **[limited conclusion]**, while **[uncertainty]** remains."*

A justified stop is a valid result. "I stopped converting S03 because its header and manifest
disagree about the source build, and converting an unknown build would produce coordinates
nobody can audit" is a better answer than any number you could have invented.
""")

co(r"""
CARD = '''\
# Evidence card - {student}
**Assigned question:** {question}

## Problem
{problem}

## Action (what I ran, and why that operation)
{action}

## Check (the specific verification, not "it ran")
{check}

## Limitation (one, stated honestly)
{limitation}

---
Inputs: {inputs}
Outputs: {outputs}
Reference/annotation: hg19 + hg38 teaching mini-chr11; RefSeq Select (UCSC hg38, 2026-10)
'''

example = CARD.format(
    student="(your name)",
    question="Student 3 - hg19-to-hg38 conversion and rejects",
    problem=("The collaborator's S01 variants are on hg19; our annotation is hg38. "
             "chr11:5,248,232 means a different base in each assembly."),
    action=("Converted the 1-based site to BED on hg19, then mapped the interval with a "
            "chain file. Variants were re-checked against the hg38 sequence rather than "
            "having POS edited."),
    check=("The lifted checkpoint is chr11:5,227,002, which matches the independently "
           "documented ClinVar position (RCV005357130). 3 of 6 records converted; 3 were "
           "rejected with their source coordinates retained."),
    limitation=("One verified locus does not validate every interval. The insertion, the "
                "symbolic <DEL> and the breakend were never converted, and the teaching "
                "chain covers only chr11:5,246,500-5,280,000."),
    inputs="S01.hg19.vcf, regions.S01.hg19.tsv, hg19ToHg38.teaching.over.chain",
    outputs="small.hg38.vcf, rejected.hg19.vcf, mapped.hg38.bed, unmapped.hg19.bed")

(OUT / "evidence_card_example.md").write_text(example)
print(example)
""")

md(r"""
---
## Appendix · the sequence from the morning activity

The 63 bases everyone pasted into BLAST were the first 63 coding bases of *HBB*. Because
*HBB* sits on the **minus** strand, they are the reverse complement of the reference. The
cell below re-derives the coordinates you recorded on slide 25 — from the actual sequence,
not from the slide.
""")

co(r"""
HBB63 = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTG"
rcomp = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))

assert rcomp(hg38.fetch(5226959, 5227021)) == HBB63
assert rcomp(hg19.fetch(5248189, 5248251)) == HBB63
print("hg38 chr11:5,226,959-5,227,021  (minus strand) -> the 63-base query  OK")
print("hg19 chr11:5,248,189-5,248,251  (minus strand) -> the 63-base query  OK")
print("local shift:", 5248189 - 5226959, "bases")

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, b1 in enumerate(bases):
    for j, b2 in enumerate(bases):
        for k, b3 in enumerate(bases):
            CODON[b1 + b2 + b3] = aas[i * 16 + j * 4 + k]
protein = "".join(CODON[HBB63[i:i + 3]] for i in range(0, 63, 3))
print("translation:", protein)
assert protein == "MVHLTPEEKSAVTALWGKVNV"
print("codon 7 (c.19-21):", HBB63[18:21], "->", CODON[HBB63[18:21]],
      "| the checkpoint base c.20 is the middle one")
""")

md(r"""
And there is the loop closing. `c.20` — the base in the middle of codon 7 — is the exact
position we have been calling `HBB_checkpoint` all session: hg19 `chr11:5,248,232`, hg38
`chr11:5,227,002`. It is also one of the most studied single bases in human genetics (the
sickle-cell β-globin variant), which is precisely *why* its coordinates are documented in a
public database on both assemblies and therefore make a trustworthy checkpoint.

We used it for its coordinates and nothing else. The synthetic genotypes in S01–S03 are
teaching fiction, and no result in this notebook is a clinical statement about anybody.

---
### Before you close the notebook
Everything you produced is in `practice_session/outputs/`:

| file | what it proves |
|---|---|
| `compatibility_table.tsv` | which inputs you trusted, and why S03 stopped |
| `mapped.hg38.bed` / `unmapped.hg19.bed` | what converted and what did not |
| `small.hg38.vcf` / `rejected.hg19.vcf` | allele-aware conversion with visible rejects |
| `read_comparison.tsv` | reads move; CIGARs need not |
| `features.hg38.tsv` | annotation with `0` and `NA` kept distinct |
| `session_hg19.xml` / `session_hg38.xml` | one session per assembly |
| `evidence_card_example.md` | the shape of your 4-minute talk |
""")

nb["cells"] = cells
nb["metadata"] = {"kernelspec": {"display_name": "Python 3", "language": "python",
                                 "name": "python3"},
                  "language_info": {"name": "python"}}
nbf.write(nb, "practice_session/genomic_data_processing_practice.ipynb")
print("cells:", len(cells))
