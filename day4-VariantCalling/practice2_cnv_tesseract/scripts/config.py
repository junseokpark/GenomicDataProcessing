"""Shared constants for the Practice 2 data generator (all coordinates hg19/GRCh37 unless noted).

Coordinates in this file follow the convention given next to each value:
  *_1based  -> 1-based, inclusive (as written in ClinVar/LOVD records and on slides)
  start0    -> 0-based, half-open (BED convention, used by every file we write)
"""
SEED = 20261015

# hg19 chromosome lengths (UCSC hg19 chromInfo)
CHROM_LENGTHS_HG19 = {"chr4": 191154276, "chr17": 81195210, "chr22": 51304566}

# Real hg19 sequence windows (1-based inclusive) that are copied from UCSC into the
# otherwise N-masked practice reference.
WINDOWS_1BASED = [
    ("chr22", 18_000_001, 23_000_000, "22q11.2 (LCR22A-F)"),
    ("chr17", 14_600_001, 15_700_000, "17p12 (PMP22 / CMT1A-REP)"),
    ("chr4", 69_000_001, 69_800_000, "4q13.2 (UGT2B17)"),
]

BIN_SIZE = 5_000          # read-depth bin size used by the caller (bp)
READ_LEN = 100            # paired-end 2 x 100 bp
FRAG_MEAN, FRAG_SD = 350, 50
FRAG_MIN, FRAG_MAX = 200, 600
DUP_RATE = 0.03           # fraction of fragments that get one PCR/optical duplicate
MEAN_DEPTH = 3.0          # target mean read depth (x) over the real-sequence windows (before duplicates)
# Lowest depth at which the Step-3 caller recovered the 35.6-kb PMP22 duplication in every simulated
# replicate (see INSTRUCTOR_NOTES.md: 1.5x 4/8, 2.0x 3/8, 2.5x 6/8, 3.0x 8/8 replicates).
# Per-sample multipliers make library sizes differ, so library-size normalisation matters.
DEPTH_FACTOR = {"P1": 1.08, "P2": 0.94, "P3": 1.12, "C1": 1.00, "C2": 0.90, "C3": 1.15, "C4": 0.88, "C5": 1.05}
QUAL_BINS = (2, 12, 23, 37)  # binned base qualities (Illumina-style 4-level binning)

# Planted CNVs. Coordinates: 1-based inclusive (as in the public records).
PLANTED = [
    dict(cnv_id="P1_22q11_LCR22A-D_DEL", sample="P1", chrom="chr22", start_1based=18_916_842,
         end_1based=21_465_662, cn=1, zygosity="heterozygous", type="DEL",
         source="ClinVar VCV000155158 (GRCh37 record)",
         public_grch38="chr22:18339130-21111373 (ClinVar GRCh38 record)",
         note="Proximal 22q11.2 deletion (LCR22A-D), the classic 22q11.2 deletion syndrome region."),
    dict(cnv_id="P2_PMP22_DUP", sample="P2", chrom="chr17", start_1based=15_133_094,
         end_1based=15_168_674, cn=3, zygosity="heterozygous", type="DUP",
         source="LOVD PMP22_000101",
         public_grch38="chr17:15229777-15265357 (LOVD GRCh38)",
         note="Whole-gene PMP22 tandem duplication (planted as a tandem copy)."),
    dict(cnv_id="P2_UGT2B17_DEL", sample="P2", chrom="chr4", start_1based=69_375_169,
         end_1based=69_483_248, cn=1, zygosity="heterozygous", type="DEL",
         source="ClinVar VCV000147015 (copy number loss, Benign); overlaps DGV Gold losses (e.g. gssvL90945, 12-16% frequency)",
         public_grch38="see liftover step",
         note="Common, benign whole-gene deletion polymorphism."),
    dict(cnv_id="P3_22q11_distal_DEL", sample="P3", chrom="chr22", start_1based=21_029_655,
         end_1based=22_481_498, cn=1, zygosity="heterozygous", type="DEL",
         source="Lecture 'case 3' (distal 22q11.2 deletion); matching ClinVar record looked up in UCSC hg19 clinvarCnv (see truth TSV)",
         public_grch38="see liftover step",
         note="Distal 22q11.2 deletion; ClinVar classifications for similar records conflict (VUS/pathogenic)."),
    dict(cnv_id="P3_UGT2B17_DEL_HOM", sample="P3", chrom="chr4", start_1based=69_375_169,
         end_1based=69_483_248, cn=0, zygosity="homozygous", type="DEL",
         source="ClinVar VCV000147015 interval (common deletion), homozygous",
         public_grch38="see liftover step",
         note="Both copies deleted (CN 0) - frequent in the population."),
    dict(cnv_id="C2_UGT2B17_DEL", sample="C2", chrom="chr4", start_1based=69_375_169,
         end_1based=69_483_248, cn=1, zygosity="heterozygous", type="DEL",
         source="ClinVar VCV000147015 interval (common deletion)",
         public_grch38="see liftover step",
         note="Population variation in a control: controls are not CNV-free."),
]

# Fictional samples (no real patient data). Sex only matters for chrY genes in the RNA matrix.
SAMPLES = [
    dict(sample="P1", role="patient", sex="M", age=6, phenotype="congenital heart defect; speech delay (fictional)"),
    dict(sample="P2", role="patient", sex="F", age=32, phenotype="slowly progressive distal weakness; slowed nerve conduction (fictional)"),
    dict(sample="P3", role="patient", sex="M", age=9, phenotype="mild learning difficulties; parents untested (fictional)"),
    dict(sample="C1", role="control", sex="F", age=35, phenotype="unaffected control (fictional)"),
    dict(sample="C2", role="control", sex="M", age=41, phenotype="unaffected control (fictional)"),
    dict(sample="C3", role="control", sex="F", age=29, phenotype="unaffected control (fictional)"),
    dict(sample="C4", role="control", sex="M", age=52, phenotype="unaffected control (fictional)"),
    dict(sample="C5", role="control", sex="F", age=38, phenotype="unaffected control (fictional)"),
]
SHIPPED_FASTQ = ["P1", "P2", "P3", "C1"]
CONTROLS = ["C1", "C2", "C3", "C4", "C5"]


def windows_bed():
    """Windows as BED rows (chrom, start0, end, name)."""
    return [(c, s - 1, e, n) for c, s, e, n in WINDOWS_1BASED]
