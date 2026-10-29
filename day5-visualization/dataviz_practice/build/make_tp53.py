"""B8 data: TP53 somatic variants in 1,200 tumours (SIMULATED counts on the real TP53 protein).

Real anchors: the TP53 protein (UniProt P04637, 393 amino acids; sequence embedded below and checked against
every named hotspot), domain boundaries in the IARC TP53 Database convention, and the real recurrent hotspot
residues (missense R175H, R248Q/W, R273H/C, ...; nonsense R196*, R213*, R306*, R342*, Q192*, E298*). Counts are
fixed per class (876 missense = 73 %, 264 truncating = 22 %, 60 splice/other = 5 %) and drawn at random within
each class (multinomial): missense concentrated at hotspots in the DNA-binding domain, truncating = recurrent
nonsense + frameshifts spread along the protein.

Writes (paths relative to dataviz_practice/)
  data/tp53_variants.tsv        tumour_id, protein_change, position (1-based residue), variant_class
  data/tp53_domains.tsv         domain, short, aa_start, aa_end (1-based, inclusive), length_aa, source
  data/truth/tp53_truth.tsv     every variant in the simulation catalogue: expected and realised counts
Seed: np.random.default_rng(20261027 + 9). Running the script twice gives byte-identical files.
"""
from __future__ import annotations

from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
rng = np.random.default_rng(20261027 + 9)
N_TUMOURS = 1200

# UniProt P04637 (cellular tumor antigen p53, canonical isoform), 393 aa
TP53 = ("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQ"
        "KTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGN"
        "LRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP"
        "PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
L = len(TP53)
assert L == 393

DOMAINS = pd.DataFrame([
    ("Transactivation", "TAD", 1, 61),
    ("Proline-rich", "PRD", 64, 92),
    ("DNA-binding", "DBD", 94, 292),
    ("Tetramerisation", "TET", 325, 356)], columns=["domain", "short", "aa_start", "aa_end"])
DOMAINS["length_aa"] = DOMAINS["aa_end"] - DOMAINS["aa_start"] + 1
DOMAINS["source"] = "IARC TP53 Database convention (other databases differ by a few residues)"
DBD = (94, 292)

# Standard genetic code -> which amino acids one nucleotide change can reach from each residue
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {"".join(c): a for c, a in zip(product(BASES, repeat=3), AA)}
REACH: dict[str, set] = {}
for codon, aa in CODE.items():
    for i, b in product(range(3), BASES):
        if b != codon[i]:
            REACH.setdefault(aa, set()).add(CODE[codon[:i] + b + codon[i + 1:]])
MISSENSE_ALT = {aa: sorted(r - {aa, "*"}) for aa, r in REACH.items() if aa != "*"}
STOP_REACHABLE = {aa for aa, r in REACH.items() if "*" in r and aa != "*"}

# expected numbers of tumours (relative weights; normalised below)
HOTSPOTS = {  # missense
    "R175H": 62, "R248Q": 38, "R248W": 32, "R273H": 34, "R273C": 30, "R282W": 32, "R282Q": 3, "G245S": 22,
    "G245D": 10, "G245C": 3, "Y220C": 22, "R249S": 18, "V157F": 12, "H179R": 10,
    "R158H": 7, "R158L": 5, "C176F": 6, "C176Y": 4, "H179Y": 4, "H193R": 6, "I195T": 7, "Y163C": 6, "M237I": 6,
    "C238Y": 4, "S241F": 5, "C242F": 4, "R280K": 5, "E285K": 6, "V272M": 6, "C275Y": 5, "P278S": 5, "P278L": 4,
    "Y234C": 5, "Y205C": 4, "H214R": 4, "V216M": 5, "K132N": 4, "C135Y": 4, "C141Y": 5, "V143A": 3, "P151S": 4,
    "A159V": 3, "R110L": 3, "N239D": 3, "G266R": 3, "R267W": 3, "E271K": 3, "D281G": 3, "R283P": 3, "R337C": 3,
    "R337H": 2, "L344P": 2, "R213Q": 3, "E286K": 3, "G244D": 3, "G244S": 3, "R273L": 3, "R273P": 2, "R248L": 2}
NONSENSE = {"R213*": 22, "R196*": 15, "R306*": 12, "R342*": 12, "Q192*": 7, "E298*": 6}
SPLICE = {125: 8, 126: 4, 187: 8, 224: 9, 225: 4, 261: 6, 262: 3, 307: 6, 332: 4, 33: 2}  # codon next to the site
BACKGROUND = {"missense_in_dbd": 274, "missense_outside": 91, "nonsense": 54, "frameshift": 132, "inframe_del": 10}

rows = []


def add(change, pos, cls, w):
    rows.append((change, int(pos), cls, float(w)))


for ch, w in HOTSPOTS.items():
    ref, pos, alt = ch[0], int(ch[1:-1]), ch[-1]
    assert TP53[pos - 1] == ref, f"{ch}: residue {pos} is {TP53[pos - 1]}"
    assert alt in MISSENSE_ALT[ref], f"{ch} needs more than one nucleotide change"
    add(ch, pos, "missense", w)
for ch, w in NONSENSE.items():
    ref, pos = ch[0], int(ch[1:-1])
    assert TP53[pos - 1] == ref and ref in STOP_REACHABLE, ch
    add(ch, pos, "truncating", w)
for pos, w in SPLICE.items():
    add(f"X{pos}_splice", pos, "splice/other", w)

# background ("long tail"): one random allele per residue, weight spread evenly over the residues
dbd_pos = np.arange(DBD[0], DBD[1] + 1)
out_pos = np.array([p for p in range(2, L + 1) if not DBD[0] <= p <= DBD[1]])      # residue 1 (Met) excluded
for pos_set, total in [(dbd_pos, BACKGROUND["missense_in_dbd"]), (out_pos, BACKGROUND["missense_outside"])]:
    for p in pos_set:
        ref = TP53[p - 1]
        add(f"{ref}{p}{rng.choice(MISSENSE_ALT[ref])}", p, "missense", total / len(pos_set))
stop_pos = np.array([p for p in range(2, L + 1) if TP53[p - 1] in STOP_REACHABLE])
for p in stop_pos:
    add(f"{TP53[p - 1]}{p}*", p, "truncating", BACKGROUND["nonsense"] / len(stop_pos))
for p in range(1, L + 1):
    add(f"{TP53[p - 1]}{p}fs", p, "truncating", BACKGROUND["frameshift"] / L)
for p in rng.choice(np.arange(2, L + 1), 10, replace=False):
    add(f"{TP53[p - 1]}{p}del", p, "splice/other", BACKGROUND["inframe_del"] / 10)

cat = (pd.DataFrame(rows, columns=["protein_change", "position", "variant_class", "weight"])
       .groupby(["protein_change", "position", "variant_class"], as_index=False)["weight"].sum())
# class totals fixed to the target spectrum (73 / 22 / 5 %); counts WITHIN a class are drawn at random
CLASS_TOTALS = {"missense": 876, "truncating": 264, "splice/other": 60}
assert sum(CLASS_TOTALS.values()) == N_TUMOURS
cat["expected_n"] = 0.0
cat["n_tumours"] = 0
for cls_name, n_cls in CLASS_TOTALS.items():
    m = (cat["variant_class"] == cls_name).to_numpy()
    w = cat.loc[m, "weight"].to_numpy()
    cat.loc[m, "expected_n"] = n_cls * w / w.sum()
    cat.loc[m, "n_tumours"] = rng.multinomial(n_cls, w / w.sum())

# one variant per tumour, tumours in random order
draws = np.repeat(np.arange(len(cat)), cat["n_tumours"].to_numpy())
rng.shuffle(draws)
variants = cat.iloc[draws][["protein_change", "position", "variant_class"]].reset_index(drop=True)
variants.insert(0, "tumour_id", [f"T{i:04d}" for i in range(1, N_TUMOURS + 1)])

truth = (cat.assign(expected_n=cat["expected_n"].round(3))
         .sort_values(["position", "variant_class", "protein_change"])
         [["protein_change", "position", "variant_class", "expected_n", "n_tumours"]])

(DATA / "truth").mkdir(parents=True, exist_ok=True)
kw = dict(sep="\t", index=False, na_rep="", lineterminator="\n")
variants.to_csv(DATA / "tp53_variants.tsv", **kw)
DOMAINS.to_csv(DATA / "tp53_domains.tsv", **kw)
truth.to_csv(DATA / "truth" / "tp53_truth.tsv", float_format="%.6g", **kw)

if __name__ == "__main__":
    cls = variants["variant_class"].value_counts()
    print((cls / N_TUMOURS * 100).round(1).to_dict(), cls.to_dict())
    inside = variants["position"].between(*DBD)
    for c in ["missense", "truncating"]:
        m = variants["variant_class"] == c
        print(f"{c}: {int((m & inside).sum())}/{int(m.sum())} in DBD = {(m & inside).sum() / m.sum():.3f} "
              f"(length fraction {(DBD[1] - DBD[0] + 1) / L:.3f})")
    per_pos = variants.groupby("position").size().sort_values(ascending=False)
    print("positions >= 30:", per_pos[per_pos >= 30].to_dict())
    print(variants.groupby("protein_change").size().sort_values(ascending=False).head(12).to_dict())
