"""Write practice.ipynb (everyone): set-up, Exercise 1, Practice A (A1-A5), Practice B core (B1, B2), defend.

Numbers in the text are filled in from build/rnaseq_build_summary.json (the provided DESeq2 files) and from
results/*.json of a previous run (pass 1), so the prose always matches what the notebook prints.
"""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from nbbuild import NB
import nb_part_a, nb_part_b

S = json.load(open(HERE / "rnaseq_build_summary.json"))
R = {}
for ex in ["A", "B1", "B2"]:
    p = ROOT / "results" / f"{ex}.json"
    if p.exists():
        R.update({k: v["value"] for k, v in json.load(open(p)).items() if isinstance(v, dict)})
g = lambda k, d="…": R.get(k, d)
T = {"TESTED_N": S["tested"], "UNTESTED_N": S["untested"], "UP_N": S["up"], "DOWN_N": S["down"],
     "SIG_N": S["sig_fdr05"], "SIG0_N": S["sig_fdr05_no_batch"],
     "TESTED": f"{S['tested']:,}", "UNTESTED": f"{S['untested']:,}", "UP": S["up"], "DOWN": S["down"],
     "SIG": S["sig_fdr05"], "SIG0": S["sig_fdr05_no_batch"],
     "PC1": f"{g('B2.pc1_pct', S['pc_var'][0]):.0f}", "PC2": f"{g('B2.pc2_pct', S['pc_var'][1]):.0f}",
     "LIBRATIO": g("B2.libsize_ratio_max_min"), "SFRANGE": str(g("B2.sf_range")).replace("-", "–"),
     "HANDTESTED": f"{g('B2.hand_tested', 0):,}", "EXPFALSE": f"{g('B2.expected_false_if_unadjusted', 0):,}",
     "RAWP05": f"{g('B2.raw_p05_with_batch', 0):,}", "RAWP05_NP": f"{g('B2.raw_p05_not_planted', 0):,}",
     "RAWP05_NP_PCT": f"{100 * g('B2.raw_p05_not_planted', 0) / max(g('B2.raw_p05_with_batch', 1), 1):.0f}",
     "HANDFDP_PCT": f"{100 * g('B2.hand_fdp_with_batch', 0):.1f}", "HANDGAIN": g("B2.hand_true_gain_from_batch"),
     "FDRONLY_PCT": g("B2.fdr_only_not_planted_pct"), "SFR": g("B2.sf_vs_libsize_r"),
     "VOLC_PLANTED_PCT": f"{100 * (1 - g('B2.volcano_genes_not_planted', 0) / (S['up'] + S['down'])):.1f}"}
nb = NB()
nb_part_a.add(nb)
nb_part_b.add(nb, T)
out = ROOT / "practice.ipynb"
nb.write(str(out))
left = sorted(set(re.findall(r"@@(\w+)@@", out.read_text())))
print(out.name, len(nb.cells), "cells; unfilled tokens:", left, "; values from previous run:", bool(R))
