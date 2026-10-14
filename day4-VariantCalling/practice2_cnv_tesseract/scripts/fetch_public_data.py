#!/usr/bin/env python
"""Download every PUBLIC input of Practice 2 and record its provenance.

Run from the package root:   python scripts/fetch_public_data.py [--stage reference|annotation|tracks|all]

reference  : real hg19 sequence for the three windows (UCSC getData/sequence), embedded in
             full-length N-masked chr4/chr17/chr22 -> data/ref/hg19_cnv_windows.fa.gz (+.fai/.gzi),
             data/ref/windows_hg19.bed and the UCSC chain file data/ref/hg19ToHg38.over.chain.gz
annotation : hg19 protein-coding genes in the windows (UCSC hg19 ncbiRefSeqSelect, NM_ only),
             hg19 segmental duplications + ClinVar CNVs in the windows,
             hg38 MANE Select genes genome-wide (UCSC hg38 mane, MANE v1.5)
tracks     : hg38 annotation tracks around each lifted planted CNV (see refresh_tracks_hg38)

Every download is logged (URL, UTC retrieval time, UCSC dataTime/downloadTime, data version
files) in data/provenance_public_data.json and data/tracks_hg38/provenance_tracks_hg38.json.
"""
import argparse
import datetime as dt
import gzip
import json
import os
import subprocess
import sys

import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C  # noqa: E402
import ucsc_api as U  # noqa: E402

CHAIN_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
VERSION_FILES = {
    "MANE (hg38 mane track)": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/mane/README_versions.txt",
    "ClinGen dosage map (hg38)": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/bbi/clinGen/clinGenDosageVersion.txt",
}
PROV_PATH = "data/provenance_public_data.json"


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def load_prov():
    if os.path.exists(PROV_PATH):
        return json.load(open(PROV_PATH))
    return {"description": "Provenance of every public file downloaded for Practice 2", "items": {}}


def save_prov(p):
    U.save_json(p, PROV_PATH)


# ----------------------------------------------------------------------------- reference
def stage_reference(prov):
    os.makedirs("data/ref", exist_ok=True)
    seqs, items = {}, []
    for chrom, s1, e1, name in C.WINDOWS_1BASED:
        dna, p = U.get_sequence("hg19", chrom, s1 - 1, e1)
        assert len(dna) == e1 - s1 + 1, (chrom, len(dna))
        p.update({"window_1based": f"{chrom}:{s1}-{e1}", "N_bases": dna.count("N")})
        items.append(p)
        seqs[chrom] = (s1 - 1, dna)
        print(f"  {chrom}:{s1:,}-{e1:,}  {len(dna):,} bp  N={dna.count('N'):,}")
    # full-length chromosomes, N everywhere except the windows
    fa_plain = "data/ref/hg19_cnv_windows.fa"
    with open(fa_plain, "w") as fh:
        for chrom in ["chr4", "chr17", "chr22"]:
            L = C.CHROM_LENGTHS_HG19[chrom]
            s0, dna = seqs[chrom]
            full = "N" * s0 + dna + "N" * (L - s0 - len(dna))
            assert len(full) == L
            fh.write(f">{chrom}\n")
            for i in range(0, L, 60):
                fh.write(full[i:i + 60] + "\n")
    subprocess.run(["bgzip", "-f", "-@", "4", "-i", fa_plain], check=True)
    subprocess.run(["samtools", "faidx", fa_plain + ".gz"], check=True)
    with open("data/ref/windows_hg19.bed", "w") as fh:
        for c, s0, e, n in C.windows_bed():
            fh.write(f"{c}\t{s0}\t{e}\t{n.replace(' ', '_')}\n")
    # chain file
    r = requests.get(CHAIN_URL, timeout=120)
    r.raise_for_status()
    open("data/ref/hg19ToHg38.over.chain.gz", "wb").write(r.content)
    prov["items"]["hg19_reference_windows"] = {
        "file": "data/ref/hg19_cnv_windows.fa.gz",
        "how": "UCSC REST API getData/sequence (hg19), uppercased, embedded in full-length N-masked chromosomes",
        "chrom_lengths_hg19": C.CHROM_LENGTHS_HG19, "downloads": items}
    prov["items"]["chain_hg19ToHg38"] = {
        "file": "data/ref/hg19ToHg38.over.chain.gz", "url": CHAIN_URL, "retrieved_utc": now(),
        "last_modified": r.headers.get("Last-Modified"), "bytes": len(r.content)}


# ----------------------------------------------------------------------------- annotation
def _genepred_rows(items, genome):
    rows = []
    for it in items:
        rows.append(dict(chrom=it["chrom"], start0=int(it["txStart"]), end=int(it["txEnd"]),
                         gene=it["name2"], strand=it["strand"], transcript=it["name"],
                         cds_start0=int(it["cdsStart"]), cds_end=int(it["cdsEnd"]),
                         exon_count=int(it["exonCount"])))
    return rows


def stage_annotation(prov):
    os.makedirs("data/annotation", exist_ok=True)
    os.makedirs("data/tracks_hg19", exist_ok=True)
    # hg19 genes in the windows: RefSeq Select + MANE subset, protein-coding (NM_) transcripts
    rows, provs = [], []
    for c, s0, e, n in C.windows_bed():
        items, p = U.get_track("hg19", "ncbiRefSeqSelect", c, s0, e)
        provs.append(p)
        rows += _genepred_rows(items, "hg19")
    g19 = pd.DataFrame(rows)
    g19 = g19[g19.transcript.str.startswith("NM_")].drop_duplicates(["gene"]).sort_values(["chrom", "start0"])
    g19.to_csv("data/annotation/genes_hg19_refseqselect_windows.tsv", sep="\t", index=False)
    prov["items"]["genes_hg19"] = {"file": "data/annotation/genes_hg19_refseqselect_windows.tsv",
                                   "track": "hg19 ncbiRefSeqSelect (NCBI RefSeq Select and MANE subset); NM_ (protein-coding) transcripts only; one row per gene",
                                   "downloads": provs, "n_genes": int(len(g19))}
    print(f"  hg19 genes in windows: {len(g19)}")
    # hg19 segmental duplications and ClinVar CNVs in the windows
    for track, out in [("genomicSuperDups", "data/tracks_hg19/genomicSuperDups_hg19_windows.tsv.gz"),
                       ("clinvarCnv", "data/tracks_hg19/clinvarCnv_hg19_windows.tsv.gz")]:
        allrows, provs = [], []
        for c, s0, e, n in C.windows_bed():
            items, p = U.get_track("hg19", track, c, s0, e)
            provs.append(p)
            allrows += items
        df = pd.DataFrame(allrows)
        df.to_csv(out, sep="\t", index=False)
        prov["items"][f"{track}_hg19"] = {"file": out, "track": f"hg19 {track}", "downloads": provs, "n_items": int(len(df))}
        print(f"  hg19 {track}: {len(df)} items")
    # hg38 MANE genome-wide (one MANE Select transcript per gene)
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    rows, provs = [], []
    for c in chroms:
        items, p = U.get_track("hg38", "mane", c)
        provs.append(p)
        for it in items:
            rows.append(dict(chrom=it["chrom"], start0=int(it["chromStart"]), end=int(it["chromEnd"]),
                             gene=it["geneName2"], strand=it["strand"], transcript=it["name"],
                             refseq=it.get("ncbiId", ""), gene_id=it.get("geneName", ""),
                             gene_type=it.get("geneType", ""), mane_status=it.get("maneStat", "")))
    mane = pd.DataFrame(rows)
    mane = mane[(mane.mane_status == "MANE Select") & (mane.gene_type == "protein_coding")]
    mane = mane.drop_duplicates("gene").sort_values(["chrom", "start0"]).reset_index(drop=True)
    mane.to_csv("data/annotation/genes_hg38_mane_select.tsv.gz", sep="\t", index=False)
    versions = {}
    for k, url in VERSION_FILES.items():
        try:
            versions[k] = {"url": url, "text": requests.get(url, timeout=60).text.strip()}
        except Exception as e:
            versions[k] = {"url": url, "error": repr(e)}
    prov["items"]["genes_hg38_mane"] = {"file": "data/annotation/genes_hg38_mane_select.tsv.gz",
                                        "track": "hg38 mane (MANE Select, protein_coding only; MANE Plus Clinical rows dropped)",
                                        "version_files": versions, "downloads": provs, "n_genes": int(len(mane))}
    print(f"  hg38 MANE Select protein-coding genes: {len(mane)}")


# ----------------------------------------------------------------------------- tracks
HG38_TRACKS = [
    # (track, why, filter)
    ("mane", "genes (MANE Select/Plus Clinical transcripts)", None),
    ("clinGenHaplo", "ClinGen haploinsufficiency (HI) scores, genes and regions", None),
    ("clinGenTriplo", "ClinGen triplosensitivity (TS) scores, genes and regions", None),
    ("clinGenGeneDisease", "ClinGen gene-disease validity", None),
    ("pliByGene", "gnomAD v2.1.1 LoF constraint (pLI, LOEUF) by gene", None),
    ("clinvarCnv", "ClinVar CNVs >= 50 bp", None),
    ("clinvarMain", "ClinVar SNVs/indels < 50 bp, kept only Pathogenic / Likely pathogenic", "plp"),
    ("phastConsElements100way", "conserved elements (100 vertebrates)", None),
    ("encodeCcreCombined", "ENCODE3 candidate cis-regulatory elements", None),
    ("genomicSuperDups", "segmental duplications (>1 kb, >90% identity)", None),
    ("gnomadStructuralVariants", "gnomAD v4.1 SVs (population frequency; NOT a TESSERACT input)", None),
    ("dgvGold", "DGV Gold Standard CNVs (population; NOT a TESSERACT input)", None),
]
SKIPPED_TRACKS = {"decipher": "DECIPHER CNVs: protectedData=True in /list/tracks (restricted; cannot be downloaded via the API)"}


def _is_plp(sig):
    s = str(sig).lower()
    return ("pathogenic" in s) and ("conflicting" not in s) and ("benign" not in s) and ("uncertain" not in s)


def refresh_tracks_hg38(regions, outdir="data/tracks_hg38", verify_names=True):
    """Live (re)download of the hg38 tracks for a list of regions.

    regions: list of dicts with keys region_id, chrom, start0, end (hg38, BED convention).
    Writes one gzipped TSV per track (all regions stacked, column region_id added) and
    outdir/provenance_tracks_hg38.json. Returns the provenance dict.
    """
    os.makedirs(outdir, exist_ok=True)
    prov = {"retrieved_utc": now(), "regions": regions, "tracks": {}, "skipped": SKIPPED_TRACKS}
    if verify_names:
        tracks, p = U.list_tracks("hg38")
        prov["track_list"] = p
        for t, _, _ in HG38_TRACKS:
            if t not in tracks:
                raise RuntimeError(f"track {t} not found in hg38 /list/tracks")
            if tracks[t].get("protectedData"):
                raise RuntimeError(f"track {t} is restricted")
        prov["track_meta"] = {t: {k: tracks[t].get(k) for k in ("shortLabel", "longLabel", "type", "dataVersion", "bigDataUrl")}
                              for t, _, _ in HG38_TRACKS}
    for track, why, filt in HG38_TRACKS:
        rows, provs = [], []
        for reg in regions:
            items, p = U.get_track("hg38", track, reg["chrom"], reg["start0"], reg["end"])
            p["region_id"] = reg["region_id"]
            provs.append(p)
            for it in items:
                it = {k: v for k, v in it.items() if not k.startswith("_mouseOver") and k not in ("_jsonHgvsTable",)}
                it["region_id"] = reg["region_id"]
                rows.append(it)
        df = pd.DataFrame(rows)
        n_raw = len(df)
        if filt == "plp" and len(df):
            df = df[df["clinSign"].map(_is_plp)]
        out = os.path.join(outdir, f"{track}.tsv.gz")
        df.to_csv(out, sep="\t", index=False)
        prov["tracks"][track] = {"why": why, "file": out, "n_items_downloaded": int(n_raw), "n_items_kept": int(len(df)),
                                 "filter": "clinSign Pathogenic/Likely pathogenic only (no conflicting/benign/uncertain)" if filt else None,
                                 "downloads": provs}
        print(f"  {track:26s} kept {len(df):6d} / {n_raw:6d}")
    U.save_json(prov, os.path.join(outdir, "provenance_tracks_hg38.json"))
    return prov


def lift_points(points, chain="data/ref/hg19ToHg38.over.chain.gz", tmp="data/tracks_hg38/_tmp_points"):
    """Lift single base positions (chrom, pos0) with CrossMap; returns list of (chrom, pos0) or None."""
    with open(tmp + ".bed", "w") as fh:
        for i, (c, p0) in enumerate(points):
            fh.write(f"{c}\t{p0}\t{p0 + 1}\tp{i}\n")
    subprocess.run(["CrossMap", "bed", chain, tmp + ".bed", tmp + ".out.bed"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = [None] * len(points)
    for line in open(tmp + ".out.bed"):
        c, s, e, name = line.split("\t")[:4]
        res[int(name.strip()[1:])] = (c, int(s))
    for f in (tmp + ".bed", tmp + ".out.bed", tmp + ".out.bed.unmap"):
        if os.path.exists(f):
            os.remove(f)
    return res


def track_regions_from_planted(flank_frac=0.20, min_flank=50_000):
    """hg38 regions = breakpoint-lifted planted CNV (union with the public GRCh38 outer span
    for the 22q11.2 proximal deletion) +/- max(20% of size, 50 kb)."""
    uniq = {}
    for ev in C.PLANTED:
        key = (ev["chrom"], ev["start_1based"], ev["end_1based"])
        uniq.setdefault(key, ev)
    regions = []
    for (chrom, s1, e1), ev in uniq.items():
        a, b = lift_points([(chrom, s1 - 1), (chrom, e1 - 1)])
        assert a and b and a[0] == b[0] == chrom, (ev["cnv_id"], a, b)
        s0, e = min(a[1], b[1]), max(a[1], b[1]) + 1
        if ev["cnv_id"].startswith("P1_"):
            s0 = min(s0, 18_339_129)   # ClinVar GRCh38 outer span chr22:18,339,130-21,111,373
            e = max(e, 21_111_373)
        fl = max(int(flank_frac * (e - s0)), min_flank)
        name = {"P1_22q11_LCR22A-D_DEL": "22q11_proximal", "P2_PMP22_DUP": "PMP22",
                "P3_22q11_distal_DEL": "22q11_distal"}.get(ev["cnv_id"], "UGT2B17")
        regions.append(dict(region_id=name, chrom=chrom, start0=max(0, s0 - fl), end=e + fl,
                            core_start0=s0, core_end=e, flank_bp=fl))
    return regions


def stage_tracks(prov):
    regions = track_regions_from_planted()
    for r in regions:
        print(f"  region {r['region_id']:15s} {r['chrom']}:{r['start0']:,}-{r['end']:,}")
    p = refresh_tracks_hg38(regions)
    prov["items"]["tracks_hg38"] = {"dir": "data/tracks_hg38", "provenance": "data/tracks_hg38/provenance_tracks_hg38.json",
                                    "regions": regions, "skipped": SKIPPED_TRACKS}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["reference", "annotation", "tracks", "all"])
    a = ap.parse_args()
    prov = load_prov()
    if a.stage in ("reference", "all"):
        print("[reference]"); stage_reference(prov); save_prov(prov)
    if a.stage in ("annotation", "all"):
        print("[annotation]"); stage_annotation(prov); save_prov(prov)
    if a.stage in ("tracks", "all"):
        print("[tracks]"); stage_tracks(prov); save_prov(prov)
