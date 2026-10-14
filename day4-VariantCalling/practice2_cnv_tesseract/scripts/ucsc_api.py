"""Small helper for the UCSC Genome Browser REST API (https://api.genome.ucsc.edu).

Used by fetch_public_data.py and by refresh_tracks_hg38() to re-download the hg38
annotation tracks live. Only public tracks are used; restricted tracks (for example
DECIPHER, which is flagged protectedData in /list/tracks) cannot be downloaded through
the API and are skipped on purpose.
"""
import datetime as _dt
import json
import time

import requests

API = "https://api.genome.ucsc.edu"
UA_TIMEOUT = 120


def _get(url, params, tries=4):
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, params=params, timeout=UA_TIMEOUT)
            if r.status_code == 200:
                return r.json()
            last = f"HTTP {r.status_code}: {r.text[:200]}"
        except Exception as e:  # network hiccup -> retry with back-off
            last = repr(e)
        time.sleep(2 * (i + 1))
    raise RuntimeError(f"UCSC API request failed: {url} {params} -> {last}")


def _query_string(params):
    return ";".join(f"{k}={v}" for k, v in params.items())


def get_sequence(genome, chrom, start0, end):
    """Return (uppercase DNA, provenance dict) for a 0-based half-open interval."""
    params = {"genome": genome, "chrom": chrom, "start": start0, "end": end}
    d = _get(f"{API}/getData/sequence", params)
    prov = {"url": f"{API}/getData/sequence?{_query_string(params)}",
            "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "downloadTime": d.get("downloadTime"), "dataTime": d.get("dataTime"),
            "length": len(d["dna"])}
    return d["dna"].upper(), prov


def get_track(genome, track, chrom, start0=None, end=None, max_items=-1):
    """Return (list of item dicts, provenance dict) for one track on one chromosome/region."""
    params = {"genome": genome, "track": track, "chrom": chrom}
    if start0 is not None:
        params.update({"start": int(start0), "end": int(end)})
    params["maxItemsOutput"] = max_items
    d = _get(f"{API}/getData/track", params)
    items = d.get(track, [])
    if isinstance(items, dict):  # some table tracks return {chrom: [...]} when no region is given
        items = [x for v in items.values() for x in v]
    prov = {"genome": genome, "track": track,
            "url": f"{API}/getData/track?{_query_string(params)}",
            "retrieved_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "downloadTime": d.get("downloadTime"), "dataTime": d.get("dataTime"),
            "bigDataUrl": d.get("bigDataUrl"), "itemsReturned": d.get("itemsReturned", len(items)),
            "maxItemsLimit": d.get("maxItemsLimit", False)}
    return items, prov


def list_tracks(genome):
    d = _get(f"{API}/list/tracks", {"genome": genome, "trackLeavesOnly": 1})
    return d[genome], {"url": f"{API}/list/tracks?genome={genome};trackLeavesOnly=1",
                       "dataTime": d.get("dataTime"), "downloadTime": d.get("downloadTime")}


def save_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
