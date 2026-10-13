#!/usr/bin/env python3
"""Execute both notebooks headless from a clean copy and check them against checkpoints/.

What is checked
  1. static: notebooks ship without outputs; code cells import only gzip, pathlib, numpy,
     pandas and matplotlib; no network calls anywhere in the notebooks or build scripts
  2. both notebooks run top to bottom with nbclient in a fresh copy of the package, with no
     errors and nothing written to stderr (no warnings), in under 60 s together
  3. every output table equals its checkpoint (text equal, or numerically equal to 1e-9)
  4. the teaching facts hold: Part A tally 6/2/3/1, Illumina spans 0 reads of the target,
     HiFi and ONT >= 3 each, the QC funnel adds up, every planted case lands where documented
  5. every figure was written

Exit status 0 when everything passes, 1 otherwise.

Usage
  python build/verify_package.py [--keep] [--kernel python3] [--executor auto|nbclient|ipython-shell]
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import nbformat
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TRANSPORT = None
NOTEBOOKS = ["practice1_sample_to_platform.ipynb", "practice2_cohort_qc.ipynb"]
TABLES = ["practice1_decisions.tsv", "practice1_read_stats.tsv",
          "practice2_qc_funnel.tsv", "practice2_excluded_samples.tsv"]
FIGURES = ["read_length_hist.png", "target_spanning_reads.png", "pca.png", "batch_ccds.png"]
ALLOWED_IMPORTS = {"gzip", "pathlib", "numpy", "pandas", "matplotlib"}
NETWORK = re.compile(r"\b(requests|urllib|urlopen|http\.client|socket|ftplib|wget|curl|"
                     r"pip install|boto3|paramiko)\b")
BUDGET_S = 60.0
EXPECTED_P1 = {"S01": "PASS", "S02": "REVIEW", "S03": "FAIL", "S04": "PASS", "S05": "PASS",
               "S06": "FAIL", "S07": "PASS", "S08": "PASS", "S09": "PASS", "S10": "FAIL",
               "S11": "REVIEW", "S12": "BLOCKED"}

problems = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        problems.append(msg)


def static_checks():
    for name in NOTEBOOKS:
        nb = nbformat.read(ROOT / name, as_version=4)
        code = [c for c in nb.cells if c.cell_type == "code"]
        check(all(not c.outputs and c.execution_count is None for c in code),
              f"{name}: shipped without outputs")
        mods = set()
        for c in code:
            for m in re.finditer(r"^\s*(?:import|from)\s+([A-Za-z_]\w*)", c.source, re.M):
                mods.add(m.group(1))
        check(mods <= ALLOWED_IMPORTS, f"{name}: imports {sorted(mods)} within {sorted(ALLOWED_IMPORTS)}")
        check(not any(NETWORK.search(c.source) for c in code), f"{name}: no network calls in code")
    for p in sorted((ROOT / "build").glob("*.py")):
        if p.name == Path(__file__).name:
            continue
        check(not NETWORK.search(p.read_text()), f"build/{p.name}: no network calls")


def clean_copy(dest):
    shutil.copytree(ROOT, dest, ignore=shutil.ignore_patterns(
        "__pycache__", ".ipynb_checkpoints", ".git", ".tmp", "*.tar.gz", "review_*.txt", "verify_log.txt"))
    out = dest / "outputs"
    for p in out.iterdir():
        if p.name != ".gitkeep":
            p.unlink() if p.is_file() else shutil.rmtree(p)
    return dest


def kernel_transport():
    """Which ZMQ transport can a Jupyter kernel use here? None = no kernel possible."""
    try:
        import zmq
        import nbclient  # noqa: F401  (the executor itself must be installed too)
        import ipykernel  # noqa: F401
    except ImportError:
        return None
    ctx = zmq.Context.instance()
    for transport in ("tcp", "ipc"):
        sock = ctx.socket(zmq.ROUTER)
        try:
            if transport == "tcp":
                sock.bind_to_random_port("tcp://127.0.0.1")
            else:
                sock.bind("ipc://" + str(Path(tempfile.mkdtemp(prefix="nbk", dir="/tmp")) / "probe"))
            return transport
        except (zmq.ZMQError, OSError):
            continue
        finally:
            sock.close(linger=0)
    return None


def run_with_nbclient(path, kernel, transport):
    from jupyter_client.manager import KernelManager
    from nbclient import NotebookClient
    from nbclient.exceptions import CellExecutionError
    nb = nbformat.read(path, as_version=4)
    km = KernelManager(kernel_name=kernel, transport=transport)
    if transport == "ipc":                       # short socket path: Unix sockets cap ~100 chars
        km.ip = str(Path(tempfile.mkdtemp(prefix="nbk", dir="/tmp")) / "k")
    error = None
    try:
        NotebookClient(nb, km=km, timeout=120, startup_timeout=30, kernel_name=kernel,
                       resources={"metadata": {"path": str(path.parent)}}).execute()
    except CellExecutionError as exc:
        error = (str(exc).strip().splitlines() or [repr(exc)])[-1]
    outputs = []
    for i, cell in enumerate(nb.cells):
        for o in cell.get("outputs", []):
            stream = o.get("name") if o.get("output_type") == "stream" else o.get("output_type")
            outputs.append((i, stream, o.get("text", "")))
    return error, outputs


SHELL_RUNNER = r"""
import contextlib, io, json, os, sys, warnings
import nbformat
from IPython.core.interactiveshell import InteractiveShell
path = os.path.abspath(sys.argv[1])
os.chdir(os.path.dirname(path))
nb = nbformat.read(path, as_version=4)
shell = InteractiveShell.instance()
report = {"error": None, "outputs": []}
for i, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err), \
            warnings.catch_warnings(record=True) as caught:
        result = shell.run_cell(cell.source, store_history=True)
    text_err = err.getvalue() + "".join(f"{w.category.__name__}: {w.message}\n" for w in caught)
    report["outputs"].append([i, "stdout", out.getvalue()])
    if text_err:
        report["outputs"].append([i, "stderr", text_err])
    if not result.success:
        exc = result.error_in_exec or result.error_before_exec
        report["error"] = f"cell {i}: {type(exc).__name__}: {exc}"
        break
print("@@REPORT@@" + json.dumps(report))
"""


def run_with_shell(path):
    """Fallback for machines that forbid kernel sockets: run the same cells, in order, in a
    fresh interpreter through IPython's InteractiveShell (the engine inside a Jupyter kernel)."""
    scratch = Path(tempfile.gettempdir())
    env = {**os.environ, "MPLBACKEND": "module://matplotlib_inline.backend_inline",
           "MPLCONFIGDIR": os.environ.get("MPLCONFIGDIR", str(scratch / "day2-mplconfig")),
           "IPYTHONDIR": os.environ.get("IPYTHONDIR", str(scratch / "day2-ipython"))}
    proc = subprocess.run([sys.executable, "-c", SHELL_RUNNER, str(path)], capture_output=True,
                          text=True, timeout=300, env=env)
    marker = proc.stdout.rfind("@@REPORT@@")
    if marker < 0:
        return f"runner crashed: {proc.stderr.strip()[-300:]}", []
    report = json.loads(proc.stdout[marker + len("@@REPORT@@"):])
    return report["error"], [tuple(o) for o in report["outputs"]]


def run_notebook(path, kernel, executor):
    t0 = time.perf_counter()
    if executor == "nbclient":
        error, outputs = run_with_nbclient(path, kernel, TRANSPORT)
    else:
        error, outputs = run_with_shell(path)
    seconds = time.perf_counter() - t0
    stderr = [f"cell {i}: {text.strip()[:200]}" for i, kind, text in outputs
              if kind == "stderr" and text.strip()]
    mismatch = [f"cell {i}" for i, kind, text in outputs if "matches checkpoint: False" in text]
    return seconds, error, stderr, mismatch


def same_table(a, b):
    x = pd.read_csv(a, sep="\t", dtype=str, keep_default_na=False)
    y = pd.read_csv(b, sep="\t", dtype=str, keep_default_na=False)
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        return False, f"shape/columns differ: {x.shape} {list(x.columns)} vs {y.shape} {list(y.columns)}"
    for col in x.columns:
        for i, (u, v) in enumerate(zip(x[col], y[col])):
            if u == v:
                continue
            try:
                if np.isclose(float(u), float(v), rtol=1e-9, atol=1e-9):
                    continue
            except ValueError:
                pass
            return False, f"row {i}, column {col}: {u!r} vs {v!r}"
    return True, f"{x.shape[0]} rows equal"


def teaching_facts(out):
    dec = pd.read_csv(out / "practice1_decisions.tsv", sep="\t")
    check(dict(zip(dec.sample_id, dec.decision)) == EXPECTED_P1, "Part A: every decision as specified")
    check(dec.decision.value_counts().to_dict() == {"PASS": 6, "FAIL": 3, "REVIEW": 2, "BLOCKED": 1},
          "Part A: tally PASS 6, REVIEW 2, FAIL 3, BLOCKED 1")
    st = pd.read_csv(out / "practice1_read_stats.tsv", sep="\t").set_index("platform")
    span = st.reads_spanning_target
    check(span["illumina_2x150"] == 0, "Part B: Illumina spans 0 reads of the 8-kb target")
    check(span["pacbio_hifi"] >= 3 and span["ont_ligation"] >= 3,
          f"Part B: HiFi spans {span['pacbio_hifi']}, ONT spans {span['ont_ligation']} (>= 3 each)")
    f = pd.read_csv(out / "practice2_qc_funnel.tsv", sep="\t")
    n = f.n_remaining.to_numpy()
    check(all(n[i] == n[i - 1] - f.n_removed[i] for i in range(1, len(n))), "Practice 2: funnel adds up")
    ex = pd.read_csv(out / "practice2_excluded_samples.tsv", sep="\t")
    check(len(ex) + n[-1] == n[0] and ex.sample_id.is_unique,
          f"Practice 2: {len(ex)} excluded + {n[-1]} in the EUR set = {n[0]} delivered")
    key = pd.read_csv(ROOT / "checkpoints" / "practice2_planted_cases.tsv", sep="\t")
    where = dict(zip(ex.sample_id, ex.step))
    bad = [(r.case, r.sample_id, where.get(r.sample_id, "kept"), r.expected_step)
           for r in key.itertuples() if where.get(r.sample_id, "kept") != r.expected_step]
    check(not bad, f"Practice 2: all {len(key)} planted cases land where documented {bad[:3]}")
    readme = (ROOT / "README.md").read_text()
    stale = [s for s in key.sample_id if s not in readme]
    check(not stale, f"README names every planted sample {stale[:3]}")
    return dict(decisions=dec.decision.value_counts().to_dict(),
                spanning={k: int(v) for k, v in span.items()},
                funnel=dict(zip(f.step, f.n_remaining.astype(int))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true", help="keep the executed copy")
    ap.add_argument("--kernel", default="python3")
    ap.add_argument("--executor", choices=["auto", "nbclient", "ipython-shell"], default="auto",
                    help="auto = nbclient when a kernel can open sockets, else the IPython fallback")
    args = ap.parse_args()

    static_checks()
    global TRANSPORT
    TRANSPORT = kernel_transport()
    executor = args.executor
    if executor == "auto":
        executor = "nbclient" if TRANSPORT else "ipython-shell"
    if executor == "nbclient" and not TRANSPORT:
        sys.exit("nbclient cannot start a kernel here (sockets forbidden or nbclient/ipykernel missing)")
    print(f"     executor: {executor}" + (f" (kernel transport {TRANSPORT})" if executor == "nbclient"
          else " (no usable Jupyter kernel here: sockets forbidden or nbclient/ipykernel missing)"))
    tmp = Path(tempfile.mkdtemp(prefix="day2-verify-")).resolve()
    if ROOT in tmp.parents:                      # never copy the package into itself
        tmp = Path(tempfile.mkdtemp(prefix="day2-verify-", dir="/tmp")).resolve()
    pkg = clean_copy(tmp / "day2-3-samples")
    timings = {}
    for name in NOTEBOOKS:
        seconds, error, stderr, mismatch = run_notebook(pkg / name, args.kernel, executor)
        timings[name] = round(seconds, 2)
        check(error is None, f"{name}: runs top to bottom in {seconds:.1f} s" + (f" ({error})" if error else ""))
        check(not stderr, f"{name}: nothing on stderr" + (f" {stderr[:2]}" if stderr else ""))
        check(not mismatch, f"{name}: no 'matches checkpoint: False'" + (f" {mismatch}" if mismatch else ""))
    total = sum(timings.values())
    check(total < BUDGET_S, f"both notebooks together: {total:.1f} s (budget {BUDGET_S:.0f} s)")

    out = pkg / "outputs"
    for name in TABLES:
        if not (out / name).exists():
            check(False, f"outputs/{name} written")
            continue
        ok, why = same_table(out / name, ROOT / "checkpoints" / name)
        check(ok, f"outputs/{name} equals checkpoints/{name}: {why}")
    for name in FIGURES:
        check((out / name).exists() and (out / name).stat().st_size > 10_000, f"outputs/{name} written")
    facts = teaching_facts(out) if all((out / t).exists() for t in TABLES) else {}

    summary = {"passed": not problems, "executor": executor, "notebook_seconds": timings, "total_seconds": round(total, 2),
               "problems": problems, **facts}
    print(json.dumps(summary, default=int))
    if args.keep:
        print("executed copy kept at", pkg)
    else:
        shutil.rmtree(tmp, ignore_errors=True)
    sys.exit(0 if not problems else 1)


if __name__ == "__main__":
    main()
