"""Execute a notebook in-process and save it with outputs (used to make the *_executed.ipynb files).

Why not Jupyter? The build machine cannot open kernel sockets, so this runner executes the code
cells one by one in a single Python namespace, like a kernel would, and stores:
printed text, the value of the last expression (tables as HTML), and every matplotlib figure.
Students should simply run the notebooks in Jupyter or Colab.

usage:  python build/run_notebook.py IN.ipynb OUT.ipynb [--timeout-cell 300]
"""
from __future__ import annotations

import ast
import base64
import contextlib
import io
import json
import os
import sys
import time
import traceback
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nbformat  # noqa: E402


def _fig_outputs():
    outs = []
    for num in plt.get_fignums():
        fig = plt.figure(num)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
        outs.append(nbformat.v4.new_output("display_data", data={
            "image/png": base64.b64encode(buf.getvalue()).decode(), "text/plain": f"<Figure {num}>"}))
    plt.close("all")
    return outs


def _rich(obj):
    data = {"text/plain": repr(obj)}
    for meth, mime in [("_repr_html_", "text/html"), ("_repr_png_", "image/png"), ("_repr_markdown_", "text/markdown")]:
        f = getattr(obj, meth, None)
        if callable(f):
            try:
                v = f()
            except Exception:
                v = None
            if v is None:
                continue
            if mime == "image/png":
                v = base64.b64encode(v).decode() if isinstance(v, (bytes, bytearray)) else v
            data[mime] = v
    return data


def run(nb_in: str, nb_out: str, quiet: bool = False) -> dict:
    nb_path = Path(nb_in).resolve()
    nb = nbformat.read(str(nb_path), as_version=4)
    old_cwd = os.getcwd()
    os.chdir(nb_path.parent)
    sys.path.insert(0, str(nb_path.parent))
    outputs_holder: list = []

    def display(*objs, **kw):
        for o in objs:
            outputs_holder.append(nbformat.v4.new_output("display_data", data=_rich(o)))

    import IPython.display as ipd
    ipd.display = display
    plt.show = lambda *a, **k: None
    ns = {"__name__": "__main__", "display": display}
    timing, n_exec, error = [], 0, None
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        n_exec += 1
        cell.execution_count = n_exec
        outputs_holder.clear()
        out, err = io.StringIO(), io.StringIO()
        src = cell.source
        t0 = time.time()
        try:
            tree = ast.parse(src)
            last = None
            if tree.body and isinstance(tree.body[-1], ast.Expr):
                last = ast.Expression(tree.body.pop().value)
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                exec(compile(tree, f"<cell {n_exec}>", "exec"), ns)
                value = eval(compile(last, f"<cell {n_exec}>", "eval"), ns) if last is not None else None
            cell_outputs = []
            if out.getvalue():
                cell_outputs.append(nbformat.v4.new_output("stream", name="stdout", text=out.getvalue()))
            if err.getvalue().strip():
                cell_outputs.append(nbformat.v4.new_output("stream", name="stderr", text=err.getvalue()))
            cell_outputs += list(outputs_holder)
            cell_outputs += _fig_outputs()
            if value is not None:
                cell_outputs.append(nbformat.v4.new_output("execute_result", data=_rich(value), execution_count=n_exec))
            cell.outputs = cell_outputs
        except Exception as e:  # record the error like Jupyter and stop
            tb = traceback.format_exc().splitlines()
            cell.outputs = ([nbformat.v4.new_output("stream", name="stdout", text=out.getvalue())] if out.getvalue() else []) + [
                nbformat.v4.new_output("error", ename=type(e).__name__, evalue=str(e), traceback=tb)]
            error = {"cell": n_exec, "error": f"{type(e).__name__}: {e}", "source_head": src[:200]}
            timing.append(round(time.time() - t0, 2))
            break
        timing.append(round(time.time() - t0, 2))
        if not quiet:
            print(f"cell {n_exec:3d}  {timing[-1]:6.2f}s  {src.splitlines()[0][:70] if src else ''}", flush=True)
    os.chdir(old_cwd)
    nb.metadata["executed_by"] = "build/run_notebook.py (in-process runner)"
    nbformat.write(nb, nb_out)
    summary = {"notebook": str(nb_in), "cells_run": n_exec, "seconds": round(sum(timing), 1), "error": error,
               "slowest": sorted(timing, reverse=True)[:3]}
    print(json.dumps(summary))
    return summary


if __name__ == "__main__":
    res = run(sys.argv[1], sys.argv[2], quiet="--quiet" in sys.argv)
    sys.exit(1 if res["error"] else 0)
