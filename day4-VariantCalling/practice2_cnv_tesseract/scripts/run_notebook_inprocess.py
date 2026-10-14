#!/usr/bin/env python
"""Fallback notebook executor (instructors only).

Normally run the notebook with Jupyter, or headlessly with
    jupyter nbconvert --to notebook --execute practice2_cnv_tesseract_student.ipynb --output executed.ipynb
On locked-down machines where a Jupyter kernel cannot open its local sockets (PermissionError when
the kernel starts), this script runs every code cell, in order, inside ONE in-process IPython shell
and writes the outputs (text, tables, images, errors) into a copy of the notebook. Stops at the first
error, like nbconvert.

Usage (from the package root):
    python scripts/run_notebook_inprocess.py practice2_cnv_tesseract_student.ipynb practice2_cnv_tesseract_executed.ipynb
"""
import base64
import io
import json
import sys
import time

import nbformat
from IPython.core.displayhook import DisplayHook
from IPython.core.displaypub import DisplayPublisher
from IPython.core.interactiveshell import InteractiveShell
from traitlets.config import Config


class _State:
    outputs = []
    count = 0


def _clean(data):
    out = {}
    for k, v in data.items():
        if isinstance(v, bytes):
            v = base64.b64encode(v).decode("ascii")
        out[k] = v
    return out


class _Stream(io.TextIOBase):
    def __init__(self, name):
        self.name = name

    def write(self, s):
        if not s:
            return 0
        o = _State.outputs
        if o and o[-1].get("output_type") == "stream" and o[-1].get("name") == self.name:
            o[-1]["text"] += s
        else:
            o.append(nbformat.v4.new_output("stream", name=self.name, text=s))
        return len(s)

    def flush(self):
        pass

    def isatty(self):
        return False


class _Pub(DisplayPublisher):
    def publish(self, data, metadata=None, source=None, *, transient=None, update=False, **kw):
        _State.outputs.append(nbformat.v4.new_output("display_data", data=_clean(data), metadata=metadata or {}))

    def clear_output(self, wait=False):
        _State.outputs.clear()


class _Hook(DisplayHook):
    def __call__(self, result=None):
        if result is None or self.quiet():
            return
        data, md = self.shell.display_formatter.format(result)
        _State.outputs.append(nbformat.v4.new_output("execute_result", data=_clean(data), metadata=md,
                                                     execution_count=_State.count))


def main(src, dst):
    nb = nbformat.read(src, as_version=4)
    cfg = Config()
    cfg.HistoryManager.enabled = False
    shell = InteractiveShell.instance(config=cfg)
    shell.display_pub = _Pub(shell=shell)
    shell.displayhook = _Hook(shell=shell)
    shell.display_trap.hook = shell.displayhook      # run_cell installs this hook as sys.displayhook

    def _showtb(etype, evalue, stb):
        _State.outputs.append(nbformat.v4.new_output("error", ename=etype.__name__, evalue=str(evalue), traceback=stb))
    shell._showtraceback = _showtb
    timings, failed = [], None
    old_out, old_err = sys.stdout, sys.stderr
    n = 0
    t_all = time.time()
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        n += 1
        _State.outputs, _State.count = [], n
        t0 = time.time()
        sys.stdout, sys.stderr = _Stream("stdout"), _Stream("stderr")
        try:
            res = shell.run_cell(cell.source, store_history=False)
        finally:
            sys.stdout, sys.stderr = old_out, old_err
        cell.outputs = _State.outputs
        cell.execution_count = n
        timings.append((n, round(time.time() - t0, 2)))
        if res.error_before_exec or res.error_in_exec:
            failed = (n, repr(res.error_in_exec or res.error_before_exec))
            break
    nb.metadata.setdefault("execution", {})["inprocess_runtime_seconds"] = round(time.time() - t_all, 1)
    nbformat.write(nb, dst)
    print(json.dumps({"cells_run": n, "failed": failed, "runtime_s": round(time.time() - t_all, 1),
                      "slowest": sorted(timings, key=lambda x: -x[1])[:5]}))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
