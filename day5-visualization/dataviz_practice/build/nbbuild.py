"""Write a Jupyter notebook from a Python script: nb = NB(); nb.md("..."); nb.code("..."); nb.write(path)."""
from __future__ import annotations
from textwrap import dedent
import nbformat


class NB:
    def __init__(self):
        self.cells = []

    def md(self, text: str):
        self.cells.append(nbformat.v4.new_markdown_cell(dedent(text).strip("\n")))
        return self

    def code(self, text: str):
        self.cells.append(nbformat.v4.new_code_cell(dedent(text).strip("\n")))
        return self

    def write(self, path: str):
        nb = nbformat.v4.new_notebook(cells=self.cells, metadata={
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"}})
        for i, c in enumerate(nb.cells):          # stable ids -> byte-identical rebuilds
            c["id"] = f"c{i:03d}"
        nbformat.write(nb, path)
        return path
