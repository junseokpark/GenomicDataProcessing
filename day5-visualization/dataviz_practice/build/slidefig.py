"""Helpers for slide figures with exact pixel sizes (288 dpi, like the original deck images)."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import dviz  # noqa: E402

DPI = 288
OUT = HERE / "slide_figs"
FRAMES = pd.read_csv(HERE / "picture_frames.tsv", sep="\t")
PRACTICE_PX = (1800, 1210)          # new practice-slide figures: 6.25 x 4.2 in


def size_of(target: str):
    """Pixel size (w, h) of an existing slide image, e.g. 'image-36-1.png', or of a practice figure."""
    if target.startswith("practice_"):
        return PRACTICE_PX
    row = FRAMES.loc[FRAMES["media"] == target]
    if row.empty:
        raise KeyError(f"{target} not in picture_frames.tsv")
    w, h = eval(row.iloc[0]["px"])
    return int(w), int(h)


def new_fig(target: str, nrows=1, ncols=1, base_size: float = 12.5, **kw):
    """Figure sized for `target` (exact pixels at 288 dpi) with the slide style."""
    dviz.use_style(base_size=base_size)
    plt.rcParams["savefig.bbox"] = "standard"      # exact canvas; never tight for slides
    w, h = size_of(target)
    fig, axs = plt.subplots(nrows, ncols, figsize=(w / DPI, h / DPI), dpi=DPI, **kw)
    return fig, axs


def save_exact(fig, target: str):
    """Save without tight bbox and assert the pixel size matches the slide frame."""
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / (target if target.endswith(".png") else target + ".png")
    with plt.rc_context({"savefig.bbox": "standard"}):
        fig.savefig(p, dpi=DPI, pad_inches=0, metadata={"Software": None})
    plt.close(fig)
    got = Image.open(p).size
    want = size_of(p.name)
    assert got == want, f"{p.name}: {got} != {want}"
    return p
