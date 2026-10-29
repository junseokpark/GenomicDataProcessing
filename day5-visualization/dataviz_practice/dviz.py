"""dviz: small helpers for the D5 practice "Data visualization in genetics and genomics".

Everything here is plain numpy / pandas / scipy / matplotlib, so you can open this file and read
how each helper works. The notebooks call these functions; nothing is hidden in compiled code.

    import dviz
    dviz.use_style()                       # fonts, colours and sizes used in the lecture slides
    dviz.check("B2.up", n_up, 150, tol=10) # a CHECK: compare your number with the reference
    dviz.save_results("B2")                # write results/B2.json (used by the instructor)
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager

# --------------------------------------------------------------------------------------------
# Paths: the package root is the folder that contains this file.
# --------------------------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGDIR = RESULTS / "figures"

# --------------------------------------------------------------------------------------------
# Colours. The first block matches the lecture slides; OKABE_ITO is a colour-blind-safe set
# for categories (Okabe & Ito, 2008).
# --------------------------------------------------------------------------------------------
INK = "#17202A"        # titles and labels
SLATE = "#5B6573"      # axes, ticks, secondary text, "control" group
GREY = "#C9CFD8"       # context (not significant, background points)
LIGHT = "#E4ECFB"      # light grid lines and reference lines
BLUE = "#1F5FEA"       # the one accent colour (highlight, "treated", "case")
DARKBLUE = "#173FAF"
SKY = "#8FB0F5"
PALEBLUE = "#D3E1F9"
HIGHLIGHT = "#FCF5F0"  # shaded region (for example a deletion)
ORANGE = "#E69F00"     # Okabe-Ito orange
VERMILION = "#D55E00"  # Okabe-Ito vermilion ("up")
GREEN = "#009E73"      # Okabe-Ito bluish green
OI_BLUE = "#0072B2"    # Okabe-Ito blue ("down")
OI_SKY = "#56B4E9"
YELLOW = "#F0E442"
PURPLE = "#CC79A7"
OKABE_ITO = [ORANGE, OI_SKY, GREEN, YELLOW, OI_BLUE, VERMILION, PURPLE, "#000000"]


def _register_fonts():
    """Use the SUIT font shipped in fonts/ (SIL Open Font License) when it is present."""
    fdir = ROOT / "fonts"
    found = False
    if fdir.exists():
        for f in sorted(fdir.glob("*.ttf")):
            try:
                font_manager.fontManager.addfont(str(f))
                found = True
            except Exception:
                pass
    return found


def use_style(base_size: float = 9.0, font: str | None = None):
    """Set matplotlib defaults that follow the lecture's figure rules.

    * one sans-serif font (SUIT if available), text >= 7 pt at the final size
    * only the left and bottom axis lines; muted axis colour so the data stand out
    * vector-friendly export: editable text in SVG and PDF
    """
    has_suit = _register_fonts()
    family = font or ("SUIT" if has_suit else "DejaVu Sans")
    mpl.rcParams.update({
        # a list of families -> matplotlib falls back glyph by glyph (SUIT has no minus sign, Greek letters)
        "font.family": [family, "DejaVu Sans"],
        "font.sans-serif": [family, "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": base_size,
        "axes.titlesize": base_size * 1.15,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlecolor": INK,
        "axes.labelsize": base_size,
        "axes.labelcolor": INK,
        "axes.edgecolor": SLATE,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.prop_cycle": mpl.cycler(color=[BLUE, SLATE, ORANGE, GREEN, VERMILION, OI_SKY, PURPLE]),
        "xtick.color": SLATE, "ytick.color": SLATE,
        "xtick.labelcolor": SLATE, "ytick.labelcolor": SLATE,
        "xtick.labelsize": base_size * 0.9, "ytick.labelsize": base_size * 0.9,
        "legend.frameon": False, "legend.fontsize": base_size * 0.9,
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "axes.unicode_minus": True,
    })
    return family


def mm(x_mm: float) -> float:
    """Millimetres to inches (journals give widths in mm: 89 mm single, 183 mm double column)."""
    return x_mm / 25.4


def savefig(fig, name: str, formats=("png",), outdir: Path | str | None = None, dpi: int = 300):
    """Save a figure as results/figures/<name>.<fmt> and return the paths."""
    out = Path(outdir) if outdir else FIGDIR
    out.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        p = out / f"{name}.{fmt}"
        fig.savefig(p, dpi=dpi, metadata={"Software": None} if fmt == "png" else None)
        paths.append(p)
    return paths


def jitter(n: int, width: float = 0.08, seed: int = 0) -> np.ndarray:
    """Small random horizontal offsets so that points with equal x do not hide each other."""
    return np.random.default_rng(seed).uniform(-width, width, n)


# --------------------------------------------------------------------------------------------
# Statistics helpers (short on purpose, so you can read the method).
# --------------------------------------------------------------------------------------------
def bh_adjust(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (FDR). NaN stays NaN (= not tested).

    Sort the m p-values, multiply the i-th smallest by m / i, then make the sequence
    non-decreasing from the largest p downwards and cap it at 1.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return out
    order = np.argsort(pv)
    scaled = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(scaled[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.minimum(adj, 1.0)
    out[ok] = res
    return out


def genomic_lambda(p) -> float:
    """Genomic inflation factor: median chi-square (1 df) of the p-values / 0.4549."""
    from scipy.stats import chi2
    p = np.asarray(p, dtype=float)
    p = p[~np.isnan(p)]
    return float(np.median(chi2.isf(p, 1)) / chi2.ppf(0.5, 1))


def qq_points(p):
    """Expected and observed -log10 p for a QQ plot (sorted, smallest p last)."""
    p = np.asarray(p, dtype=float)
    p = np.sort(p[~np.isnan(p)])[::-1]
    n = p.size
    expected = -np.log10((np.arange(n, 0, -1) - 0.5) / n)
    return expected, -np.log10(p)


def manhattan_x(chrom, pos, gap: float = 2e7):
    """Cumulative genome position for a Manhattan plot, with a gap between chromosomes.

    Returns (x, ticks) where ticks is a DataFrame with the centre of every chromosome.
    """
    df = pd.DataFrame({"chrom": np.asarray(chrom), "pos": np.asarray(pos, dtype=float)})
    offsets, centres, off = {}, [], 0.0
    for c in sorted(df["chrom"].unique(), key=_chrom_key):
        mx = df.loc[df["chrom"] == c, "pos"].max()
        offsets[c] = off
        centres.append((c, off + mx / 2))
        off += mx + gap
    x = df["pos"].to_numpy() + df["chrom"].map(offsets).to_numpy(dtype=float)
    return x, pd.DataFrame(centres, columns=["chrom", "centre"])


def _chrom_key(c):
    s = str(c).replace("chr", "")
    return (0, int(s)) if s.isdigit() else (1, s)


def segment(values, min_size: int = 5, z: float = 6.0, max_segments: int = 200) -> pd.DataFrame:
    """Binary segmentation of a 1-D signal (a short, readable version of the CBS idea).

    Repeatedly split a segment at the position where the means on the two sides differ most
    (a t-like statistic). Keep the split only if the statistic is larger than `z`.
    NaN values (masked bins) are skipped and do not break segments.
    Returns one row per segment: index of first and last value, number of values and mean.
    """
    v = np.asarray(values, dtype=float)
    idx = np.flatnonzero(~np.isnan(v))
    x = v[idx]
    sd = np.median(np.abs(np.diff(x))) / (0.6745 * np.sqrt(2)) if x.size > 2 else 1.0  # robust noise
    todo, done = [(0, x.size)], []
    while todo and len(done) + len(todo) < max_segments:
        a, b = todo.pop()
        seg = x[a:b]
        n = seg.size
        if n < 2 * min_size:
            done.append((a, b))
            continue
        cs = np.cumsum(seg)
        k = np.arange(min_size, n - min_size + 1)
        left = cs[k - 1] / k
        right = (cs[-1] - cs[k - 1]) / (n - k)
        t = np.abs(left - right) / (sd * np.sqrt(1 / k + 1 / (n - k)))
        j = int(np.argmax(t))
        if t[j] > z:
            todo += [(a, a + k[j]), (a + k[j], b)]
        else:
            done.append((a, b))
    done += todo
    rows = [(int(idx[a]), int(idx[b - 1]), b - a, float(np.mean(x[a:b]))) for a, b in sorted(done)]
    return pd.DataFrame(rows, columns=["first", "last", "n_bins", "mean"])


def calibration_table(y_true, prob, n_bins: int = 10, min_count: int = 20) -> pd.DataFrame:
    """Observed fraction of positives per bin of predicted probability (bins with < min_count dropped)."""
    y = np.asarray(y_true, dtype=float)
    p = np.asarray(prob, dtype=float)
    edges = np.linspace(0, 1, n_bins + 1)
    b = np.clip(np.digitize(p, edges[1:-1]), 0, n_bins - 1)
    df = pd.DataFrame({"bin": b, "p": p, "y": y}).groupby("bin").agg(
        n=("y", "size"), mean_predicted=("p", "mean"), observed=("y", "mean")).reset_index()
    return df[df["n"] >= min_count].reset_index(drop=True)


# --------------------------------------------------------------------------------------------
# Colour-vision deficiency simulation (Machado, Oliveira & Fernandes 2009, severity 1.0).
# --------------------------------------------------------------------------------------------
_CVD = {
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968],
                              [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
    "protanopia": np.array([[0.152286, 1.052583, -0.204868],
                            [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779],
                            [-0.078411, 0.930809, 0.147602],
                            [0.004733, 0.691367, 0.303900]]),
}


def _to_linear(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def _to_srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def simulate_cvd(image, kind: str = "deuteranopia") -> np.ndarray:
    """Return how an RGB image looks with a colour-vision deficiency, or in greyscale.

    image: a path to a PNG/JPG or an array (H x W x 3 or 4, values 0-1 or 0-255).
    kind: 'deuteranopia', 'protanopia', 'tritanopia' or 'greyscale'.
    """
    if isinstance(image, (str, os.PathLike)):
        image = plt.imread(str(image))
    a = np.asarray(image, dtype=float)
    if a.max() > 1.0:
        a = a / 255.0
    rgb = a[..., :3]
    lin = _to_linear(rgb)
    if kind == "greyscale":
        y = lin @ np.array([0.2126, 0.7152, 0.0722])
        out = np.repeat(_to_srgb(y)[..., None], 3, axis=-1)
    else:
        out = _to_srgb(lin @ _CVD[kind].T)
    return out


def simulate_deuteranopia(image) -> np.ndarray:
    """Shortcut used in the slides: simulate_cvd(image, 'deuteranopia')."""
    return simulate_cvd(image, "deuteranopia")


def show_cvd(image_path, kinds=("deuteranopia", "greyscale"), width: float = 9.0):
    """Show a saved figure next to its colour-vision and greyscale simulations."""
    img = plt.imread(str(image_path))
    n = 1 + len(kinds)
    fig, axs = plt.subplots(1, n, figsize=(width, width / n * img.shape[0] / img.shape[1] + 0.4))
    panels = [("original", img[..., :3])] + [(k, simulate_cvd(img, k)) for k in kinds]
    for ax, (title, im) in zip(axs, panels):
        ax.imshow(im)
        ax.set_title(title, fontsize=9)
        ax.axis("off")
    return fig


# --------------------------------------------------------------------------------------------
# UpSet plot (Lex et al., 2014) drawn with matplotlib only.
# --------------------------------------------------------------------------------------------
def upset_counts(df: pd.DataFrame, sets) -> pd.DataFrame:
    """Count rows for every combination of True/False across the boolean columns `sets`."""
    sets = list(sets)
    b = df[sets].astype(bool)
    cnt = b.value_counts()                      # one row per observed combination
    out = cnt.reset_index()
    out.columns = sets + ["count"]
    out = out[out[sets].any(axis=1)]            # drop "in no set"
    return out.sort_values(["count"] + sets, ascending=False).reset_index(drop=True)


def upset(df: pd.DataFrame, sets, max_bars: int = 10, highlight_all: bool = True, fig=None,
          bar_color=SLATE, accent=BLUE, title: str | None = None):
    """Draw an UpSet plot: one bar per intersection (sorted by size) and a dot matrix below.

    Returns (fig, counts) so you can print the numbers behind the bars.
    """
    sets = list(sets)
    counts = upset_counts(df, sets).head(max_bars)
    sizes = df[sets].astype(bool).sum()
    if fig is None:
        fig = plt.figure(figsize=(6.5, 3.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 4], height_ratios=[3, 2], wspace=0.6, hspace=0.05)   # room for the set names
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_dot)
    x = np.arange(len(counts))
    colors = [accent if (highlight_all and bool(r[sets].all())) else bar_color for _, r in counts.iterrows()]
    ax_bar.bar(x, counts["count"], color=colors, width=0.6)
    for xi, c in zip(x, counts["count"]):
        ax_bar.text(xi, c, f"{c:,}", ha="center", va="bottom", fontsize=mpl.rcParams["font.size"] * 0.8, color=INK)
    ax_bar.set_ylabel("Variants")
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    if title:
        ax_bar.set_title(title)
    yy = np.arange(len(sets))
    for xi, (_, r) in zip(x, counts.iterrows()):
        on = [j for j, s in enumerate(sets) if bool(r[s])]
        ax_dot.scatter([xi] * len(sets), yy, s=28, color=GREY, zorder=1)
        ax_dot.scatter([xi] * len(on), on, s=28, color=INK, zorder=2)
        if len(on) > 1:
            ax_dot.plot([xi, xi], [min(on), max(on)], color=INK, lw=1.5, zorder=1)
    ax_dot.set_yticks(yy, [f"{s} ({sizes[s]:,})" for s in sets])
    ax_dot.tick_params(axis="x", bottom=False, labelbottom=False)
    for sp in ["left", "bottom"]:
        ax_dot.spines[sp].set_visible(False)
    ax_dot.tick_params(axis="y", left=False)
    ax_set.barh(yy, [sizes[s] for s in sets], color=GREY, height=0.5)
    ax_set.invert_xaxis()
    ax_set.tick_params(axis="y", left=False, labelleft=False)
    ax_set.spines["left"].set_visible(False)
    ax_set.set_xlabel("Set size", fontsize=mpl.rcParams["font.size"] * 0.85)
    ax_dot.set_ylim(-0.6, len(sets) - 0.4)
    return fig, counts


# --------------------------------------------------------------------------------------------
# CHECK cells: compare a number with its reference value and keep a record.
# --------------------------------------------------------------------------------------------
_RESULTS: dict = {}


def check(name: str, value, expected=None, tol=None, rel: float | None = None, between=None, note: str = ""):
    """Print a CHECK line and record the value.

    check("B2.up", 152, 150, tol=10)      -> passes if |value - expected| <= tol
    check("B3.lambda", 1.04, between=(0.98, 1.10))
    check("A1.r", 0.816)                    -> just records the value
    Returns True / False (or None when there is nothing to compare).
    """
    ok = None
    if between is not None:
        ok = bool(between[0] <= value <= between[1])
        ref = f"in [{between[0]}, {between[1]}]"
    elif expected is not None:
        if isinstance(expected, (int, float, np.integer, np.floating)) and not isinstance(expected, bool):
            t = tol if tol is not None else (abs(expected) * rel if rel is not None else 0)
            ok = bool(abs(float(value) - float(expected)) <= t + 1e-12)
            ref = f"{expected}" + (f" ± {t:g}" if t else "")
        else:
            ok = bool(value == expected)
            ref = f"{expected}"
    else:
        ref = "recorded"
    v = value.item() if hasattr(value, "item") else value
    _RESULTS[name] = {"value": v, "reference": ref, "ok": ok, "note": note}
    mark = {True: "✓", False: "✗", None: "·"}[ok]
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        shown = str(v)
    elif isinstance(v, int):
        shown = f"{v:,}"
    else:
        shown = f"{v:,.4g}"
    print(f"CHECK {mark} {name} = {shown}   (reference: {ref}){'  ' + note if note else ''}")
    return ok


def save_results(exercise: str, extra: dict | None = None) -> Path:
    """Write all CHECK values whose name starts with `exercise` to results/<exercise>.json."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    data = {k: v for k, v in _RESULTS.items() if k.split(".")[0] == exercise}
    if extra:
        data.update(extra)
    n_fail = sum(1 for v in data.values() if isinstance(v, dict) and v.get("ok") is False)
    p = RESULTS / f"{exercise}.json"
    p.write_text(json.dumps(data, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))
    print(f"saved results/{exercise}.json  ({len(data)} values, {n_fail} failed checks)")
    return p
