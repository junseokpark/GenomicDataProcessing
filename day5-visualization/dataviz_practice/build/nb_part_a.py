"""Cells for practice.ipynb: title, Step 0 set-up, Step 1 (Exercise 1, diagnose) and Practice A (A1-A5)."""


def add(nb):
    nb.md(r"""
    # Practice · From question to a figure you can defend

    **Diagnose → simple data, simple figures (Practice A) → realistic data: test, check, visualize (Practice B) → defend and present**

    *Genomic Data Processing · Day 5 (10/27) · Data visualization in genetics and genomics · about 90 minutes, in pairs*

    > **All data in this practice are simulated (computer-generated). No real patient data are used.**
    > The simulations use real anchors: real human gene names and hg38 coordinates (MANE Select), and published biology
    > such as genes that respond to the steroid drug dexamethasone. So the figures look like real research figures, and
    > every realistic dataset comes with a **truth file**: you can test whether a method finds the answer that was planted.

    **How to use this notebook.** Run the cells from top to bottom (Shift + Enter). Read the short text before each cell:
    it says why the step exists and which new words appear.

    * `# CHECK ·` cells compare your number with a reference value and print **✓** or **✗**. If you see ✗, stop and find
      out why before you change anything else.
    * **What changed?** says what the data or the figure look like after the step, and which number proves it.
    * **Your turn** is a small task: change one thing and look again. The answers are folded (click *Answer*).
    * Practice B has two core exercises for everyone (B1, B2). After B2, each pair opens **one elective notebook** in
      `electives/` (B3 GWAS, B4 single cell, B5 copy number, B6 variant callers, B7 classifiers, B8 TP53 lollipop).

    | Part | Steps | Minutes |
    |---|---|---|
    | Set-up | Step 0 | 5 |
    | Exercise 1 · Diagnose | Step 1 | 8 |
    | Practice A · simple data, simple figures | A1–A5 | 20 |
    | Practice B · core | B1 (two groups), B2 (RNA-seq) | 27 |
    | Practice B · elective | one of B3–B8 (own notebook) | 20 |
    | Defend and present | D1–D3, then one figure per pair | 10 + presentations |
    """)

    nb.md(r"""
    ## Step 0 · Set-up  *(5 min)*

    **Google Colab:** upload `dataviz_practice_dataset.zip`, then run in a cell `!unzip -q dataviz_practice_dataset.zip`
    and `%cd dataviz_practice`, and install the two extra packages with `!pip install -q umap-learn adjustText`.
    **Your laptop (Jupyter):** unzip the folder, open a terminal in it, run `pip install -r requirements.txt`, then open
    this notebook.

    The cell below loads the packages and `dviz.py`, a small file of helpers written for this course (figure style,
    colours, and the CHECK function). Open `dviz.py` in a text editor at any time: every helper is a few readable lines.
    """)

    nb.code(r"""
    # Step 0 · Set-up: imports, folders and figure style
    import os, sys, time
    from pathlib import Path
    ROOT = Path.cwd() if (Path.cwd() / "dviz.py").exists() else Path.cwd().parent
    sys.path.insert(0, str(ROOT))
    os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".numba_cache"))   # needed by umap-learn on some systems
    import numpy as np
    import pandas as pd
    import matplotlib
    import matplotlib.pyplot as plt
    import scipy
    from scipy import stats
    import sklearn
    import dviz
    font = dviz.use_style()
    DATA, FIG = ROOT / "data", ROOT / "results" / "figures"
    FIG.mkdir(parents=True, exist_ok=True)
    T_START = time.time()
    print(f"python {sys.version.split()[0]} · numpy {np.__version__} · pandas {pd.__version__} · "
          f"matplotlib {matplotlib.__version__} · scipy {scipy.__version__} · scikit-learn {sklearn.__version__}")
    print("figure font:", font, "· data folder:", DATA)
    """)

    nb.code(r"""
    # CHECK · every input file of this notebook is present
    needed = ["anscombe.tsv", "mapping_rates.tsv", "expression_by_group.tsv", "timecourse.tsv",
              "celltype_composition.tsv", "methylation_cpg.tsv.gz", "rnaseq_counts.tsv.gz", "rnaseq_metadata.tsv",
              "de_results.tsv", "truth/rnaseq_truth.tsv"]
    missing = [f for f in needed if not (DATA / f).exists()]
    dviz.check("S0.files_missing", len(missing), 0, note=", ".join(missing))
    """)

    # ---------------------------------------------------------------------------------------------- Ex 1
    nb.md(r"""
    ## Step 1 · Exercise 1 · Diagnose: spot the misleader  *(8 min)*

    **Why this step.** Before we make figures, we practise *reading* them critically. The eight figures below were made
    from the same data you will use today, and each one breaks at least one principle from Part 2 of the lecture.

    **Your turn (in pairs).** For each figure M1–M8 (M for *misleader*) write down (1) the principle it breaks and (2) the fix, in one line
    each. Then, if you use an AI assistant, give it **one** figure together with the critique prompt in
    `prompts/critique_prompt.md`: what did it find that you missed, and did it invent anything?
    """)

    nb.code(r"""
    # Step 1 · show the eight misleading figures from bad_figures/
    from matplotlib.image import imread
    fig, axs = plt.subplots(2, 4, figsize=(13, 6.2))
    for k, ax in enumerate(axs.flat, start=1):
        p = ROOT / "bad_figures" / f"M{k}.png"
        if p.exists():
            ax.imshow(imread(p))
        ax.set_title(f"M{k}", loc="left")
        ax.axis("off")
    fig.tight_layout()
    plt.show()
    """)

    nb.md(r"""
    | Figure | Principle it breaks | Fix |
    |---|---|---|
    | M1 | | |
    | M2 | | |
    | M3 | | |
    | M4 | | |
    | M5 | | |
    | M6 | | |
    | M7 | | |
    | M8 | | |

    <details><summary>Answer</summary>

    M1 bar of means with SE hides n = 4 and an outlier (show every point, name the test) · M2 bars do not start at zero
    (start at 0, or use dots with a labelled zoomed axis) · M3 rainbow colours, no colour bar, rows not scaled (z-score,
    diverging map centred at 0, colour bar) · M4 zero-expressing cells drawn on top (draw high values last) · M5 pie with
    many slices, 3-D shadow and a legend (sorted bars with direct labels) · M6 red-green, raw p, no thresholds, untested
    genes shown (adjusted p, thresholds, colour-blind-safe colours, counts) · M7 a claim based on distances between UMAP
    clusters (UMAP keeps neighbours, not distances) · M8 two y-axes make two series look related (both are log2 FC vs 0 h: one shared axis, direct
    labels). The full answer key is in `solutions/ANSWER_KEY_bad_figures.md`.
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- A1
    nb.md(r"""
    # Practice A · Simple data, simple figures

    Practice A uses small tables that you can read by eye. Each step makes one basic figure and asks you to **change one
    thing and look**: a setting, an axis or a colour. The goal is to see how a single choice changes what a reader
    concludes.

    ## Step A1 · Plot the data before you trust a summary  *(4 min)*

    **Why this step.** Summary numbers such as the mean or the correlation are useful, but very different data can give
    the same numbers. Anscombe built four small datasets in 1973 to show this. These are the real published values.

    **New terms**
    * **Correlation r** – how close the points are to a straight line (−1 to 1).
    * **Linear regression** – the straight line that best predicts y from x (`y = intercept + slope × x`).
    """)

    nb.code(r"""
    # Step A1a · the same summary numbers for four different datasets
    ans = pd.read_csv(DATA / "anscombe.tsv", sep="\t")
    rows = []
    for name, d in ans.groupby("dataset"):
        slope, intercept = np.polyfit(d["x"], d["y"], 1)
        rows.append({"dataset": name, "n": len(d), "mean y": d["y"].mean(), "variance y": d["y"].var(),
                     "r": np.corrcoef(d["x"], d["y"])[0, 1], "intercept": intercept, "slope": slope})
    summary_a1 = pd.DataFrame(rows).round(2)
    summary_a1
    """)

    nb.code(r"""
    # Step A1b · now plot them
    fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharex=True, sharey=True)
    xs = np.array([2, 20])
    for ax, (name, d) in zip(axs, ans.groupby("dataset")):
        ax.scatter(d["x"], d["y"], s=22, color=dviz.BLUE, zorder=3)
        ax.plot(xs, 3 + 0.5 * xs, color=dviz.SLATE, lw=1)                # the shared regression line
        ax.set_title(f"Set {name} · r = {np.corrcoef(d['x'], d['y'])[0, 1]:.2f}")
        ax.set_xlabel("x")
    axs[0].set_ylabel("y")
    fig.tight_layout()
    dviz.savefig(fig, "A1_anscombe")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · all four sets have r = 0.82 and the line y = 3 + 0.5 x
    for _, r in summary_a1.iterrows():
        dviz.check(f"A.A1_r_{r['dataset']}", round(r["r"], 2), 0.82, tol=0.005)
    dviz.check("A.A1_slope_range", round(summary_a1["slope"].max() - summary_a1["slope"].min(), 3), 0.0, tol=0.002)
    """)

    nb.md(r"""
    **What changed?** The table says the four sets are the same (r = 0.82, y = 3 + 0.5x for all four). The plots say they
    are not: Set I is a noisy line, II is a curve, III is a line with one outlier, and IV has no relationship except
    one point far to the right.

    **Your turn.** Which set would you trust a correlation for? Remove the outlier in Set III
    (`d[d["y"] < 12]`) and compute r again.

    <details><summary>Answer</summary>

    Only Set I. Without its outlier, Set III is an almost perfect line (r ≈ 1.00), so one point lowered r to 0.82.
    Rule: plot the data before you report a summary statistic.
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- A2
    nb.md(r"""
    ## Step A2 · One gene, two groups: how to show a distribution  *(5 min)*

    **Why this step.** A typical first question is "is gene X higher in patients (Case) than in healthy people
    (Control)?". The table `expression_by_group.tsv` has the expression of one gene, in TPM (transcripts per million),
    for two studies: a small one (4 + 4 people) and a larger one (40 + 40 people).

    **New terms**
    * **Histogram** – counts of values in bins; the number of bins is a choice.
    * **Kernel density estimate (KDE)** – a smooth curve made by adding a small bump for every value; the *bandwidth*
      sets how smooth it is.
    * **Box plot** – median (line), quartiles (box, the middle 50%) and whiskers to 1.5 × the box height.
    * **Violin** – a mirrored KDE; **strip plot** – every value as a point.
    """)

    nb.code(r"""
    # Step A2a · the Case group of the larger study, three histograms and three smoothing bandwidths
    ebg = pd.read_csv(DATA / "expression_by_group.tsv", sep="\t")
    big = ebg[ebg["dataset"] == "large_n"]
    case = big.loc[big["group"] == "Case", "expression"].to_numpy()
    grid = np.linspace(30, 140, 300)
    fig, axs = plt.subplots(1, 4, figsize=(13, 2.8), sharex=True)
    for ax, bins in zip(axs[:3], [3, 12, 60]):
        ax.hist(case, bins=bins, range=(40, 130), color=dviz.SKY, edgecolor="white")
        ax.set_title(f"{bins} bins")
    for bw, col in [(1.2, dviz.SLATE), (0.4, dviz.BLUE), (0.15, dviz.ORANGE)]:
        axs[3].plot(grid, stats.gaussian_kde(case, bw_method=bw)(grid), color=col, label=f"bandwidth {bw}")
    axs[3].set_title("Kernel density")
    axs[3].legend(fontsize=7)
    for ax in axs:
        ax.set_xlabel("Expression (TPM)")
    axs[0].set_ylabel("Cases (n = 40)")
    fig.tight_layout()
    dviz.savefig(fig, "A2_bins_bandwidth")
    plt.show()
    """)

    nb.code(r"""
    # Step A2b · box, violin, strip and all three together (larger study)
    groups = ["Control", "Case"]
    vals = [big.loc[big["group"] == g, "expression"].to_numpy() for g in groups]
    cols = [dviz.SLATE, dviz.BLUE]
    fig, axs = plt.subplots(1, 4, figsize=(13, 3), sharey=True)
    axs[0].boxplot(vals, positions=[0, 1], widths=0.5)
    axs[1].violinplot(vals, positions=[0, 1], showmedians=True)
    for i, (v, c) in enumerate(zip(vals, cols)):
        axs[2].scatter(i + dviz.jitter(len(v), 0.12, seed=i), v, s=10, color=c)
        axs[3].violinplot([v], positions=[i], showextrema=False)
        axs[3].boxplot([v], positions=[i], widths=0.12, showfliers=False)
        axs[3].scatter(i + 0.3 + dviz.jitter(len(v), 0.05, seed=i), v, s=8, color=c)
    for ax, t in zip(axs, ["Box", "Violin", "Strip (every point)", "Violin + box + points"]):
        ax.set_title(t)
        ax.set_xticks([0, 1], [f"{g}\nn = {len(v)}" for g, v in zip(groups, vals)])
    axs[0].set_ylabel("Expression (TPM)")
    fig.tight_layout()
    dviz.savefig(fig, "A2_box_violin_points")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · the numbers behind the box plots (medians and quartiles)
    for g, v in zip(groups, vals):
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        print(f"{g:8s} n = {len(v)}  median {med:.1f}  quartiles {q1:.1f}-{q3:.1f}  range {v.min():.1f}-{v.max():.1f}")
    dviz.check("A.A2_case_median", round(float(np.median(vals[1])), 2))
    dviz.check("A.A2_control_median", round(float(np.median(vals[0])), 2))
    dviz.check("A.A2_case_above_90", int((vals[1] > 90).sum()), note="cases in the upper subgroup (> 90 TPM)")
    """)

    nb.md(r"""
    **What changed?** With 3 bins, or a wide bandwidth, the cases look like one broad hump. With about 12 bins, or a
    bandwidth around 0.4, two groups of cases appear: most cases look like controls, and a subgroup is much higher.
    60 bins (or bandwidth 0.15) show mostly noise. The box plot hides the two groups (it shows only five numbers); the
    strip plot and the combined plot show them.

    **Your turn.** Change `bins=bins` to `bins=6` for the first histogram. Can you still see two groups? What would you
    write in a figure legend so that a reader can judge your choice?

    <details><summary>Answer</summary>

    With 6 bins the two groups are barely visible. Write the number of bins (or the KDE bandwidth) in the legend, and
    show the points when n is small enough.
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- A3
    nb.md(r"""
    ## Step A3 · Bars start at zero; dots may zoom  *(3 min)*

    **Why this step.** `mapping_rates.tsv` has the percentage of sequencing reads that mapped to the genome for three
    versions (A, B, C) of a pipeline, on the same four libraries. A bar is read by its **length**, so its axis must start
    at zero. A dot is read by its **position**, so a zoomed axis is fine if it is labelled.
    """)

    nb.code(r"""
    # Step A3 · the same three means drawn three ways
    mr = pd.read_csv(DATA / "mapping_rates.tsv", sep="\t")
    means = mr.groupby("pipeline")["pct_mapped"].mean()
    colors = [dviz.GREY, dviz.GREY, dviz.BLUE]
    fig, axs = plt.subplots(1, 3, figsize=(11, 3))
    for ax, lo, title in [(axs[0], 95.5, "Bars, axis from 95.5% (wrong)"), (axs[1], 0, "Bars, axis from 0")]:
        ax.bar(means.index, means.values, color=colors, width=0.6)
        ax.set_ylim(lo, 100)
        ax.set_title(title)
        for i, v in enumerate(means.values):
            ax.text(i, v + 0.05 * (100 - lo) / 4.5, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    for i, (p, d) in enumerate(mr.groupby("pipeline")):
        axs[2].scatter(np.full(len(d), i) + dviz.jitter(len(d), 0.06, seed=i), d["pct_mapped"], s=18, color=colors[i] if i == 2 else dviz.SLATE)
        axs[2].hlines(d["pct_mapped"].mean(), i - 0.2, i + 0.2, color=dviz.INK, lw=1.5)
    axs[2].set_xticks([0, 1, 2], list(means.index))
    axs[2].set_ylim(95.5, 98.5)
    axs[2].set_title("Dots + mean, zoomed axis (labelled)")
    for ax in axs:
        ax.set_ylabel("Reads mapped (%)")
        ax.set_xlabel("Pipeline version (4 libraries each)")
    fig.tight_layout()
    dviz.savefig(fig, "A3_bars_from_zero")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · how much bigger does C look than A when the axis starts at 95.5%?
    looks = (means["C"] - 95.5) / (means["A"] - 95.5)
    true = means["C"] / means["A"]
    print(f"means: A {means['A']:.1f} · B {means['B']:.1f} · C {means['C']:.1f}")
    dviz.check("A.A3_bar_ratio_looks", round(float(looks), 1), note="bar C height / bar A height with the axis at 95.5%")
    dviz.check("A.A3_true_ratio", round(float(true), 3), 1.02, tol=0.001)
    """)

    nb.md(r"""
    **What changed?** Same numbers, different impressions. From 95.5%, bar C is about four times taller than bar A,
    although C maps only 1.9 percentage points (2%) more reads. From zero, the three bars look almost equal, which is
    the truth. The dot plot zooms in honestly: the axis is labelled, and you also see that the four libraries agree.

    **Your turn.** In the first panel change `95.5` to `96.0`. How many times taller does C look now?

    <details><summary>Answer</summary>

    (98.0 − 96.0) / (96.1 − 96.0) = 20 times. The closer the baseline is to the smallest bar, the bigger the distortion.
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- A4
    nb.md(r"""
    ## Step A4 · Lines for time, bars for composition  *(4 min)*

    **Why this step.** Two common tables: (1) a **time course** – expression of three genes at 0, 1, 2, 4, 8 and 24 hours
    after adding dexamethasone (a steroid drug) to cells, three replicate wells per time; (2) a **composition** – the
    percentage of each blood cell type in four donors (from the single-cell data of elective B4).
    A line connects ordered values, so it suits time. Composition is about parts of a whole.

    **New terms**
    * **log2 fold change (log2 FC)** – log2 of (expression at time t / expression at 0 h). +1 = doubled, −1 = halved.
    * **SD (standard deviation)** – the spread of the replicates.
    """)

    nb.code(r"""
    # Step A4a · time course: mean ± SD of 3 replicates, the replicates as points, direct labels
    tc = pd.read_csv(DATA / "timecourse.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(6, 3.3))
    for g, col in zip(["TSC22D3", "FKBP5", "PER1"], [dviz.BLUE, dviz.ORANGE, dviz.SLATE]):
        d = tc[tc["gene"] == g]
        m = d.groupby("time_h")["log2fc_vs_0h"].agg(["mean", "std"])
        ax.errorbar(m.index, m["mean"], yerr=m["std"], color=col, marker="o", ms=4, capsize=2, lw=1.6)
        ax.scatter(d["time_h"], d["log2fc_vs_0h"], color=col, s=8, alpha=0.5)
        ax.text(24.6, m["mean"].iloc[-1], g, color=col, va="center", fontsize=9)
    ax.axhline(0, color=dviz.LIGHT, lw=1, zorder=0)
    ax.set_xticks([0, 1, 2, 4, 8, 24])
    ax.set_xlim(-0.5, 28)
    ax.set_xlabel("Hours after dexamethasone")
    ax.set_ylabel("log2 FC vs 0 h")
    ax.set_title("Mean ± SD, n = 3 wells per time")
    fig.tight_layout()
    dviz.savefig(fig, "A4_timecourse")
    plt.show()
    """)

    nb.code(r"""
    # Step A4b · composition: a pie vs sorted bars (all donors), and per donor
    comp = pd.read_csv(DATA / "celltype_composition.tsv", sep="\t")
    overall = comp.groupby("cell_type")["n_cells"].sum().sort_values()
    pct = 100 * overall / overall.sum()
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.4), gridspec_kw={"width_ratios": [1, 1.2, 1.4]})
    colour = dict(zip(pct.index[::-1], dviz.OKABE_ITO))         # one colour per cell type, used in both panels
    axs[0].pie(pct.values[::-1], colors=[colour[c] for c in pct.index[::-1]], startangle=90, counterclock=False)
    axs[0].set_title("Pie: which is bigger?")
    axs[1].barh(pct.index, pct.values, color=dviz.SKY)
    for i, v in enumerate(pct.values):
        axs[1].text(v + 0.4, i, f"{v:.1f}%", va="center", fontsize=8)
    axs[1].set_xlabel("Cells (%)")
    axs[1].set_title(f"Sorted bars, n = {overall.sum():,} cells")
    wide = comp.pivot(index="donor", columns="cell_type", values="percent")[list(pct.index[::-1])]
    left = np.zeros(len(wide))
    for ct in wide.columns:
        axs[2].barh(wide.index, wide[ct], left=left, color=colour[ct], label=ct)
        left += wide[ct].to_numpy()
    axs[2].set_xlabel("Cells (%)")
    axs[2].set_title("Stacked, per donor")
    axs[2].legend(ncol=2, fontsize=7, bbox_to_anchor=(1.0, -0.25), loc="upper right")
    fig.tight_layout()
    dviz.savefig(fig, "A4_composition")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · percentages add up to 100 in every donor; the two cell types the lecture compared
    dviz.check("A.A4_donor_sum_max_error", round(float((comp.groupby("donor")["percent"].sum() - 100).abs().max()), 1), 0.0, tol=0.15)
    dviz.check("A.A4_cd14_mono_pct", round(float(pct["CD14 Mono"]), 1), 19.0, tol=0.05)
    dviz.check("A.A4_memory_cd4_pct", round(float(pct["Memory CD4 T"]), 1), 16.3, tol=0.05)
    dviz.check("A.A4_fkbp5_24h", round(float(tc.query("gene == 'FKBP5' and time_h == 24")["log2fc_vs_0h"].mean()), 2))
    """)

    nb.md(r"""
    **What changed?** The time course shows three different shapes: TSC22D3 rises fast and stays up, FKBP5 rises slowly
    and is still rising at 24 h, PER1 peaks at 1 h and falls back. Direct labels replace the legend. In the pie you cannot
    tell whether CD14 monocytes (19.0%) are more frequent than memory CD4 T cells (16.3%); the sorted bars answer at a
    glance. In the stacked bars only the first segment shares a baseline, so the middle cell types are hard to compare
    across donors.

    **Your turn.** Remove `marker="o"` and the scatter line from A4a. What information is lost?

    <details><summary>Answer</summary>

    You no longer see when the samples were taken (the uneven spacing 0, 1, 2, 4, 8, 24 h) or how much the replicates
    vary. A line between sparse time points also suggests values that were never measured.
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- A5
    nb.md(r"""
    ## Step A5 · Many points: transparency and 2-D density  *(4 min)*

    **Why this step.** `methylation_cpg.tsv.gz` has 60,000 CpG sites measured in the same sample by two methods:
    bisulfite sequencing (WGBS, the reference) and nanopore sequencing. The **β value** is the fraction of reads that are
    methylated at a CpG (0 = never, 1 = always). With 60,000 points a normal scatter plot becomes a solid blob.

    **New terms**
    * **Pearson r** – correlation of the values (sensitive to outliers); **Spearman ρ** – correlation of the ranks.
    * **Hexbin** – a 2-D histogram with hexagonal bins; colour = how many points fall in each bin.
    """)

    nb.code(r"""
    # Step A5 · opaque points, transparent points, and a hexbin density
    me = pd.read_csv(DATA / "methylation_cpg.tsv.gz", sep="\t")
    x, y = me["beta_wgbs"], me["beta_nanopore"]
    r, rho = x.corr(y), x.corr(y, method="spearman")
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    axs[0].scatter(x, y, s=2, color=dviz.BLUE, rasterized=True)
    axs[0].set_title(f"{len(me):,} points, opaque")
    axs[1].scatter(x, y, s=2, color=dviz.BLUE, alpha=0.04, rasterized=True)
    axs[1].set_title("Transparency (alpha 0.04)")
    hb = axs[2].hexbin(x, y, gridsize=45, bins="log", cmap="viridis", mincnt=1)
    fig.colorbar(hb, ax=axs[2], label="CpGs (log scale)")
    axs[2].set_title(f"Hexbin · r = {r:.3f}, ρ = {rho:.3f}")
    for ax in axs:
        ax.set_xlabel("Bisulfite (WGBS) β")
    axs[0].set_ylabel("Nanopore β")
    fig.tight_layout()
    dviz.savefig(fig, "A5_methylation_hexbin")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · correlation, and the disagreement at intermediate methylation vs near 0 or 1
    mid = me[(x > 0.3) & (x < 0.7)]
    ext = me[(x < 0.1) | (x > 0.9)]
    d_mid = (mid["beta_wgbs"] - mid["beta_nanopore"]).abs().median()
    d_ext = (ext["beta_wgbs"] - ext["beta_nanopore"]).abs().median()
    dviz.check("A.A5_pearson_r", round(float(r), 3))
    dviz.check("A.A5_spearman_rho", round(float(rho), 3))
    dviz.check("A.A5_median_abs_diff_mid", round(float(d_mid), 3), note=f"β 0.3-0.7, n = {len(mid):,}")
    dviz.check("A.A5_median_abs_diff_extreme", round(float(d_ext), 3), note=f"β < 0.1 or > 0.9, n = {len(ext):,}")
    dviz.check("A.A5_mid_worse_than_extreme", bool(d_mid > 2 * d_ext), True)
    """)

    nb.md(r"""
    **What changed?** The opaque plot hides where most points are. Transparency helps; the hexbin with a log colour scale
    shows that most CpGs sit near 0 or 1, and that the two methods disagree more at intermediate methylation (β 0.3–0.7),
    where sampling a limited number of reads matters most. The correlation (r ≈ 0.98) summarises the overall trend but not
    this. Spearman's ρ is lower because thousands of CpGs are tied near 0 or 1.

    **Your turn.** Remove `bins="log"` from the hexbin. What happens to the intermediate region?

    <details><summary>Answer</summary>

    With a linear colour scale the two dense corners take almost the whole colour range, so the intermediate region
    looks empty. A log scale is the right choice when counts span several orders of magnitude; say so on the colour bar.
    </details>
    """)
