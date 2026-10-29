"""Cells for practice.ipynb: Practice B core (B1 two groups, B2 bulk RNA-seq) and the Defend / present steps."""


def add(nb, tokens=None):
    _md, _code = nb.md, nb.code
    T = tokens or {}
    def sub(t):
        for k, v in T.items():
            t = t.replace("@@" + k + "@@", str(v))
        return t
    nb.md = lambda t: _md(sub(t))
    nb.code = lambda t: _code(sub(t))
    # ---------------------------------------------------------------------------------------------- intro B
    nb.md(r"""
    # Practice B · Realistic data: test, check, then visualize

    In Practice B a figure is the **last** step of an analysis. Every exercise has the same four beats:

    1. **Question** – the biological question in one sentence.
    2. **Test** – run the statistic or method that a real analysis would use, and print the numbers.
    3. **Check** – compare the numbers with the planted truth or a reference value (`# CHECK ·` cells).
    4. **Visualize** – make the figure that shows the result honestly: n, the test, thresholds, units.

    You may write the plotting code yourself, use the cells provided, or ask an AI assistant with a prompt card from
    `prompts/prompt_cards.md`. Whichever you choose, the CHECK cells decide whether the figure is right.

    ## Step B1 · Two groups: which test, which error bar?  *(7 min)*

    ### Question
    Is the gene from Step A2 higher in Case than in Control, in the small study (4 + 4 people) and in the larger one
    (40 + 40)?

    **New terms**
    * **Welch t-test** – compares the two **means**; assumes roughly normal values; allows different spreads.
    * **Mann–Whitney U test** – compares the **ranks** (who is larger than whom); no normality assumption.
    * **SD** – spread of the data. **SEM** = SD / √n – how precisely the mean is known.
      **95% CI** ≈ mean ± t × SEM – the range that would contain the true mean in 95% of repeated studies.

    ### Test
    """)

    nb.code(r"""
    # Step B1a · two tests on two studies
    rows = []
    for ds in ["small_n", "large_n"]:
        d = ebg[ebg["dataset"] == ds]
        ctrl = d.loc[d["group"] == "Control", "expression"].to_numpy()
        case = d.loc[d["group"] == "Case", "expression"].to_numpy()
        welch = stats.ttest_ind(case, ctrl, equal_var=False)
        mwu = stats.mannwhitneyu(case, ctrl, alternative="two-sided")       # exact p-value for small n
        rows.append({"study": ds, "n per group": len(ctrl), "median Control": np.median(ctrl),
                     "median Case": np.median(case), "mean Control": ctrl.mean(), "mean Case": case.mean(),
                     "Welch p": welch.pvalue, "Mann-Whitney p": mwu.pvalue})
    b1 = pd.DataFrame(rows)
    b1.round(4)
    """)

    nb.md(r"""
    ### Check
    """)

    nb.code(r"""
    # CHECK · the two tests disagree for n = 4; both agree for n = 40
    s, l = b1.iloc[0], b1.iloc[1]
    dviz.check("B1.small_median_control", round(s["median Control"], 1), 60.8, tol=0.05)
    dviz.check("B1.small_median_case", round(s["median Case"], 1), 89.8, tol=0.05)
    dviz.check("B1.small_mwu_p", round(s["Mann-Whitney p"], 4), round(2 / 70, 4), tol=0.0001,
               note="2/70: the smallest two-sided p possible with 4 + 4")
    dviz.check("B1.small_welch_p", round(s["Welch p"], 3), 0.107, tol=0.002)
    dviz.check("B1.large_mwu_p", float(f"{l['Mann-Whitney p']:.2g}"), between=(0, 0.001))
    dviz.check("B1.large_welch_p", float(f"{l['Welch p']:.2g}"), between=(0, 0.001))
    """)

    nb.code(r"""
    # Step B1b · three error bars for the same 4 Case values: SD, SEM, 95% CI
    v = ebg.query("dataset == 'small_n' and group == 'Case'")["expression"].to_numpy()
    n, m, sd = len(v), v.mean(), v.std(ddof=1)
    sem = sd / np.sqrt(n)
    tcrit = stats.t.ppf(0.975, n - 1)
    bars = pd.Series({"SD": sd, "SEM": sem, "95% CI": tcrit * sem}).round(1)
    print(f"n = {n}, mean = {m:.1f} TPM, t(0.975, df = {n - 1}) = {tcrit:.2f}")
    bars
    """)

    nb.code(r"""
    # CHECK · the 95% CI is t times wider than the SEM bar (t = 3.18 for n = 4)
    dviz.check("B1.ci_over_sem", round(float(bars["95% CI"] / bars["SEM"]), 2), 3.18, tol=0.01)
    dviz.check("B1.small_case_sd", round(float(sd), 1))
    """)

    nb.md(r"""
    ### Visualize
    """)

    nb.code(r"""
    # Step B1c · the figure: every point, n under each label, the test named (no stars)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, ds, title in [(axs[0], "small_n", "Small study"), (axs[1], "large_n", "Larger study")]:
        d = ebg[ebg["dataset"] == ds]
        vv = [d.loc[d["group"] == g, "expression"].to_numpy() for g in ["Control", "Case"]]
        ax.boxplot(vv, positions=[0, 1], widths=0.35, showfliers=False, medianprops={"color": dviz.INK})
        for i, (x, c) in enumerate(zip(vv, [dviz.SLATE, dviz.BLUE])):
            ax.scatter(i + dviz.jitter(len(x), 0.09, seed=i), x, s=14, color=c, zorder=3)
        p = b1.set_index("study").loc[ds, "Mann-Whitney p"]
        ax.set_title(f"{title} · Mann-Whitney p = {p:.2g}")
        ax.set_xticks([0, 1], [f"Control\nn = {len(vv[0])}", f"Case\nn = {len(vv[1])}"])
        ax.set_ylim(0, 160)
        ax.set_ylabel("Expression (TPM)")
    labels = ["data", "mean ± SD", "mean ± SEM", "mean ± 95% CI"]
    axs[2].scatter(np.zeros(n) + dviz.jitter(n, 0.05), v, color=dviz.GREY, s=16)
    for k, (lab, half, c) in enumerate(zip(labels[1:], bars.values, [dviz.ORANGE, dviz.BLUE, dviz.GREEN]), start=1):
        axs[2].errorbar(k, m, yerr=half, fmt="o", color=c, capsize=4, lw=1.6)
        axs[2].text(k + 0.12, m + half, f"±{half:.1f}", ha="left", va="center", fontsize=8, color=c)
    axs[2].set_xticks(range(4), labels, fontsize=8)
    axs[2].set_title(f"Same 4 Cases, three error bars")
    axs[2].set_ylim(0, 200)
    fig.tight_layout()
    dviz.savefig(fig, "B1_two_groups")
    plt.show()
    """)

    nb.md(r"""
    **What changed?** In the small study the Mann–Whitney test gives p = 0.029 (all four cases are above all four
    controls), while the Welch t-test gives p = 0.11, because one case at 150 TPM makes the spread large. With n = 4,
    0.029 is the smallest p-value the Mann–Whitney test can produce, so neither result is convincing: show the points.
    In the larger study both tests give p < 0.001, and the figure shows what the p-value cannot: only a subgroup of cases
    is high. The SEM bar is the shortest, but the 95% CI, the honest statement of uncertainty, is 3.2 times wider.

    **Your turn.** A figure says "mean ± error, n = 4, *". Which two pieces of information are missing?

    <details><summary>Answer</summary>

    What the error bar is (SD, SEM or CI) and which test produced the star (and its p-value). Also whether n counts
    people (biological replicates) or wells from one sample (technical replicates).
    </details>
    """)

    # ---------------------------------------------------------------------------------------------- B2
    nb.md(r"""
    ## Step B2 · Bulk RNA-seq: from counts to a volcano plot you can defend  *(20 min)*

    ### Question
    **Which genes change when airway smooth-muscle cells are treated with dexamethasone?** Dexamethasone is a steroid
    (glucocorticoid) drug used for asthma. Six cultures were treated and six were left untreated (controls). The twelve
    RNA libraries were sequenced in **two batches** (A and B), with 3 controls and 3 treated samples in each batch.

    The data are simulated in the style of a real study (Himes et al., 2014): the gene names are real (8,000 human
    protein-coding genes), and the strongest induced genes include known glucocorticoid-response genes such as *FKBP5*,
    *TSC22D3* and *PER1*. The file `data/truth/rnaseq_truth.tsv` says which genes were planted as changed: use it only in
    the CHECK cells.

    ### B2.1 · Look at the data first  *(3 min)*

    **New terms**
    * **Count matrix** – one row per gene, one column per sample; each number is the reads assigned to that gene.
    * **Library size** – total reads of a sample. Samples are sequenced to different depths.
    * **Batch** – samples processed together (same day, kit or run). Batches can shift the measurements.
    """)

    nb.code(r"""
    # Step B2.1 · read the counts and the sample table
    counts = pd.read_csv(DATA / "rnaseq_counts.tsv.gz", sep="\t", index_col="gene")
    meta = pd.read_csv(DATA / "rnaseq_metadata.tsv", sep="\t")
    annot = counts[["gene_id", "chrom"]]
    counts = counts.drop(columns=["gene_id", "chrom"])[meta["sample"]]     # genes x samples, same order as meta
    print(counts.shape, "(genes, samples)")
    print(pd.crosstab(meta["condition"], meta["batch"]))
    counts.iloc[:5, :6]
    """)

    nb.code(r"""
    # Step B2.1b · library sizes (QC): do the samples have similar sequencing depth?
    libsize = counts.sum(axis=0) / 1e6
    fig, ax = plt.subplots(figsize=(8, 2.6))
    cols = [dviz.BLUE if c == "treated" else dviz.SLATE for c in meta["condition"]]
    ax.bar(meta["sample"], libsize.values, color=cols, width=0.65)
    ax.axvline(5.5, color=dviz.GREY, lw=1, ls=":")
    ax.text(2.5, libsize.max() * 1.02, "batch A", ha="center", fontsize=8)
    ax.text(8.5, libsize.max() * 1.02, "batch B", ha="center", fontsize=8)
    ax.set_ylabel("Reads (millions)")
    ax.set_title("Library size per sample (grey = control, blue = treated)")
    ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout()
    dviz.savefig(fig, "B2_library_sizes")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · 8,000 genes x 12 samples; balanced design (3 + 3 in each batch)
    dviz.check("B2.n_genes", counts.shape[0], 8000)
    dviz.check("B2.n_samples", counts.shape[1], 12)
    dviz.check("B2.balanced", bool((pd.crosstab(meta["condition"], meta["batch"]) == 3).all().all()), True)
    dviz.check("B2.libsize_ratio_max_min", round(float(libsize.max() / libsize.min()), 2), note="deepest / shallowest library")
    """)

    nb.md(r"""
    **What changed?** Nothing yet: we only looked. The deepest library has @@LIBRATIO@@ times more reads than the
    shallowest, so raw counts cannot be compared between samples. The design is **balanced**: each batch has three
    controls and three treated samples, so batch and treatment can be separated.

    ### B2.2 · Normalise: make the samples comparable  *(3 min)*

    **Why this step.** A gene with 200 reads in a deep library may be *less* expressed than a gene with 150 reads in a
    shallow one. DESeq2 divides each sample by a **size factor** computed with the *median-of-ratios* method: compare
    every gene with its average across samples, and take the median ratio per sample. The median ignores the few genes
    that really change.

    **New terms**
    * **Geometric mean** – the average on the log scale; robust to one very large value.
    * **log2(x + 1)** – a log scale on which a difference of 1 means a two-fold change; + 1 avoids log(0).
    """)

    nb.code(r"""
    # Step B2.2 · median-of-ratios size factors (the DESeq2 method) in four lines
    expressed = (counts > 0).all(axis=1)                         # genes with a count in every sample
    log_geo_mean = np.log(counts[expressed]).mean(axis=1)        # log of the "average sample"
    size_factors = np.exp((np.log(counts[expressed]).sub(log_geo_mean, axis=0)).median(axis=0))
    norm = counts / size_factors                                 # normalised counts
    log2n = np.log2(norm + 1)
    pd.DataFrame({"library size (M)": libsize.round(2), "size factor": size_factors.round(3)}).T
    """)

    nb.code(r"""
    # CHECK · size factors are centred near 1 and follow the library sizes (closely, not exactly)
    dviz.check("B2.sf_geomean", round(float(np.exp(np.log(size_factors).mean())), 3), between=(0.97, 1.03))
    dviz.check("B2.sf_vs_libsize_r", round(float(np.corrcoef(size_factors, libsize)[0, 1]), 2), between=(0.85, 1.0))
    dviz.check("B2.sf_range", f"{size_factors.min():.2f}-{size_factors.max():.2f}")
    """)

    nb.md(r"""
    **What changed?** Each sample was divided by its size factor (@@SFRANGE@@ here), so a normalised count of 200 now means
    the same in every sample. The size factors follow the library sizes closely (r = @@SFR@@) but not exactly: the
    median-of-ratios method also corrects for genes that change strongly and take up a larger share of the reads.

    ### B2.3 · Quality control with PCA: biology or batch?  *(4 min)*

    **Why this step.** Before testing any gene, check what makes the samples different overall. **PCA** (principal
    component analysis) finds the directions in which the samples differ most. If samples split by batch more than by
    treatment, every later result is suspect. We use the 1,000 most variable genes among genes with a mean normalised
    count of at least 10 (very low counts are mostly noise).

    **New terms**
    * **PC1, PC2** – the first and second directions of largest variance; **% variance** says how much of the total
      variation each one captures.
    """)

    nb.code(r"""
    # Step B2.3 · PCA on the 1,000 most variable genes (after removing low-count genes)
    from sklearn.decomposition import PCA
    keep = norm.mean(axis=1) >= 10
    top = log2n[keep].var(axis=1).nlargest(1000).index
    X = log2n.loc[top].T                                         # samples x genes
    pca = PCA(n_components=5).fit(X - X.mean())
    S = pca.transform(X - X.mean())
    vr = pca.explained_variance_ratio_ * 100
    fig, ax = plt.subplots(figsize=(4.8, 3.8))
    for cond, col in [("control", dviz.SLATE), ("treated", dviz.BLUE)]:
        for b, mk in [("A", "o"), ("B", "^")]:
            msk = ((meta["condition"] == cond) & (meta["batch"] == b)).to_numpy()
            ax.scatter(S[msk, 0], S[msk, 1], color=col, marker=mk, s=60, edgecolor="white", label=f"{cond}, batch {b}")
    ax.axhline(0, color=dviz.LIGHT, lw=1, zorder=0)
    ax.axvline(0, color=dviz.LIGHT, lw=1, zorder=0)
    ax.set_xlabel(f"PC1 ({vr[0]:.0f}% variance)")
    ax.set_ylabel(f"PC2 ({vr[1]:.0f}% variance)")
    ax.set_title("1,000 most variable genes")
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fig.tight_layout()
    dviz.savefig(fig, "B2_pca")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · what do PC1 and PC2 correspond to?
    treated = (meta["condition"] == "treated").to_numpy().astype(float)
    batch_b = (meta["batch"] == "B").to_numpy().astype(float)
    dviz.check("B2.pc1_pct", round(float(vr[0]), 1))
    dviz.check("B2.pc2_pct", round(float(vr[1]), 1))
    dviz.check("B2.pc1_vs_treatment_r", round(abs(float(np.corrcoef(S[:, 0], treated)[0, 1])), 2), between=(0.9, 1.0))
    dviz.check("B2.pc2_vs_batch_r", round(abs(float(np.corrcoef(S[:, 1], batch_b)[0, 1])), 2), between=(0.9, 1.0))
    """)

    nb.md(r"""
    **What changed?** Twelve samples with 1,000 genes each became twelve points in two dimensions. PC1 (@@PC1@@% of the
    variation) separates treated from control: treatment is the biggest effect, as we hope. PC2 (@@PC2@@%)
    separates batch A from batch B. That is a warning, not a result: the batch effect is real and must be in the model.

    **Your turn.** Colour by batch and use the marker shape for condition. Which version of the figure would you show in
    a paper about the drug, and which in a quality-control report?

    <details><summary>Answer</summary>

    For the drug paper, colour = condition (the message) and shape = batch (the check). For a QC report, colour = batch,
    because the question is "is there a batch effect?". Either way, two variables need two channels, and the shape keeps
    working in greyscale.
    </details>

    ### B2.4 · Test every gene, with and without batch in the model  *(5 min)*

    **Why this step.** For each gene we fit a small **linear model** on the log2 normalised expression:
    `expression = baseline + batch effect + treatment effect + noise`, and test whether the treatment effect is zero.
    Without the batch term, the batch differences end up in the "noise" and real changes look less certain.
    This is a simplified version of what limma or DESeq2 do; we write it out so that you can see each part.

    **New terms**
    * **p-value histogram** – if most genes do not change, their p-values are spread evenly between 0 and 1 (a flat
      floor); real changes add a spike near 0.
    * **Multiple testing** – with @@HANDTESTED@@ tests at p < 0.05, about @@EXPFALSE@@ genes pass by chance alone.
    * **Benjamini–Hochberg (BH)** – adjusts p-values so that among the genes you call, the expected fraction of false
      ones (the **false discovery rate, FDR**) is at most 5%.
    """)

    nb.code(r"""
    # Step B2.4a · one linear model per gene (all genes at once with matrix algebra)
    Y = log2n[keep].to_numpy().T                                  # samples x genes
    def treatment_p(design):
        beta, *_ = np.linalg.lstsq(design, Y, rcond=None)         # fit every gene at once
        resid = Y - design @ beta
        df = design.shape[0] - design.shape[1]                    # residual degrees of freedom
        se = np.sqrt((resid ** 2).sum(axis=0) / df * np.linalg.inv(design.T @ design)[-1, -1])
        t = beta[-1] / se                                         # the last column is "treated"
        return 2 * stats.t.sf(np.abs(t), df), beta[-1]
    one = np.ones(len(meta))
    p_with, lfc_with = treatment_p(np.column_stack([one, batch_b, treated]))
    p_without, _ = treatment_p(np.column_stack([one, treated]))
    print(f"{keep.sum():,} genes tested by hand; raw p < 0.05: with batch {np.sum(p_with < 0.05):,}, "
          f"without batch {np.sum(p_without < 0.05):,}")
    """)

    nb.code(r"""
    # Step B2.4b · p-value histograms: check that the test behaved before reading any result
    fig, axs = plt.subplots(1, 2, figsize=(9, 2.8), sharey=True)
    for ax, p, t in [(axs[0], p_with, "Model with batch"), (axs[1], p_without, "Model without batch")]:
        ax.hist(p, bins=40, range=(0, 1), color=dviz.SKY, edgecolor="white")
        ax.axhline(len(p) / 40 * 0.93, color=dviz.SLATE, ls="--", lw=1)
        ax.set_title(t)
        ax.set_xlabel("raw p-value")
    axs[0].set_ylabel("Genes")
    axs[0].text(0.55, len(p_with) / 40 * 1.25, "flat floor = genes that do not change", fontsize=7, color=dviz.SLATE)
    fig.tight_layout()
    dviz.savefig(fig, "B2_pvalue_histograms")
    plt.show()
    """)

    nb.code(r"""
    # Step B2.4c · Benjamini-Hochberg, then compare with the truth (allowed in a CHECK)
    truth = pd.read_csv(DATA / "truth" / "rnaseq_truth.tsv", sep="\t", index_col="gene")
    is_de = truth.loc[log2n.index[keep], "is_de"].to_numpy() == 1
    rows = []
    for name, p in [("with batch", p_with), ("without batch", p_without)]:
        called = dviz.bh_adjust(p) < 0.05
        rows.append({"model": name, "called (FDR 5%)": int(called.sum()),
                     "truly changed": int((called & is_de).sum()), "not planted": int((called & ~is_de).sum()),
                     "false discovery proportion": round(float((called & ~is_de).sum() / max(called.sum(), 1)), 3)})
    b2_tests = pd.DataFrame(rows)
    b2_tests
    """)

    nb.code(r"""
    # CHECK · BH kept the false fraction near 5%; adding batch found more real changes
    w, wo = b2_tests.iloc[0], b2_tests.iloc[1]
    dviz.check("B2.hand_called_with_batch", int(w["called (FDR 5%)"]))
    dviz.check("B2.hand_called_without_batch", int(wo["called (FDR 5%)"]))
    dviz.check("B2.hand_fdp_with_batch", float(w["false discovery proportion"]), between=(0.0, 0.08))
    dviz.check("B2.hand_true_gain_from_batch", int(w["truly changed"] - wo["truly changed"]), between=(1, 100))
    dviz.check("B2.expected_false_if_unadjusted", int(round(keep.sum() * 0.05)), note="tests x 0.05 if nothing changed")
    dviz.check("B2.hand_tested", int(keep.sum()))
    dviz.check("B2.raw_p05_with_batch", int((p_with < 0.05).sum()))
    dviz.check("B2.raw_p05_not_planted", int(((p_with < 0.05) & ~is_de).sum()))
    """)

    nb.md(r"""
    **What changed?** We now have one p-value per gene. Both histograms have a flat floor and a spike at 0, so the
    tests behaved. After BH, about 5% of the called genes were not planted, as BH promises. With batch in the model we
    find more truly changed genes, because the batch differences no longer inflate the noise. In a balanced design
    ignoring batch costs **power**. In an unbalanced design, for example all treated samples in one batch, it would also
    create **false** discoveries, and no model could rescue it.

    **Your turn.** How many genes have a raw p < 0.05 in the model with batch, and how many of those were not planted?
    Compare with `keep.sum() * 0.05`.

    <details><summary>Answer</summary>

    `np.sum(p_with < 0.05)` is @@RAWP05@@, of which @@RAWP05_NP@@ were not planted, close to the @@EXPFALSE@@ expected by chance
    (@@HANDTESTED@@ tests × 0.05). That is why we adjust: with raw p < 0.05, @@RAWP05_NP_PCT@@% of the "hits" would be false.
    </details>

    ### B2.5 · What the standard tool says: DESeq2  *(4 min)*

    **Why this step.** Counts are not normally distributed: low counts are noisy, and the spread grows with the mean.
    **DESeq2** models counts with a *negative binomial* distribution, normalises with the size factors from B2.2, borrows
    information across genes to estimate their variability, and **shrinks** the fold changes of uncertain genes toward 0.
    The instructor ran PyDESeq2 (the Python version of DESeq2) with the design `~ batch + condition`; the result is
    `data/de_results.tsv`. The optional cell below re-runs it if `pydeseq2` is installed (about 20 s).

    **New terms**
    * **log2FC_shrunk** – fold change after shrinkage (apeGLM): reliable genes keep their value, noisy ones move to 0.
    * **padj = NA** – the gene was **not tested** (all zeros, too few reads, or an outlier). Never plot NA as 0 or 1.
    """)

    nb.code(r"""
    # Step B2.5a · read the DESeq2 results and count up / down genes with BOTH thresholds
    de = pd.read_csv(DATA / "de_results.tsv", sep="\t")
    tested = de["padj"].notna()
    up = tested & (de["padj"] < 0.05) & (de["log2FC_shrunk"] > 1)
    dn = tested & (de["padj"] < 0.05) & (de["log2FC_shrunk"] < -1)
    fdr_only = tested & (de["padj"] < 0.05)
    de0 = pd.read_csv(DATA / "de_results_no_batch.tsv", sep="\t")
    print(f"tested {tested.sum():,} · not tested (padj NA) {(~tested).sum():,}")
    print(f"FDR < 5%: {fdr_only.sum():,} genes with batch in the design, {(de0['padj'] < 0.05).sum():,} without")
    print(f"FDR < 5% and |log2FC_shrunk| > 1: {up.sum():,} up, {dn.sum():,} down")
    """)

    nb.code(r"""
    # Step B2.5b (optional) · re-run DESeq2 yourself (needs: pip install pydeseq2)
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        RERUN = False                    # set to True to re-run (about 20 s)
        if RERUN:
            dds = DeseqDataSet(counts=counts.T, metadata=meta.set_index("sample")[["condition", "batch"]],
                               design="~batch + condition", quiet=True)
            dds.deseq2()
            ds = DeseqStats(dds, contrast=["condition", "treated", "control"], quiet=True)
            ds.summary()
            print((ds.results_df["padj"] < 0.05).sum(), "genes at FDR < 5% in your run")
    except ImportError:
        print("pydeseq2 is not installed: the provided de_results.tsv is used (that is fine).")
    """)

    nb.code(r"""
    # CHECK · counts, and how many highlighted genes were really planted
    planted = truth.loc[de["gene"], "is_de"].to_numpy() == 1
    dviz.check("B2.tested", int(tested.sum()), @@TESTED_N@@)
    dviz.check("B2.untested", int((~tested).sum()), @@UNTESTED_N@@)
    dviz.check("B2.up", int(up.sum()), @@UP_N@@)
    dviz.check("B2.down", int(dn.sum()), @@DOWN_N@@)
    dviz.check("B2.fdr05_with_batch", int(fdr_only.sum()), @@SIG_N@@)
    dviz.check("B2.fdr05_without_batch", int((de0["padj"] < 0.05).sum()), @@SIG0_N@@)
    dviz.check("B2.volcano_genes_not_planted", int(((up | dn) & ~planted).sum()), 0, tol=3)
    dviz.check("B2.fdr_only_not_planted_pct", round(100 * float((fdr_only & ~planted).sum() / fdr_only.sum()), 1),
               note="% of FDR-only genes that were not planted")
    dviz.check("B2.any_gene_up_and_down", bool((up & dn).any()), False)
    """)

    nb.md(r"""
    **What changed?** DESeq2 tested @@TESTED@@ genes; @@UNTESTED@@ were not tested (padj NA) and must be left out of the plot and
    counted in the legend. With batch in the design DESeq2 calls @@SIG@@ genes at FDR < 5%, against @@SIG0@@ without batch.
    With **both** thresholds (FDR < 5% and |log2FC| > 1) there are @@UP@@ up and @@DOWN@@ down, and @@VOLC_PLANTED_PCT@@% of them were planted.
    Among the FDR-only genes, @@FDRONLY_PCT@@% were not planted: most have small fold changes in highly expressed genes, where a
    tiny error in the estimated variability is enough to give a small p-value. Statistical significance is not biological
    importance; the fold-change threshold protects you here.

    ### B2.6 · Visualize: volcano and MA plots  *(4 min)*
    """)

    nb.code(r"""
    # Step B2.6a · volcano: grey first (draw order!), thresholds drawn, top genes labelled, counts in the legend
    x = de["log2FC_shrunk"]
    y = -np.log10(de["padj"])
    rest = tested & ~up & ~dn
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    ax.scatter(x[rest], y[rest], s=4, color=dviz.GREY, rasterized=True, label=f"not changed ({rest.sum():,})")
    ax.scatter(x[up], y[up], s=9, color=dviz.VERMILION, label=f"up ({up.sum()})")
    ax.scatter(x[dn], y[dn], s=9, color=dviz.OI_BLUE, label=f"down ({dn.sum()})")
    ax.axhline(-np.log10(0.05), color=dviz.SLATE, ls="--", lw=0.8)
    for v in (-1, 1):
        ax.axvline(v, color=dviz.SLATE, ls="--", lw=0.8)
    lim = np.ceil(np.nanmax(np.abs(x[tested])))
    ax.set_xlim(-lim, lim)
    top8 = de[tested].nsmallest(8, "padj")
    texts = [ax.text(r["log2FC_shrunk"], -np.log10(r["padj"]), r["gene"], fontsize=7, style="italic") for _, r in top8.iterrows()]
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, expand=(1.2, 1.6), arrowprops=dict(arrowstyle="-", color=dviz.SLATE, lw=0.5))
    except ImportError:
        pass
    ax.set_xlabel("Shrunken log2 fold change (treated vs control)")
    ax.set_ylabel("−log10 adjusted p (BH)")
    ax.set_title(f"FDR < 0.05 and |log2FC| > 1 · n = 6 vs 6")
    ax.legend(fontsize=7, loc="upper center", markerscale=1.5, frameon=True, facecolor="white", edgecolor="none", framealpha=1)
    fig.text(0.01, 0.005, f"{(~tested).sum()} genes not tested (padj NA) are not shown. DESeq2, design ~ batch + condition.",
             fontsize=6.5, color=dviz.SLATE)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    dviz.savefig(fig, "B2_volcano")
    plt.show()
    """)

    nb.code(r"""
    # Step B2.6b · MA plot: raw vs shrunken fold changes against mean expression
    A = np.log2(de["baseMean"] + 1)
    sig = up | dn
    fig, axs = plt.subplots(1, 2, figsize=(10, 3), sharey=True)
    for ax, col, t in [(axs[0], "log2FC", "Raw log2 fold change"), (axs[1], "log2FC_shrunk", "Shrunken log2 fold change")]:
        ax.scatter(A[tested & ~sig], de.loc[tested & ~sig, col], s=3, color=dviz.GREY, rasterized=True)
        ax.scatter(A[sig], de.loc[sig, col], s=5, color=dviz.BLUE)
        ax.axhline(0, color=dviz.SLATE, lw=0.8)
        ax.set_title(t)
        ax.set_xlabel("A = log2 mean normalised count")
    axs[0].set_ylabel("M = log2 fold change")
    axs[0].set_ylim(-6, 6)
    fig.tight_layout()
    dviz.savefig(fig, "B2_ma_plot")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · the labelled genes, and the known glucocorticoid-response genes among the "up" genes
    known = truth.index[truth["known_gr_target"] == 1]
    dviz.check("B2.top8_genes", ", ".join(top8["gene"]))
    dviz.check("B2.known_gr_up", int(de.loc[up, "gene"].isin(known).sum()), between=(12, len(known)), note=f"of {len(known)} known targets")
    low_a = tested & (A < 4)
    dviz.check("B2.shrinkage_low_counts", round(float(de.loc[low_a, "log2FC"].abs().median() / max(de.loc[low_a, "log2FC_shrunk"].abs().median(), 1e-6)), 1),
               note="raw / shrunken |log2FC| for genes with A < 4")
    """)

    nb.md(r"""
    **What changed?** The volcano now carries the evidence a reader needs to check the claim: adjusted p-values, shrunken
    fold changes, both thresholds drawn, colour-blind-safe colours, the counts, and a note about the untested genes.
    The most significant genes include the planted glucocorticoid-response genes (for example *FKBP5* and *TSC22D3*).
    In the MA plot, raw fold changes scatter widely at low counts (left), and shrinkage pulls those uncertain genes
    toward 0 (right).

    **Your turn.** Replace `padj` by `pvalue` on the y-axis and in the threshold. How many genes would you call up now?
    Are the extra genes planted?

    <details><summary>Answer</summary>

    More genes pass, and most of the extra ones are not planted. The y-axis must use the adjusted p-value when the
    threshold is an FDR.
    </details>

    ### B2.7 · Heatmap of the top 40 genes  *(optional, 3 min)*

    **Why this step.** A heatmap shows a pattern across genes and samples at once. Three choices must be stated:
    scaling (z-score each gene, so the colour shows the pattern rather than the expression level), the colour map
    (diverging, centred at 0, with a colour bar), and clustering (distance and linkage). Annotation bars show whether the
    clusters follow treatment or batch.
    """)

    nb.code(r"""
    # Step B2.7 · z-scored rows, correlation distance, average linkage, annotation bars
    from scipy.cluster.hierarchy import linkage, leaves_list
    from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
    top40 = de[tested].nsmallest(40, "padj")["gene"]
    Z = log2n.loc[top40]
    Z = Z.sub(Z.mean(axis=1), axis=0).div(Z.std(axis=1), axis=0)          # z-score each gene (row)
    r_order = leaves_list(linkage(Z.to_numpy(), "average", metric="correlation"))
    c_order = leaves_list(linkage(Z.T.to_numpy(), "average", metric="correlation"))
    Zo = Z.iloc[r_order, c_order]
    bwo = LinearSegmentedColormap.from_list("bwo", [dviz.DARKBLUE, dviz.SKY, "white", "#F2C185", dviz.VERMILION])
    fig = plt.figure(figsize=(7, 6.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.25, 0.25, 9], width_ratios=[20, 0.6], hspace=0.05, wspace=0.05)
    m2 = meta.set_index("sample").loc[Zo.columns]
    for k, (col, cmap) in enumerate([("condition", {"control": dviz.SLATE, "treated": dviz.BLUE}),
                                     ("batch", {"A": dviz.GREY, "B": dviz.ORANGE})]):
        a = fig.add_subplot(gs[k, 0])
        a.imshow([[matplotlib.colors.to_rgb(cmap[v]) for v in m2[col]]], aspect="auto")
        a.set_yticks([0], [col], fontsize=7)
        a.set_xticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
    ax = fig.add_subplot(gs[2, 0])
    im = ax.imshow(Zo, aspect="auto", cmap=bwo, norm=TwoSlopeNorm(0, -2.5, 2.5), interpolation="nearest")
    ax.set_yticks(range(len(Zo)), Zo.index, fontsize=6, style="italic")
    ax.set_xticks(range(Zo.shape[1]), Zo.columns, rotation=90, fontsize=7)
    fig.colorbar(im, cax=fig.add_subplot(gs[2, 1]), label="row z-score")
    fig.suptitle("Top 40 DE genes · z-scored rows · 1 − Pearson r, average linkage", x=0.08, ha="left", fontsize=9, weight="bold")
    dviz.savefig(fig, "B2_heatmap")
    plt.show()
    """)

    nb.code(r"""
    # CHECK · rows are scaled; the columns split first by treatment
    dviz.check("B2.heatmap_row_mean_max", round(float(Z.mean(axis=1).abs().max()), 6), 0.0, tol=1e-6)
    dviz.check("B2.heatmap_row_sd", round(float(Z.std(axis=1).mean()), 3), 1.0, tol=0.001)
    first_half = m2["condition"].iloc[:6]
    dviz.check("B2.heatmap_columns_split_by_condition", bool(first_half.nunique() == 1), True)
    """)

    nb.md(r"""
    **What changed?** The 40 genes form two blocks, up and down with treatment, and the column dendrogram separates
    treated from control first. The batch annotation bar is mixed within each block: for these genes, treatment
    dominates.

    ### B2 · Summary

    | Step | What you did | Number that proves it |
    |---|---|---|
    | B2.1 | looked at counts and design | 8,000 genes × 12 samples, 3 + 3 per batch |
    | B2.2 | median-of-ratios normalisation | size factors @@SFRANGE@@ |
    | B2.3 | PCA for QC | PC1 ≈ treatment (@@PC1@@%), PC2 ≈ batch (@@PC2@@%) |
    | B2.4 | @@HANDTESTED@@ tests, BH | @@HANDFDP_PCT@@% of called genes not planted; batch in the model finds @@HANDGAIN@@ more real changes |
    | B2.5 | DESeq2 results | @@TESTED@@ tested, @@UNTESTED@@ not tested; @@UP@@ up, @@DOWN@@ down |
    | B2.6 | volcano + MA | thresholds drawn, counts in the legend |

    **Presentation sentence (example).** *"Which genes respond to dexamethasone? @@UP@@ genes go up and @@DOWN@@ go down at
    FDR 5% and |log2FC| > 1 (DESeq2, 6 vs 6, batch in the model). Limitation: a batch effect explains @@PC2@@% of the
    variation; the design is balanced, so we could model it. The AI wrote the first plotting code; we checked every count
    against the notebook."*
    """)

    # ---------------------------------------------------------------------------------------------- Defend
    nb.md(r"""
    # Defend and present  *(10 min)*

    ## Step D1 · Verify: do the figure, the legend and the numbers agree?

    Every number in a figure or a legend must come from printed output. The cell below lists all CHECK values from this
    notebook. Compare them with your figure legends and with the reference table on the slide.
    """)

    nb.code(r"""
    # Step D1 · all CHECK values of this notebook
    res = pd.DataFrame(dviz._RESULTS).T[["value", "reference", "ok"]]
    print(f"{(res['ok'] == True).sum()} passed, {(res['ok'] == False).sum()} failed, {res['ok'].isna().sum()} recorded only")
    res
    """)

    nb.md(r"""
    ## Step D2 · Accessibility: colour vision and greyscale

    About 1 in 12 men (and 1 in 200 women) of Northern European ancestry have red–green colour-vision deficiency.
    `dviz.show_cvd()` shows a saved figure as it looks with deuteranopia (no green-sensitive cones) and in greyscale.
    Try it on your volcano and on the misleading volcano M6.
    """)

    nb.code(r"""
    # Step D2 · simulate deuteranopia and greyscale for your volcano and for M6
    for path in [FIG / "B2_volcano.png", ROOT / "bad_figures" / "M6.png"]:
        if path.exists():
            f = dviz.show_cvd(path, width=10)
            f.suptitle(path.name, x=0.01, ha="left", fontsize=9)
            plt.show()
    """)

    nb.md(r"""
    **What to look for.** In M6 the red and green groups become the same olive tone, so up and down can no longer be told
    apart. In your volcano, vermilion and blue stay distinct, and the threshold lines and counts still carry the message.
    If two groups look the same, switch to the Okabe–Ito colours **and** add a second channel (shape or direct labels).

    ## Step D3 · Words, AI and the log  *(optional)*

    1. Ask an AI assistant for **alt text (≤ 60 words)** and a **figure legend** using only your code and the printed
       numbers (not the image). Correct anything it could not know (biology, sample origin) and anything that does not
       match the printed numbers.
    2. Record what you did in `ai_use_log.md`: tool, date, prompt, what you accepted, what you changed, how you checked.

    ## Present one figure (2 minutes per pair)

    **Question** (one sentence) → **Figure** (explain the axes first) → **Takeaway** with a number → **Limitation** →
    **AI and checks** (what the AI did and how you verified it). Use the review checklist in `checklist.md`.
    """)

    nb.code(r"""
    # Save all CHECK values (used by the instructor) and the run time
    for ex in ["S0", "A", "B1", "B2"]:
        dviz.save_results(ex)
    print(f"Notebook time: {time.time() - T_START:.0f} s")
    """)
    nb.md, nb.code = _md, _code
