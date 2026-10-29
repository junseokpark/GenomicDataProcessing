# Prompts for review, accessibility and refinement

These are optional. Use them with any assistant. Give the assistant the **code and the printed numbers**, not only the
image: models misread values and colours from pictures (lecture Part 5).

## 1. Critique prompt (Exercise 1 and your own figure)

Why the definitions? Models detected misleading charts only moderately well when asked "is this figure good?", and
considerably better when the prompt defined the misleaders (Alexander et al., 2024).

```text
Review this figure as a strict reviewer. Check for these misleaders and quote the evidence for each:
- truncated axis (bars not starting at zero) or dual y-axes
- bar of means hiding n or the individual values
- rainbow / jet or red–green colour maps
- overplotting or a draw order that hides points
- missing units, n, thresholds or multiple-testing correction
- distance or size claims from a UMAP / t-SNE
- angle or area used for a precise comparison (pie, bubble)
List the problems by severity, each with a one-line fix.
Say "not visible" rather than guessing any value you cannot read.
```

After it answers, compare with `solutions/ANSWER_KEY_bad_figures.md`: what did it find that you missed, and what did it
**invent**?

## 2. Alt text and legend (Defend · D2)

```text
Here is the code that made my figure and the numbers it printed:
<paste the code cell>
<paste the printed CHECK output>
Write (1) alt text of at most 60 words: chart type, what is on each axis, the main pattern with one number;
(2) a figure legend: what each mark, colour and line means, n and its unit, the test and correction, thresholds.
Use only information in the code and the printed numbers. Write "unknown" for anything you cannot tell from them.
```

Then correct anything it could not know (biology, sample origin) and anything that does not match the printed numbers.

## 3. One change at a time (refinement loop)

Render → look → ask for **one** change → render again. One change at a time shows you what each change did and lets you
catch the assistant breaking something else.

```text
The figure is almost right. Make only this change and keep everything else identical: <one change>.
```

Examples: "labels overlap; use adjustText"; "start the y-axis at 0"; "draw cells with expression > 0 last";
"add n under each group label"; "use marker shape for batch as well as colour".

Do **not** ask the assistant to make a result "look stronger" or "more separated": that changes the picture, not the
analysis.

## 4. Agent-skill prompts (if you use a coding agent with skills installed)

Read any skill's `SKILL.md` before installing it; install only what you need.

* **scientific-visualization** — "Show every sample, mean ± SD, n = 3 independent experiments, name the test, 89 mm SVG."
* **scientific-critical-thinking** — "List the questions my PI will ask about this figure."
* **diagram-design** — "RNA-seq workflow, sample → counts → DESeq2 → figures, left to right, one accent on DESeq2, SVG."

Never use an image-generation model for figures of data or results.
