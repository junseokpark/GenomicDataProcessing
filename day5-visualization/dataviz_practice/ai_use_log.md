# AI-use log

AI is **optional** in this practice. If you use an assistant for any figure — code, a critique, alt text, a legend —
add one entry per figure. A journal may ask for this information (Springer Nature: document LLM use in the Methods).

Rules (lecture Part 5):
1. **Verify** — every number on the figure comes from printed, checked output (`CHECK ·` cells), not from the assistant's text.
2. **Disclose** — say what the assistant did.
3. **No generative-AI images of data or results.**
4. **Protect data** — never paste patient data or controlled-access genomic data into public tools. (The practice data are simulated, so they are safe to share.)

## Entries

| Figure (file) | Tool + version | Date | Prompt (card # or text) | Accepted | Changed by you | How you verified |
|---|---|---|---|---|---|---|
| *Example:* `results/figures/B2_volcano.png` | *Assistant X, version* | *2026-10-27* | *Prompt card 2, then "labels overlap; use adjustText"* | *Plotting code; colours; label placement* | *It used `pvalue`; changed to `padj`. It filled NA with 0; dropped untested genes and counted them* | *Up/down counts match CHECK B2.volcano (177 up, 106 down); CVD simulation; checklist* |
| | | | | | | |
| | | | | | | |

## Alt text and legend (if the assistant drafted them)

* **Alt text (≤ 60 words):**
* **Legend:**
* **What you corrected** (things the assistant could not know from the code — biology, sample origin — or that did not match the printed numbers):
