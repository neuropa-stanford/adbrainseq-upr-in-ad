#!/usr/bin/env python3
"""R1.1 rebuttal-ready figure — UPR/ERAD-core GO terms only, UP vs DOWN only.

Drops (i) the three generic parent terms (protein localization / transport / peptide
transport), which track the global transcriptional shift rather than UPR biology, and
(ii) the BOTH column, which is dominated by query size. Terms ranked by P value.

Panel A: heatmap, cell number = query genes in term, ns = not enriched.
Panel B: response to ER stress (GO:0034976) across the same 12 conditions.
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
rows = list(csv.DictReader(open(os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv"))))

NAMES = {
    "GO:0034976": "Response to ER stress",
    "GO:0030968": "ER unfolded protein response",
    "GO:0006986": "Response to unfolded protein",
    "GO:0034620": "Cellular response to unfolded protein",
    "GO:0035966": "Response to topologically incorrect protein",
    "GO:0035967": "Cellular response to topol. incorrect protein",
    "GO:0036498": "IRE1-mediated UPR",
    "GO:1905897": "Regulation of response to ER stress",
    "GO:1903573": "Neg. regulation of response to ER stress",
    "GO:0036503": "ERAD pathway",
    "GO:0006457": "Protein folding",
    "GO:0034975": "Protein folding in ER",
    "GO:0006888": "ER to Golgi vesicle-mediated transport",
}
DATASETS = ["Thapsigargin", "Mizuno", "Nativio"]
THRESH = [("no cutoff (sign only)", "no cutoff"), ("p<0.05", "p<0.05")]
DIRS = ["UP", "DOWN"]
cols = [(d, t[0], t[1], dr) for d in DATASETS for t in THRESH for dr in DIRS]
idx = {(r["GO"], r["dataset"], r["threshold"], r["direction"]): r for r in rows}

best = {}
for r in rows:
    if r["GO"] in NAMES and r["adjP"] and r["significant"] == "True" and r["direction"] != "BOTH":
        p = float(r["adjP"])
        best[r["GO"]] = min(best.get(r["GO"], 1.0), p if p > 0 else 1e-330)
ORDER = sorted(NAMES, key=lambda g: best.get(g, 2.0))

cmap = LinearSegmentedColormap.from_list(
    "upr", ["#fdf6ec", "#ffd9a0", "#f79f4a", "#e0562b", "#a01f28", "#5c0f1c"])
norm = Normalize(0, 16)

fig = plt.figure(figsize=(13.6, 9.4))
gs = fig.add_gridspec(2, 1, height_ratios=[2.9, 1.15], hspace=0.26,
                      left=0.375, right=0.955, top=0.865, bottom=0.135)

ax = fig.add_subplot(gs[0])
for i, go in enumerate(ORDER):
    for j, (ds, thr, _, dr) in enumerate(cols):
        r = idx.get((go, ds, thr, dr))
        if not r or not r["adjP"]:
            face, txt = "#ffffff", ""
        elif r["significant"] != "True":
            face, txt = "#eceae6", "ns"
        else:
            v = -math.log10(float(r["adjP"]))
            face, txt = cmap(norm(v)), r["genes_in_term"]
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", lw=1.2))
        if txt:
            dark = txt != "ns" and norm(-math.log10(float(r["adjP"]))) > .58
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=8.4,
                    color="#9a968f" if txt == "ns" else ("white" if dark else "#2b2b2b"))

ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
ax.set_yticklabels([f"{i+1}.  {NAMES[g]}  ({g})" for i, g in enumerate(ORDER)], fontsize=9.4)
ax.set_xticks(np.arange(len(cols)) + .5)
ax.set_xticklabels([c[3] for c in cols], fontsize=9,
                   fontweight="bold")
for t, (_, _, _, dr) in zip(ax.get_xticklabels(), cols):
    t.set_color("#c0392b" if dr == "UP" else "#2c5f8a")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
for j in range(2, len(cols), 2):
    ax.plot([j, j], [0, len(ORDER)], color="#4a4a4a" if j % 4 == 0 else "#c4c0b9",
            lw=2.4 if j % 4 == 0 else 1.1, zorder=5)
for j in range(0, len(cols), 2):
    ax.text(j + 1, -0.28, cols[j][2], ha="center", va="bottom", fontsize=9.4,
            style="italic", color="#333")
for k, ds in enumerate(DATASETS):
    ax.text(k * 4 + 2, -1.05, ds, ha="center", va="bottom", fontsize=13, fontweight="bold")
ax.text(2, -1.72, "positive control", ha="center", fontsize=9, color="#7a7a7a", style="italic")
ax.text(8, -1.72, "AD brain cohorts", ha="center", fontsize=9, color="#7a7a7a", style="italic")

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                  fraction=0.021, pad=0.013, extend="max")
cb.set_label("$-\\log_{10}$ adjusted $P$", fontsize=9.5)
cb.ax.tick_params(labelsize=8.5)
ax.set_title("a   UPR/ERAD gene sets are enriched among UP-regulated genes in the positive control\n"
             "     and exclusively among DOWN-regulated genes in both AD cohorts",
             loc="left", fontsize=12.8, fontweight="bold", pad=64)

# ------------------------------------------------------------------ panel b
ax2 = fig.add_subplot(gs[1])
vals, sigs, labs, colors = [], [], [], []
for (ds, thr, thr_lab, dr) in cols:
    r = idx.get(("GO:0034976", ds, thr, dr))
    sig = bool(r and r["significant"] == "True")
    vals.append(-math.log10(float(r["adjP"])) if sig else 0)
    sigs.append(sig); labs.append(f"{dr}\n{thr_lab}")
    colors.append("#c0392b" if dr == "UP" else "#2c5f8a")
x = np.arange(len(cols))
ax2.bar(x, vals, color=colors, width=.66, edgecolor="white")
ax2.set_ylim(0, max(vals) * 1.30)
for i, (v, s) in enumerate(zip(vals, sigs)):
    ax2.text(i, v + max(vals) * .02, ("%.1f" % v) if s else "ns", ha="center",
             va="bottom", fontsize=8.4, color="#333" if s else "#a8a29a")
ax2.axhline(-math.log10(0.05), color="#999", ls="--", lw=1)
ax2.text(len(cols) - .35, -math.log10(0.05) + .3, "P = 0.05", fontsize=7.6,
         color="#8a8a8a", ha="right")
ax2.set_xticks(x); ax2.set_xticklabels(labs, fontsize=8.4)
ax2.set_ylabel("$-\\log_{10}$ adjusted $P$", fontsize=9.5)
ax2.set_xlim(-.65, len(cols) - .35)
for j in range(4, len(cols), 4):
    ax2.axvline(j - .5, color="#4a4a4a", lw=1.6)
for k, ds in enumerate(DATASETS):
    ax2.text(k * 4 + 1.5, max(vals) * 1.17, ds, ha="center", fontsize=10.5, fontweight="bold")
for s in ["top", "right"]:
    ax2.spines[s].set_visible(False)
ax2.set_title("b   Response to endoplasmic reticulum stress (GO:0034976): the cutoff strengthens the "
              "positive control and abolishes the AD signal",
              loc="left", fontsize=12, fontweight="bold", pad=8)

fig.text(0.375, 0.018,
         "GO:BP over-representation (g:Profiler, g:SCS-adjusted, whole-genome background). Terms ranked by the "
         "strongest adjusted P reached in any condition; cell number = query genes in the term; ns = not enriched.\n"
         "Generic parent terms (protein localization/transport, peptide transport) and the direction-agnostic "
         "union are omitted here and reported in full in Supplementary Table SX. Inputs: SuppD1 (thapsigargin), "
         "SuppD2 (Mizuno, GSE173955), SuppD3 (Nativio, GSE159699).",
         fontsize=8, color="#666", va="bottom")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q1_UPR_GO_rebuttal_clean.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q1_UPR_GO_rebuttal_clean.{png,svg}")
