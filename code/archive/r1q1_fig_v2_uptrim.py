#!/usr/bin/env python3
"""R1.1 figure — trimmed variant. UP block = Thapsigargin only (the positive control, the only
UP-enriched dataset); DOWN block = the three AD cohorts (Nativio, Mizuno, Combined(M∩N)).
Panel b = representative conserved genes, sorted most-down first. Outputs *_v2_UPcontrol.*"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "axes.unicode_minus": False})
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
rows = list(csv.DictReader(open(os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv"))))
rep = list(csv.DictReader(open(os.path.join(OUT, "FIG_representative_conserved_UPR_genes.csv"))))

NAMES = {"GO:0034976": "Response to ER stress", "GO:0030968": "ER unfolded protein response",
    "GO:0006986": "Response to unfolded protein", "GO:0034620": "Cellular response to unfolded protein",
    "GO:0035966": "Response to topologically incorrect protein",
    "GO:0035967": "Cellular response to topol. incorrect protein", "GO:0036498": "IRE1-mediated UPR",
    "GO:1905897": "Regulation of response to ER stress", "GO:1903573": "Neg. regulation of response to ER stress",
    "GO:0036503": "ERAD pathway", "GO:0006457": "Protein folding", "GO:0034975": "Protein folding in ER",
    "GO:0006888": "ER to Golgi vesicle-mediated transport"}
THRESH = [("no cutoff (sign only)", "no cut"), ("p<0.05", "p<0.05")]
DS_ABBR = {"Thapsigargin": "Thaps", "Nativio": "Nativio", "Mizuno": "Mizuno", "Combined(M∩N)": "M∩N"}
# UP: Thapsigargin only  |  DOWN: the three AD cohorts
cols = [("Thapsigargin", t[0], t[1], "UP") for t in THRESH] \
     + [(d, t[0], t[1], "DOWN") for d in ["Nativio", "Mizuno", "Combined(M∩N)"] for t in THRESH]
SPLIT = 2                                                              # UP block ends at col 2
idx = {(r["GO"], r["dataset"], r["threshold"], r["direction"]): r for r in rows}

best = {}
for r in rows:
    if r["GO"] in NAMES and r["adjP"] and r["significant"] == "True" and r["direction"] != "BOTH":
        best[r["GO"]] = min(best.get(r["GO"], 1.0), float(r["adjP"]))
ORDER = sorted(NAMES, key=lambda g: best.get(g, 2.0))[:10]

cmapA = LinearSegmentedColormap.from_list("upr", ["#fdf6ec", "#ffd9a0", "#f79f4a", "#e0562b", "#a01f28", "#5c0f1c"])
normA = Normalize(0, 16)

fig = plt.figure(figsize=(12.0, 7.8))
gs = fig.add_gridspec(1, 2, width_ratios=[1.85, 1.0], left=0.30, right=0.955, top=0.74, bottom=0.09, wspace=0.42)

ax = fig.add_subplot(gs[0, 0])
for i, go in enumerate(ORDER):
    for j, (ds, thr, _, dr) in enumerate(cols):
        r = idx.get((go, ds, thr, dr))
        if not r or not r["adjP"]:
            face, txt = "#ffffff", ""
        elif r["significant"] != "True":
            face, txt = "#eceae6", "ns"
        else:
            v = -math.log10(float(r["adjP"])); face, txt = cmapA(normA(v)), r["genes_in_term"]
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", lw=1.1))
        if txt:
            dark = txt != "ns" and normA(-math.log10(float(r["adjP"]))) > .58
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=8.4,
                    color="#9a968f" if txt == "ns" else ("white" if dark else "#2b2b2b"))
ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
ax.set_yticklabels([f"{i+1}.  {NAMES[g]}  ({g})" for i, g in enumerate(ORDER)], fontsize=9)
ax.set_xticks(np.arange(len(cols)) + .5); ax.set_xticklabels([c[2] for c in cols], fontsize=7.4, rotation=90)
for t, (_, _, _, dr) in zip(ax.get_xticklabels(), cols):
    t.set_color("#c0392b" if dr == "UP" else "#2c5f8a")
ax.tick_params(length=0)
for s in ax.spines.values(): s.set_visible(False)
for j in range(2, len(cols), 2):
    split = (j == SPLIT)
    ax.plot([j, j], [0, len(ORDER)], color="#222" if split else "#c9c5be", lw=3.4 if split else 1.0, zorder=5)
for g in range(0, len(cols), 2):
    ax.text(g + 1, -0.42, DS_ABBR[cols[g][0]], ha="center", va="bottom", fontsize=9, fontweight="bold", color="#333")
ax.text(SPLIT / 2, -1.6, "UP", ha="center", va="bottom", fontsize=15, fontweight="bold", color="#c0392b")
ax.text((SPLIT + len(cols)) / 2, -1.6, "DOWN", ha="center", va="bottom", fontsize=15, fontweight="bold", color="#2c5f8a")
cb = fig.colorbar(plt.cm.ScalarMappable(norm=normA, cmap=cmapA), ax=ax, fraction=0.025, pad=0.015, extend="max")
cb.set_label("$-\\log_{10}$ adj. $P$", fontsize=9); cb.ax.tick_params(labelsize=8)
fig.text(0.30, 0.955, "a   UPR/ERAD GO terms — enriched only among UP genes in the positive control, "
         "only among DOWN genes in the AD cohorts", fontsize=11, fontweight="bold", ha="left", va="top")
fig.text(0.30, 0.925, "cell number = genes in term · ns = not enriched", fontsize=9, color="#666", ha="left", va="top")

axc = fig.add_subplot(gs[0, 1])
for r in rep: r["_mean"] = (float(r["Mizuno_log2FC"]) + float(r["Nativio_log2FC"])) / 2
rep_sorted = sorted(rep, key=lambda r: r["_mean"])
M = np.array([[float(r["Mizuno_log2FC"]), float(r["Nativio_log2FC"])] for r in rep_sorted])
cmapC = LinearSegmentedColormap.from_list("dn", ["#1b4f72", "#2c7fb8", "#a6cee3", "#f2f0ec"])
normC = Normalize(-0.7, 0)
for i, r in enumerate(rep_sorted):
    for j in range(2):
        axc.add_patch(Rectangle((j, i), 1, 1, facecolor=cmapC(normC(M[i, j])), edgecolor="white", lw=1))
        axc.text(j + .5, i + .5, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6.8,
                 color="white" if M[i, j] < -0.30 else "#2b2b2b")
axc.set_xlim(0, 2); axc.set_ylim(len(rep_sorted), 0)
axc.set_yticks(np.arange(len(rep_sorted)) + .5); axc.set_yticklabels([r["gene"] for r in rep_sorted], fontsize=7.6)
axc.set_xticks([.5, 1.5]); axc.set_xticklabels(["Mizuno", "Nativio"], fontsize=8.5)
axc.tick_params(length=0)
for s in axc.spines.values(): s.set_visible(False)
for i, r in enumerate(rep_sorted):
    axc.text(2.08, i + .5, r["branch_function"], va="center", fontsize=6.3, color="#777")
axc.set_title("b   Representative conserved genes\n     (log$_2$FC, sorted most-down first)", loc="left", fontsize=10.4, fontweight="bold", pad=8)
cbc = fig.colorbar(plt.cm.ScalarMappable(norm=normC, cmap=cmapC), ax=axc, orientation="horizontal", fraction=0.05, pad=0.10, aspect=24)
cbc.set_label("log$_2$FC (negative = down)", fontsize=7.5); cbc.ax.tick_params(labelsize=7)

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q1_Figure_GO_consolidated_v2_UPcontrol.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q1_Figure_GO_consolidated_v2_UPcontrol.{png,svg}")
