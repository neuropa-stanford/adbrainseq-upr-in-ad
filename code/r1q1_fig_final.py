#!/usr/bin/env python3
"""R1.1 figure (team-meeting revision).
Panel a: UPR/ERAD GO terms (top 10 by significance) x [Thapsigargin | Nativio | Mizuno | Combined(M∩N)],
         UP columns grouped then DOWN columns, no-cutoff & p<0.05. Cell = genes in term; colour = -log10 adjP.
Panel b: representative conserved genes, concordantly down in both cohorts, SORTED by fold-change.
(old Panel B — the ER-stress bars — removed per team feedback.)
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "axes.unicode_minus": False,
    "pdf.fonttype": 42, "ps.fonttype": 42,   # embed TrueType (crisp), not Type-3 bitmap
})
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
rows = list(csv.DictReader(open(os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv"))))
rep = list(csv.DictReader(open(os.path.join(OUT, "FIG_representative_conserved_UPR_genes.csv"))))

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
DATASETS = ["Thapsigargin", "Nativio", "Mizuno", "Combined(M∩N)"]      # Nativio before Mizuno
DS_SHORT = {"Thapsigargin": "Thapsigargin", "Nativio": "Nativio", "Mizuno": "Mizuno",
            "Combined(M∩N)": "Combined (M∩N)"}
THRESH = [("no cutoff (sign only)", "no cut"), ("p<0.05", "p<0.05")]
DIRS = ["UP", "DOWN"]
# ALL UP columns first (across datasets), then ALL DOWN columns — split the table into UP | DOWN blocks
cols = [(d, t[0], t[1], dr) for dr in DIRS for d in DATASETS for t in THRESH]
idx = {(r["GO"], r["dataset"], r["threshold"], r["direction"]): r for r in rows}

best = {}
for r in rows:
    if r["GO"] in NAMES and r["adjP"] and r["significant"] == "True" and r["direction"] != "BOTH":
        best[r["GO"]] = min(best.get(r["GO"], 1.0), float(r["adjP"]))
ORDER = sorted(NAMES, key=lambda g: best.get(g, 2.0))[:9]               # top 9 (term 10 = Protein folding in ER: not enriched in any DOWN cohort — dropped)

cmapA_up = LinearSegmentedColormap.from_list(
    "up", ["#fdf6ec", "#ffd9a0", "#f79f4a", "#e0562b", "#a01f28", "#5c0f1c"])   # UP = warm/red
cmapA_dn = LinearSegmentedColormap.from_list(
    "dn", ["#f0f5fb", "#a6cee3", "#4a90c2", "#2c5f8a", "#173a5c"])              # DOWN = blue (matches panel b)
normA = Normalize(0, 16)
PANEL_A_LETTER, PANEL_B_LETTER = "A", "B"   # panel labels (clean, matching the manuscript style)

fig = plt.figure(figsize=(8.0, 4.5))
gs = fig.add_gridspec(1, 2, width_ratios=[3.05, 1.05],
                      left=0.235, right=0.925, top=0.82, bottom=0.11, wspace=0.95)

# =========================================================== Panel a
ax = fig.add_subplot(gs[0, 0])
for i, go in enumerate(ORDER):
    for j, (ds, thr, _, dr) in enumerate(cols):
        r = idx.get((go, ds, thr, dr))
        if not r or not r["adjP"]:
            face, txt = "#ffffff", ""
        elif r["significant"] != "True":
            face, txt = "#eceae6", "ns"
        else:
            v = -math.log10(float(r["adjP"]))
            face, txt = (cmapA_up if dr == "UP" else cmapA_dn)(normA(v)), r["genes_in_term"]
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", lw=1.1))
        if txt:
            dark = txt != "ns" and normA(-math.log10(float(r["adjP"]))) > .58
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=6.6,
                    color="#9a968f" if txt == "ns" else ("white" if dark else "#2b2b2b"))
ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
# two-line GO term labels: name on top, GO id underneath
ax.set_yticklabels([f"{NAMES[g]}\n({g})" for g in ORDER], fontsize=8.3)
# bottom x-tick = threshold, coloured by direction
ax.set_xticks(np.arange(len(cols)) + .5)
ax.set_xticklabels([c[2] for c in cols], fontsize=7.3, rotation=90)
for t, (_, _, _, dr) in zip(ax.get_xticklabels(), cols):
    t.set_color("#c0392b" if dr == "UP" else "#2c5f8a")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
# dividers: medium between datasets (every 2 cols); THICK at the central UP | DOWN split
half = len(cols) // 2
for j in range(2, len(cols), 2):
    split = (j == half)
    ax.plot([j, j], [0, len(ORDER)], color="#222" if split else "#c9c5be",
            lw=3.4 if split else 1.0, zorder=5)
# dataset labels (each spans 2 cols) — abbreviated, repeated under UP and DOWN blocks
DS_ABBR = {"Thapsigargin": "Thaps", "Nativio": "Nativio", "Mizuno": "Mizuno", "Combined(M∩N)": "M∩N"}
# dataset-identity colours matching the bottom violins: Mizuno = GSE173955 green, Nativio = GSE159699 purple,
# M∩N (Mizuno∩Nativio) = blend of the two. Thaps (in-vitro reference) keeps its direction colour.
DS_COLOR = {"Mizuno": "#82D06F", "Nativio": "#838EBB", "Combined(M∩N)": "#82AF95"}
for g in range(0, len(cols), 2):
    ds = cols[g][0]
    col = DS_COLOR.get(ds, "#c0392b" if cols[g][3] == "UP" else "#2c5f8a")
    ax.text(g + 0.25, -0.28, DS_ABBR[ds], ha="left", va="bottom", rotation=30, fontsize=8,
            fontweight="bold", color=col)
# UP / DOWN big block labels (each spans half the table)
ax.text(half / 2, -1.6, "UP", ha="center", va="bottom", fontsize=15, fontweight="bold", color="#c0392b")
ax.text(half + half / 2, -1.6, "DOWN", ha="center", va="bottom", fontsize=15, fontweight="bold", color="#2c5f8a")
# two same-size colorbars side by side in the gap: DOWN(blue) has NO ticks/label; UP(red) carries the numbers + label
# order matches panel a: UP(red) LEFT, DOWN(blue) RIGHT. Numbers only on UP, placed on its OUTER (left) side.
cax_up = fig.add_axes([0.628, 0.30, 0.013, 0.24])   # red UP moved right, close to DOWN
cb1 = fig.colorbar(plt.cm.ScalarMappable(norm=normA, cmap=cmapA_up), cax=cax_up, extend="max")
cb1.ax.yaxis.set_ticks_position("left"); cb1.ax.tick_params(labelsize=7)
cax_dn = fig.add_axes([0.650, 0.30, 0.013, 0.24])
cb2 = fig.colorbar(plt.cm.ScalarMappable(norm=normA, cmap=cmapA_dn), cax=cax_dn, extend="max")
cb2.ax.tick_params(left=False, right=False, labelleft=False, labelright=False)  # blue: colour only
# angled UP/DOWN labels so the two close bars don't collide
fig.text(0.632, 0.552, "UP", fontsize=8.4, color="#c0392b", fontweight="bold", rotation=30, ha="left", va="bottom")
fig.text(0.653, 0.552, "DOWN", fontsize=8.4, color="#2c5f8a", fontweight="bold", rotation=30, ha="left", va="bottom")
fig.text(0.640, 0.285, "$-$Log$_{10}$ adj. $P$", fontsize=7.5, ha="center", va="top")
fig.text(0.012, 0.988, PANEL_A_LETTER, fontsize=18.3, fontweight="bold", ha="left", va="top")
# NOTE: "cell number = genes in term; ns = not enriched" moved to the FIGURE LEGEND (removed from panel).

# =========================================================== Panel b (conserved genes: grouped by UPR gene group, gradient within group)
axc = fig.add_subplot(gs[0, 1])
def group_of(bf):
    b = bf.lower()
    if b.startswith("perk") or b.startswith("isr"): return "PERK / ISR"
    if b.startswith("ire1"): return "IRE1"
    if b.startswith("atf6"): return "ATF6"
    if b.startswith("erad"): return "ERAD"
    if b.startswith("chaperone"): return "Chaperone"
    if b.startswith("pdi"): return "PDI"
    return "Other"
for r in rep:
    r["_mean"] = (float(r["Mizuno_log2FC"]) + float(r["Nativio_log2FC"])) / 2
    r["_grp"] = group_of(r["branch_function"])
groups = {}
for r in rep: groups.setdefault(r["_grp"], []).append(r)
# groups ordered most-down first (by group mean); genes within a group most-down first
gorder = sorted(groups, key=lambda g: np.mean([x["_mean"] for x in groups[g]]))
cmapC = LinearSegmentedColormap.from_list("dn", ["#1b4f72", "#2c7fb8", "#a6cee3", "#f2f0ec"])
normC = Normalize(-0.7, 0)
GAP = 0.85
layout, spans, y = [], [], 0.0
for g in gorder:
    members = sorted(groups[g], key=lambda r: r["_mean"])
    y0 = y
    for r in members:
        layout.append((y, r)); y += 1
    spans.append((g, y0, y)); y += GAP
total = y - GAP
for yy, r in layout:
    vals = [float(r["Mizuno_log2FC"]), float(r["Nativio_log2FC"])]
    for j, v in enumerate(vals):
        axc.add_patch(Rectangle((j, yy), 1, 1, facecolor=cmapC(normC(v)), edgecolor="white", lw=1))
        axc.text(j + .5, yy + .5, f"{v:.2f}", ha="center", va="center", fontsize=6.8,
                 color="white" if v < -0.30 else "#2b2b2b")
axc.set_xlim(0, 2); axc.set_ylim(total, 0)
axc.set_yticks([yy + .5 for yy, _ in layout])
axc.set_yticklabels([r["gene"] for _, r in layout], fontsize=7.6)
axc.set_xticks([.5, 1.5]); axc.set_xticklabels(["Mizuno", "Nativio"], fontsize=8.5)
axc.tick_params(length=0)
for s in axc.spines.values():
    s.set_visible(False)
# group brackets + labels on the right; faint sub-function under each gene name is dropped for clarity
for g, y0, y1 in spans:
    axc.plot([2.10, 2.10], [y0 + 0.08, y1 - 0.08], color="#8a8a8a", lw=1.4, clip_on=False)
    axc.text(2.20, (y0 + y1) / 2, g, va="center", ha="left", fontsize=7.4,
             color="#444", fontweight="bold", clip_on=False)
axc.set_title(PANEL_B_LETTER, loc="left", fontsize=18.3, fontweight="bold", pad=8)
cbc = fig.colorbar(plt.cm.ScalarMappable(norm=normC, cmap=cmapC), ax=axc,
                   orientation="horizontal", fraction=0.05, pad=0.10, aspect=24)
cbc.set_label("Log$_2$fc", fontsize=7.5)
cbc.ax.tick_params(labelsize=7)

for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUT, f"R1Q1_Figure_GO_consolidated.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q1_Figure_GO_consolidated.{png,svg}")
print("Panel a terms (1-10):")
for i, g in enumerate(ORDER, 1): print(f"  {i}. {NAMES[g]} ({g})")
