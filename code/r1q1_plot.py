#!/usr/bin/env python3
"""Figure for R1.1 / Jonathan's request: 18 UPR-related GO terms x
(Thapsigargin | Mizuno | Nativio) x (no cutoff | p<0.05) x (UP | DOWN | BOTH).

Panel A: full heatmap, colour = -log10(adjusted P), cell text = genes in term.
Panel B: the headline contrast — "response to ER stress" across all 18 conditions.
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

OUT = "/data/adbrainseq/Stanford U/ClaudeAgentwithWIKILLM/subprojects/adbrainseq-manuscript/R1Q1_figure_data"
CSV = os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv")

rows = list(csv.DictReader(open(CSV)))

# ---- row order: UPR-core first, then folding/ERAD, then generic trafficking
ORDER = [
    ("GO:0034976", "Response to ER stress"),
    ("GO:0030968", "ER unfolded protein response"),
    ("GO:0006986", "Response to unfolded protein"),
    ("GO:0034620", "Cellular response to unfolded protein"),
    ("GO:0035966", "Response to topologically incorrect protein"),
    ("GO:0035967", "Cellular response to topol. incorrect protein"),
    ("GO:0036498", "IRE1-mediated UPR"),
    ("GO:1905897", "Regulation of response to ER stress"),
    ("GO:1903573", "Neg. regulation of response to ER stress"),
    ("GO:0036503", "ERAD pathway"),
    ("GO:0030433", "Ubiquitin-dependent ERAD pathway"),
    ("GO:0006457", "Protein folding"),
    ("GO:0034975", "Protein folding in ER"),
    ("GO:0006888", "ER to Golgi vesicle-mediated transport"),
    ("GO:0042886", "Amide transport"),
    ("GO:0015833", "Peptide transport"),
    ("GO:0015031", "Protein transport"),
    ("GO:0008104", "Protein localization"),
]
GROUP_SPLIT = 11
OBSOLETE = {"GO:0030433", "GO:0042886"}  # not present in the current g:Profiler GO:BP release          # after "Ubiquitin-dependent ERAD" = UPR/ERAD core
DATASETS = ["Thapsigargin", "Mizuno", "Nativio"]
THRESH = [("no cutoff (sign only)", "no cutoff"), ("p<0.05", "p<0.05")]
DIRS = ["UP", "DOWN", "BOTH"]

cols = [(d, t[0], t[1], dr) for d in DATASETS for t in THRESH for dr in DIRS]
idx = {(r["GO"], r["dataset"], r["threshold"], r["direction"]): r for r in rows}

P = np.full((len(ORDER), len(cols)), np.nan)
N = np.zeros_like(P)
SIG = np.zeros_like(P, dtype=bool)
for i, (go, _) in enumerate(ORDER):
    for j, (ds, thr, _, dr) in enumerate(cols):
        r = idx.get((go, ds, thr, dr))
        if r and r["adjP"]:
            p = float(r["adjP"])
            P[i, j] = min(-math.log10(p), 90) if p > 0 else 90
            N[i, j] = int(r["genes_in_term"] or 0)
            SIG[i, j] = r["significant"] == "True"

cmap = LinearSegmentedColormap.from_list(
    "upr", ["#f7f4ef", "#ffd9a0", "#f79f4a", "#e0562b", "#a01f28", "#5c0f1c"])
norm = Normalize(0, 60)

fig = plt.figure(figsize=(17.5, 11.2))
gs = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.15], hspace=0.20,
                      left=0.335, right=0.965, top=0.905, bottom=0.075)

# ============================================================ Panel A heatmap
ax = fig.add_subplot(gs[0])
ax.set_facecolor("white")
for i in range(len(ORDER)):
    for j in range(len(cols)):
        v = P[i, j]
        if np.isnan(v):
            face, txt = ("#f4f1ec" if ORDER[i][0] in OBSOLETE else "#ffffff"), ""
        elif not SIG[i, j]:
            face, txt = "#eceae6", "ns"
        else:
            face, txt = cmap(norm(v)), f"{int(N[i,j])}"
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face,
                               edgecolor="white", linewidth=1.1))
        if txt:
            lum = 0 if not SIG[i, j] else norm(v)
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=7.4,
                    color="#9a968f" if txt == "ns" else ("white" if lum > .55 else "#2b2b2b"))

ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
ax.set_yticklabels([f"{n}  ({g})" for g, n in ORDER], fontsize=9)
ax.set_xticks(np.arange(len(cols)) + .5)
ax.set_xticklabels([c[3] for c in cols], fontsize=8.2)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)

# group separators / headers
for j in range(3, len(cols), 3):
    lw = 2.6 if j % 6 == 0 else 1.2
    ax.plot([j, j], [0, len(ORDER)], color="#4a4a4a" if j % 6 == 0 else "#b9b5ae",
            lw=lw, zorder=5)
ax.plot([0, len(cols)], [GROUP_SPLIT, GROUP_SPLIT], color="#4a4a4a", lw=2.0, zorder=5)

for k, (thr_key, thr_lab) in enumerate(THRESH * 3):
    pass
for j in range(0, len(cols), 3):
    thr_lab = cols[j][2]
    ax.text(j + 1.5, -0.32, thr_lab, ha="center", va="bottom", fontsize=9.2,
            style="italic", color="#333")
for k, ds in enumerate(DATASETS):
    ax.text(k * 6 + 3, -1.15, ds, ha="center", va="bottom", fontsize=12.5,
            fontweight="bold")
ax.text(3, -1.85, "positive control (acute UPR activation)", ha="center",
        fontsize=8.6, color="#7a7a7a", style="italic")
ax.text(12, -1.85, "AD brain cohorts", ha="center", fontsize=8.6,
        color="#7a7a7a", style="italic")


ax.annotate("", xy=(-0.345, 0.02), xytext=(-0.345, GROUP_SPLIT - 0.02),
            xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
            arrowprops=dict(arrowstyle="-", lw=2.4, color="#8c8c8c"))
ax.text(-0.365, GROUP_SPLIT / 2, "UPR / ERAD core", rotation=90, ha="center",
        va="center", fontsize=10, fontweight="bold", color="#555",
        transform=ax.get_yaxis_transform())
ax.annotate("", xy=(-0.345, GROUP_SPLIT + 0.02), xytext=(-0.345, len(ORDER) - 0.02),
            xycoords=("axes fraction", "data"), textcoords=("axes fraction", "data"),
            arrowprops=dict(arrowstyle="-", lw=2.4, color="#8c8c8c"))
ax.text(-0.365, (GROUP_SPLIT + len(ORDER)) / 2, "folding /\ntrafficking", rotation=90,
        ha="center", va="center", fontsize=10, fontweight="bold", color="#555",
        transform=ax.get_yaxis_transform())
for i, (g, n) in enumerate(ORDER):
    if g in OBSOLETE:
        ax.text(len(cols) / 2, i + .5, "term obsolete / merged in current GO release "
                "(not testable)", ha="center", va="center", fontsize=8.4,
                style="italic", color="#a09a92")

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                  fraction=0.018, pad=0.012, extend="max")
cb.set_label("$-\\log_{10}$ adjusted $P$ (g:SCS)", fontsize=9.5)
cb.ax.tick_params(labelsize=8.5)
ax.set_title("A   UPR-related GO terms across datasets, DEG thresholds and directions"
             "   —   cell number = genes in term, ns = not enriched",
             loc="left", fontsize=12.5, fontweight="bold", pad=52)

# ============================================================ Panel B bars
ax2 = fig.add_subplot(gs[1])
go_head = "GO:0034976"
vals, labs, colors, sigs = [], [], [], []
for (ds, thr, thr_lab, dr) in cols:
    r = idx.get((go_head, ds, thr, dr))
    p = float(r["adjP"]) if r and r["adjP"] else np.nan
    s = bool(r and r["significant"] == "True")
    vals.append(min(-math.log10(p), 90) if p and p > 0 and s else 0)
    sigs.append(s)
    labs.append(f"{dr}\n{thr_lab}")
    colors.append({"UP": "#c0392b", "DOWN": "#2c5f8a", "BOTH": "#6b6b6b"}[dr])
x = np.arange(len(cols))
ax2.bar(x, vals, color=colors, width=.72, edgecolor="white")
for i, (v, s) in enumerate(zip(vals, sigs)):
    ax2.text(i, v + .5, ("%.0f" % v) if s else "ns", ha="center", va="bottom",
             fontsize=8, color="#333" if s else "#999")
ax2.axhline(-math.log10(0.05), color="#999", ls="--", lw=1)
ax2.text(len(cols) - .3, -math.log10(0.05) + .4, "P = 0.05", fontsize=7.5,
         color="#888", ha="right")
ax2.set_xticks(x); ax2.set_xticklabels(labs, fontsize=7.6)
ax2.set_ylabel("$-\\log_{10}$ adj. $P$", fontsize=9.5)
ax2.set_xlim(-.7, len(cols) - .3)
for j in range(6, len(cols), 6):
    ax2.axvline(j - .5, color="#4a4a4a", lw=1.6)
ax2.set_ylim(0, max(vals) * 1.22)
for k, ds in enumerate(DATASETS):
    ax2.text(k * 6 + 2.5, max(vals) * 1.12, ds, ha="center",
             fontsize=10.5, fontweight="bold")
for s in ["top", "right"]:
    ax2.spines[s].set_visible(False)
ax2.set_title("B   Headline term: response to endoplasmic reticulum stress (GO:0034976)",
              loc="left", fontsize=12.5, fontweight="bold", pad=8)

fig.text(0.335, 0.014,
         "g:Profiler (GO:BP, g:SCS-adjusted, whole-genome background, all_results). "
         "Input gene lists from SuppD1 (thapsigargin, DESeq2), SuppD2 (Mizuno, GSE173955), "
         "SuppD3 (Nativio, GSE159699; Welch t-test on counts).\n"
         "UP/DOWN = genes changed in that direction; BOTH = union (direction-agnostic) and is "
         "therefore dominated by query size — interpret directional columns, not BOTH.",
         fontsize=7.6, color="#666", va="bottom")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q1_UPR_GO_matrix_18terms.{ext}"),
                dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved to", OUT)
