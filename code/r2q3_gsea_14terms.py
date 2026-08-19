#!/usr/bin/env python3
"""R2.3 figure — GSEA NES heatmap, single panel, terms 1-14 (clean; no leading numbers).
Same 18 GO terms are tested (BH across all 18); only the 14 testable UPR/ERAD-related
terms are displayed. Panel b (ER-stress bar chart) removed per author request.
Rank metric = sign(log2FC) x -log10(P) (Nativio: Welch t). Output: R2Q3_GSEA_14terms.{png,svg}
"""
import csv, os, textwrap
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
rows = list(csv.DictReader(open(os.path.join(OUT, "R1Q1_UPR_GSEA_18terms.csv"))))

CORE = {"GO:0034976", "GO:0030968", "GO:0006986", "GO:0034620", "GO:0035966",
        "GO:0035967", "GO:0036498", "GO:1905897", "GO:1903573", "GO:0036503",
        "GO:0030433"}
DATASETS = ["Thapsigargin", "Mizuno", "Nativio"]
VARIANTS = [("full ranked list", "full list"), ("p<0.05 subset", "p<0.05 subset")]
cols = [(d, v[0], v[1]) for d in DATASETS for v in VARIANTS]
idx = {(r["GO"], r["dataset"], r["variant"]): r for r in rows}
NAMES = {r["GO"]: r["term"] for r in rows}

best = {}
for r in rows:
    if r["padj"] and r["significant"] == "True":
        best[r["GO"]] = min(best.get(r["GO"], 1.0), float(r["padj"]))
ORDER = sorted(NAMES, key=lambda g: (float(idx[(g, "Thapsigargin", "full ranked list")]["K"] or 0) < 10,
                                     best.get(g, 2.0)))
ORDER = ORDER[:14]   # show terms 1-14 only

cmap = LinearSegmentedColormap.from_list(
    "nes", ["#1b4f72", "#2c7fb8", "#a6cee3", "#f2f0ec", "#f8b878", "#e0562b", "#8c1c13"])
norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)

fig = plt.figure(figsize=(12.8, 6.6))
ax = fig.add_axes([0.40, 0.175, 0.55, 0.66])   # leave room: left for labels, top for headers, bottom for footnote

for i, go in enumerate(ORDER):
    for j, (ds, var, _) in enumerate(cols):
        r = idx.get((go, ds, var))
        untestable = (not r) or (not r["NES"]) or r["NES"] == "nan" or int(r["K"] or 0) < 10
        if untestable:
            face, txt, col = "#f4f1ec", "untestable", "#a09a92"
        else:
            nes = float(r["NES"])
            sig = r["significant"] == "True"
            face = cmap(norm(nes)) if sig else "#eceae6"
            txt = f"{nes:+.2f}" if sig else "ns"
            col = "white" if sig and abs(nes) > 1.75 else ("#2b2b2b" if sig else "#9a968f")
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", lw=1.2))
        if txt and txt != "untestable":
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=8.6, color=col)

ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
lab = []
for i, g in enumerate(ORDER):
    k = idx.get((g, "Nativio", "full ranked list"), {}).get("K", "")
    lab.append(f"{NAMES[g]}  (n={k})" if k and k != "0" else f"{NAMES[g]}")
ax.set_yticklabels(lab, fontsize=9.2)
for t in ax.get_yticklabels():          # all labels black + bold (core/non-core distinction removed)
    t.set_color("#111"); t.set_fontweight("bold")
ax.set_xticks(np.arange(len(cols)) + .5)
ax.set_xticklabels([c[2] for c in cols], fontsize=8.8, style="italic")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
for j in range(2, len(cols), 2):
    ax.plot([j, j], [0, len(ORDER)], color="#4a4a4a", lw=2.4, zorder=5)
for k, ds in enumerate(DATASETS):
    ax.text(k * 2 + 1, -0.75, ds, ha="center", va="bottom", fontsize=13, fontweight="bold")

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.022, pad=0.014)
cb.set_label("NES  (+ enriched among up-regulated,\n−  enriched among down-regulated)", fontsize=9)
cb.ax.tick_params(labelsize=8.5)
ax.set_title("Preranked GSEA on the UPR-related GO terms —\n"
             "direction comes from the sign of NES, so no\n"
             "UP/DOWN split and no cutoff are needed",
             loc="left", fontsize=12.6, fontweight="bold", pad=58)

# footnote wrapped to the plot width so it does not overhang the heatmap (removes right margin)
_foot = ("Weighted (p=1) preranked GSEA; NES and P from 10,000 size-matched gene-set permutations, "
         "BH-adjusted across all 18 candidate terms within each condition (the 14 testable UPR/ERAD-related "
         "terms are shown). Rank metric = sign(log2FC) × −log10(P) (Nativio: Welch t statistic). "
         "Gene sets = QuickGO human protein annotations for each GO term plus descendants; "
         "n = set genes present in the Nativio ranked list.")
fig.text(0.055, 0.028, "\n".join(textwrap.wrap(_foot, 118)),
         fontsize=8, color="#666", va="bottom")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q3_GSEA_14terms.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R2Q3_GSEA_14terms.{png,svg}")
