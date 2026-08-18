#!/usr/bin/env python3
"""R1.5 as a heatmap. Rows = curated UPR gene sets (same sets as the UPR analysis: ER-stress 260,
PERK 31, IRE1 32, ATF6 74, ERAD 75). DOWN | UP blocks. Count cells shaded per-column; the
'Reproduced in snRNA (%)' columns share a 0-100% green scale so the DOWN-vs-UP contrast is visible."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

SETS = ["ER-stress (260)", "PERK (31)", "IRE1 (32)", "ATF6 (74)", "ERAD (75)"]
# DOWN: Mizuno↓, Nativio↓, Bulk∩↓, snRNA neuron↓, reproduced_count, reproduced_pct
DOWN = {"ER-stress (260)": (109, 165, 100, 177, 91, 91), "PERK (31)": (12, 21, 12, 22, 11, 91),
        "IRE1 (32)": (8, 22, 7, 21, 5, 71), "ATF6 (74)": (32, 52, 28, 59, 26, 92),
        "ERAD (75)": (39, 55, 37, 57, 34, 91)}
UP = {"ER-stress (260)": (95, 39, 30, 95, 13, 43), "PERK (31)": (15, 6, 6, 15, 5, 83),
      "IRE1 (32)": (20, 6, 5, 20, 2, 40), "ATF6 (74)": (35, 15, 11, 36, 7, 63),
      "ERAD (75)": (23, 7, 5, 37, 4, 80)}
CNTLAB = ["Mizuno", "Nativio", "Bulk\noverlapped", "snRNA\nneuron/glia", "Reproduced\nin snRNA (%)"]

cmap_dn = LinearSegmentedColormap.from_list("dn", ["#eef4fb", "#9dc3e6", "#2c5f8a"])
cmap_up = LinearSegmentedColormap.from_list("up", ["#fdeeea", "#efab93", "#c0392b"])
cmap_rep = LinearSegmentedColormap.from_list("rep", ["#f1f7ee", "#a6d08a", "#2e7d32"])

fig, ax = plt.subplots(figsize=(12.6, 4.6))
fig.subplots_adjust(top=0.72, bottom=0.13, left=0.105, right=0.99)
ncol = 10  # 5 down + 5 up
nrow = len(SETS)

def col_norm(vals):
    v = np.array(vals, float); lo, hi = v.min(), v.max()
    return (v - lo) / (hi - lo) if hi > lo else np.zeros_like(v)

# per-column normalized shading for the 4 count columns in each block; reproduced on 0-100 scale
for block, data, cmap_c, xoff in [("DOWN", DOWN, cmap_dn, 0), ("UP", UP, cmap_up, 5)]:
    counts = np.array([data[s][:4] for s in SETS], float)          # rows x 4
    for c in range(4):
        norm_c = col_norm(counts[:, c])
        for r in range(nrow):
            j = xoff + c
            ax.add_patch(Rectangle((j, r), 1, 1, facecolor=cmap_c(0.15 + 0.85 * norm_c[r]),
                                   edgecolor="white", lw=1.4))
            ax.text(j + .5, r + .5, f"{int(counts[r, c])}", ha="center", va="center",
                    fontsize=9, color="#1b1b1b")
    # reproduced % column (5th of block), shared 0-100 green scale
    for r, s in enumerate(SETS):
        pct = data[s][5]; cnt = data[s][4]
        j = xoff + 4
        ax.add_patch(Rectangle((j, r), 1, 1, facecolor=cmap_rep(pct / 100.0),
                               edgecolor="white", lw=1.4))
        ax.text(j + .5, r + .5, f"{cnt} ({pct}%)", ha="center", va="center",
                fontsize=8.6, fontweight="bold", color="white" if pct > 55 else "#1b1b1b")

ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0)
ax.set_yticks(np.arange(nrow) + .5); ax.set_yticklabels(SETS, fontsize=9.5, fontweight="bold")
ax.set_xticks(np.arange(ncol) + .5)
ax.set_xticklabels(CNTLAB + CNTLAB, fontsize=7.8)
ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.plot([5, 5], [0, nrow], color="#222", lw=3.2, zorder=6)                 # DOWN | UP divider
ax.text(2.5, -0.55, "DOWN-regulated", ha="center", fontsize=12.5, fontweight="bold", color="#2c5f8a")
ax.text(7.5, -0.55, "UP-regulated", ha="center", fontsize=12.5, fontweight="bold", color="#c0392b")
fig.text(0.105, 0.95, "Cross-modality directional reproducibility of UPR gene sets in AD brain "
         "transcriptomics", fontsize=11.5, fontweight="bold", ha="left", va="top")
fig.text(0.105, 0.905, "Same curated UPR gene sets (ER-stress 260, PERK 31, IRE1 32, ATF6 74, ERAD 75). "
         "Cell = gene count; green column = % reproduced in snRNA (of bulk-overlapped).",
         fontsize=8.6, color="#666", ha="left", va="top")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_reproducibility_heatmap.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_reproducibility_heatmap.{png,svg}")
