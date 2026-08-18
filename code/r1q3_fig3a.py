#!/usr/bin/env python3
"""Reproduce manuscript Figure 3a: combined tSNE of all nuclei (n=48 donors), coloured by
broad.cell.type with on-plot cell-type labels. Same palette + Pericytes (Per) added. Local; no pandas."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

MATH = ("/data/adbrainseq/Stanford U/Collaboration support/"
        "Prof. Eun-hye Joe Ajou/2019 Mathys/Gene Expression (RNA seq)/filtered_column_metadata.txt")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
       "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
COL = {"Ex": "#8CBF43", "In": "#33A45C", "Ast": "#F0806A", "Mic": "#28C2D4", "Oli": "#4F86C6",
       "Opc": "#9C77B4", "End": "#E3A43A", "Per": "#B94FA0"}
PLOT_ORDER = ["Ex", "Oli", "In", "Ast", "Opc", "Mic", "End", "Per"]
# on-plot labels: (text_x, text_y, text, leader_target_or_None)
LAB = {
    "Ex":  (16, 3, "Excitatory\nneurons (Ex)", None),
    "In":  (-46, 42, "Inhibitory\nneurons (In)", None),
    "Ast": (47, 20, "Astrocytes\n(Ast)", None),
    "Mic": (12, -44, "Microglia\n(Mic)", None),
    "Oli": (-24, -28, "Oligodendro-\ncytes (OligD)", None),
    "Opc": (-53, -29, "Oligo-\ndendrocyte\nprogenitor\ncells (OPC)", (-55, -3)),
    "End": (30, 56, "Endothelial\ncells (End)", (9, 54)),
    "Per": (-16, 57, "Pericytes\n(Per)", (5, 50)),
}

x, y, ct = [], [], []
with open(MATH) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        try: xx, yy = float(r["tsne1"]), float(r["tsne2"])
        except (TypeError, ValueError): continue
        x.append(xx); y.append(yy); ct.append(r["broad.cell.type"])
x = np.array(x); y = np.array(y); ct = np.array(ct)
ndon = 48

fig, ax = plt.subplots(figsize=(6.6, 6.2))
for c in PLOT_ORDER:
    m = ct == c
    ax.scatter(x[m], y[m], s=2.2, c=COL[c], linewidths=0, rasterized=True)
for c, (tx, ty, txt, tgt) in LAB.items():
    if tgt is not None:
        ax.annotate("", xy=tgt, xytext=(tx, ty), arrowprops=dict(arrowstyle="-", color="#333", lw=0.7))
    ax.text(tx, ty, txt, ha="center", va="center", fontsize=9.6, color="#111",
            linespacing=1.0, fontweight="normal")
ax.set_xlim(-64, 60); ax.set_ylim(-62, 64)
ax.set_xticks([-60, -30, 0, 30, 60]); ax.set_yticks([-60, -30, 0, 30, 60])
ax.set_xlabel("tsne1", fontsize=14.3); ax.set_ylabel("tsne2", fontsize=14.3)
ax.tick_params(labelsize=14.3)
ax.set_title(f"syn18485175   (n={ndon})", fontsize=14.3, loc="center")
ax.text(-0.14, 1.03, "A", transform=ax.transAxes, fontsize=26, fontweight="bold")
fig.subplots_adjust(left=0.12, right=0.97, top=0.92, bottom=0.10)
for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUT, f"R1Q3_Figure3a_reproduced.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q3_Figure3a_reproduced.{png,svg}  cells:", len(x))
