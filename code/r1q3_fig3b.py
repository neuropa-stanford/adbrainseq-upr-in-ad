#!/usr/bin/env python3
"""Reproduce manuscript Figure 3b: stacked-violin of canonical cell-type marker genes across
broad.cell.type, from the local Mathys count matrix. Same palette + Pericytes (Per, marker AMBP).
Rows = markers, cols = cell types; each violin = log-normalised expression, filled by cell-type
colour; per-row shared scale so the on-diagonal marker stands out. Local; sparse; no pandas."""
import os, csv
import numpy as np
from scipy.io import mmread
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False, "pdf.fonttype": 42, "ps.fonttype": 42})

MDIR = ("/data/adbrainseq/Stanford U/Collaboration support/"
        "Prof. Eun-hye Joe Ajou/2019 Mathys/Gene Expression (RNA seq)")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
       "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
COL = {"Ex": "#8CBF43", "In": "#33A45C", "Ast": "#F0806A", "Mic": "#28C2D4", "Oli": "#4F86C6",
       "Opc": "#9C77B4", "End": "#E3A43A", "Per": "#B94FA0"}
CT = ["Ex", "In", "Ast", "Mic", "Oli", "Opc", "End", "Per"]
CTLAB = {"Ex": "Excitatory neurons (Ex)", "In": "Inhibitory neurons (In)", "Ast": "Astrocytes (Ast)",
         "Mic": "Microglia (Mic)", "Oli": "Oligodendrocytes (OligD)", "Opc": "Oligodendrocyte progenitor cells (OPC)",
         "End": "Endothelial cells (End)", "Per": "Pericytes (Per)"}
# marker gene per cell type (row), diagonal
MARK = [("NRGN", "Ex"), ("GAD1", "In"), ("AQP4", "Ast"), ("CSF1R", "Mic"),
        ("MBP", "Oli"), ("VCAN", "Opc"), ("FLT1", "End"), ("PDGFRB", "Per")]
# NB: Mathys used AMBP for pericytes but AMBP is undetected in this filtered matrix (as the team
# also observed); PDGFRB is the canonical pericyte marker and is highest in Per here.

genes = [l.strip() for l in open(os.path.join(MDIR, "filtered_gene_row_names.txt"))]
gidx = {g: i for i, g in enumerate(genes)}
cell_ct = []
with open(os.path.join(MDIR, "filtered_column_metadata.txt")) as f:
    rd = csv.reader(f, delimiter="\t"); hdr = next(rd); ci = hdr.index("broad.cell.type")
    for r in rd: cell_ct.append(r[ci])
cell_ct = np.array(cell_ct)

print("reading matrix ...", flush=True)
M = mmread(os.path.join(MDIR, "filtered_count_matrix.mtx")).tocsr()
colsum = np.asarray(M.sum(axis=0)).ravel().astype(float); colsum[colsum == 0] = 1.0
def norm_log(g):
    v = np.asarray(M.getrow(gidx[g]).todense()).ravel().astype(float)
    return np.log1p(v / colsum * 1e4)
# expr[marker][celltype] = array of log-norm expression
expr = {}
for mk, _ in MARK:
    le = norm_log(mk)
    expr[mk] = {c: le[cell_ct == c] for c in CT}

nR, nC = len(MARK), len(CT)
fig, axes = plt.subplots(nR, nC, figsize=(11.5, 9.2), sharex="col")
for i, (mk, home) in enumerate(MARK):
    rowmax = max(np.percentile(expr[mk][c], 99.5) if len(expr[mk][c]) else 0 for c in CT) or 1.0
    for j, c in enumerate(CT):
        ax = axes[i, j]; v = expr[mk][c]
        if len(v) > 2 and v.max() > 0:
            parts = ax.violinplot([v], positions=[0], widths=0.95, showmeans=False, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(COL[c]); pc.set_alpha(0.9); pc.set_edgecolor("#333"); pc.set_linewidth(0.5)
        ax.set_ylim(-0.05 * rowmax, rowmax * 1.05); ax.set_xlim(-0.7, 0.7)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values(): sp.set_visible(False)
        ax.axhline(0, color="#bbb", lw=0.5)
        if j == 0:
            ax.set_ylabel(mk, rotation=0, ha="right", va="center", fontsize=20.4,
                          fontweight="bold", color=COL[home], labelpad=12)
        if i == 0:
            ax.text(0.0, 1.08, CTLAB[c], transform=ax.transAxes, rotation=30, rotation_mode="anchor",
                    ha="left", va="bottom", fontsize=20.4, fontweight="bold", color=COL[c])
fig.text(0.02, 0.985, "B", fontsize=37, fontweight="bold", va="top")
fig.subplots_adjust(left=0.11, right=0.985, top=0.74, bottom=0.02, hspace=0.15, wspace=0.15)
for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUT, f"R1Q3_Figure3b_reproduced.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q3_Figure3b_reproduced.{png,svg}")
for mk, home in MARK:
    print(f"  {mk:6s} home={home:3s} mean-in-home={expr[mk][home].mean():.2f}  "
          f"max-other={max(expr[mk][c].mean() for c in CT if c!=home):.2f}")
