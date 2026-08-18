#!/usr/bin/env python3
"""R1.5 representative-gene heatmap, clean version. DOWN | UP panels, DEDUPLICATED genes, sorted by
magnitude of change (most-changed at top). DOWN = reproduced in neurons (Ex,In); UP = reproduced in
glia (Mic,Oli). Branch shown as a right-side annotation, no group separators."""
import sys, os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import openpyxl
sys.path.insert(0, ".")
from r1q1_gomatrix import load_mizuno, load_nativio

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
OUT = os.path.join(BASE, "Major Revision", "processed raw data")
def rd(f): return set(l.strip() for l in open(os.path.join(OUT, f)) if l.strip())
PERK, IRE1, ATF6, ERAD, ERS = rd("geneset_PERK.txt"), rd("geneset_IRE1.txt"), rd("geneset_ATF6.txt"), rd("geneset_ERAD.txt"), rd("ERstress_260_geneset.txt")
def branch(g):
    if g in PERK: return "PERK/ISR"
    if g in IRE1: return "IRE1"
    if g in ATF6: return "ATF6"
    if g in ERAD: return "ERAD"
    return "chaperone/ER-stress"
UNIV = PERK | IRE1 | ATF6 | ERAD | ERS

miz = {g: v[0] for g, v in load_mizuno().items()}
nat = {g: v[0] for g, v in load_nativio().items()}
wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD6_snRNSeqDB.xlsx"), read_only=True); ws = wb.active
CT = {"Ex": 2, "In": 4, "Mic": 8, "Oli": 10}
sn = {}
for r in ws.iter_rows(min_row=2, values_only=True):
    if not r[0]: continue
    d = {}
    for ct, ci in CT.items():
        try:
            fc = float(r[ci])
            if fc > 0 and math.isfinite(fc): d[ct] = math.log2(fc)
        except (TypeError, ValueError): pass
    sn[str(r[0]).strip()] = d

N = 12
def collect(direction):
    rows = []
    for g in UNIV:
        if g not in miz or g not in nat or g not in sn: continue
        s = sn[g]
        if direction == "down":
            if not (all(k in s for k in ("Ex", "In")) and miz[g] < 0 and nat[g] < 0 and s["Ex"] < 0 and s["In"] < 0): continue
            rows.append((g, miz[g], nat[g], s["Ex"], s["In"], branch(g)))
        else:
            if not (all(k in s for k in ("Mic", "Oli")) and miz[g] > 0 and nat[g] > 0 and s["Mic"] > 0 and s["Oli"] > 0): continue
            rows.append((g, miz[g], nat[g], s["Mic"], s["Oli"], branch(g)))
    rows.sort(key=lambda r: (r[1] + r[2] + r[3] + r[4]) / 4, reverse=(direction == "up"))
    return rows[:N]
DN, UP = collect("down"), collect("up")

cmap = LinearSegmentedColormap.from_list("bwr2", ["#1b4f72", "#2c7fb8", "#cfe3f2", "#ffffff", "#f6d5cd", "#d1584a", "#8b1a10"])
norm = Normalize(-1.5, 1.5)
maxrows = max(len(DN), len(UP))
fig = plt.figure(figsize=(11.8, 0.34 * maxrows + 1.8))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], left=0.085, right=0.905, top=0.80, bottom=0.15, wspace=0.62)

def panel(pos, rows, collab, title):
    ax = fig.add_subplot(pos)
    for i, r in enumerate(rows):
        for j, v in enumerate(r[1:5]):
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(v)), edgecolor="white", lw=0.8))
            ax.text(j + .5, i + .5, f"{v:.2f}", ha="center", va="center", fontsize=6.6,
                    color="white" if abs(v) > 0.8 else "#2b2b2b")
        ax.text(4.12, i + .5, r[5], ha="left", va="center", fontsize=6.6, color="#888")
    ax.set_xlim(0, 4); ax.set_ylim(maxrows, 0)
    ax.set_yticks(np.arange(len(rows)) + .5); ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xticks(np.arange(4) + .5); ax.set_xticklabels(collab, fontsize=8.5, rotation=45, ha="left")
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left", pad=14)
    return ax
panel(gs[0, 0], DN, ["Mizuno", "Nativio", "Ex", "In"], "a  DOWN in neurons")
panel(gs[0, 1], UP, ["Mizuno", "Nativio", "Mic", "Oli"], "b  UP in glia")
cax = fig.add_axes([0.37, 0.055, 0.27, 0.017])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label("log$_2$FC (blue = down · red = up)", fontsize=8); cb.ax.tick_params(labelsize=7)
fig.text(0.085, 0.965, "Representative overlapped UPR genes — top " + str(N) + " by magnitude, "
         "deduplicated · sorted most-changed at top (DOWN in neurons vs UP in glia)",
         fontsize=10.8, fontweight="bold", ha="left", va="top")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_representative_gene_heatmap_clean.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_representative_gene_heatmap_clean.{png,svg}")
print("DOWN:", [r[0] for r in DN]); print("UP:", [r[0] for r in UP])
