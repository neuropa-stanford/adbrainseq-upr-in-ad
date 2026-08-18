#!/usr/bin/env python3
"""R1.5 representative-gene heatmap, DOWN | UP.
DOWN panel: up to 5 reproduced-DOWN genes per UPR set (down in both bulk + both NEURONS), cols =
Mizuno, Nativio, snRNA Ex, snRNA In. UP panel: up to 5 reproduced-UP genes per set (up in both bulk +
both GLIA), cols = Mizuno, Nativio, snRNA Mic, snRNA Oli. Diverging log2FC colormap (blue down / red up)."""
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
def rd(f): return [l.strip() for l in open(os.path.join(OUT, f)) if l.strip()]
SETS = [("ER-stress", "ERstress_260_geneset.txt"), ("PERK", "geneset_PERK.txt"),
        ("IRE1", "geneset_IRE1.txt"), ("ATF6", "geneset_ATF6.txt"), ("ERAD", "geneset_ERAD.txt")]
gsets = {n: rd(f) for n, f in SETS}
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

def select(direction):
    """returns list of (set_name, [ (gene, v_miz, v_nat, v_c1, v_c2) ... ]) ; c1,c2 = Ex/In or Mic/Oli"""
    out = []
    for name, _ in SETS:
        rows = []
        for g in gsets[name]:
            if g not in miz or g not in nat or g not in sn: continue
            s = sn[g]
            if direction == "down":
                if not (all(k in s for k in ("Ex", "In")) and miz[g] < 0 and nat[g] < 0 and s["Ex"] < 0 and s["In"] < 0): continue
                rows.append((g, miz[g], nat[g], s["Ex"], s["In"]))
            else:
                if not (all(k in s for k in ("Mic", "Oli")) and miz[g] > 0 and nat[g] > 0 and s["Mic"] > 0 and s["Oli"] > 0): continue
                rows.append((g, miz[g], nat[g], s["Mic"], s["Oli"]))
        rows.sort(key=lambda r: (r[1] + r[2]) / 2, reverse=(direction == "up"))   # sort by bulk (Miz+Nat) magnitude
        out.append((name, rows[:5]))
    return out

DN = select("down"); UP = select("up")
cmap = LinearSegmentedColormap.from_list("bwr2", ["#1b4f72", "#2c7fb8", "#cfe3f2", "#ffffff", "#f6d5cd", "#d1584a", "#8b1a10"])
norm = Normalize(-1.5, 1.5)

def flat(sel):
    genes, bounds = [], []
    for name, rows in sel:
        a = len(genes)
        for r in rows: genes.append(r)
        bounds.append((name, a, len(genes)))
    return genes, bounds
gdn, bdn = flat(DN); gup, bup = flat(UP)
maxrows = max(len(gdn), len(gup))

fig = plt.figure(figsize=(11.5, 0.30 * maxrows + 2.1))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], left=0.10, right=0.955, top=0.80, bottom=0.13, wspace=0.55)

def panel(axpos, genes, bounds, collab, title):
    ax = fig.add_subplot(axpos)
    for i, r in enumerate(genes):
        vals = r[1:5]
        for j, v in enumerate(vals):
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(v)), edgecolor="white", lw=0.8))
            ax.text(j + .5, i + .5, f"{v:.2f}", ha="center", va="center", fontsize=6.2,
                    color="white" if abs(v) > 0.8 else "#2b2b2b")
    ax.set_xlim(0, 4); ax.set_ylim(maxrows, 0)
    ax.set_yticks(np.arange(len(genes)) + .5); ax.set_yticklabels([r[0] for r in genes], fontsize=7.2)
    ax.set_xticks(np.arange(4) + .5); ax.set_xticklabels(collab, fontsize=8, rotation=45, ha="left")
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    for name, a, b in bounds:
        if a > 0: ax.plot([0, 4], [a, a], color="#333", lw=1.3)
        ax.text(-1.35, (a + b) / 2, name, ha="right", va="center", fontsize=8, fontweight="bold", color="#1b2a3a")
    ax.set_title(title, fontsize=10.5, fontweight="bold", loc="left", pad=14)
    return ax

panel(gs[0, 0], gdn, bdn, ["Mizuno", "Nativio", "Ex", "In"], "a  DOWN — reproduced in neurons")
panel(gs[0, 1], gup, bup, ["Mizuno", "Nativio", "Mic", "Oli"], "b  UP — reproduced in glia")
cax = fig.add_axes([0.37, 0.05, 0.27, 0.015])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
cb.set_label("log$_2$FC (blue = down · red = up)", fontsize=8); cb.ax.tick_params(labelsize=7)
fig.text(0.10, 0.965, "Representative overlapped UPR genes — 5 per UPR branch (not deduplicated) · "
         "DOWN in neurons vs UP in glia · sorted most-changed at top",
         fontsize=10.8, fontweight="bold", ha="left", va="top")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_representative_gene_heatmap_DownUp.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_representative_gene_heatmap_DownUp.{png,svg}")
for n, r in DN: print(f"  DOWN {n}: {len(r)}")
for n, r in UP: print(f"  UP   {n}: {len(r)}")
