#!/usr/bin/env python3
"""R1.5 representative-gene heatmap. For each curated UPR gene set (ER-stress, PERK, IRE1, ATF6, ERAD),
show up to 10 representative genes = the strongest reproduced-down genes (down in both bulk cohorts AND
both neuron types), as a log2FC heatmap across Mizuno, Nativio (bulk) and snRNA Ex, In (neurons)."""
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
SETS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
        ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"), ("ERAD (75)", "geneset_ERAD.txt")]
gsets = {n: rd(f) for n, f in SETS}

miz = {g: v[0] for g, v in load_mizuno().items()}
nat = {g: v[0] for g, v in load_nativio().items()}
wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD6_snRNSeqDB.xlsx"), read_only=True); ws = wb.active
CT = {"Ex": 2, "In": 4}
sn = {}
for r in ws.iter_rows(min_row=2, values_only=True):
    if not r[0]: continue
    d = {}
    for ct, ci in CT.items():
        try:
            fc = float(r[ci])
            if fc > 0 and math.isfinite(fc): d[ct] = math.log2(fc)
        except (TypeError, ValueError): pass
    if len(d) == 2: sn[str(r[0]).strip()] = d
wb.close()

# per set: reproduced-down genes (down in both bulk AND both neurons), top 10 by mean log2FC
rows_data, group_bounds, labels = [], [], []
seen_used = set()
for name, _ in SETS:
    pool = [g for g in gsets[name] if g in miz and g in nat and g in sn
            and miz[g] < 0 and nat[g] < 0 and sn[g]["Ex"] < 0 and sn[g]["In"] < 0]
    pool = sorted(pool, key=lambda g: (miz[g] + nat[g] + sn[g]["Ex"] + sn[g]["In"]) / 4)[:10]
    start = len(rows_data)
    for g in pool:
        rows_data.append((g, miz[g], nat[g], sn[g]["Ex"], sn[g]["In"]))
    group_bounds.append((name, start, len(rows_data)))

genes = [r[0] for r in rows_data]
M = np.array([[r[1], r[2], r[3], r[4]] for r in rows_data], float)
COLS = ["Mizuno", "Nativio", "snRNA Ex", "snRNA In"]
cmap = LinearSegmentedColormap.from_list("bwr2", ["#1b4f72", "#2c7fb8", "#cfe3f2", "#ffffff", "#f6d5cd", "#d1584a", "#8b1a10"])
vmax = 1.5; norm = Normalize(-vmax, vmax)

nrow = len(genes); ncol = len(COLS)
fig, ax = plt.subplots(figsize=(6.2, 0.28 * nrow + 1.8))
fig.subplots_adjust(left=0.42, right=0.86, top=0.90, bottom=0.06)
for i in range(nrow):
    for j in range(ncol):
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cmap(norm(M[i, j])), edgecolor="white", lw=0.8))
        ax.text(j + .5, i + .5, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=6.2,
                color="white" if abs(M[i, j]) > 0.75 else "#2b2b2b")
ax.set_xlim(0, ncol); ax.set_ylim(nrow, 0)
ax.set_yticks(np.arange(nrow) + .5); ax.set_yticklabels(genes, fontsize=7)
ax.set_xticks(np.arange(ncol) + .5); ax.set_xticklabels(COLS, fontsize=8.5, rotation=45, ha="left")
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
# set-group separators + labels on the far left
for name, a, b in group_bounds:
    if a > 0: ax.plot([0, ncol], [a, a], color="#333", lw=1.6)
    ax.text(-1.65, (a + b) / 2, name.split(" (")[0], ha="left", va="center", rotation=90,
            fontsize=8.5, fontweight="bold", color="#1b2a3a")
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, fraction=0.04, pad=0.03)
cb.set_label("log$_2$FC (blue = down)", fontsize=8); cb.ax.tick_params(labelsize=7)
fig.suptitle("Representative UPR genes — strongest reproduced-down per gene set\n"
             "(log$_2$FC across bulk cohorts and snRNA neurons)", fontsize=10, fontweight="bold", x=0.42, ha="left")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_representative_gene_heatmap.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_representative_gene_heatmap.{png,svg}")
for name, a, b in group_bounds: print(f"  {name}: {b-a} genes")
