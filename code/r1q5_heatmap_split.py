#!/usr/bin/env python3
"""R1.5 reproducibility heatmap with snRNA split: neuron -> ExN, InN ; glia -> Mic, OligD.
DOWN block: Mizuno↓, Nativio↓, Bulk∩↓, ExN↓, InN↓, Reproduced(%).
UP block:   Mizuno↑, Nativio↑, Bulk∩↑, Mic↑, OligD↑, Reproduced(%).
Reproduced-DOWN = bulk-overlapped-down AND down in BOTH neuron types; Reproduced-UP = bulk-overlapped-up
AND up in BOTH glia types (unchanged definition, only the snRNA display column is split)."""
import sys, os, math, csv
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
GS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
      ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"), ("ERAD (75)", "geneset_ERAD.txt")]
sets = {n: set(l.strip() for l in open(os.path.join(OUT, f)) if l.strip()) for n, f in GS}
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
    if len(d) == 4: sn[str(r[0]).strip()] = d
wb.close()

DOWN, UP, tabrows = {}, {}, []
for name, _ in GS:
    U = [g for g in sets[name] if g in miz and g in nat and g in sn]
    md = sum(1 for g in U if miz[g] < 0); mu = sum(1 for g in U if miz[g] > 0)
    nd = sum(1 for g in U if nat[g] < 0); nu = sum(1 for g in U if nat[g] > 0)
    bd = set(g for g in U if miz[g] < 0 and nat[g] < 0); bu = set(g for g in U if miz[g] > 0 and nat[g] > 0)
    exn = sum(1 for g in U if sn[g]["Ex"] < 0); inn = sum(1 for g in U if sn[g]["In"] < 0)
    micu = sum(1 for g in U if sn[g]["Mic"] > 0); oliu = sum(1 for g in U if sn[g]["Oli"] > 0)
    repd = bd & set(g for g in U if sn[g]["Ex"] < 0 and sn[g]["In"] < 0)
    repu = bu & set(g for g in U if sn[g]["Mic"] > 0 and sn[g]["Oli"] > 0)
    pd = 100 * len(repd) // max(len(bd), 1); pu = 100 * len(repu) // max(len(bu), 1)
    DOWN[name] = (md, nd, len(bd), exn, inn, len(repd), pd)
    UP[name] = (mu, nu, len(bu), micu, oliu, len(repu), pu)
    tabrows.append([name, md, nd, len(bd), exn, inn, f"{len(repd)} ({pd}%)",
                    mu, nu, len(bu), micu, oliu, f"{len(repu)} ({pu}%)"])

hdr = ["UPR gene set", "Miz↓", "Nat↓", "Bulk∩↓", "ExN↓", "InN↓", "Reproduced↓ (%)",
       "Miz↑", "Nat↑", "Bulk∩↑", "Mic↑", "OligD↑", "Reproduced↑ (%)"]
with open(os.path.join(OUT, "R1Q5_reproducibility_table_split.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(hdr); w.writerows(tabrows)
print("=== split table ===")
print("  ".join(hdr))
for r in tabrows: print(r)

# ---- heatmap ----
SETS = [n for n, _ in GS]
CNTLAB_D = ["Mizuno", "Nativio", "Bulk\noverlapped", "ExN", "InN", "Reproduced\nin snRNA (%)"]
CNTLAB_U = ["Mizuno", "Nativio", "Bulk\noverlapped", "Mic", "OligD", "Reproduced\nin snRNA (%)"]
# flat, unified fills — cell colour is no longer a (redundant) count gradient; the number in each
# cell carries the value. DOWN = light sky-blue, UP = light pink, Reproduced = mild green.
FILL_DN, FILL_UP, FILL_REP = "#d3e7f7", "#f8ddd6", "#a9d69a"

fig, ax = plt.subplots(figsize=(15.5, 4.9))
fig.subplots_adjust(top=0.66, bottom=0.04, left=0.09, right=0.995)
nrow = len(SETS)
for block, data, fill_c, xoff in [("DOWN", DOWN, FILL_DN, 0), ("UP", UP, FILL_UP, 6)]:
    counts = np.array([data[s][:5] for s in SETS], float)
    for c in range(5):
        for r in range(nrow):
            j = xoff + c
            ax.add_patch(Rectangle((j, r), 1, 1, facecolor=fill_c, edgecolor="white", lw=1.4))
            ax.text(j + .5, r + .5, f"{int(counts[r, c])}", ha="center", va="center", fontsize=8.6, color="#1b1b1b")
    for r, s in enumerate(SETS):
        cnt, pct = data[s][5], data[s][6]; j = xoff + 5
        ax.add_patch(Rectangle((j, r), 1, 1, facecolor=FILL_REP, edgecolor="white", lw=1.4))
        ax.text(j + .5, r + .5, f"{cnt} ({pct}%)", ha="center", va="center", fontsize=8, fontweight="bold",
                color="#1b1b1b")
ax.set_xlim(0, 12); ax.set_ylim(nrow, 0)
ax.set_yticks(np.arange(nrow) + .5); ax.set_yticklabels(SETS, fontsize=9.5, fontweight="bold")
ax.set_xticks(np.arange(12) + .5); ax.set_xticklabels(CNTLAB_D + CNTLAB_U, fontsize=7.4)
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
ax.tick_params(length=0)
for sp in ax.spines.values(): sp.set_visible(False)
ax.plot([6, 6], [0, nrow], color="#222", lw=3.2, zorder=6)
ax.text(3, -1.55, "DOWN-regulated", ha="center", fontsize=12.5, fontweight="bold", color="#2c5f8a")
ax.text(9, -1.55, "UP-regulated", ha="center", fontsize=12.5, fontweight="bold", color="#c0392b")
fig.text(0.09, 0.985, "Cross-modality directional reproducibility of UPR gene sets in AD brain "
         "transcriptomics", fontsize=11.5, fontweight="bold", ha="left", va="top")
fig.text(0.09, 0.945, "Same curated UPR gene sets. snRNA neurons split into ExN, InN; glia split into Mic, OligD. "
         "Cell number = gene count; green column = number (%) reproduced in snRNA of the bulk-overlapped set.",
         fontsize=8.4, color="#666", ha="left", va="top")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_reproducibility_heatmap_split.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("\nsaved R1Q5_reproducibility_heatmap_split.{png,svg} + table_split.csv")
