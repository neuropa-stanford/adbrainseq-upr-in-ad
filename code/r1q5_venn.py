#!/usr/bin/env python3
"""R1.5 Venn diagrams: bulk (Mizuno, Nativio) vs snRNA directional concordance,
per gene set (ER-stress 260, PERK 31, IRE1 32, ATF6 74).
Row 1 = DOWN (Mizuno-down / Nativio-down / snRNA-neuron-down) — reproduced in bulk.
Row 2 = UP   (Mizuno-up / Nativio-up / snRNA-glia-up) — glial up-regulation masked in bulk.
"""
import sys, math, os
sys.path.insert(0, '.')
from r1q1_gomatrix import load_mizuno, load_nativio
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from matplotlib_venn import venn3, venn3_circles

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
OUT = os.path.join(BASE, "Major Revision", "processed raw data")

GS = {}
GS["ER-stress (260)"] = set(l.strip() for l in open(os.path.join(OUT, "ERstress_260_geneset.txt")) if l.strip())
for b in ["PERK", "IRE1", "ATF6"]:
    GS[f"{b} ({{}})"] = set(l.strip() for l in open(os.path.join(OUT, f"geneset_{b}.txt")) if l.strip())

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

names = list(GS.keys())
fig, axes = plt.subplots(2, 4, figsize=(19, 9.5))
fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.04, hspace=0.40, wspace=0.14)
COLD = ("#2c7fb8", "#7bccc4", "#c0392b")   # Mizuno, Nativio, snRNA colors

for ci, name in enumerate(names):
    genes = GS[name]
    U = [g for g in genes if g in miz and g in nat and g in sn]
    title = name.format(len(genes)) if "{}" in name else name
    # ---- DOWN row ----
    mizD = set(g for g in U if miz[g] < 0)
    natD = set(g for g in U if nat[g] < 0)
    snND = set(g for g in U if sn[g]["Ex"] < 0 and sn[g]["In"] < 0)
    ndc = len(mizD & natD & snND)
    axD = axes[0][ci]
    venn3([mizD, natD, snND], set_labels=("Mizuno\nbulk↓", "Nativio\nbulk↓", "snRNA\nneuron↓"),
          set_colors=COLD, alpha=0.55, ax=axD)
    for t in axD.texts:
        t.set_fontsize(8.5)
    axD.set_title(f"{title}   ({len(U)} genes)\n{ndc} concordantly DOWN in all three",
                  fontsize=10.5, fontweight="bold", color="#1b4f72")
    # ---- UP row ----
    mizU = set(g for g in U if miz[g] > 0)
    natU = set(g for g in U if nat[g] > 0)
    snGU = set(g for g in U if sn[g]["Mic"] > 0 and sn[g]["Oli"] > 0)
    guc = len(mizU & natU & snGU)
    axU = axes[1][ci]
    venn3([mizU, natU, snGU], set_labels=("Mizuno\nbulk↑", "Nativio\nbulk↑", "snRNA\nglia↑"),
          set_colors=COLD, alpha=0.55, ax=axU)
    for t in axU.texts:
        t.set_fontsize(8.5)
    axU.set_title(f"only {guc} concordantly UP in all three", fontsize=10.5,
                  fontweight="bold", color="#8c1c13")

fig.text(0.5, 0.955, "Cross-modality directional concordance of UPR gene sets — bulk (Mizuno, Nativio) vs single-nucleus RNA-seq",
         ha="center", fontsize=14, fontweight="bold")
fig.text(0.5, 0.925, "Top: neuronal downregulation is reproduced in bulk.   "
         "Bottom: glial upregulation is largely absent from the neuron-dominated bulk signal (masked).",
         ha="center", fontsize=10.5, color="#555")
fig.text(0.012, 0.66, "DOWN\n(neurons)", fontsize=12, fontweight="bold", color="#1b4f72", rotation=90, va="center")
fig.text(0.012, 0.24, "UP\n(glia)", fontsize=12, fontweight="bold", color="#8c1c13", rotation=90, va="center")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_venn_bulk_vs_snRNA.{ext}"),
                dpi=200 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_venn_bulk_vs_snRNA.{png,svg}")
