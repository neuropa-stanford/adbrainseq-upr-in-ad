#!/usr/bin/env python3
"""Mathys 2024 (Nature; ROSMAP · PFC · 48 donors · Braak I-VI) independent-cohort replication in the
manuscript Figure-4 style. Per broad cell type: per-cell-type DE logFC (logFC_nb) of the
Response-to-ER-stress set (260 genes) vs NFT-tangle burden (path=nft; continuous tau pathology, NOT a
binary AD/non-AD contrast). Green violins, red mean bar, blue quartiles, black points, sign test."""
import os, glob, gzip, csv, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from matplotlib.patches import Patch
from scipy.stats import binomtest

OUT = os.path.dirname(os.path.abspath(__file__))
SCR = ("/tmp/"
       "adbrainseq_work/mathys24/dereg")
UPR = set(l.strip() for l in open(os.path.join(OUT, "ERstress_260_geneset.txt")) if l.strip())
SENS = {"EIF2AK3", "ERN1", "ATF6"}
MCL = {"Ast": "Ast", "Exc": "Ex", "Inh": "In", "Mic": "Mic", "Oli": "Oli", "Opc": "OPC"}
CELLS = [("Ex", "Excitatory\nneurons"), ("In", "Inhibitory\nneurons"), ("Ast", "Astrocytes"),
         ("Mic", "Microglia"), ("Oli", "Oligodendrocytes"), ("OPC", "Oligodendrocyte\nprogenitor cells")]
GREEN = "#4a9e3e"; RED = "#e2231a"; BLUE = "#2b3a8f"

raw = {v: {} for v in MCL.values()}
for fp in glob.glob(os.path.join(SCR, "aggregated_fullset.*.tsv.gz")):
    pre = os.path.basename(fp).split(".")[1].split("_")[0]
    if pre not in MCL: continue
    with gzip.open(fp, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("path") != "nft" or row.get("region") != "PFC": continue
            g = row["gene"]
            if g not in UPR or g in SENS: continue
            try: raw[MCL[pre]].setdefault(g, []).append(float(row["logFC_nb"]))
            except (TypeError, ValueError): continue
data = {c: {g: st.mean(v) for g, v in d.items()} for c, d in raw.items()}

def pstar(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

# scale for display: logFC_nb (per-unit-tangle) is small; standardize per cell type is NOT done -
# keep raw but set y-range to the data; direction/sign is the message.
allv = [v for c, _ in CELLS for v in data[c].values()]
lim = np.percentile(np.abs(allv), 99) * 1.15
fig, ax = plt.subplots(figsize=(11, 5.8))
ax.axhline(0, color="#000", lw=0.8, ls=":", zorder=1)
YL = (-lim, lim)
for i, (c, lab) in enumerate(CELLS):
    vals = np.array(list(data[c].values()))
    dn = int((vals < 0).sum()); n = len(vals); p = binomtest(dn, n, 0.5).pvalue
    parts = ax.violinplot([np.clip(vals, *YL)], positions=[i], widths=0.8, showmeans=False, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(GREEN); pc.set_alpha(0.9); pc.set_edgecolor("#1f1f1f"); pc.set_linewidth(0.8)
    jit = (np.random.RandomState(i).rand(n) - 0.5) * 0.30
    ax.scatter(i + jit, np.clip(vals, *YL), s=4, color="#111", alpha=0.5, edgecolor="none", zorder=3)
    q1, q3 = np.percentile(vals, [25, 75])
    ax.hlines([q1, q3], i - 0.34, i + 0.34, color=BLUE, lw=1.3, zorder=5)
    ax.hlines(vals.mean(), i - 0.36, i + 0.36, color=RED, lw=3.0, zorder=6)
    ax.text(i, YL[1] * 0.78, pstar(p), ha="center", fontsize=13, fontweight="bold")
ax.set_xlim(-0.6, len(CELLS) - 0.4); ax.set_ylim(*YL)
ax.set_xticks(range(len(CELLS))); ax.set_xticklabels([l for _, l in CELLS], fontsize=9.5)
ax.set_ylabel("logFC vs NFT tangle burden\n(Response to ER stress, genes=260)", fontsize=10.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(handles=[Patch(facecolor=GREEN, edgecolor="#1f1f1f", label="Mathys 2024 (logFC vs NFT tangles)"),
                   plt.Line2D([0], [0], color=RED, lw=3, label="set mean"),
                   plt.Line2D([0], [0], color=BLUE, lw=1.3, label="quartiles")],
          frameon=False, fontsize=8.5, loc="upper right")
fig.suptitle("Mathys 2024 (ROSMAP) — Response to ER stress declines in neurons, rises in oligodendrocytes with tau pathology",
             fontsize=11.5, fontweight="bold", y=0.99)
fig.text(0.5, 0.925, "ROSMAP · prefrontal cortex · 48 donors · Braak I–VI · per-cell-type DE vs NFT-tangle burden (nft, "
         "continuous tau — not AD/non-AD) · each point = one UPR gene · sign test vs 0.",
         ha="center", fontsize=8.3, color="#555")
fig.subplots_adjust(top=0.86, bottom=0.11, left=0.10, right=0.98)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q9_Mathys2024_figure.{ext}"), dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved R2Q9_Mathys2024_figure.{png,svg}")
for c, _ in CELLS:
    vals = list(data[c].values()); dn = sum(v < 0 for v in vals)
    print(f"  {c:4s} n={len(vals)} mean={st.mean(vals):+.4f} {100*dn//len(vals)}% down p={binomtest(dn,len(vals),0.5).pvalue:.1e}")
