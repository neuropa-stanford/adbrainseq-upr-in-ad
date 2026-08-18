#!/usr/bin/env python3
"""Donor-level snRNA sensor-gene log2FC violins (from R2Q1_snRNA_sensor_donormeans.npy).
Rows = EIF2AK3/PERK, ERN1/IRE1, ATF6 ; cols = Ex, In, Ast, Mic, Oli, OPC.
Y = per-donor log2FC vs the Braak I-II detected-donor mean (common scale). DETECTED donors are shown
as filled points in the violin; donors where the gene is UNDETECTED (0 counts in that cell type) are
shown as hollow markers on the 'not detected' floor strip, so every donor is displayed. Kruskal-Wallis
across the three Braak groups on the linear donor means (all donors, incl. undetected)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from scipy.stats import kruskal

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
       "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
rows = np.load(os.path.join(OUT, "R2Q1_snRNA_sensor_donormeans.npy"), allow_pickle=True)

SENSORS = [("EIF2AK3", "PERK"), ("ERN1", "IRE1"), ("ATF6", "ATF6")]
CELLS = [("Ex", "Excitatory"), ("In", "Inhibitory"), ("Ast", "Astrocyte"),
         ("Mic", "Microglia"), ("Oli", "Oligodendrocyte"), ("Opc", "OPC")]
GRPS = [("low", "Braak\nI–II", "#b0b0b0"), ("int", "Braak\nIII–IV", "#e8a33d"), ("late", "Braak\nV–VI", "#c0392b")]
YLIM = (-2.6, 2.2)
FLOOR = YLIM[0] + 0.18

fc = {s: {ct: {g: [] for g, _, _ in GRPS} for ct, _ in CELLS} for s, _ in SENSORS}   # detected log2FC
lin = {s: {ct: {g: [] for g, _, _ in GRPS} for ct, _ in CELLS} for s, _ in SENSORS}  # linear (all)
und = {s: {ct: {g: 0 for g, _, _ in GRPS} for ct, _ in CELLS} for s, _ in SENSORS}   # undetected count
for r in rows:
    g, ct, grp, v, f = r[0], r[1], r[4], float(r[5]), r[6]
    if g not in fc or ct not in fc[g]: continue
    lin[g][ct][grp].append(v)
    if f is not None and np.isfinite(float(f)): fc[g][ct][grp].append(float(f))
    else: und[g][ct][grp] += 1

def pstar(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

fig, axes = plt.subplots(len(SENSORS), len(CELLS), figsize=(16.5, 8.8), sharey="row")
for gi, (gene, alias) in enumerate(SENSORS):
    for ci, (ct, ctlab) in enumerate(CELLS):
        ax = axes[gi, ci]
        ax.axhspan(-0.585, 0.585, color="#f4f4f4", zorder=0)
        ax.axhline(0, color="#999", lw=0.9, ls="--", zorder=1)
        ax.axhline(FLOOR + 0.12, color="#ddd", lw=0.6, zorder=0)          # floor strip divider
        vlists = [np.array(fc[gene][ct][g]) for g, _, _ in GRPS]
        good = [i for i, v in enumerate(vlists) if len(v) >= 2]
        if good:
            parts = ax.violinplot([vlists[i] for i in good], positions=good, widths=0.78,
                                  showmeans=False, showextrema=False)
            for pc, i in zip(parts["bodies"], good):
                col = GRPS[i][2]; pc.set_facecolor(col); pc.set_alpha(0.32); pc.set_edgecolor(col); pc.set_linewidth(1.1)
        for pos, (g, _, col) in enumerate(GRPS):
            v = vlists[pos]
            if len(v):
                jit = (np.random.RandomState(gi*20+ci*3+pos).rand(len(v)) - 0.5) * 0.20
                ax.scatter(pos + jit, v, s=12, color=col, edgecolor="white", linewidth=0.4, zorder=3)
                ax.hlines(v.mean(), pos-0.22, pos+0.22, color="#222", lw=1.5, zorder=4)
            nu = und[gene][ct][g]                                          # undetected -> hollow floor markers
            if nu:
                jit = (np.random.RandomState(99+gi*20+ci*3+pos).rand(nu) - 0.5) * 0.20
                ax.scatter(pos + jit, np.full(nu, FLOOR), s=13, facecolors="none",
                           edgecolors=col, linewidth=0.8, zorder=3)
        # KW on linear donor means (all donors, incl undetected/0)
        gl = [np.array(lin[gene][ct][g]) for g, _, _ in GRPS]
        if all(len(x) >= 3 for x in gl):
            p = kruskal(*gl).pvalue
            y = YLIM[1] - 0.30
            ax.plot([0, 0, 2, 2], [y, y+.10, y+.10, y], lw=1.0, color="#555")
            ax.plot([1, 1], [y+.10, y+.05], lw=1.0, color="#555")
            ax.text(1, y+.12, f"KW {pstar(p)}  p={p:.2g}", ha="center", va="bottom", fontsize=6.8)
        ax.set_ylim(*YLIM); ax.set_xlim(-.6, 2.6); ax.set_xticks([0, 1, 2])
        ax.set_xticklabels([f"{lab}\n(n={len(fc[gene][ct][g])+und[gene][ct][g]})"
                            for g, lab, _ in GRPS], fontsize=6.6)
        ax.tick_params(labelsize=7, length=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        if ci == 0: ax.set_ylabel(f"{gene} / {alias}\n\nlog$_2$FC vs Braak I–II", fontsize=8.6, fontweight="bold")
        if gi == 0: ax.set_title(ctlab, fontsize=10.5, fontweight="bold", pad=6)

fig.suptitle("snRNA UPR-sensor transcripts are unchanged across Braak stage (donor-level; Mathys 2019, n = donors)",
             fontsize=12.5, fontweight="bold", y=1.02)
fig.text(0.5, 0.965, "Per-donor log$_2$FC vs the Braak I–II detected-donor mean of each cell type (common scale). "
         "Filled points = detected donors · hollow points on floor = gene not detected in that donor's cell type · "
         "bar = mean of detected · dashed line = no change · grey band = ±1.5-fold · Kruskal–Wallis (linear means, all donors).\n"
         "Sensor mRNA level is not a readout of sensor activation (Reviewer 2, comment 1); as internal reference genes they show no consistent shift.",
         ha="center", va="top", fontsize=8, color="#555")
fig.subplots_adjust(top=0.885, bottom=0.05, left=0.07, right=0.995, hspace=0.40, wspace=0.12)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q1_sensor_expression_snRNA_donorlevel_violins.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R2Q1_sensor_expression_snRNA_donorlevel_violins.{png,svg}")
