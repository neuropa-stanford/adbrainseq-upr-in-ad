#!/usr/bin/env python3
"""snRNA sensor-gene DONOR-LEVEL log2FC by BRANCH colour, manuscript style (Braak III/IV and V/VI
each vs Braak I/II). EIF2AK3/PERK = cyan, ERN1/IRE1 = maroon, ATF6 = blue (III/IV lighter, V/VI
darker shade of the branch colour). Rows = branch, cols = Ex, In, Ast, Mic, Oli, OPC.
Each point = one donor; filled = detected, hollow floor markers = gene not detected in that donor's
cell type. Bar = mean of detected. Two-sided Mann-Whitney U vs Braak I/II (donor linear means)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from scipy.stats import mannwhitneyu

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
       "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
rows = np.load(os.path.join(OUT, "R2Q1_snRNA_sensor_donormeans.npy"), allow_pickle=True)

# branch: gene, alias, III/IV colour, V/VI colour, edge
BRANCH = [("EIF2AK3", "PERK", "#B6E4EE", "#4FB8CC", "#1F8EA0"),
          ("ERN1", "IRE1", "#DBA6B7", "#9C2A4E", "#611026"),
          ("ATF6", "ATF6", "#AEB7DF", "#3E52A0", "#232E63")]
CELLS = [("Ex", "Excitatory"), ("In", "Inhibitory"), ("Ast", "Astrocyte"),
         ("Mic", "Microglia"), ("Oli", "Oligodendrocyte"), ("Opc", "OPC")]
COND = [("int", "III/IV"), ("late", "V/VI")]
YLIM = (-2.6, 2.2); FLOOR = YLIM[0] + 0.18

fc = {b[0]: {ct: {c: [] for c, _ in COND} for ct, _ in CELLS} for b in BRANCH}
lin = {b[0]: {ct: {g: [] for g in ("low", "int", "late")} for ct, _ in CELLS} for b in BRANCH}
und = {b[0]: {ct: {c: 0 for c, _ in COND} for ct, _ in CELLS} for b in BRANCH}
for r in rows:
    g, ct, grp, v, f = r[0], r[1], r[4], float(r[5]), r[6]
    if g not in fc or ct not in fc[g]: continue
    lin[g][ct][grp].append(v)
    if grp in ("int", "late"):
        if f is not None and np.isfinite(float(f)): fc[g][ct][grp].append(float(f))
        else: und[g][ct][grp] += 1

def pstar(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

fig, axes = plt.subplots(len(BRANCH), len(CELLS), figsize=(16.5, 8.8), sharey="row")
for gi, (gene, alias, c_int, c_late, edge) in enumerate(BRANCH):
    cols = {"int": c_int, "late": c_late}
    for ci, (ct, ctlab) in enumerate(CELLS):
        ax = axes[gi, ci]
        ax.axhspan(-0.585, 0.585, color="#f4f4f4", zorder=0)
        ax.axhline(0, color="#999", lw=0.9, ls="--", zorder=1)
        ax.axhline(FLOOR + 0.12, color="#e0e0e0", lw=0.6, zorder=0)
        for pos, (cond, clab) in enumerate(COND):
            d = np.array(fc[gene][ct][cond]); col = cols[cond]
            if len(d) >= 2:
                parts = ax.violinplot([d], positions=[pos], widths=0.75, showmeans=False, showextrema=False)
                for pc in parts["bodies"]:
                    pc.set_facecolor(col); pc.set_alpha(0.55); pc.set_edgecolor(edge); pc.set_linewidth(1.1)
            if len(d):
                jit = (np.random.RandomState(gi*20+ci*3+pos).rand(len(d)) - 0.5) * 0.20
                ax.scatter(pos + jit, d, s=13, color=col, edgecolor="white", linewidth=0.4, zorder=3)
                ax.hlines(d.mean(), pos-0.24, pos+0.24, color=edge, lw=1.7, zorder=4)
            nu = und[gene][ct][cond]
            if nu:
                jit = (np.random.RandomState(88+gi*20+ci*3+pos).rand(nu) - 0.5) * 0.20
                ax.scatter(pos + jit, np.full(nu, FLOOR), s=13, facecolors="none", edgecolors=edge, linewidth=0.8, zorder=3)
            # MWU vs Braak I/II (linear donor means, all donors)
            a = np.array(lin[gene][ct][cond]); b = np.array(lin[gene][ct]["low"])
            if len(a) >= 2 and len(b) >= 2:
                p = mannwhitneyu(a, b, alternative="two-sided").pvalue
                ax.text(pos, YLIM[1]-0.22, pstar(p), ha="center", va="top", fontsize=7.5, color="#333")
        ax.set_ylim(*YLIM); ax.set_xlim(-.6, 1.6); ax.set_xticks([0, 1])
        ax.set_xticklabels([f"{clab}\nvs I/II\n(n={len(fc[gene][ct][c])+und[gene][ct][c]})"
                            for c, clab in COND], fontsize=6.8)
        ax.tick_params(labelsize=7, length=0)
        for sp in ("top", "right"): ax.spines[sp].set_visible(False)
        if ci == 0: ax.set_ylabel(f"{gene} / {alias}\n\nlog$_2$FC vs Braak I/II", fontsize=8.6, fontweight="bold", color=edge)
        if gi == 0: ax.set_title(ctlab, fontsize=10.5, fontweight="bold", pad=6)

fig.suptitle("snRNA UPR-sensor transcripts by branch — Braak III/IV and V/VI vs I/II (donor-level; Mathys 2019)",
             fontsize=12.5, fontweight="bold", y=1.02)
fig.text(0.5, 0.965, "Per-donor log$_2$FC vs the Braak I/II detected-donor mean of each cell type. Filled = detected donors · "
         "hollow floor markers = gene not detected in that donor's cell type · bar = mean of detected · dashed line = no change · "
         "grey band = ±1.5-fold · two-sided Mann–Whitney U vs Braak I/II (linear means, all donors).",
         ha="center", va="top", fontsize=8, color="#555")
fig.subplots_adjust(top=0.885, bottom=0.055, left=0.07, right=0.995, hspace=0.42, wspace=0.12)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q1_sensor_snRNA_branch.{ext}"), dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved R2Q1_sensor_snRNA_branch.{png,svg}")
