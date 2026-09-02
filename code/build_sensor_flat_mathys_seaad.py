#!/usr/bin/env python3
"""UPR-sensor 'flat gene' donor-level log2FC violins — Mathys (top) + SEA-AD (bottom), one PDF.
Sensors are internal reference genes (mRNA level != activation): they stay flat across Braak stage.
Rows: EIF2AK3/PERK, ERN1/IRE1, ATF6/ATF6 ; cols: Excitatory, Inhibitory, Microglia, Oligodendrocyte.
Y = per-donor log2FC vs the Braak I-II donor mean of that cell type. Point = donor, bar = mean,
grey band = +-1.5-fold, dashed = no change. Same cell-type colours/labels/style as Figures 4-5.
Mathys per-donor values from R2Q1_snRNA_sensor_donormeans.npy; SEA-AD from SEAAD_donor_pseudobulk_UPRunion.csv."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
from matplotlib.gridspec import GridSpec
from scipy.stats import f_oneway

DRD = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(DRD); SP = DRD
SENSORS = [("EIF2AK3", "PERK"), ("ERN1", "IRE1"), ("ATF6", "ATF6")]
CELLS = [("Ex", "Excitatory", "#746fbd"), ("In", "Inhibitory", "#c2b7de"),
         ("Mic", "Microglia", "#cdea9a"), ("Oli", "Oligodendrocyte", "#9ad4e7")]
GRPS = ["low", "mid", "high"]; GLAB = ["Braak\nI-II", "Braak\nIII-IV", "Braak\nV-VI"]
YLIM = (-2.6, 2.2); FLOOR = YLIM[0] + 0.18; BAND = 0.585
def pstar(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

# ---- Mathys: per-donor log2FC from the .npy (col5 linear, col6 log2FC); grp int->mid, late->high ----
MAP = {"low": "low", "int": "mid", "late": "high"}
rows = np.load(os.path.join(DRD, "R2Q1_snRNA_sensor_donormeans.npy"), allow_pickle=True)
def blank(): return {g: {c: {gr: [] for gr in GRPS} for c, _, _ in CELLS} for g, _ in SENSORS}
Mfc, Mlin, Mund = blank(), blank(), {g: {c: {gr: 0 for gr in GRPS} for c, _, _ in CELLS} for g, _ in SENSORS}
for r in rows:
    g, ct, grp = r[0], r[1], MAP.get(r[4])
    if g not in Mfc or ct not in Mfc[g] or grp is None: continue
    Mlin[g][ct][grp].append(float(r[5]))
    f = r[6]
    if f is not None and np.isfinite(float(f)): Mfc[g][ct][grp].append(float(f))
    else: Mund[g][ct][grp] += 1

# ---- SEA-AD: per-donor log2FC vs Braak I-II mean, from UPRunion pseudobulk (linearise via expm1) ----
sr = list(csv.DictReader(open(os.path.join(DRD, "SEAAD_donor_pseudobulk_UPRunion.csv"))))
Sfc, Slin, Sund = blank(), blank(), {g: {c: {gr: 0 for gr in GRPS} for c, _, _ in CELLS} for g, _ in SENSORS}
for g, _ in SENSORS:
    for c, _, _ in CELLS:
        byg = {gr: [np.expm1(float(r[g])) for r in sr if r["cell_type"] == c and r["braak_group"] == gr] for gr in GRPS}
        low_mean = np.mean(byg["low"]) if byg["low"] else np.nan
        for gr in GRPS:
            for lin in byg[gr]:
                Slin[g][c][gr].append(lin)
                if lin > 0 and low_mean > 0: Sfc[g][c][gr].append(np.log2(lin / low_mean))
                else: Sund[g][c][gr] += 1

DATA = [("Mathys 2019", Mfc, Mlin, Mund),
        ("SEA-AD", Sfc, Slin, Sund)]

fig = plt.figure(figsize=(12.0, 15.6))
gs = GridSpec(6, 4, figure=fig, left=0.085, right=0.99, top=0.935, bottom=0.035, hspace=0.55, wspace=0.14)
panel_axes = {}
for bk, (hdr, fc, lin, und) in enumerate(DATA):
    for si, (gene, alias) in enumerate(SENSORS):
        rr = bk*3 + si
        for ci, (ct, ctlab, col) in enumerate(CELLS):
            ax = fig.add_subplot(gs[rr, ci]); panel_axes[(rr, ci)] = ax
            ax.axhspan(-BAND, BAND, color="#f4f4f4", zorder=0)
            ax.axhline(0, color="#999", lw=0.9, ls="--", zorder=1)
            vlists = [np.array(fc[gene][ct][gr]) for gr in GRPS]
            good = [i for i, v in enumerate(vlists) if len(v) >= 2]
            if good:
                for pc, i in zip(ax.violinplot([vlists[i] for i in good], positions=good, widths=0.78,
                                                showmeans=False, showextrema=False)["bodies"], good):
                    pc.set_facecolor(col); pc.set_alpha(0.42); pc.set_edgecolor(col); pc.set_linewidth(1.1)
            for pos, gr in enumerate(GRPS):
                v = vlists[pos]
                if len(v):
                    jit = (np.random.RandomState(bk*80+si*20+ci*3+pos).rand(len(v)) - 0.5) * 0.22
                    vc = np.clip(v, YLIM[0]+0.06, YLIM[1]-0.06)   # keep every dot inside the panel (some log2FC exceed ±2)
                    ax.scatter(pos + jit, vc, s=12, color=col, edgecolor="white", linewidth=0.4, zorder=3, clip_on=True, rasterized=True)
                    ax.hlines(v.mean(), pos-0.24, pos+0.24, color="#222", lw=1.6, zorder=4)
                nu = und[gene][ct][gr]
                if nu:
                    jit = (np.random.RandomState(9+bk*80+si*20+ci*3+pos).rand(nu) - 0.5) * 0.22
                    ax.scatter(pos + jit, np.full(nu, FLOOR), s=13, facecolors="none", edgecolors=col, linewidth=0.8, zorder=3, clip_on=True, rasterized=True)
            gl = [np.array(lin[gene][ct][gr]) for gr in GRPS]
            if all(len(x) >= 3 for x in gl):
                st = pstar(f_oneway(*gl).pvalue)
                if st != "ns":   # one-way ANOVA across the 3 Braak groups; show only significant stars
                    ax.text(1, YLIM[1]-0.08, st, ha="center", va="top", fontsize=11, fontweight="bold")
            ax.set_ylim(*YLIM); ax.set_xlim(-.6, 2.6); ax.set_xticks([0, 1, 2])
            ax.set_xticklabels([f"{GLAB[i]}\n(n={len(fc[gene][ct][gr])+und[gene][ct][gr]})" for i, gr in enumerate(GRPS)], fontsize=6.6)
            ax.tick_params(labelsize=7, length=0)
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            if ci == 0: ax.set_ylabel(f"{gene} / {alias}\n\nlog$_2$FC vs Braak I-II", fontsize=8.4, fontweight="bold")
            else: ax.tick_params(labelleft=False)
            if si == 0: ax.set_title(ctlab, fontsize=10.5, fontweight="bold", pad=6)
    y = 0.952 if bk == 0 else 0.475
    fig.text(0.085, y, hdr, fontsize=11.5, fontweight="bold", ha="left")

# panel letters A..X at each panel's top-left corner (row-major: A-D row0 ... U-X row5)
import string
fig.canvas.draw()
_LET = string.ascii_uppercase
_k = 0
for _rr in range(6):
    for _ci in range(4):
        _p = panel_axes[(_rr, _ci)].get_position()
        fig.text(max(_p.x0 - 0.020, 0.002), min(_p.y1 + 0.004, 0.999), _LET[_k],
                 fontsize=12, fontweight="bold", va="bottom", ha="left")
        _k += 1

pdf = os.path.join(OUTDIR, "Figure_UPRsensor_flat_Mathys_SEAAD_20260827.pdf")
fig.savefig(pdf, facecolor="white", dpi=600)   # dot layer rasterized at 600 dpi -> clips cleanly in every viewer
import shutil
shutil.copyfile(pdf, os.path.join(OUTDIR, "Figure_UPRsensor_flat_Mathys_SEAAD_20260827.ai"))
print("wrote", pdf, "| panel letters A-X")
