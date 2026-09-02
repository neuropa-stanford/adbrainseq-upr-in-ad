#!/usr/bin/env python3
"""Figure 5 I-P: SEA-AD reproduction of the Mathys A-H UPR-branch panels, matched to the original .ai.
Each branch (PERK/IRE1/ATF6/ERAD) = neuron panel (Ex,In) + glia panel (Mic,OligD) sharing the branch axis
(glia hides y-numbers). Columns positioned to overlie the original A-H columns exactly (neuron plots wider
than glia, like the .ai). Metric = gene-level log2(mean_group / mean_low) → shared -2..2 INTEGER axis (data
clipped to +-2; stars/brackets in the 2.2-2.7 headroom). Per-violin star hugs each violin, bracket sits
clearly above. Sizes measured off the .ai: x-labels 7.5, y-ticks 8.2, y-title 9.0, letters 18, stars 7."""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import FuncFormatter
from scipy.stats import wilcoxon, gaussian_kde

DRD = os.path.dirname(os.path.abspath(__file__)); SP = DRD
PB = os.path.join(DRD, "SEAAD_donor_pseudobulk_UPRunion.csv")
BR = [("PERK", "geneset_PERK.txt", "I", "J"), ("IRE1", "geneset_IRE1.txt", "K", "L"),
      ("ATF6", "geneset_ATF6.txt", "M", "N"), ("ERAD", "geneset_ERAD.txt", "O", "P")]
NEUR = [("Ex", "Ex", "#746fbd"), ("In", "In", "#c2b7de")]
GLIA = [("Mic", "Mic", "#cdea9a"), ("Oli", "OligD", "#9ad4e7")]
RED, BLUE = "#e2231a", "#2b3a8f"
FS_X, FS_YTICK, FS_YTITLE, FS_STAR, FS_LET, FS_HDR = 7.5, 8.2, 9.0, 7.0, 18.0, 9.5
BW = 0.42
def pstar(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

rows = list(csv.DictReader(open(PB)))
allg = set(k for k in rows[0].keys() if k not in ("cell_type", "donor", "braak_group", "braak_num"))
byct = {c: {"low": [], "mid": [], "high": []} for c in ("Ex", "In", "Mic", "Oli")}
for r in rows:
    if r["cell_type"] in byct and r["braak_group"] in byct[r["cell_type"]]:
        byct[r["cell_type"]][r["braak_group"]].append(r)
def load_set(fn): return [g.strip() for g in open(os.path.join(DRD, fn)) if g.strip() and g.strip() in allg]
def gene_fc(c, grp, genes):   # true log2 fold-change = log2(mean_group / mean_low)
    lo = np.array([[float(rr[g]) for g in genes] for rr in byct[c]["low"]]).mean(0)
    gg = np.array([[float(rr[g]) for g in genes] for rr in byct[c][grp]]).mean(0)
    return np.log2((gg + 1e-9) / (lo + 1e-9))

FC, GN = {}, {}
for sn, fn, _, _ in BR:
    genes = load_set(fn); GN[sn] = len(genes)
    for c in ("Ex", "In", "Mic", "Oli"):
        for g in ("mid", "high"): FC[(sn, c, g)] = gene_fc(c, g, genes)
YL = (-2.0, 2.0); TICKS = np.array([-2, -1, 0, 1, 2])
fmt = FuncFormatter(lambda v, _: ("%g" % v))

# columns overlie the original A-H columns (composed at page x-offset 11, scale 1.0 → frac=(page_x-11)/590)
COL = [(0.0692, 0.189), (0.263, 0.169), (0.527, 0.189), (0.719, 0.170)]   # (left, width) col0..col3
LETX = [0.088, 0.273, 0.546, 0.732]
PH = 0.275; ROWB = [0.600, 0.150]; LET_DY = 0.018

fig = plt.figure(figsize=(590/72, 380/72))
def panel(col_i, y, cells, sn, show_y, letter, rng):
    left, pw = COL[col_i]
    ax = fig.add_axes([left, y, pw, PH]); ax.axhline(0, color="#000", lw=0.7, ls=":")
    for i, (c, clab, col) in enumerate(cells):
        tops = []
        for grp in ("mid", "high"):
            xp = i + (0.5 if grp == "high" else -0.5) * BW
            raw = FC[(sn, c, grp)]; v = np.clip(raw, *YL)
            for pc in ax.violinplot([v], positions=[xp], widths=BW*0.9, showextrema=False)["bodies"]:
                pc.set_facecolor(col); pc.set_alpha(0.95); pc.set_edgecolor("#1f1f1f"); pc.set_linewidth(0.5)
            if len(np.unique(v)) > 1:
                dens = gaussian_kde(v)(v); hw = (BW*0.44) * dens / dens.max()
            else: hw = np.zeros(len(v))
            ax.scatter(xp + (rng.rand(len(v))*2-1)*hw*0.85, v, s=1.3, color="#222", alpha=0.4, edgecolor="none", zorder=4)
            q1, q3 = np.percentile(v, [25, 75]); ax.hlines([q1, q3], xp-BW*.46, xp+BW*.46, color=BLUE, lw=0.9, zorder=5)
            ax.hlines(v.mean(), xp-BW*.5, xp+BW*.5, color=RED, lw=1.9, zorder=6)
            sy = min(v.max() + 0.16, 1.55)
            st = pstar(wilcoxon(raw).pvalue)
            if st != "ns": ax.text(xp, sy, st, ha="center", va="bottom", fontsize=FS_STAR)
            tops.append(sy)
        pb = wilcoxon(FC[(sn, c, "mid")], FC[(sn, c, "high")]).pvalue
        if pstar(pb) != "ns":
            yb = max(tops) + 0.62
            ax.plot([i-.5*BW, i-.5*BW, i+.5*BW, i+.5*BW], [yb-0.14, yb, yb, yb-0.14], color="#000", lw=0.8)
            ax.text(i, yb+0.02, pstar(pb), ha="center", va="bottom", fontsize=FS_STAR)
    ax.set_xlim(-0.6, len(cells)-0.4); ax.set_ylim(-2.2, 2.7); ax.set_xticks([])
    ax.set_yticks(TICKS); ax.yaxis.set_major_formatter(fmt)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for i, (c, clab, col) in enumerate(cells):
        for grp, dx in [("mid", -0.5*BW), ("high", 0.5*BW)]:
            ax.text(i+dx, -0.035, f"{clab} {'III,IV' if grp=='mid' else 'V,VI'}", ha="right", va="top",
                    fontsize=FS_X, rotation=45, transform=tr, rotation_mode="anchor")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    if show_y: ax.tick_params(axis="y", labelsize=FS_YTICK)
    else: ax.tick_params(axis="y", labelleft=False, length=0)
    fig.text(LETX[col_i], y + PH + LET_DY, letter, fontsize=FS_LET, fontweight="bold", va="bottom", ha="left")
    return ax

for bi, (sn, fn, Ln, Lg) in enumerate(BR):
    row, cp = bi // 2, bi % 2; rowb = ROWB[row]; ncol, gcol = cp*2, cp*2+1; rng = np.random.RandomState(3)
    axN = panel(ncol, rowb, NEUR, sn, True, Ln, rng)
    panel(gcol, rowb, GLIA, sn, False, Lg, rng)
    axN.set_ylabel(f"Log$_2$fc Braak III/IV, V/VI / Braak I/II\n({sn} gene set, genes={GN[sn]})", fontsize=FS_YTITLE)

fig.text(0.015, 0.975, "SEA-AD (Allen, middle temporal gyrus) — independent reproduction of Figure 5 (A-H)",
         fontsize=FS_HDR, fontweight="bold", ha="left", va="top")
pdf = os.path.join(SP, "fig5_IJKL_seaad.pdf")
fig.savefig(pdf, facecolor="white"); print("wrote", pdf, "| genes", GN)
