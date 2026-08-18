#!/usr/bin/env python3
"""SEA-AD (independent) Braak V,VI vs 0,I,II — manuscript Figure-4 style.
Per cell type: distribution of per-gene log2FC of the Response-to-ER-stress set (green violin, red mean,
blue quartiles, black points) — SAME contrast/units as our Fig 4. Significance is DONOR-LEVEL
(Mann-Whitney on the per-donor UPR-set score, high-Braak n=49 vs low-Braak n=6), not a gene sign test."""
import os, csv, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from scipy.stats import mannwhitneyu

OUT = os.path.dirname(os.path.abspath(__file__))
FC = os.path.join(OUT, "SEAAD_Braak56vs012_log2FC_perCellType.csv")
PB = os.path.join(OUT, "SEAAD_Braak56vs012_donor_pseudobulk.csv")
SENS = {"EIF2AK3", "ERN1", "ATF6"}
CELLS = [("Ex", "Excitatory\nneurons"), ("In", "Inhibitory\nneurons"), ("Ast", "Astrocytes"),
         ("Mic", "Microglia"), ("Oli", "Oligodendrocytes"), ("OPC", "Oligodendrocyte\nprogenitor cells")]
GREEN = "#4a9e3e"; RED = "#e2231a"; BLUE = "#2b3a8f"

# per-gene log2FC (UPR-target)
fc = {c: [] for c, _ in CELLS}
for r in csv.DictReader(open(FC)):
    if r["gene_set"] == "UPR-target" and r["cell_type"] in fc:
        try: fc[r["cell_type"]].append(float(r["log2FC_Braak56_vs_012"]))
        except ValueError: pass

# donor-level UPR-set score per donor (mean over UPR-target genes), high vs low
rows = list(csv.DictReader(open(PB)))
genecols = [g for g in rows[0].keys() if g not in ("cell_type", "donor", "braak_group") and g not in SENS]
donor = {c: {"high": [], "low": []} for c, _ in CELLS}
for r in rows:
    if r["cell_type"] not in donor: continue
    vals = [float(r[g]) for g in genecols if r[g] not in ("", "nan")]
    if vals: donor[r["cell_type"]][r["braak_group"]].append(st.mean(vals))

def pstar(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"

lim = max(abs(v) for c, _ in CELLS for v in fc[c]) * 1.1
fig, ax = plt.subplots(figsize=(11, 5.8))
ax.axhline(0, color="#000", lw=0.8, ls=":", zorder=1)
YL = (-lim, lim)
for i, (c, lab) in enumerate(CELLS):
    vals = np.array(fc[c])
    parts = ax.violinplot([np.clip(vals, *YL)], positions=[i], widths=0.8, showmeans=False, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(GREEN); pc.set_alpha(0.9); pc.set_edgecolor("#1f1f1f"); pc.set_linewidth(0.8)
    jit = (np.random.RandomState(i).rand(len(vals)) - 0.5) * 0.30
    ax.scatter(i + jit, np.clip(vals, *YL), s=4, color="#111", alpha=0.5, edgecolor="none", zorder=3)
    q1, q3 = np.percentile(vals, [25, 75])
    ax.hlines([q1, q3], i - 0.34, i + 0.34, color=BLUE, lw=1.3, zorder=5)
    ax.hlines(vals.mean(), i - 0.36, i + 0.36, color=RED, lw=3.0, zorder=6)
    # DONOR-LEVEL significance
    hi, lo = donor[c]["high"], donor[c]["low"]
    p = mannwhitneyu(hi, lo, alternative="two-sided").pvalue if len(hi) >= 2 and len(lo) >= 2 else float("nan")
    ax.text(i, YL[1] * 0.80, pstar(p), ha="center", fontsize=13, fontweight="bold")
    ax.text(i, YL[1] * 0.62, f"p={p:.2g}", ha="center", fontsize=6.5, color="#555")
ax.set_xlim(-0.6, len(CELLS) - 0.4); ax.set_ylim(*YL)
ax.set_xticks(range(len(CELLS))); ax.set_xticklabels([l for _, l in CELLS], fontsize=9.5)
ax.set_ylabel("log$_2$FC  Braak V,VI / 0,I,II\n(Response to ER stress, genes=257)", fontsize=10.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(facecolor=GREEN, edgecolor="#1f1f1f", label="SEA-AD per-gene log$_2$FC"),
                   plt.Line2D([0], [0], color=RED, lw=3, label="set mean"),
                   plt.Line2D([0], [0], color=BLUE, lw=1.3, label="quartiles")],
          frameon=False, fontsize=8.5, loc="lower right")
fig.suptitle("SEA-AD (independent) — ER-stress set trends down in neurons with Braak stage (donor-level: n.s.)",
             fontsize=12, fontweight="bold", y=0.99)
fig.text(0.5, 0.925, "Allen · MTG · Braak V,VI (n=49 donors) vs 0,I,II (n=6) · per-gene log2FC direction is neuron-down / glia-up, "
         "but the DONOR-LEVEL test (Mann–Whitney on per-donor UPR-set score) is not significant — n=6 low-Braak limits power.",
         ha="center", fontsize=8.0, color="#555")
fig.subplots_adjust(top=0.86, bottom=0.11, left=0.10, right=0.98)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q9_SEAAD_Braak_figure.{ext}"), dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved R2Q9_SEAAD_Braak_figure.{png,svg}")
for c, _ in CELLS:
    hi, lo = donor[c]["high"], donor[c]["low"]
    p = mannwhitneyu(hi, lo, alternative="two-sided").pvalue if len(hi) >= 2 and len(lo) >= 2 else float("nan")
    print(f"  {c:4s} gene-log2FC mean={np.mean(fc[c]):+.4f} {100*sum(v<0 for v in fc[c])//len(fc[c])}% down | "
          f"donor-level: high {np.mean(hi):+.3f} (n={len(hi)}) vs low {np.mean(lo):+.3f} (n={len(lo)}) p={p:.2g}")
