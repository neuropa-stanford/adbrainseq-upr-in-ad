#!/usr/bin/env python3
"""Figure 4 middle row (revision): DONOR-LEVEL ER-stress module score per cell type, by Braak stage.
Ex,In (neurons) + Ast,Mic,Oli,OPC (glia). Per-donor points + group mean + 95% CI; donor-level Cohen's d
(V-VI vs I-II) and BH-adjusted p. Goes UNDER the original gene-level violins (b,c). Vector, Arial.
Panel letters drawn later in assembly overlay."""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

OUT = os.path.dirname(os.path.abspath(__file__))
DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
        "Human AD brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
info = {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    info.setdefault(r["Subject"], r)
braak = {s: int(float(info[s]["braaksc"])) for s in info}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")
GLAB = {"low": "I–II", "int": "III–IV", "late": "V–VI"}
scores = {}
for key, donors, sc in np.load(os.path.join(OUT, "R2Q2_scores.npy"), allow_pickle=True):
    ct, sn = key.split("|"); scores[(ct, sn)] = (donors.split(","), np.array([float(x) for x in sc.split(",")]))
summ = {}
for r in csv.DictReader(open(os.path.join(OUT, "R2Q2_donorlevel_module_scores.csv"))):
    summ[(r["cell_type"], r["gene_set"], r["comparison"])] = r
CELL = {"Ex": "Excitatory\nneurons", "In": "Inhibitory\nneurons", "Ast": "Astrocytes",
        "Mic": "Microglia", "Oli": "Oligodendro-\ncytes", "Opc": "OPCs"}
COL = {"Ex": "#8CBF43", "In": "#33A45C", "Ast": "#F0806A", "Mic": "#28C2D4", "Oli": "#4F86C6", "Opc": "#9C77B4"}
ORDER = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]; ES = "ER-stress (260)"

fig = plt.figure(figsize=(8.1, 3.25))
# spacer column (index 2) puts a gap between the neuron pair (Ex,In) and the glia group; neurons hug the left (near D)
gs = fig.add_gridspec(1, 7, width_ratios=[1, 1, 0.55, 1, 1, 1, 1], wspace=0.5,
                      left=0.115, right=0.995, top=0.90, bottom=0.16)
axes = [fig.add_subplot(gs[0, c]) for c in [0, 1, 3, 4, 5, 6]]
for ax, ct in zip(axes, ORDER):
    donors, sc = scores[(ct, ES)]
    g = {k: np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == k]) for k in ("low", "int", "late")}
    ax.axhline(0, color="#bbb", lw=0.6, ls=":")
    for xi, k in enumerate(("low", "int", "late")):
        v = g[k]; jit = (np.random.RandomState(xi).rand(len(v)) - 0.5) * 0.26
        ax.scatter(np.full(len(v), xi) + jit, v, s=11, color=COL[ct], alpha=0.75, edgecolor="#333", linewidth=0.25, zorder=3)
        m = v.mean(); ci = 1.96 * v.std(ddof=1) / math.sqrt(len(v))
        ax.hlines(m, xi - 0.24, xi + 0.24, color="#111", lw=1.8, zorder=4)
        ax.vlines(xi, m - ci, m + ci, color="#111", lw=1.0, zorder=4)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels([GLAB[k] for k in ("low", "int", "late")], fontsize=9.5, rotation=30, ha="right")
    ax.set_xlim(-0.5, 2.5); ax.tick_params(axis="y", labelsize=9.5)
    ax.set_title(CELL[ct], fontsize=10.5, fontweight="normal", color="black", pad=4)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("ER-stress module score\n(per donor, z)", fontsize=9.5, color="black", labelpad=2)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(os.path.dirname(OUT), f"f4_donor_row.{ext}"),
                dpi=200 if ext == "png" else None, facecolor="white")
print("saved f4_donor_row.{pdf,png}")
