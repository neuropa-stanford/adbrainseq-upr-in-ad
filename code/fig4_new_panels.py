#!/usr/bin/env python3
"""Figure 4 NEW panels (revision, R2.2/R2.3/R3): donor-level as the headline inference.
  4b  neuron ER-stress module score per donor by Braak stage (Ex, In)
  4c  glia   ER-stress module score per donor by Braak stage (Ast, Mic, Oli, OPC)
       -> per-donor points + group mean + 95% CI; donor-level Cohen's d (V-VI vs I-II), adj-p, n=donors.
  4g  specificity: Cohen's d (V-VI vs I-II) per cell type x gene set incl. an expression-matched
       internal control (mRNA transport). Neurons: UPR sets AND control move together (global shift);
       only oligodendrocyte IRE1 exceeds the control (d=1.45 vs 0.76).
Data: R2Q2_scores.npy (per-donor module scores) + R2Q2_donorlevel_module_scores.csv (d/CI/p/q).
Vector, Arial. Standalone preview of the new content before assembling the full Figure 4."""
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

# ---- donors -> Braak group ----
info = {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    info.setdefault(r["Subject"], r)
braak = {s: int(float(info[s]["braaksc"])) for s in info}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")
GLAB = {"low": "I–II", "int": "III–IV", "late": "V–VI"}

# ---- per-donor ER-stress module scores ----
scores = {}
for key, donors, sc in np.load(os.path.join(OUT, "R2Q2_scores.npy"), allow_pickle=True):
    ct, sn = key.split("|")
    scores[(ct, sn)] = (donors.split(","), np.array([float(x) for x in sc.split(",")]))

# ---- summary stats (d, CI, p_adj) ----
summ = {}
for r in csv.DictReader(open(os.path.join(OUT, "R2Q2_donorlevel_module_scores.csv"))):
    summ[(r["cell_type"], r["gene_set"], r["comparison"])] = r

CELL = {"Ex": "Excitatory\nneurons", "In": "Inhibitory\nneurons", "Ast": "Astrocytes",
        "Mic": "Microglia", "Oli": "Oligodendrocytes", "Opc": "OPCs"}
COL = {"Ex": "#8CBF43", "In": "#33A45C", "Ast": "#F0806A", "Mic": "#28C2D4", "Oli": "#4F86C6", "Opc": "#9C77B4"}
NEUR, GLIA = ["Ex", "In"], ["Ast", "Mic", "Oli", "Opc"]
ES = "ER-stress (260)"

def draw_modulescore(ax, ct):
    donors, sc = scores[(ct, ES)]
    g = {k: np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == k]) for k in ("low", "int", "late")}
    ax.axhline(0, color="#bbb", lw=0.7, ls=":")
    for xi, k in enumerate(("low", "int", "late")):
        v = g[k]; jit = (np.random.RandomState(xi).rand(len(v)) - 0.5) * 0.28
        ax.scatter(np.full(len(v), xi) + jit, v, s=14, color=COL[ct], alpha=0.75, edgecolor="#333", linewidth=0.3, zorder=3)
        m = v.mean(); se = v.std(ddof=1) / math.sqrt(len(v)); ci = 1.96 * se
        ax.hlines(m, xi - 0.22, xi + 0.22, color="#111", lw=2.2, zorder=4)
        ax.vlines(xi, m - ci, m + ci, color="#111", lw=1.2, zorder=4)
        ax.text(xi, ax.get_ylim()[0], f"n={len(v)}", ha="center", va="bottom", fontsize=6.5, color="#666")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels([GLAB[k] for k in ("low", "int", "late")], fontsize=8.5)
    ax.set_xlim(-0.5, 2.5); ax.set_title(CELL[ct].replace("\n", " "), fontsize=9.5, fontweight="bold", color=COL[ct])
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    # donor-level headline: late vs low
    s = summ.get((ct, ES, "late vs low"))
    if s:
        d = float(s["cohens_d"]); q = float(s["p_adj_BH"])
        star = "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "ns"
        ax.text(0.5, 0.97, f"d={d:+.2f}, q={q:.2g} {star}", transform=ax.transAxes,
                ha="center", va="top", fontsize=7.4, fontweight="bold")

# ===== figure: 4b + 4c (row of 6) on top, 4g (full width) below =====
fig = plt.figure(figsize=(11.5, 7.4))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.85], hspace=0.42, wspace=0.42,
                      left=0.075, right=0.985, top=0.90, bottom=0.10)
axes = [fig.add_subplot(gs[0, j]) for j in range(6)]
for ax, ct in zip(axes, NEUR + GLIA):
    draw_modulescore(ax, ct)
axes[0].set_ylabel("ER-stress module score\n(per donor, z)", fontsize=9)
# panel letters + neuron/glia brackets
fig.text(0.02, 0.965, "b", fontsize=15, fontweight="bold")
fig.text(0.36, 0.965, "c", fontsize=15, fontweight="bold")
fig.text(0.145, 0.925, "Neurons (down)", ha="center", fontsize=9, fontweight="bold", color="#555")
fig.text(0.66, 0.925, "Glia (up)", ha="center", fontsize=9, fontweight="bold", color="#555")
fig.text(0.5, 0.945, "Donor-level module score (donor = unit of inference; n = 10 / 21 / 17 donors)  — gene-level Wilcoxon shown in 4b,c violins is secondary",
         ha="center", fontsize=8, color="#666")

# ---- 4g specificity ----
axg = fig.add_subplot(gs[1, :])
SETS = [("ER-stress (260)", "#3b78b5"), ("PERK (31)", "#e8862e"), ("IRE1 (32)", "#4aa24a"),
        ("ATF6 (74)", "#c0392b"), ("ERAD (75)", "#8e6bb0"), ("Internal control: mRNA transport (91)", "#333333")]
SLAB = {"ER-stress (260)": "ER-stress", "PERK (31)": "PERK", "IRE1 (32)": "IRE1", "ATF6 (74)": "ATF6",
        "ERAD (75)": "ERAD", "Internal control: mRNA transport (91)": "control (mRNA transp.)"}
order = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
nS = len(SETS); bw = 0.13
for si, (sn, c) in enumerate(SETS):
    ds = [float(summ[(ct, sn, "late vs low")]["cohens_d"]) if (ct, sn, "late vs low") in summ else 0 for ct in order]
    xpos = np.arange(len(order)) + (si - (nS - 1) / 2) * bw
    axg.bar(xpos, ds, bw, color=c, label=SLAB[sn], edgecolor="none")
axg.axhline(0, color="#111", lw=0.9)
axg.set_xticks(range(len(order))); axg.set_xticklabels([CELL[c].replace("\n", " ") for c in order], fontsize=9)
axg.set_ylabel("Cohen's d\n(Braak V–VI vs I–II)", fontsize=9)
axg.legend(frameon=False, fontsize=7.5, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.14), handletextpad=0.4, columnspacing=1.0)
for sp in ("top", "right"): axg.spines[sp].set_visible(False)
fig.text(0.02, 0.44, "g", fontsize=15, fontweight="bold")
axg.text(0.5, -0.24, "Neurons: UPR sets and the internal control move together (global shift).   "
         "Oligodendrocytes: IRE1 exceeds the matched control (d = 1.45 vs 0.76).",
         transform=axg.transAxes, ha="center", fontsize=8, color="#555")

for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(os.path.dirname(OUT), f"Figure4_newpanels_bcg.{ext}"),
                dpi=200 if ext == "png" else None, facecolor="white")
print("saved Figure4_newpanels_bcg.{png,svg,pdf}")
