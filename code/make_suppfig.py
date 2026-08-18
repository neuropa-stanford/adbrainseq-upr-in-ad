#!/usr/bin/env python3
"""R1#1 Supplementary Figure: concordant ER-stress genes across Mizuno & Nativio cohorts.
Panel A = GO:BP enrichment (-log10 p) of the 118 concordant genes (g:Profiler).
Panel B = log2FC heatmap of representative conserved UPR genes in both cohorts.
Inputs computed from submitted SuppD2/D3/D4. Outputs: TIFF (ANC) + PNG + EPS, 300 dpi.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.linewidth": 0.8, "axes.edgecolor": "#333333",
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})

# ---- Panel A data ----
goA = []
with open(os.path.join(HERE, "FIG_GO_118concordant_bar.csv")) as f:
    for r in csv.DictReader(f):
        goA.append((r["term"], float(r["neg_log10_p"]), int(r["intersection"]), int(r["term_size"])))
goA.sort(key=lambda x: x[1])  # ascending so largest on top after barh
labelsA = [t[0].replace("endoplasmic reticulum", "ER") for t in goA]
valsA = [t[1] for t in goA]
annA = [f"{t[2]}/{t[3]}" for t in goA]

# ---- Panel B data ----
genes, branch, miz, nat = [], [], [], []
with open(os.path.join(HERE, "FIG_representative_conserved_UPR_genes.csv")) as f:
    for r in csv.DictReader(f):
        genes.append(r["gene"]); branch.append(r["branch_function"])
        miz.append(float(r["Mizuno_log2FC"])); nat.append(float(r["Nativio_log2FC"]))
# order rows by functional group then by mean magnitude
order = sorted(range(len(genes)), key=lambda i: (branch[i], (miz[i]+nat[i])/2))
genes = [genes[i] for i in order]; branch = [branch[i] for i in order]
miz = [miz[i] for i in order]; nat = [nat[i] for i in order]
M = np.array([miz, nat]).T  # rows=genes, cols=[Mizuno,Nativio]

fig = plt.figure(figsize=(7.2, 8.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.55], hspace=0.28)

# ===== Panel A =====
axA = fig.add_subplot(gs[0])
yA = np.arange(len(labelsA))
axA.barh(yA, valsA, color="#2c7fb8", height=0.62, zorder=3)
for i, (v, a) in enumerate(zip(valsA, annA)):
    axA.text(v + max(valsA)*0.012, i, a, va="center", ha="left", fontsize=7.2, color="#444")
axA.set_yticks(yA); axA.set_yticklabels(labelsA, fontsize=8.2)
axA.set_xlabel("GO:BP enrichment  $-\\log_{10}(p_{adj})$", fontsize=8.6)
axA.set_xlim(0, max(valsA)*1.16)
axA.grid(axis="x", color="#e5e5e5", lw=0.7, zorder=0)
axA.set_axisbelow(True)
for s in ("top", "right"): axA.spines[s].set_visible(False)
axA.set_title("A   Concordant ER-stress genes (n=118) recapitulate the UPR/ERAD program",
              loc="left", fontsize=9.6, fontweight="bold", pad=8)

# ===== Panel B =====
axB = fig.add_subplot(gs[1])
# convention: downregulation (negative log2FC) = BLUE, up = red -> use coolwarm (low=blue)
norm = TwoSlopeNorm(vmin=-0.75, vcenter=0.0, vmax=0.05)
im = axB.imshow(M, cmap="coolwarm", norm=norm, aspect="auto")
axB.set_xticks([0, 1]); axB.set_xticklabels(["Mizuno\n(Japanese)", "Nativio\n(American)"], fontsize=8.4)
axB.set_yticks(np.arange(len(genes)))
# gene + functional annotation folded into the y-tick label (no right-side collision)
axB.set_yticklabels([f"{g}  ·  {b}" for g, b in zip(genes, branch)], fontsize=7.6)
for i in range(len(genes)):
    for j in range(2):
        axB.text(j, i, f"{M[i,j]:.2f}", ha="center", va="center", fontsize=6.8,
                 color="white" if M[i, j] < -0.30 else "#222")
axB.set_xlim(-0.5, 1.5)
axB.set_title("B   Representative conserved UPR genes — $\\log_2$FC (AD vs low-Braak)",
              loc="left", fontsize=9.6, fontweight="bold", pad=8)
axB.tick_params(length=0)
for s in axB.spines.values(): s.set_visible(False)
cbar = fig.colorbar(im, ax=axB, fraction=0.04, pad=0.03, aspect=24)
cbar.set_label("$\\log_2$FC\n(blue = down)", fontsize=7.2); cbar.ax.tick_params(labelsize=6.8)

fig.suptitle("Supplementary Figure — Cross-cohort concordance of ER-stress transcriptional downregulation",
             fontsize=10.2, fontweight="bold", y=0.985)
fig.text(0.5, 0.008,
         "118/195 (Nativio) and 118/181 (Mizuno) ER-stress GO:0034976 genes concordantly downregulated in both cohorts. "
         "Enrichment: g:Profiler GO:BP, g:SCS<0.05.",
         ha="center", fontsize=6.6, color="#666")

for ext in ("tif", "png", "eps"):
    out = os.path.join(HERE, f"SuppFig_R1Q1_concordant_ERstress.{ext}")
    kw = {"dpi": 300, "bbox_inches": "tight"}
    if ext == "tif":
        kw["pil_kwargs"] = {"compression": "tiff_lzw"}
    fig.savefig(out, **kw)
    print("wrote", out)
print("done")
