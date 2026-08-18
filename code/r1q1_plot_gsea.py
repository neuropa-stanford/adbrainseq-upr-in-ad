#!/usr/bin/env python3
"""GSEA counterpart of the ORA figure — same 18 GO terms, same datasets.

Panel a: NES heatmap (colour = NES, red = up, blue = down; cell number = genes of the
         set present in the ranked list; ns = not significant after BH).
Panel b: response to ER stress (GO:0034976) NES across the six conditions.
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
rows = list(csv.DictReader(open(os.path.join(OUT, "R1Q1_UPR_GSEA_18terms.csv"))))

CORE = {"GO:0034976", "GO:0030968", "GO:0006986", "GO:0034620", "GO:0035966",
        "GO:0035967", "GO:0036498", "GO:1905897", "GO:1903573", "GO:0036503",
        "GO:0030433"}
DATASETS = ["Thapsigargin", "Mizuno", "Nativio"]
VARIANTS = [("full ranked list", "full list"), ("p<0.05 subset", "p<0.05 subset")]
cols = [(d, v[0], v[1]) for d in DATASETS for v in VARIANTS]
idx = {(r["GO"], r["dataset"], r["variant"]): r for r in rows}
NAMES = {r["GO"]: r["term"] for r in rows}

best = {}
for r in rows:
    if r["padj"] and r["significant"] == "True":
        best[r["GO"]] = min(best.get(r["GO"], 1.0), float(r["padj"]))
ORDER = sorted(NAMES, key=lambda g: (float(idx[(g, "Thapsigargin", "full ranked list")]["K"] or 0) < 10,
                                     best.get(g, 2.0)))

cmap = LinearSegmentedColormap.from_list(
    "nes", ["#1b4f72", "#2c7fb8", "#a6cee3", "#f2f0ec", "#f8b878", "#e0562b", "#8c1c13"])
norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)

fig = plt.figure(figsize=(12.8, 9.6))
gs = fig.add_gridspec(2, 1, height_ratios=[2.85, 1.1], hspace=0.28,
                      left=0.40, right=0.95, top=0.855, bottom=0.135)

ax = fig.add_subplot(gs[0])
for i, go in enumerate(ORDER):
    for j, (ds, var, _) in enumerate(cols):
        r = idx.get((go, ds, var))
        untestable = (not r) or (not r["NES"]) or r["NES"] == "nan" or int(r["K"] or 0) < 10
        if untestable:
            face, txt, col = "#f4f1ec", "untestable", "#a09a92"
        else:
            nes = float(r["NES"])
            sig = r["significant"] == "True"
            face = cmap(norm(nes)) if sig else "#eceae6"
            txt = f"{nes:+.2f}" if sig else "ns"
            col = "white" if sig and abs(nes) > 1.75 else ("#2b2b2b" if sig else "#9a968f")
        ax.add_patch(Rectangle((j, i), 1, 1, facecolor=face, edgecolor="white", lw=1.2))
        if txt and txt != "untestable":
            ax.text(j + .5, i + .5, txt, ha="center", va="center", fontsize=8.6, color=col)
        elif txt and j == 0:
            k = int(r["K"] or 0) if r and r["K"] else 0
            msg = (f"fewer than 10 set genes in the ranked list (n = {k}) — not tested"
                   if k else "term obsolete in the current GO release — not testable")
            ax.text(len(cols) / 2, i + .5, msg, ha="center", va="center",
                    fontsize=8.4, style="italic", color="#a09a92")

ax.set_xlim(0, len(cols)); ax.set_ylim(len(ORDER), 0)
ax.set_yticks(np.arange(len(ORDER)) + .5)
lab = []
for i, g in enumerate(ORDER):
    k = idx.get((g, "Nativio", "full ranked list"), {}).get("K", "")
    lab.append(f"{i+1}.  {NAMES[g]}  ({g}, n={k})" if k and k != "0" else f"{i+1}.  {NAMES[g]}  ({g})")
ax.set_yticklabels(lab, fontsize=9.2)
for t, g in zip(ax.get_yticklabels(), ORDER):
    t.set_color("#111" if g in CORE else "#8a8a8a")
    t.set_fontweight("bold" if g in CORE else "normal")
ax.set_xticks(np.arange(len(cols)) + .5)
ax.set_xticklabels([c[2] for c in cols], fontsize=8.8, style="italic")
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
for j in range(2, len(cols), 2):
    ax.plot([j, j], [0, len(ORDER)], color="#4a4a4a", lw=2.4, zorder=5)
for k, ds in enumerate(DATASETS):
    ax.text(k * 2 + 1, -0.95, ds, ha="center", va="bottom", fontsize=13, fontweight="bold")
ax.text(1, -1.62, "positive control", ha="center", fontsize=9, color="#7a7a7a", style="italic")
ax.text(4, -1.62, "AD brain cohorts", ha="center", fontsize=9, color="#7a7a7a", style="italic")

cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                  fraction=0.022, pad=0.014)
cb.set_label("NES  (+ enriched among up-regulated,\n−  enriched among down-regulated)", fontsize=9)
cb.ax.tick_params(labelsize=8.5)
ax.set_title("a   Preranked GSEA on the same GO terms — direction comes from the sign of NES,\n"
             "     so no UP/DOWN split and no cutoff are needed",
             loc="left", fontsize=12.6, fontweight="bold", pad=62)

# ------------------------------------------------------------------ panel b
ax2 = fig.add_subplot(gs[1])
vals, labs, colors, texts = [], [], [], []
for (ds, var, vlab) in cols:
    r = idx.get(("GO:0034976", ds, var))
    nes = float(r["NES"]) if r and r["NES"] else 0.0
    sig = r and r["significant"] == "True"
    vals.append(nes)
    colors.append(("#c0392b" if nes > 0 else "#2c5f8a") if sig else "#cfcac3")
    labs.append(f"{ds}\n{vlab}")
    q = float(r["padj"]) if r and r["padj"] else float("nan")
    qs = ("%.1e" % q) if q < 1e-3 else ("%.3f" % q)
    texts.append(f"NES {nes:+.2f}\nBH p = {qs}" if sig else
                 f"NES {nes:+.2f}\nns (BH p = {q:.2f})")
x = np.arange(len(cols))
ax2.bar(x, vals, color=colors, width=.6, edgecolor="white")
for i, (v, t) in enumerate(zip(vals, texts)):
    ax2.text(i, v + (0.18 if v > 0 else -0.18), t, ha="center",
             va="bottom" if v > 0 else "top", fontsize=7.8, color="#333")
ax2.axhline(0, color="#333", lw=1)
ax2.set_ylim(min(vals) * 1.75, max(vals) * 1.55)
ax2.set_xticks(x); ax2.set_xticklabels(labs, fontsize=8.4)
ax2.set_ylabel("NES", fontsize=10)
ax2.set_xlim(-.6, len(cols) - .4)
for j in range(2, len(cols), 2):
    ax2.axvline(j - .5, color="#4a4a4a", lw=1.5)
for s in ["top", "right"]:
    ax2.spines[s].set_visible(False)
ax2.set_title("b   Response to endoplasmic reticulum stress (GO:0034976): significant and "
              "oppositely signed in the control and in both AD cohorts",
              loc="left", fontsize=11.8, fontweight="bold", pad=8)

fig.text(0.40, 0.016,
         "Weighted (p=1) preranked GSEA; NES and P from 10,000 size-matched gene-set permutations, "
         "BH-adjusted across the 18 terms within each condition. Rank metric = sign(log2FC) × −log10(P) "
         "(Nativio: Welch t statistic).\n"
         "Gene sets = QuickGO human protein annotations for each GO term plus descendants; n = set genes "
         "present in the Nativio ranked list. Bold labels = UPR/ERAD-core terms. "
         "GO:0030433 and GO:0042886 are obsolete in the current GO release.",
         fontsize=8, color="#666", va="bottom")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q1_UPR_GSEA_18terms.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q1_UPR_GSEA_18terms.{png,svg}")
