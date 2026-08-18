#!/usr/bin/env python3
"""R2.6 specificity figure: UPR sets vs expression-matched random sets and vs the
genome-wide background (cameraPR), in the two bulk AD cohorts."""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
rows = [r for r in csv.DictReader(open(os.path.join(
    OUT, "R2q6_specificity_expression_matched_and_camera.csv"))) if r["z_vs_matched_null"]]

GROUPS = ["UPR/ERAD core", "folding/trafficking", "generic parent", "negative control"]
GC = {"UPR/ERAD core": "#c0392b", "folding/trafficking": "#e08a3c",
      "generic parent": "#6b6b6b", "negative control": "#4a7fa5"}
DATASETS = ["Mizuno", "Nativio"]

fig, axes = plt.subplots(1, 2, figsize=(18.5, 8.6), sharex=False)
fig.subplots_adjust(left=0.235, right=0.975, top=0.845, bottom=0.145, wspace=1.02)

for ax, ds in zip(axes, DATASETS):
    sub = [r for r in rows if r["dataset"] == ds]
    sub.sort(key=lambda r: (GROUPS.index(r["group"]), float(r["z_vs_matched_null"])))
    y = np.arange(len(sub))[::-1]
    z = [float(r["z_vs_matched_null"]) for r in sub]
    cols = [GC[r["group"]] for r in sub]
    alph = [1.0 if float(r["padj_matched_null"]) < 0.05 else 0.42 for r in sub]
    for yy, zz, c, a in zip(y, z, cols, alph):
        ax.barh(yy, zz, color=c, alpha=a, height=.72,
                edgecolor="#333" if a == 1.0 else "none", lw=.8)
    for yy, r in zip(y, sub):
        q1, q2 = float(r["padj_matched_null"]), float(r["padj_camera"])
        zz = float(r["z_vs_matched_null"])
        txt = f"BH {q1:.2f} | camera {q2:.3f}"
        ax.text(zz + (0.22 if zz >= 0 else -0.22), yy, txt, va="center",
                ha="left" if zz >= 0 else "right", fontsize=6.6,
                color="#222" if q1 < 0.05 else "#9a9a9a")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['term'][:42]}  (n={r['K_in_data']})" for r in sub], fontsize=8.6)
    for t, r in zip(ax.get_yticklabels(), sub):
        t.set_color(GC[r["group"]])
        t.set_fontweight("bold" if r["group"] == "UPR/ERAD core" else "normal")
    ax.axvline(0, color="#333", lw=1)
    for v in (-1.96, 1.96):
        ax.axvline(v, color="#aaa", ls="--", lw=.9)
    lim = max(6.5, max(abs(np.array(z))) * 1.55)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("z of the set's mean statistic vs 10,000 expression-matched random sets\n"
                  "(negative = shifted down more than equally expressed genes)", fontsize=8.6)
    ax.tick_params(length=0)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    gm = float(sub[0]["global_mean_stat"])
    ax.set_title(f"{ds}   (global mean statistic = {gm:+.2f})", fontsize=12.5,
                 fontweight="bold", loc="left", pad=10)

handles = [plt.Rectangle((0, 0), 1, 1, color=GC[g]) for g in GROUPS]
axes[0].legend(handles, GROUPS, fontsize=8.4, loc="lower left", frameon=False,
               bbox_to_anchor=(0, -0.255), ncol=4, columnspacing=1.2, handlelength=1.2)

fig.suptitle("R2.6 — is the UPR change specific, or part of the global transcriptional shift?",
             fontsize=14.5, fontweight="bold", x=0.235, ha="left", y=0.965)
fig.text(0.235, 0.912,
         "Solid bars = significant after BH (expression-matched null). Dashed lines = ±1.96. "
         "Each bar is also annotated with the cameraPR competitive-test BH p (set vs all other genes, "
         "VIF-corrected).",
         fontsize=8.8, color="#666", ha="left")
fig.text(0.235, 0.022,
         "Rank metric: sign(log2FC) × −log10(P) (Nativio: Welch t). Random sets matched to each real set for "
         "size and expression-decile composition, 10,000 iterations.\n"
         "READ-OUT: in neither cohort do the UPR/ERAD-core sets move significantly beyond "
         "expression-matched background — the bulk UPR signature tracks the global proteostasis/trafficking "
         "shift rather than exceeding it.",
         fontsize=8.2, color="#555", ha="left", va="bottom")

for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2q6_specificity_test.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R2q6_specificity_test.{png,svg}")
