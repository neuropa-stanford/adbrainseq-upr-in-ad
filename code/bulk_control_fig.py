#!/usr/bin/env python3
"""Figure: internal-control GO gene sets are unchanged (log2FC~0) in both bulk cohorts,
in contrast to the broader genome; validates them as stable internal controls."""
import sys, os, statistics
import numpy as np
sys.path.insert(0, '.')
from r1q1_gomatrix import load_mizuno
from nativio_cpm import load_nativio_cpm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
miz = {g: v[0] for g, v in load_mizuno().items()}
nat = load_nativio_cpm()

CTRL = [("mRNA transport", "controlset_MRNA_TRANSPORT.txt"),
        ("RNA export from nucleus", "controlset_RNA_EXPORT_FROM_NUCLEUS.txt"),
        ("Transcription elongation", "controlset_DNA_TEMPLATED_TRANSCRIPTION_ELONGATION.txt"),
        ("Ribosome assembly", "controlset_RIBOSOME_ASSEMBLY.txt"),
        ("RNA destabilization", "controlset_RNA_DESTABILIZATION.txt")]
UPR = [("ER-stress", "ERstress_260_geneset.txt"), ("ERAD", "geneset_ERAD.txt"),
       ("PERK", "geneset_PERK.txt")]

def vals(fn, d):
    g = [l.strip() for l in open(os.path.join(OUT, fn)) if l.strip()]
    return [d[x] for x in g if x in d]

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2), sharey=True)
fig.subplots_adjust(left=0.20, right=0.98, top=0.88, bottom=0.10, wspace=0.06)
labels = [c[0] for c in CTRL] + ["—"] + [u[0] for u in UPR]
for ax, (dn, d) in zip(axes, [("Mizuno", miz), ("Nativio (CPM)", nat)]):
    data, cols = [], []
    for name, fn in CTRL:
        data.append(vals(fn, d)); cols.append("#2c7fb8")
    data.append([]); cols.append("white")           # spacer
    for name, fn in UPR:
        data.append(vals(fn, d)); cols.append("#c0392b")
    y = np.arange(len(data))[::-1]
    for yy, v, c in zip(y, data, cols):
        if not v: continue
        med = statistics.median(v)
        ax.scatter(v, np.full(len(v), yy) + np.random.uniform(-0.13, 0.13, len(v)),
                   s=7, color=c, alpha=0.35, edgecolor="none")
        ax.plot([np.percentile(v, 25), np.percentile(v, 75)], [yy, yy], color=c, lw=6, alpha=0.5,
                solid_capstyle="butt")
        ax.plot([med, med], [yy - 0.3, yy + 0.3], color="#111", lw=2)
        ax.text(1.55, yy, f"{med:+.2f}", va="center", fontsize=8,
                color="#12507e" if c == "#2c7fb8" else "#8c1c13", fontweight="bold")
    ax.axvline(0, color="#333", lw=1.2)
    ax.set_xlim(-1.6, 1.9)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
    for tk, c in zip(ax.get_yticklabels(), cols):
        tk.set_color("#12507e" if c == "#2c7fb8" else ("#8c1c13" if c == "#c0392b" else "white"))
        tk.set_fontweight("bold")
    ax.set_xlabel("log$_2$FC (AD vs control)", fontsize=10)
    ax.set_title(dn, fontsize=12.5, fontweight="bold")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
axes[0].text(-1.55, len(labels) - 0.4, "internal control\n(GO, 40–100 genes)", fontsize=9,
             color="#12507e", fontweight="bold", va="top")
axes[0].text(-1.55, 2.4, "UPR gene sets", fontsize=9, color="#8c1c13", fontweight="bold", va="top")
fig.suptitle("Internal-control GO gene sets are unchanged (log$_2$FC ≈ 0) in both bulk AD cohorts",
             fontsize=13.5, fontweight="bold", y=0.965)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"bulk_internal_control.{ext}"),
                dpi=250 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved bulk_internal_control.{png,svg}")
