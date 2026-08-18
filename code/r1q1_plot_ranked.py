#!/usr/bin/env python3
"""R1.1 ranked version — GO terms ranked by P value within every condition,
with the number of genes in each term annotated on the bar (reference-figure style).

Layout: 6 rows (dataset x DEG threshold) x 3 columns (UP | DOWN | BOTH).
Within each panel the 18 GO terms are sorted by adjusted P (most significant on top).
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "axes.unicode_minus": False,
})
from matplotlib.colors import LinearSegmentedColormap, Normalize

OUT = "/data/adbrainseq/anc/Major Revision/processed raw data"
CSV = os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv")
rows = list(csv.DictReader(open(CSV)))

SHORT = {
    "GO:0034976": "Response to ER stress",
    "GO:0030968": "ER unfolded protein response",
    "GO:0006986": "Response to unfolded protein",
    "GO:0034620": "Cellular resp. to unfolded protein",
    "GO:0035966": "Resp. to topol. incorrect protein",
    "GO:0035967": "Cell. resp. to topol. incorrect prot.",
    "GO:0036498": "IRE1-mediated UPR",
    "GO:1905897": "Regulation of resp. to ER stress",
    "GO:1903573": "Neg. regulation of resp. to ER stress",
    "GO:0036503": "ERAD pathway",
    "GO:0030433": "Ubiquitin-dependent ERAD pathway",
    "GO:0006457": "Protein folding",
    "GO:0034975": "Protein folding in ER",
    "GO:0006888": "ER to Golgi vesicle-mediated transport",
    "GO:0042886": "Amide transport",
    "GO:0015833": "Peptide transport",
    "GO:0015031": "Protein transport",
    "GO:0008104": "Protein localization",
}
CORE = {"GO:0034976", "GO:0030968", "GO:0006986", "GO:0034620", "GO:0035966",
        "GO:0035967", "GO:0036498", "GO:1905897", "GO:1903573", "GO:0036503",
        "GO:0030433"}
# Restrict to the SAME 9 UPR/ERAD/folding terms shown in the consolidated Figure 1a, so the two
# figures use an identical GO list (drop generic parent terms, GO-obsolete terms, and the
# lowest-ranked regulation/IRE1 terms that fall outside Figure 1a's top list).
DROP = {"GO:0008104", "GO:0015031", "GO:0015833",   # generic parent terms
        "GO:0030433", "GO:0042886",                 # obsolete in current GO release
        "GO:1905897", "GO:0036498", "GO:1903573",   # ranks 11-13, not in Figure 1a
        "GO:0034975"}                               # Protein folding in ER — not enriched in any DOWN cohort (dropped from Fig 1a)
DATASETS = ["Thapsigargin", "Nativio", "Mizuno"]
THRESH = [("no cutoff (sign only)", "no cutoff"), ("p<0.05", "p<0.05")]
# 3 columns per row: [this dataset · UP] | [this dataset · DOWN] | [Combined (M∩N) · DOWN].
# The 3rd column is NOT the union of that dataset's own up+down genes. It is the cross-cohort
# INTERSECTION: genes present in BOTH Mizuno and Nativio AND moving in the same (down) direction —
# the reproduced AD-brain signal. Since it is a Mizuno∩Nativio quantity it is shown beside each
# AD cohort (identical panel at a given threshold) and omitted for the Thapsigargin control, which
# is an in-vitro positive control with no cross-cohort intersection.
COMBINED = "Combined(M∩N)"
COL_COLOR = {"UP": "#c0392b", "DOWN": "#2c5f8a", "COMBINED": "#6a3d9a"}

idx = {}
for r in rows:
    idx.setdefault((r["dataset"], r["threshold"], r["direction"]), []).append(r)

cmap = LinearSegmentedColormap.from_list(
    "upr", ["#fbe7c6", "#f7b05b", "#e0562b", "#a01f28", "#5c0f1c"])
norm = Normalize(0, 60)
NS_COL = "#d9d5cf"

fig, axes = plt.subplots(6, 3, figsize=(22, 25.5))
fig.subplots_adjust(left=0.125, right=0.99, top=0.945, bottom=0.035,
                    hspace=0.42, wspace=0.62)

ranked_out = []
seen_out = set()
for ri, ds in enumerate(DATASETS):
    for ti, (thr_key, thr_lab) in enumerate(THRESH):
        row = ri * 2 + ti
        # column spec: (lookup_dataset, direction, is_combined_column)
        colspec = [(ds, "UP", False), (ds, "DOWN", False), (COMBINED, "DOWN", True)]
        for ci, (look_ds, dr, is_comb) in enumerate(colspec):
            ax = axes[row][ci]
            # Combined (M∩N) column: omit for the Thapsigargin positive control — it is an
            # in-vitro dataset, not one of the two AD cohorts, so no cross-cohort intersection exists.
            if is_comb and ds == "Thapsigargin":
                ax.axis("off")
                ax.text(.5, .5, "Combined (M∩N)\nnot applicable\n(in-vitro positive control —\n"
                        "no cross-cohort intersection)",
                        ha="center", va="center", fontsize=8.6, color="#b0aca4",
                        style="italic", transform=ax.transAxes)
                continue
            recs = []
            for r in idx.get((look_ds, thr_key, dr), []):
                if r["GO"] in DROP:
                    continue
                if not r["adjP"]:
                    continue
                p = float(r["adjP"])
                recs.append((r["GO"], -math.log10(p) if p > 0 else 330, p,
                             int(r["genes_in_term"] or 0),
                             int(r["term_size"] or 0), r["significant"] == "True"))
            recs.sort(key=lambda x: -x[1])
            n_q = idx.get((look_ds, thr_key, dr), [{}])[0].get("n_query_genes", "")

            for rank, (go, nl, p, ng, ts, sig) in enumerate(recs, 1):
                # de-dup the Combined column, which is shown beside both AD cohorts
                key = (look_ds, thr_lab, dr, go)
                if key in seen_out:
                    continue
                seen_out.add(key)
                ranked_out.append({"dataset": look_ds, "threshold": thr_lab,
                                   "direction": dr, "rank": rank, "GO": go,
                                   "term": SHORT[go], "adjP": p,
                                   "genes_in_term": ng, "term_size": ts,
                                   "significant": sig})

            if not recs:
                ax.axis("off")
                ax.text(.5, .5, "no genes in query", ha="center", va="center",
                        fontsize=9, color="#999", transform=ax.transAxes)
                continue

            y = np.arange(len(recs))[::-1]
            vals = [min(r[1], 90) for r in recs]
            cols_ = [cmap(norm(min(r[1], 90))) if r[5] else NS_COL for r in recs]
            ax.barh(y, vals, color=cols_, height=.78,
                    edgecolor=["#7a2418" if r[5] else "#bfbab3" for r in recs], lw=.5)
            xmax = max(vals) * 1.30 + 1
            for yy, r in zip(y, recs):
                lab = f"{r[3]}/{r[4]} genes" if r[5] else f"ns  ({r[3]}/{r[4]})"
                ax.text(min(r[1], 90) + xmax * .012, yy, lab, va="center",
                        fontsize=7.2, color="#333" if r[5] else "#a8a29a")
            ax.set_yticks(y)
            ax.set_yticklabels(
                [f"{SHORT[r[0]]}  ({r[0]})" for r in recs], fontsize=7.3)
            for tick, r in zip(ax.get_yticklabels(), recs):
                tick.set_color("#111" if r[0] in CORE else "#8a8a8a")
                if r[0] in CORE:
                    tick.set_fontweight("bold")
            ax.axvline(-math.log10(0.05), color="#aaa", ls="--", lw=.9, zorder=0)
            ax.set_xlim(0, xmax)
            ax.set_xlabel("$-\\log_{10}$ adjusted $P$", fontsize=8.5)
            ax.tick_params(axis="x", labelsize=8)
            ax.tick_params(length=0)
            for s in ["top", "right", "left"]:
                ax.spines[s].set_visible(False)
            if is_comb:
                title = f"Combined (M∩N) · {thr_lab} · DOWN   (n = {int(n_q):,} genes)"
                tcol = COL_COLOR["COMBINED"]
            else:
                title = f"{ds} · {thr_lab} · {dr}   (n = {int(n_q):,} genes)"
                tcol = COL_COLOR[dr]
            ax.set_title(title, fontsize=10.2, fontweight="bold", loc="left", pad=6, color=tcol)

fig.suptitle("UPR-related GO terms ranked by adjusted P value within each condition — "
             "columns: [dataset · UP] | [dataset · DOWN] | [Combined (M∩N) · DOWN]",
             fontsize=15, fontweight="bold", x=0.125, ha="left", y=0.982)
fig.text(0.125, 0.966,
         "Bar label = query genes in term / total term size. Bold = UPR/ERAD-core terms; grey = folding/"
         "trafficking terms. Grey bars = not significant (ns). Dashed line = P 0.05. The same 9 UPR-related "
         "terms as consolidated Figure 1a are shown. "
         "Combined (M∩N) = the cross-cohort INTERSECTION (genes in both Mizuno and Nativio, concordantly "
         "DOWN) — NOT a union; identical beside each AD cohort, not applicable to the Thapsigargin control.",
         fontsize=9, color="#666", ha="left")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q1_UPR_GO_ranked_10terms.{ext}"),
                dpi=250 if ext == "png" else None, facecolor="white",
                bbox_inches="tight")

fn = os.path.join(OUT, "R1Q1_UPR_GO_ranked_by_pvalue.csv")
with open(fn, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ranked_out[0].keys()))
    w.writeheader(); w.writerows(ranked_out)
print("saved figure + ", fn, len(ranked_out), "rows")
