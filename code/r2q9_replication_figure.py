#!/usr/bin/env python3
"""R2.9 replication figure (two independent cohorts, both large & pathology-gradient-based).
  Mathys 2024 (Nature; ROSMAP, PFC, vs NFT tangles)   -- logFC_nb
  SEA-AD / Gabitto 2024 (Allen; MTG, vs Continuous Pseudo-progression Score; INDEPENDENT of ROSMAP)
Per broad cell type: DIRECTION of the UPR-associated set (ER-stress 260) vs the UPR-sensor mRNAs
(EIF2AK3/ERN1/ATF6, internal control). Direction index = frac(up) - frac(down). The neuronal decrease
reproduces in both cohorts (highly significant); the sensor mRNAs decline far less (internal reference)."""
import os, glob, gzip, csv, statistics as st
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
from scipy.stats import binomtest

OUT = os.path.dirname(os.path.abspath(__file__))
UPR = set(l.strip() for l in open(os.path.join(OUT, "ERstress_260_geneset.txt")) if l.strip())
SENS = {"EIF2AK3", "ERN1", "ATF6"}
SCR = ("/tmp/"
       "adbrainseq_work")

def didx(vals):
    vals = [v for v in vals if v is not None]
    if len(vals) < 2: return None
    up = sum(v > 0 for v in vals); dn = sum(v < 0 for v in vals); n = len(vals)
    return (up - dn) / n, binomtest(dn, n, 0.5).pvalue, n

# ---- Mathys 2024 (PFC, nft) ----
M = os.path.join(SCR, "mathys24/dereg")
MCL = {"Ast": "Ast", "Exc": "Ex", "Inh": "In", "Mic": "Mic", "Oli": "Oli", "Opc": "OPC"}
m24 = {v: {} for v in MCL.values()}
for fp in glob.glob(os.path.join(M, "aggregated_fullset.*.tsv.gz")):
    pre = os.path.basename(fp).split(".")[1].split("_")[0]
    if pre not in MCL: continue
    with gzip.open(fp, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("path") != "nft" or row.get("region") != "PFC": continue
            try: m24[MCL[pre]].setdefault(row["gene"], []).append(float(row["logFC_nb"]))
            except (TypeError, ValueError, KeyError): pass
m24 = {k: {g: st.mean(v) for g, v in d.items()} for k, d in m24.items()}

# ---- SEA-AD (MTG, CPS) from extracted UPR-gene rows ----
sea = {c: {} for c in ["Ex", "In", "Ast", "Mic", "Oli", "OPC"]}
for fp in glob.glob(os.path.join(SCR, "seaad", "*.csv")):
    c = os.path.basename(fp).split("__")[0]
    if c not in sea: continue
    for row in csv.reader(open(fp)):
        if len(row) < 8: continue
        try: sea[c].setdefault(row[0], []).append(float(row[7]))
        except (TypeError, ValueError): continue
sea = {k: {g: st.mean(v) for g, v in d.items()} for k, d in sea.items()}

COHORTS = [("Mathys 2024  ·  ROSMAP · PFC · 48 donors · Braak I–VI  (vs NFT tangles)", m24, ["Ex", "In", "Ast", "Mic", "Oli", "OPC"]),
           ("SEA-AD 2024  ·  Allen · MTG · 84 donors · Braak 0–VI  (independent of ROSMAP)", sea, ["Ex", "In", "Ast", "Mic", "Oli", "OPC"])]
def pstar(p): return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""

fig, axes = plt.subplots(1, 2, figsize=(13, 4.9))
for ax, (title, data, order) in zip(axes, COHORTS):
    xs = np.arange(len(order))
    vals, stars, cols = [], [], []
    for c in order:
        r = didx([data[c][g] for g in (UPR - SENS) if g in data[c]])
        v = r[0] if r else 0; vals.append(v); stars.append(pstar(r[1]) if r else "")
        cols.append("#3a6ea5" if c in ("Ex", "In") else "#c0563a" if v > 0 else "#9aa7b3")
    ax.bar(xs, vals, 0.62, color=cols, edgecolor="white")
    for xi, (v, sstar) in enumerate(zip(vals, stars)):
        if sstar: ax.text(xi, v + (0.02 if v >= 0 else -0.08), sstar, ha="center", fontsize=12)
    ax.axhline(0, color="#444", lw=0.9); ax.set_ylim(-0.66, 0.62)
    ax.axvspan(-0.5, 1.5, color="#eef3f8", zorder=0)                     # neuron block shading
    ax.set_xticks(xs); ax.set_xticklabels([f"{c}\n{'neuron' if c in ('Ex','In') else 'glia'}" for c in order], fontsize=9)
    ax.set_ylabel("UPR-set direction index\n(+up / –down)", fontsize=9.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.set_title(title, fontsize=9.6, fontweight="bold")
    ax.text(0.5, 0.57, "neurons", ha="center", fontsize=8.5, color="#3a6ea5", fontweight="bold")
fig.suptitle("Neuronal down-regulation of UPR-associated transcripts replicates in two independent AD snRNA cohorts",
             fontsize=12, fontweight="bold", y=1.0)
fig.text(0.5, 0.925, "Direction of the ER-stress/UPR set (~230 genes) vs pathology progression per broad cell type · "
         "sign test · SEA-AD is independent of ROSMAP (different consortium and brain region).",
         ha="center", fontsize=8.5, color="#555")
fig.subplots_adjust(top=0.85, bottom=0.13, left=0.07, right=0.99, wspace=0.22)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q9_replication_figure.{ext}"), dpi=300 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved R2Q9_replication_figure.{png,svg}")
for title, data, order in COHORTS:
    print(f"\n{title}")
    for c in order:
        u = didx([data[c][g] for g in (UPR - SENS) if g in data[c]])
        s = didx([data[c][g] for g in SENS if g in data[c]])
        print(f"  {c:4s} UPR dir={u[0]:+.2f} p={u[1]:.1e} (n={u[2]}) | sensor dir={s[0]:+.2f} (n={s[2]})" if u and s else f"  {c}: n/a")
