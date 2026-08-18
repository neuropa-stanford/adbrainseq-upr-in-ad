#!/usr/bin/env python3
"""R2.2 figure: donor-level module scores (donor = unit of inference), with individual donor points,
plus the effect-size summary including the internal control."""
import csv, os, sys, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})

DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
        "Human AD brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

seen, info = [], {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    s = r["Subject"]
    if s not in info: seen.append(s); info[s] = r
braak = {s: int(float(info[s]["braaksc"])) for s in seen}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")

SETS = [("ER-stress", "ERstress_260_geneset.txt"), ("PERK", "geneset_PERK.txt"),
        ("IRE1", "geneset_IRE1.txt"), ("ATF6", "geneset_ATF6.txt"), ("ERAD", "geneset_ERAD.txt"),
        ("control\n(mRNA transp.)", "controlset_MRNA_TRANSPORT.txt")]
gsets = {n: [l.strip() for l in open(os.path.join(OUT, f)) if l.strip()] for n, f in SETS}
CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
LAB = {"Ex": "Excitatory neurons", "In": "Inhibitory neurons", "Ast": "Astrocytes",
       "Mic": "Microglia", "Oli": "Oligodendrocytes", "Opc": "OPCs"}

def module_scores(ct):
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dcol = hdr.index("donor"); ccol = hdr.index("celltype")
    gcols = [j for j in range(2, len(hdr)) if j not in (dcol, ccol)]
    genes = [hdr[j] for j in gcols]; donors = [r[dcol] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gcols] for r in tab[1:]], float)
    mu = np.nanmean(X, axis=0); sd = np.nanstd(X, axis=0, ddof=1)
    ok = (sd > 0) & np.isfinite(sd)
    Z = np.full_like(X, np.nan); Z[:, ok] = (X[:, ok] - mu[ok]) / sd[ok]
    gi = {g: j for j, g in enumerate(genes)}
    out = {}
    for sn, gl in gsets.items():
        cols = [gi[g] for g in gl if g in gi and ok[gi[g]]]
        if len(cols) >= 10: out[sn] = (donors, np.nanmean(Z[:, cols], axis=1))
    return out

def cohens_d(a, b):
    na, nb = len(a), len(b)
    sp = math.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1)) / (na+nb-2))
    return (np.mean(a)-np.mean(b))/sp if sp > 0 else np.nan

allsc = {ct: module_scores(ct) for ct in CELLS}

fig = plt.figure(figsize=(17, 9.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.42, left=0.06, right=0.985,
                      top=0.88, bottom=0.08)

# ---------- Row 1: donor-level ER-stress module score per cell type ----------
gs1 = gs[0].subgridspec(1, 6, wspace=0.32)
GC = {"low": "#7f8c8d", "int": "#e08a3c", "late": "#c0392b"}
for k, ct in enumerate(CELLS):
    ax = fig.add_subplot(gs1[k])
    donors, sc = allsc[ct]["ER-stress"]
    for gi_, g in enumerate(["low", "int", "late"]):
        v = np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == g])
        v = v[np.isfinite(v)]
        ax.scatter(np.full(len(v), gi_) + np.random.uniform(-.13, .13, len(v)), v,
                   s=26, color=GC[g], alpha=.75, edgecolor="#333", lw=.4, zorder=3)
        ax.plot([gi_-.26, gi_+.26], [v.mean()]*2, color="#111", lw=2.2, zorder=4)
        se = v.std(ddof=1)/math.sqrt(len(v))
        ax.plot([gi_, gi_], [v.mean()-1.96*se, v.mean()+1.96*se], color="#111", lw=1.2, zorder=4)
        ax.text(gi_, ax.get_ylim()[0], f"n={len(v)}", ha="center", va="bottom", fontsize=7.5, color="#555")
    ax.axhline(0, color="#bbb", lw=.8, ls="--")
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["I–II", "III–IV", "V–VI"], fontsize=9)
    ax.set_title(LAB[ct], fontsize=10, fontweight="bold")
    if k == 0: ax.set_ylabel("ER-stress module score\n(per donor)", fontsize=9.5)
    ax.set_xlabel("Braak stage", fontsize=9)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)

# ---------- Row 2: Cohen's d (late vs low), all sets incl. control ----------
ax2 = fig.add_subplot(gs[1])
setnames = [s[0] for s in SETS]
W = 0.13
for si, sn in enumerate(setnames):
    ds = []
    for ct in CELLS:
        donors, sc = allsc[ct][sn]
        a = np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == "late"])
        b = np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == "low"])
        ds.append(cohens_d(a[np.isfinite(a)], b[np.isfinite(b)]))
    x = np.arange(len(CELLS)) + (si - 2.5) * W
    isctrl = "control" in sn
    ax2.bar(x, ds, width=W, label=sn.replace("\n", " "),
            color="#333333" if isctrl else plt.cm.tab10(si), alpha=.95 if isctrl else .85,
            edgecolor="#111" if isctrl else "none", lw=1.4 if isctrl else 0, zorder=3)
ax2.axhline(0, color="#333", lw=1)
ax2.set_xticks(np.arange(len(CELLS))); ax2.set_xticklabels([LAB[c] for c in CELLS], fontsize=10)
ax2.set_ylabel("Cohen's d  (Braak V–VI vs I–II)", fontsize=10)
ax2.legend(fontsize=8.5, ncol=6, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.16))
for s in ["top", "right"]: ax2.spines[s].set_visible(False)
ax2.text(0.5, -0.62, "neurons: UPR sets and the internal control move together (global shift)   |   "
         "oligodendrocytes: IRE1 exceeds the control (d = 1.45 vs 0.76)",
         transform=ax2.transAxes, ha="center", fontsize=9, color="#555")

fig.suptitle("Donor-level analysis (donor = unit of inference, n = 10/21/17 donors): "
             "UPR module scores by Braak stage, with an internal control gene set",
             fontsize=13.5, fontweight="bold", y=0.965)
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q2_donorlevel.{ext}"), dpi=250 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
print("saved R2Q2_donorlevel.{png,svg}")
