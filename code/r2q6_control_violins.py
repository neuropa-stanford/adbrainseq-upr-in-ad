#!/usr/bin/env python3
"""R2.6 violins — for two neutral control gene sets (Heike control; mRNA-transport internal control),
the per-gene log2 fold-change distribution in bulk (Mizuno, Nativio) and in Mathys snRNA per cell type
(late Braak V-VI vs low I-II). Shows: controls are flat in bulk but track the cell-type shift in snRNA."""
import csv, os, sys, glob, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})
sys.path.insert(0, ".")
from r1q1_gomatrix import load_mizuno, load_nativio

OUT = os.path.dirname(os.path.abspath(__file__))
DATA = glob.glob("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
                 "**/data_extraction", recursive=True)[0]

seen, info = [], {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    s = r["Subject"]
    if s not in info: seen.append(s); info[s] = r
braak = {s: int(float(info[s]["braaksc"])) for s in seen}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")

miz = {g: v[0] for g, v in load_mizuno().items()}

# Nativio: CPM-normalized per-gene log2FC (library-size corrected; matches our internal-control finding
# and the R2.4 normalization decision), from SuppD3 per-sample raw counts (12 AD, 10 aged controls).
import openpyxl
BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD3_NativioSeq.xlsx"), read_only=True); ws = wb.active
ng, cAD, cOLD = [], [], []
for r in ws.iter_rows(min_row=3, values_only=True):
    sym = r[1]
    if not sym or str(sym).strip() in ("", "None"): continue
    try:
        ad = [float(x) for x in r[6:18]]; ol = [float(x) for x in r[18:28]]
    except (TypeError, ValueError): continue
    if len(ad) != 12 or len(ol) != 10: continue
    ng.append(str(sym).strip()); cAD.append(ad); cOLD.append(ol)
wb.close()
cnt = np.hstack([np.array(cAD, float), np.array(cOLD, float)])
cpm = cnt / cnt.sum(0) * 1e6
mAD = cpm[:, :12].mean(1); mOLD = cpm[:, 12:].mean(1)
nat = {}
for j, g in enumerate(ng):
    lfc = math.log2((mAD[j] + 1) / (mOLD[j] + 1))
    if g not in nat: nat[g] = lfc

CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
CLAB = {"Ex": "Excit.\nneuron", "In": "Inhib.\nneuron", "Ast": "Astro", "Mic": "Micro", "Oli": "Oligo", "Opc": "OPC"}
# snRNA per-gene log2FC late/low per cell type
sn = {}
for ct in CELLS:
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dc = hdr.index("donor"); cc = hdr.index("celltype")
    gc = [j for j in range(2, len(hdr)) if j not in (dc, cc)]
    genes = [hdr[j] for j in gc]; don = [r[dc] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gc] for r in tab[1:]], float)
    late = np.array([grp(braak[d]) == "late" for d in don]); low = np.array([grp(braak[d]) == "low" for d in don])
    ml = np.nanmean(X[late], 0); mo = np.nanmean(X[low], 0)
    eps = 1e-3
    lfc = np.log2((ml + eps) / (mo + eps))
    sn[ct] = {g: lfc[j] for j, g in enumerate(genes)}

def series(geneset):
    """return dict group -> array of per-gene log2FC"""
    out = {}
    out["Mizuno"] = np.array([miz[g] for g in geneset if g in miz], float)
    out["Nativio"] = np.array([nat[g] for g in geneset if g in nat], float)
    for ct in CELLS:
        out[ct] = np.array([sn[ct][g] for g in geneset if g in sn[ct] and np.isfinite(sn[ct][g])], float)
    return out

SETS = [("Heike control gene set", "Heike_control_geneset_45.txt"),
        ("mRNA-transport internal control", "controlset_MRNA_TRANSPORT.txt")]

fig, axes = plt.subplots(2, 1, figsize=(11, 9))
groups = ["Mizuno", "Nativio"] + CELLS
glabels = ["Mizuno\n(bulk)", "Nativio\n(bulk)"] + [CLAB[c] for c in CELLS]
# colors: bulk grey, neurons blue, glia/other red-ish
gcol = {"Mizuno": "#8a8f98", "Nativio": "#8a8f98", "Ex": "#2c5f8a", "In": "#3a78ad",
        "Ast": "#c0392b", "Mic": "#c0392b", "Oli": "#a01f28", "Opc": "#d98a80"}

summary = []
for ax, (title, fn) in zip(axes, SETS):
    gs = [l.strip() for l in open(os.path.join(OUT, fn)) if l.strip()]
    ser = series(gs)
    data = [ser[g] for g in groups]
    parts = ax.violinplot(data, showmedians=True, widths=0.85)
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(gcol[groups[i]]); pc.set_alpha(0.55); pc.set_edgecolor("#444"); pc.set_linewidth(0.6)
    for key in ["cmedians", "cmaxes", "cmins", "cbars"]:
        parts[key].set_color("#333"); parts[key].set_linewidth(1.0)
    # overlay medians as text
    for i, d in enumerate(data):
        med = np.median(d) if len(d) else np.nan
        ax.text(i + 1, ax.get_ylim()[1] if False else 0, "", ha="center")
        summary.append((title, groups[i], len(d), round(float(med), 3)))
        ax.text(i + 1, np.median(d), f"{np.median(d):+.2f}", ha="center", va="bottom",
                fontsize=7.5, color="#111", fontweight="bold")
    ax.axhline(0, color="#999", ls="--", lw=1)
    ax.axvline(2.5, color="#333", lw=1.4)  # bulk | snRNA divider
    ax.set_xticks(range(1, len(groups) + 1)); ax.set_xticklabels(glabels, fontsize=9)
    ax.set_ylabel("per-gene log$_2$FC\n(AD vs ctrl / late vs low Braak)", fontsize=9.5)
    ax.set_title(f"{title}  (n = {len(gs)} genes)   —   bulk: flat   |   snRNA: tracks the cell-type shift",
                 fontsize=11, fontweight="bold", loc="left")
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    ax.text(1.5, ax.get_ylim()[0], "BULK", ha="center", va="bottom", fontsize=8, color="#666", style="italic")
    ax.text(5.5, ax.get_ylim()[0], "snRNA (Mathys, per cell type)", ha="center", va="bottom",
            fontsize=8, color="#666", style="italic")

fig.suptitle("Neutral control gene sets: flat in bulk RNA-seq, but track the neuron-down / glia-up "
             "cell-type shift in snRNA", fontsize=12.5, fontweight="bold", y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.97])
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R2Q6_control_violins.{ext}"), dpi=250 if ext == "png" else None,
                facecolor="white", bbox_inches="tight")
# save summary
with open(os.path.join(OUT, "R2Q6_control_violins_medians.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["gene_set", "group", "n_genes", "median_log2FC"]); w.writerows(summary)
print("saved R2Q6_control_violins.{png,svg} + medians.csv")
for s in summary: print(s)
