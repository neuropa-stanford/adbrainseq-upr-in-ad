#!/usr/bin/env python3
"""Figure 4 v2.0 — COHORT-STACKED single figure (Mathys on top, SEA-AD directly below, per cell type).
A = the two GO panels merged into ONE grouped graph (hue=direction, shade=cohort). Then violins: B,C
(Mathys neurons | glia) over D,E (SEA-AD neurons | glia). Then donor module scores: F,G (Mathys) over
H,I (SEA-AD). Same data/fonts/dots/stars/colours as v1.0 (build_figure4_combined.py). Fully vector.
NOTE: distinct from the older crop-based build_figure4_v2.py; this is the matplotlib cohort-stacked v2.0."""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from matplotlib.ticker import FuncFormatter
from scipy.stats import wilcoxon, spearmanr, gaussian_kde

DRD = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(DRD)
SENS = {"EIF2AK3", "ERN1", "ATF6"}
NEUR = [("Ex", "Excitatory\nneurons"), ("In", "Inhibitory\nneurons")]
GLIA = [("Ast", "Astrocytes"), ("Mic", "Microglia"), ("Oli", "Oligodendro-\ncytes"), ("OPC", "OPCs")]
CELLS = NEUR + GLIA
C_MID, C_HIGH = "#3c8f34", "#c8d94a"; RED, BLUE = "#e2231a", "#2b3a8f"
CTCOL = {"Ex": "#8cbf43", "In": "#33a45c", "Ast": "#e07b54", "Mic": "#2ca0a8", "Oli": "#4f86c6", "OPC": "#9c77b4"}
FS_LAB, FS_TICK, FS_YT, FS_STAR, FS_LET, FS_LEG, FS_RHO = 8.0, 7.0, 8.0, 9.5, 13.0, 7.0, 6.8
def pstar(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
def lead0(v, _): return ("%g" % v) if v == int(v) else ("%.1f" % v)

def load(fname):
    rows = list(csv.DictReader(open(os.path.join(DRD, fname))))
    meta = {"cell_type", "donor", "braak_group", "braak_num"}
    genes = [g for g in rows[0].keys() if g not in meta and g not in SENS]
    byct = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for r in rows:
        c = r["cell_type"]
        if c in byct and r["braak_group"] in byct[c]:
            byct[c][r["braak_group"]].append(np.array([float(r[g]) for g in genes]))
    def gene_fc(c, grp):
        lo = np.vstack(byct[c]["low"]).mean(0); g = np.vstack(byct[c][grp]).mean(0); return g - lo
    score = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for c, _ in CELLS:
        order = [(v, "low") for v in byct[c]["low"]] + [(v, "mid") for v in byct[c]["mid"]] + [(v, "high") for v in byct[c]["high"]]
        M = np.vstack([v for v, _ in order]); mu = M.mean(0); sd = M.std(0); sd[sd == 0] = 1; Z = (M - mu) / sd
        for k, (_, grp) in enumerate(order): score[c][grp].append(float(Z[k].mean()))
    return dict(gene_fc=gene_fc, score=score, ngene=len(genes))

def load_mathys():
    # per-donor mean expression (mean_DGE_by_donor_cell_type_*) + donor_info Braak; reproduces the
    # original Fig-4 B/C violins + D/E donor module scores and fixes the panel-G flat artefact.
    DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/Human AD "
            "brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
    ES = [l.strip() for l in open(os.path.join(DRD, "ERstress_260_geneset.txt")) if l.strip() and l.strip() not in SENS]
    di = list(csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))))
    grpname = ["low" if int(float(r["braaksc"])) <= 2 else ("mid" if int(float(r["braaksc"])) <= 4 else "high") for r in di]
    FILE = {"Ex": "Ex", "In": "In", "Ast": "Ast", "Mic": "Mic", "Oli": "Oli", "OPC": "Opc"}
    raw, present = {}, {}
    for c, _ in CELLS:
        rows = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{FILE[c]}.csv"))))
        gidx = {g: j for j, g in enumerate(rows[0])}; present[c] = set(g for g in ES if g in gidx); raw[c] = (gidx, rows[1:])
    genes = [g for g in ES if all(g in present[c] for c, _ in CELLS)]
    byct = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for c, _ in CELLS:
        gidx, data = raw[c]; cols = [gidx[g] for g in genes]
        for i, r in enumerate(data): byct[c][grpname[i]].append(np.array([float(r[j]) for j in cols]))
    def gene_fc(c, grp):
        lo = np.vstack(byct[c]["low"]).mean(0); g = np.vstack(byct[c][grp]).mean(0); return g - lo
    score = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for c, _ in CELLS:
        order = [(v, "low") for v in byct[c]["low"]] + [(v, "mid") for v in byct[c]["mid"]] + [(v, "high") for v in byct[c]["high"]]
        M = np.vstack([v for v, _ in order]); mu = M.mean(0); sd = M.std(0, ddof=1); sd[sd == 0] = 1; Z = (M - mu) / sd
        for k, (_, grp) in enumerate(order): score[c][grp].append(float(Z[k].mean()))
    return dict(gene_fc=gene_fc, score=score, ngene=len(genes))

MA = load_mathys()
SA = load("SEAAD_3group_donor_pseudobulk.csv")
GX = {"low": 0, "mid": 1, "high": 2}
GO_VALS = {
    "MATHYS": {"Ex": (12.54, True), "In": (12.62, True), "Ast": (5.75, False), "Mic": (11.86, False), "Oli": (8.99, False), "OPC": (4.83, False)},
    "SEAAD":  {"Ex": (14.87, True), "In": (21.47, True), "Ast": (0.79, True),  "Mic": (2.15, False), "Oli": (2.39, False), "OPC": (1.59, False)},
}
RED_D, RED_L = "#d1201f", "#f2a19b"
BLU_D, BLU_L = "#23348c", "#7387cf"
def go_merged(ax):
    yc = np.arange(len(CELLS))[::-1]; h = 0.36
    for i, (c, lab) in enumerate(CELLS):
        mv, md = GO_VALS["MATHYS"][c]; sv, sd = GO_VALS["SEAAD"][c]
        ax.barh(yc[i] + h*0.58, mv, color=(BLU_D if md else RED_D), height=h)
        ax.barh(yc[i] - h*0.58, sv, color=(BLU_L if sd else RED_L), height=h)
    ax.axvline(1.3, color="#2b3a8f", lw=1.0)
    ax.set_xlim(0, max(max(GO_VALS["MATHYS"][c][0], GO_VALS["SEAAD"][c][0]) for c, _ in CELLS) * 1.10)
    ax.set_yticks(yc); ax.set_yticklabels([l for _, l in CELLS], fontsize=FS_LAB, linespacing=0.85)
    ax.set_ylim(-0.6, len(CELLS) - 0.4)
    ax.set_xlabel(r"$-$log$_{10}$($p$-value)  ·  GO Response to ER stress", fontsize=FS_LAB); ax.tick_params(labelsize=FS_TICK)
    H = [Patch(fc=RED_D, label="Mathys Up"), Patch(fc=BLU_D, label="Mathys Down"),
         Patch(fc=RED_L, label="SEA-AD Up"), Patch(fc=BLU_L, label="SEA-AD Down")]
    order = [0, 2, 1, 3]
    ax.legend(handles=[H[i] for i in order], fontsize=FS_LEG, frameon=False, loc="lower right",
              ncol=2, handlelength=1.1, columnspacing=1.1, labelspacing=0.35, borderpad=0.3)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

def violin_panel(ax, D, cells, ylab, ylim, legend, show_y=True):
    bw = 0.40; rng = np.random.RandomState(4); star_pos = []
    ax.axhline(0, color="#000", lw=0.7, ls=":")
    for i, (c, lab) in enumerate(cells):
        for grp, col in [("mid", C_MID), ("high", C_HIGH)]:
            xp = i + (0.5 if grp == "high" else -0.5) * bw
            v = np.clip(D["gene_fc"](c, grp), *ylim)
            for pc in ax.violinplot([v], positions=[xp], widths=bw*0.92, showextrema=False)["bodies"]:
                pc.set_facecolor(col); pc.set_alpha(0.92); pc.set_edgecolor("#1f1f1f"); pc.set_linewidth(0.6)
            if len(np.unique(v)) > 1:
                dens = gaussian_kde(v)(v); hw = (bw*0.46) * dens / dens.max()
            else: hw = np.zeros(len(v))
            jit = (rng.rand(len(v))*2-1) * hw * 0.85
            ax.scatter(xp + jit, v, s=1.4, color="#222", alpha=0.40, edgecolor="none", zorder=3)
            q1, q3 = np.percentile(v, [25, 75]); ax.hlines([q1, q3], xp-bw*.5, xp+bw*.5, color=BLUE, lw=1.0, zorder=5)
            ax.hlines(v.mean(), xp-bw*.52, xp+bw*.52, color=RED, lw=2.0, zorder=6)
            st = pstar(wilcoxon(D["gene_fc"](c, grp)).pvalue)
            if st != "ns": star_pos.append((xp, st))
    ax.set_xlim(-0.6, len(cells)-0.4); ax.set_ylim(*ylim); ax.set_xticks([])
    trs = blended_transform_factory(ax.transData, ax.transAxes)
    for xp, st in star_pos: ax.text(xp, 1.015, st, transform=trs, ha="center", va="bottom", fontsize=FS_STAR)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for i, (c, lab) in enumerate(cells): ax.text(i, -0.03, lab, ha="center", va="top", fontsize=FS_LAB, transform=tr, linespacing=0.95)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    if show_y:
        ax.set_ylabel(ylab, fontsize=FS_YT); ax.tick_params(axis="y", labelsize=FS_TICK)
    else:
        ax.tick_params(axis="y", labelleft=False, length=0)
    if legend:
        ax.legend(handles=[Patch(fc=C_MID, ec="#1f1f1f", label="Braak III,IV / low"),
                           Patch(fc=C_HIGH, ec="#1f1f1f", label="Braak V,VI / low")],
                  fontsize=FS_LEG, frameon=True, loc="lower left", bbox_to_anchor=(0.0, 1.11), ncol=2, columnspacing=1.0, handlelength=1.0)

def donor_panel(fig, sub, D, cells, ylab, show_y=True):
    axes = [fig.add_subplot(sub[j]) for j in range(len(cells))]
    for idx, (c, lab) in enumerate(cells):
        ax = axes[idx]; ax.axhline(0, color="#bbb", lw=0.7, ls=":")
        for grp in ("low", "mid", "high"):
            v = np.array(D["score"][c][grp]); x = GX[grp]
            jit = (np.random.RandomState(idx*3+x).rand(len(v))-0.5)*0.5
            ax.scatter(np.full(len(v), x)+jit, v, s=10, color=CTCOL[c], alpha=0.75, edgecolor="#333", lw=0.3, zorder=3)
            m = v.mean(); se = v.std(ddof=1)/math.sqrt(len(v)) if len(v) > 1 else 0
            ax.plot([x-0.28, x+0.28], [m, m], color="#111", lw=1.8); ax.plot([x, x], [m-1.96*se, m+1.96*se], color="#111", lw=1.0)
        allv = np.array(D["score"][c]["low"]+D["score"][c]["mid"]+D["score"][c]["high"])
        bkn = np.array([1]*len(D["score"][c]["low"])+[3.5]*len(D["score"][c]["mid"])+[5.5]*len(D["score"][c]["high"]))
        rho, spv = spearmanr(allv, bkn)
        ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["I–II", "III–IV", "V–VI"], fontsize=FS_TICK, rotation=45, ha="right")
        ax.set_title(lab, fontsize=FS_LAB, linespacing=0.95, pad=10)
        if spv < 0.05:
            ax.text(0.5, 1.02, f"ρ={rho:+.2f} {pstar(spv)}", transform=ax.transAxes, ha="center", va="bottom", fontsize=FS_RHO, color="#b00")
        if idx == 0 and show_y: ax.set_ylabel(ylab, fontsize=FS_YT)
        for s2 in ("top", "right"): ax.spines[s2].set_visible(False)
        ax.tick_params(axis="y", labelsize=FS_TICK); ax.yaxis.set_major_formatter(FuncFormatter(lead0))
    return axes

def vlim(D):
    allv = np.concatenate([np.abs(D["gene_fc"](c, g)) for c, _ in CELLS for g in ("mid", "high")])
    m = np.percentile(allv, 96.5); return (-m*1.28, m*1.28)
YLM, YLS = vlim(MA), vlim(SA)

fig = plt.figure(figsize=(7.4, 10.4))
outer = fig.add_gridspec(5, 1, height_ratios=[0.66, 1.0, 1.0, 0.96, 0.96], hspace=0.62,
                         left=0.11, right=0.985, top=0.955, bottom=0.03)
WR = [0.335, 0.665]
YT = "Log$_2$FC  III,IV or V,VI / low\n(Response to ER stress, genes=%d)"
DYT = "Response to ER stress score\n(per donor, z; genes=%d)"

goGS = outer[0].subgridspec(1, 2, width_ratios=[0.60, 0.40])
axGO = fig.add_subplot(goGS[0]); go_merged(axGO)

vblocks = []
for ri, (D, yl, lg) in enumerate([(MA, YLM, True), (SA, YLS, False)]):
    sub = outer[1+ri].subgridspec(1, 2, width_ratios=WR, wspace=0.10)
    axn = fig.add_subplot(sub[0]); axg = fig.add_subplot(sub[1])
    violin_panel(axn, D, NEUR, YT % D["ngene"], yl, legend=lg, show_y=True)
    violin_panel(axg, D, GLIA, "", yl, legend=False, show_y=False)
    vblocks.append((axn, axg))

dblocks = []
for ri, D in enumerate([MA, SA]):
    sub = outer[3+ri].subgridspec(1, 2, width_ratios=WR, wspace=0.22)
    an = donor_panel(fig, sub[0].subgridspec(1, 2, wspace=0.55), D, NEUR, DYT % D["ngene"], show_y=True)
    ag = donor_panel(fig, sub[1].subgridspec(1, 4, wspace=0.55), D, GLIA, "", show_y=False)
    dblocks.append((an[0], ag[0]))

fig.canvas.draw()
def letter(ax, L, dx=-0.052, dy=0.012):
    p = ax.get_position(); fig.text(max(p.x0 + dx, 0.004), p.y1 + dy, L, fontsize=FS_LET, fontweight="bold", va="bottom", ha="left")
letter(axGO, "A", dx=-0.052, dy=0.006)
(cn, ce), (dn, df) = vblocks
letter(cn, "B"); letter(ce, "C", dx=-0.030)
letter(dn, "D"); letter(df, "E", dx=-0.030)
(gn, gi), (hn, hj) = dblocks
letter(gn, "F"); letter(gi, "G", dx=-0.040)
letter(hn, "H"); letter(hj, "I", dx=-0.040)

import shutil
pdf = os.path.join(OUTDIR, "Figure4_REVISION_v2.0_20260826.pdf")
fig.savefig(pdf, facecolor="white"); shutil.copyfile(pdf, os.path.join(OUTDIR, "Figure4_REVISION_v2.0_20260826.ai"))
print("wrote", pdf, "| Mathys", MA["ngene"], "SEA-AD", SA["ngene"])
