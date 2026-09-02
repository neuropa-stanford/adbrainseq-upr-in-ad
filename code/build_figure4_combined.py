#!/usr/bin/env python3
"""Figure 4 — Mathys (left) vs SEA-AD (right), fully REGENERATED from data so every panel shares the
same style: gene DOTS in all violins, identical fonts / star sizes / violin widths / cell labels,
tight 2-column layout. Rows: GO | neuron violin | glia violin | neuron donor | glia donor. F/G/H
correlation scatter dropped. Letters A-J."""
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

DRD = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
       "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")
OUTDIR = os.path.dirname(DRD)
SENS = {"EIF2AK3", "ERN1", "ATF6"}
NEUR = [("Ex", "Excitatory\nneurons"), ("In", "Inhibitory\nneurons")]
GLIA = [("Ast", "Astrocytes"), ("Mic", "Microglia"), ("Oli", "Oligodendro-\ncytes"), ("OPC", "OPCs")]
CELLS = NEUR + GLIA
C_MID, C_HIGH = "#3c8f34", "#c8d94a"; RED, BLUE = "#e2231a", "#2b3a8f"
# GO bar palette: hue = direction (red=up / blue=down), shade = cohort (dark=Mathys / light=SEA-AD)
RED_D, RED_L = "#d1201f", "#f2a19b"; BLU_D, BLU_L = "#23348c", "#7387cf"
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
    ngene = len(genes)
    go = {}
    for c, _ in CELLS:
        v = gene_fc(c, "high"); p = wilcoxon(v).pvalue if len(v) > 5 else 1.0
        go[c] = (-math.log10(max(p, 1e-300)), np.median(v) < 0)
    return dict(gene_fc=gene_fc, score=score, go=go, ngene=ngene)

def load_mathys():
    # Published Mathys pipeline: per-donor mean expression (mean_DGE_by_donor_cell_type_*) with the
    # Braak annotation from donor_info.csv (braaksc). Rows in donor_info order (ids 1..48). This
    # reproduces the original Figure-4 B/C violins and D/E donor module scores (verified); it fixes the
    # panel-G flat artefact that the MATHYS24 pseudobulk (different Braak annotation) produced.
    DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/Human AD "
            "brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
    ES = [l.strip() for l in open(os.path.join(DRD, "ERstress_260_geneset.txt")) if l.strip() and l.strip() not in SENS]
    di = list(csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))))
    braak = [int(float(r["braaksc"])) for r in di]
    grpname = ["low" if b <= 2 else ("mid" if b <= 4 else "high") for b in braak]
    FILE = {"Ex": "Ex", "In": "In", "Ast": "Ast", "Mic": "Mic", "Oli": "Oli", "OPC": "Opc"}
    raw, present = {}, {}
    for c, _ in CELLS:
        rows = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{FILE[c]}.csv"))))
        hdr = rows[0]; gidx = {g: j for j, g in enumerate(hdr)}
        gp = [g for g in ES if g in gidx]; present[c] = set(gp)
        raw[c] = (gidx, rows[1:])
    genes = [g for g in ES if all(g in present[c] for c, _ in CELLS)]   # common set across all cell types
    byct = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for c, _ in CELLS:
        gidx, data = raw[c]; cols = [gidx[g] for g in genes]
        for i, r in enumerate(data):
            byct[c][grpname[i]].append(np.array([float(r[j]) for j in cols]))
    def gene_fc(c, grp):
        lo = np.vstack(byct[c]["low"]).mean(0); g = np.vstack(byct[c][grp]).mean(0); return g - lo
    score = {c: {"low": [], "mid": [], "high": []} for c, _ in CELLS}
    for c, _ in CELLS:
        order = [(v, "low") for v in byct[c]["low"]] + [(v, "mid") for v in byct[c]["mid"]] + [(v, "high") for v in byct[c]["high"]]
        M = np.vstack([v for v, _ in order]); mu = M.mean(0); sd = M.std(0, ddof=1); sd[sd == 0] = 1; Z = (M - mu) / sd
        for k, (_, grp) in enumerate(order): score[c][grp].append(float(Z[k].mean()))
    return dict(gene_fc=gene_fc, score=score, go=None, ngene=len(genes))

MA = load_mathys()
SA = load("SEAAD_3group_donor_pseudobulk.csv")
GX = {"low": 0, "mid": 1, "high": 2}

# published GO enrichment (-log10 p, is_down) per cell type — measured from the original Figure-4 A / SuppFig4 A
GO_VALS = {
    "MATHYS": {"Ex": (12.54, True), "In": (12.62, True), "Ast": (5.75, False), "Mic": (11.86, False), "Oli": (8.99, False), "OPC": (4.83, False)},
    "SEAAD":  {"Ex": (14.87, True), "In": (21.47, True), "Ast": (0.79, True),  "Mic": (2.15, False), "Oli": (2.39, False), "OPC": (1.59, False)},
}
def go_panel(ax, key):
    gv = GO_VALS[key]; yA = np.arange(len(CELLS))[::-1]
    rc, bc = (RED_D, BLU_D) if key == "MATHYS" else (RED_L, BLU_L)   # dark=Mathys / light=SEA-AD
    for i, (c, lab) in enumerate(CELLS):
        lp, down = gv[c]; ax.barh(yA[i], lp, color=(bc if down else rc), height=0.62)
    ax.axvline(1.3, color="#2b3a8f", lw=1.0)
    ax.set_xlim(0, max(v for v, _ in gv.values()) * 1.10)
    ax.set_yticks(yA); ax.set_yticklabels([l for _, l in CELLS], fontsize=FS_LAB, linespacing=0.85)
    ax.set_ylim(-0.6, len(CELLS) - 0.4)
    ax.set_xlabel(r"$-$log$_{10}$($p$-value)", fontsize=FS_LAB); ax.tick_params(labelsize=FS_TICK)
    ax.legend(handles=[Patch(fc=rc, label="Up"), Patch(fc=bc, label="Down")], fontsize=FS_LEG, frameon=False, loc="lower right", handlelength=1.0, borderpad=0.3)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)

def violin_panel(ax, D, cells, ylab, ylim, legend):
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
            if st != "ns": star_pos.append((xp, st))   # omit non-significant markers
    ax.set_xlim(-0.6, len(cells)-0.4); ax.set_ylim(*ylim); ax.set_xticks([])
    trs = blended_transform_factory(ax.transData, ax.transAxes)   # stars ABOVE the plot (clear of violins)
    for xp, st in star_pos:
        ax.text(xp, 1.015, st, transform=trs, ha="center", va="bottom", fontsize=FS_STAR)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for i, (c, lab) in enumerate(cells): ax.text(i, -0.03, lab, ha="center", va="top", fontsize=FS_LAB, transform=tr, linespacing=0.95)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.set_ylabel(ylab, fontsize=FS_YT); ax.tick_params(axis="y", labelsize=FS_TICK)
    if legend:
        ax.legend(handles=[Patch(fc=C_MID, ec="#1f1f1f", label="Braak III,IV / low"),
                           Patch(fc=C_HIGH, ec="#1f1f1f", label="Braak V,VI / low")],
                  fontsize=FS_LEG, frameon=True, loc="lower left", bbox_to_anchor=(0.0, 1.11), ncol=2, columnspacing=1.0, handlelength=1.0)

def donor_panel(fig, sub, D, cells, ylab):
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
        if spv < 0.05:   # show ρ only for significant donor-level correlations
            ax.text(0.5, 1.02, f"ρ={rho:+.2f} {pstar(spv)}", transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=FS_RHO, color="#b00")
        if idx == 0: ax.set_ylabel(ylab, fontsize=FS_YT)
        for s2 in ("top", "right"): ax.spines[s2].set_visible(False)
        ax.tick_params(axis="y", labelsize=FS_TICK); ax.yaxis.set_major_formatter(FuncFormatter(lead0))
    return axes

# violin y-limits per cohort — use a high percentile (not max) so a few outlier genes don't
# stretch the axis and squash the violins (esp. SEA-AD, whose effects are small).
def vlim(D):
    allv = np.concatenate([np.abs(D["gene_fc"](c, g)) for c, _ in CELLS for g in ("mid", "high")])
    m = np.percentile(allv, 96.5); return (-m*1.28, m*1.28)
YLM, YLS = vlim(MA), vlim(SA)

fig = plt.figure(figsize=(7.4, 9.5))
gs = fig.add_gridspec(5, 2, height_ratios=[0.85, 1.0, 1.0, 0.92, 0.92], hspace=0.80, wspace=0.42,
                      left=0.115, right=0.985, top=0.945, bottom=0.035)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])
axE = fig.add_subplot(gs[2, 0]); axF = fig.add_subplot(gs[2, 1])
go_panel(axA, "MATHYS"); go_panel(axB, "SEAAD")   # regenerated GO bars with the published enrichment values
YT = "Log$_2$FC  III,IV or V,VI / low\n(Response to ER stress, genes=%d)"
violin_panel(axC, MA, NEUR, YT % MA["ngene"], YLM, True)
violin_panel(axD, SA, NEUR, YT % SA["ngene"], YLS, True)
violin_panel(axE, MA, GLIA, YT % MA["ngene"], YLM, False)
violin_panel(axF, SA, GLIA, YT % SA["ngene"], YLS, False)
DYT = "Response to ER stress score\n(per donor, z; genes=%d)"
donor_panel(fig, gs[3, 0].subgridspec(1, 2, wspace=0.5), MA, NEUR, DYT % MA["ngene"])
donor_panel(fig, gs[3, 1].subgridspec(1, 2, wspace=0.5), SA, NEUR, DYT % SA["ngene"])
donor_panel(fig, gs[4, 0].subgridspec(1, 4, wspace=0.55), MA, GLIA, DYT % MA["ngene"])
donor_panel(fig, gs[4, 1].subgridspec(1, 4, wspace=0.55), SA, GLIA, DYT % SA["ngene"])

# C/D are only 2 cell types (neurons) — shrink their width 40% so the cells match E/F's width
for ax in (axC, axD):
    p = ax.get_position(); ax.set_position([p.x0, p.y0, p.width * 0.60, p.height])

# letters — placed at the figure's LEFT edge of each column, clear of y-titles/content
fig.canvas.draw()
LX = {0: 0.012, 1: 0.512}   # x for left / right column letters
# A/B (GO, no rotated y-title) sit just above their panel; the rest need extra lift to clear the y-title top
OFF = {"A": 0.003, "B": 0.003}
for ax, L, col in [(axA, "A", 0), (axB, "B", 1), (axC, "C", 0), (axD, "D", 1), (axE, "E", 0), (axF, "F", 1)]:
    p = ax.get_position(); fig.text(LX[col], p.y1 + OFF.get(L, 0.020), L, fontsize=FS_LET, fontweight="bold", va="bottom", ha="left")
for cell, L, col in [(gs[3, 0], "G", 0), (gs[3, 1], "H", 1), (gs[4, 0], "I", 0), (gs[4, 1], "J", 1)]:
    p = cell.get_position(fig); fig.text(LX[col], p.y1 + 0.0155, L, fontsize=FS_LET, fontweight="bold", va="bottom", ha="left")

import shutil
pdf = os.path.join(OUTDIR, "Figure4_REVISION_20260826.pdf")
fig.savefig(pdf, facecolor="white"); shutil.copyfile(pdf, os.path.join(OUTDIR, "Figure4_REVISION_20260826.ai"))
print("wrote", pdf, "| Mathys genes", MA["ngene"], "SEA-AD", SA["ngene"])
