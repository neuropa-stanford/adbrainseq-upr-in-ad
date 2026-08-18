#!/usr/bin/env python3
"""Assemble Figure 4 revision (author feedback):
  A  cell-type GO (original)  -> scaled to 75%, centred
  B,C  original gene-level ER-stress violins (neurons / glia)
  D,E  NEW donor-level ER-stress module score (neurons / glia) -- black titles, enlarged, no d/q text
  F,G,H  original TRIB3/TMED2 correlations (re-lettered from D,E,F)
Continuous letters A..H (Arial Bold 18.3 pt). All vector."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "pdf.fonttype": 42, "ps.fonttype": 42})
import pypdf
from pypdf import Transformation

SCR = "/tmp/adbrainseq_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(HERE)
A = pypdf.PdfReader(os.path.join(SCR, "f4_A.pdf")).pages[0]
BC = pypdf.PdfReader(os.path.join(SCR, "f4_BC.pdf")).pages[0]
DN = pypdf.PdfReader(os.path.join(OUTDIR, "f4_donor_row.pdf")).pages[0]
DEF = pypdf.PdfReader(os.path.join(SCR, "f4_bot.pdf")).pages[0]
def wh(p): return float(p.mediabox.width), float(p.mediabox.height)
wa, ha = wh(A); wbc, hbc = wh(BC); wd, hd = wh(DN); wf, hf = wh(DEF)
Wref = wbc
M, GAP = 14.0, 10.0
sA = 0.75                       # panel A -> 75%
sBC = Wref / wbc; sD = Wref / wd; sF = Wref / wf
haS = ha * sA; waS = wa * sA
hbcS = hbc * sBC; hdS = hd * sD; hfS = hf * sF
Wp = Wref + 2 * M
Hp = M + hfS + GAP + hdS + GAP + hbcS + GAP + haS + M
yF = M
yD = M + hfS + GAP
yBC = yD + hdS + GAP
yA = yBC + hbcS + GAP
xA = M + (Wref - waS) / 2       # centre A

w = pypdf.PdfWriter(); pg = w.add_blank_page(width=Wp, height=Hp)
pg.merge_transformed_page(A,   Transformation().scale(sA).translate(xA, yA))
pg.merge_transformed_page(BC,  Transformation().scale(sBC).translate(M, yBC))
pg.merge_transformed_page(DN,  Transformation().scale(sD).translate(M, yD))
pg.merge_transformed_page(DEF, Transformation().scale(sF).translate(M, yF))

# ---- lettering overlay ----
LP = 18.3
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
def letter(x, y, s): ax.text(x, y, s, fontsize=LP, fontweight="bold", family="Arial", color="black", ha="left", va="top", zorder=6)
def cover(x, y, wdt, hgt): ax.add_patch(Rectangle((x, y), wdt, hgt, facecolor="white", edgecolor="none", zorder=5))
# A: cover the shrunk original 'A' (measured page pos ~x102,y835) and redraw at 18.3 pt in place
cover(97, 831, 20, 28); letter(99, 857, "A")
# B,C kept from original (in BC band)
# D,E on donor row (D = neurons top-left, E = glia / Astrocytes ~ 3rd of 6 panels)
dTop = yD + hdS
letter(M, dTop + 1, "D")
letter(M + Wref * 0.39, dTop + 1, "E")
# F,G,H : cover original D,E,F letters (page x 32/188/341) in the DEF band, redraw
fTop = yF + hfS
for xl, cw, new in [(32, 24, "F"), (188, 28, "G"), (341, 26, "H")]:
    cover(xl - 5, fTop - 26, cw, 34); letter(xl, fTop + 3, new)
OVL = os.path.join(SCR, "f4_letters.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
pg.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "Figure4_REVISION_20260814.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| page %.0f x %.0f | A@%.0f%%" % (Wp, Hp, sA * 100))
