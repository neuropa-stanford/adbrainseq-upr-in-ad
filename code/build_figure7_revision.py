#!/usr/bin/env python3
"""Assemble Figure 7 (was Figure 6, IHC) per author (corrected):
  a           H&E / phospho-Tau (AT8), control vs tauopathy  -> KEPT from old Figure 6 panel a
  b..g        07222026_Human brain staining.ai  (representative IF + quantification):
              A neurons img -> b ; B neuron quant -> c ; C microglia img -> d ; D microglia quant -> e ;
              E oligo img -> f ; F oligo quant -> g
Panel a keeps its native 'a'. Human-brain-staining 'A..F' are covered and re-lettered b..g in the
established style (Arial Bold 18.3 pt). All vector.
(The Neuron_MicroG_OligD multi-case montage is a separate SUPPLEMENTARY figure, built elsewhere.)"""
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
A = pypdf.PdfReader(os.path.join(SCR, "f7_a.pdf")).pages[0]      # 590 x 190 (AT8, keeps its own 'a')
HS = pypdf.PdfReader(os.path.join(SCR, "f7_hs.pdf")).pages[0]    # 601 x 755 (Human brain staining, A..F)
def wh(p): return float(p.mediabox.width), float(p.mediabox.height)
wa, ha = wh(A); whs, hhs = wh(HS)

Wt = 595.0
sa = Wt / wa; shs = Wt / whs
haS = ha * sa; hhsS = hhs * shs
MARG, GAP = 14.0, 16.0
Wp = 612.0
Hp = MARG + hhsS + GAP + haS + MARG
xa = xhs = (Wp - Wt) / 2
yHS = MARG                      # Human-brain-staining body (bottom)
yA = MARG + hhsS + GAP          # AT8 panel a (top)

w = pypdf.PdfWriter(); pg = w.add_blank_page(width=Wp, height=Hp)
pg.merge_transformed_page(HS, Transformation().scale(shs).translate(xhs, yHS))
pg.merge_transformed_page(A,  Transformation().scale(sa).translate(xa, yA))

# ---- re-letter A..F -> b..g (cover source letter, draw new) ; panel a keeps native 'a' ----
# HS baked page height 755 covers original pdf-y [17,772]; original y-from-top Yt -> baked-y = 775 - Yt
# assembled-y = yHS + baked-y*shs ; assembled-x = xhs + (orig_x-7)*shs
def place(Yt_top, Yt_bot, x0, x1):
    ay_top = yHS + (775 - Yt_top) * shs
    ay_bot = yHS + (775 - Yt_bot) * shs
    ax0 = xhs + (x0 - 7) * shs
    ax1 = xhs + (x1 - 7) * shs
    return ax0, ax1, ay_bot, ay_top
# source A..F: (Yt_top, Yt_bot, x0, x1) -> new letter
SRC = [
    (21.2, 35.4, 7.5, 15.0, "b"),   # A neurons img
    (69.1, 83.3, 518.8, 526.3, "c"),# B neuron quant
    (274.0, 288.2, 7.5, 15.0, "d"), # C microglia img
    (315.6, 329.9, 518.8, 526.3, "e"),# D microglia quant
    (530.6, 544.8, 7.5, 14.5, "f"), # E oligo img
    (565.9, 580.1, 518.8, 525.2, "g"),# F oligo quant
]
LP = 18.3
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
for yt, yb, x0, x1, new in SRC:
    ax0, ax1, ay_bot, ay_top = place(yt, yb, x0, x1)
    ax.add_patch(Rectangle((ax0 - 3, ay_bot - 3), (ax1 - ax0) + 7, (ay_top - ay_bot) + 6,
                           facecolor="white", edgecolor="none", zorder=5))
    ax.text(ax0 - 1, ay_top + 1.5, new, fontsize=LP, fontweight="bold", family="Arial",
            color="black", ha="left", va="top", zorder=6)
OVL = os.path.join(SCR, "f7_letters.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
pg.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "Figure7_REVISION_20260814.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| page %.0f x %.0f" % (Wp, Hp))
