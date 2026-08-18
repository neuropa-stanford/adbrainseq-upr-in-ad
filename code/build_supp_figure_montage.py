#!/usr/bin/env python3
"""Supplementary Figure — reviewer-requested per-donor IHC montage (07222026_Neuron_MicroG_OligD):
all cases (#1-#3, triplicate) of Neuron/TRIB3, Microglia CD45/TMED2, OligD MOG/TMED2, Braak 0 vs VI.
Relabel section heads 'A:/B:/C:' -> a/b/c in the SAME style as the main figures (Arial Bold 18.3 pt),
keeping the 'Neurons/Microglia/OligD' descriptors. All vector."""
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
M_ = pypdf.PdfReader(os.path.join(SCR, "f7_m.pdf")).pages[0]     # 564 x 576 montage (crop of new file, baked)
wm, hm = float(M_.mediabox.width), float(M_.mediabox.height)

LM = 22.0                       # left margin holds the panel letters clear of the section titles
Wp = wm + LM + 6; Hp = hm + 12
xm = LM; ym = 6
w = pypdf.PdfWriter(); pg = w.add_blank_page(width=Wp, height=Hp)
pg.merge_transformed_page(M_, Transformation().translate(xm, ym))

# section heads on f7_m (baked from crop [7,211,571,787], height 576): baked-y = 581 - Yt ; baked-x = orig_x-7
# assembled = + (xm, ym)
SEC = [(5.9, 16.6, 9.2, 17.4, "a"), (199.4, 210.0, 9.2, 17.4, "b"), (398.1, 408.8, 9.2, 17.4, "c")]
LP = 18.3
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
for yt, yb, x0, x1, new in SEC:
    ay_top = ym + (581 - yt); ay_bot = ym + (581 - yb)
    ax0 = xm + (x0 - 7); ax1 = xm + (x1 - 7)
    # cover only the 'A:/B:/C:' prefix (keep the Neurons/Microglia/OligD descriptor)
    ax.add_patch(Rectangle((ax0 - 3, ay_bot - 3), (ax1 - ax0) + 5, (ay_top - ay_bot) + 6,
                           facecolor="white", edgecolor="none", zorder=5))
    # panel letter in the left margin, clear of the descriptor
    ax.text(4, ay_top + 2.0, new, fontsize=LP, fontweight="bold", family="Arial",
            color="black", ha="left", va="top", zorder=6)
OVL = os.path.join(SCR, "suppfig_letters.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
pg.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "SupplementaryFigure_IHC_percase_montage_20260814.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| page %.0f x %.0f" % (Wp, Hp))
