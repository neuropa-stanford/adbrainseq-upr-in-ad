#!/usr/bin/env python3
"""Build merged Figure 1 (revision):
  Top   = GO consolidated vector panel (A matrix + B conserved genes)  [R1Q1_Figure_GO_consolidated.pdf]
  Bottom= original Figure-1 lower panels (violins + top-10 bars + SRPX), rasterised at 600 dpi from the
          published Figure1.pdf, relettered C,D,E,F,G, with the LEFT Y-AXIS TITLES re-typeset uniformly
          (violin panels C/F used a larger title than bar panels D/G in the original — equalise them).
Everything else in the bottom panels (bars, violins, axis numbers, stars, x-labels) is the untouched
original raster. Panel letters = Arial Bold 18.3 pt (style manual). Output = Figure1_REVISION_v{N}.pdf.
Run from this directory; expects f1_600.png + R1Q1_Figure_GO_consolidated.pdf in the scratchpad paths below.
"""
import os
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pypdf import PdfReader, PdfWriter, Transformation

plt.rcParams.update({
    "font.family": "Arial", "font.sans-serif": ["Arial"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "mathtext.default": "regular", "axes.unicode_minus": False,
})

SCR = "/tmp/adbrainseq_work"
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)  # Major Revision/
GO_PDF = os.path.join(HERE, "R1Q1_Figure_GO_consolidated.pdf")
F1_600 = os.path.join(SCR, "f1_600.png")

DPI = 600
LP = 18.3          # panel-letter size (pt)  -- style manual
TS = 7.8           # UNIFORM y-axis-title size (pt) for C,D,F,G (fits the two-line violin titles)

# ---- reconstruct the bottom-panel crop from the 600-dpi render of the original Figure1 ----
im = np.asarray(Image.open(F1_600).convert("RGB"))
g = np.asarray(Image.open(F1_600).convert("L"))
cut = 1272
bot = int(1227 * DPI / 150) + 2
xs = np.where((g < 245).any(axis=0))[0]
x0, x1 = xs.min(), xs.max() + 1
crop = im[cut:bot, x0:x1]
ch, cw = crop.shape[:2]

# ---- panel letters (Arial Bold 18.3 pt) -- verified positions ----
relabel = [(69, 169, "C"), (1371, 199, "D"), (4071, 199, "E"), (76, 2044, "F"), (1396, 2088, "G")]

# ---- uniform re-typeset Y-axis titles (crop pixel coords) ----
# Robust cover: white-box the whole title STRIP (left of the axis numbers, within the panel's y-range;
# nothing but the old title lives there), then redraw one uniform title. Boxes bounded to clear the axis
# numbers (right) and the green GSE dataset labels (below the violins).
# each: (white-box x0,x1,y0,y1), list of (cx,cy,text)
L2 = "Log$_2$fc"
RESP = "(Response to ER stress, genes = 260)"
titles = [
    # C violin (GSE173955): axis numbers ~x398, GSE label lives at x>500 (not in strip) -> strip x[0,393] y[100,1560]; two lines
    (dict(bx=(0, 393, 100, 1560)),
     [(122, 875, f"{L2} Braak V,VI / Braak 0,I,II,III"), (56, 875, RESP)]),
    # D bars: numbers ~x1611 -> strip x[1400,1600]; x-labels y>1550 -> y[200,1540]; one line
    (dict(bx=(1400, 1600, 200, 1540)),
     [(1505, 990, f"{L2} Braak V,VI / Braak 0,I,II,III")]),
    # F violin (GSE159699): old title y~1981..3620, GSE label at x>500 (not in strip) -> strip x[0,393] y[1990,3635]; two lines
    (dict(bx=(0, 393, 1990, 3635)),
     [(122, 2820, f"{L2} Braak V,VI / Braak 0,I,II"), (56, 2820, RESP)]),
    # G bars: old title top y~1983, numbers ~x1556 -> strip x[1400,1554] y[1960,3345]; one line
    (dict(bx=(1400, 1554, 1960, 3345)),
     [(1505, 2820, f"{L2} Braak V,VI / Braak 0,I,II")]),
]

fig = plt.figure(figsize=(cw / DPI, ch / DPI), dpi=DPI)
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.imshow(crop, zorder=0); ax.set_xlim(0, cw); ax.set_ylim(ch, 0)

# re-typeset titles
for spec, lines in titles:
    bx = spec["bx"]
    ax.add_patch(Rectangle((bx[0], bx[2]), bx[1] - bx[0], bx[3] - bx[2],
                           facecolor="white", edgecolor="none", zorder=5))
    for cx, cy, txt in lines:
        ax.text(cx, cy, txt, rotation=90, rotation_mode="anchor",
                ha="center", va="center", fontsize=TS, family="Arial",
                color="black", zorder=6)

# panel letters
for lx, ly, L in relabel:
    ax.text(lx, ly, L, fontsize=LP, fontweight="bold", family="Arial", color="black",
            ha="center", va="center",
            bbox=dict(boxstyle="square,pad=0.42", facecolor="white", edgecolor="white"), zorder=10)

BF = os.path.join(SCR, "f1_BF_relab.pdf")
fig.savefig(BF, dpi=DPI); plt.close(fig)

# ---- merge: GO (vector) on top, relabelled B-F raster on bottom ----
go = PdfReader(GO_PDF).pages[0]; gw, gh = float(go.mediabox.width), float(go.mediabox.height)
bf = PdfReader(BF).pages[0];     fw, fh = float(bf.mediabox.width), float(bf.mediabox.height)
M, GAP = 18, 12
Wp = max(fw, gw) + 2 * M
Hp = gh + GAP + fh + 2 * M
w = PdfWriter(); page = w.add_blank_page(width=Wp, height=Hp)
page.merge_transformed_page(go, Transformation().translate(tx=M + (Wp - 2 * M - gw) / 2, ty=Hp - M - gh))
page.merge_transformed_page(bf, Transformation().translate(tx=M + (Wp - 2 * M - fw) / 2, ty=M))

V = 19
final = os.path.join(OUTDIR, f"Figure1_REVISION_v{V}_20260813.pdf")
with open(final, "wb") as fh_:
    w.write(fh_)
# also stable name
with open(os.path.join(OUTDIR, "Figure1_MERGED_final.pdf"), "wb") as fh_:
    w.write(fh_)
print("wrote", final)
print("page pt", round(Wp, 1), round(Hp, 1), "| GO", round(gw, 1), round(gh, 1), "| BF", round(fw, 1), round(fh, 1))
