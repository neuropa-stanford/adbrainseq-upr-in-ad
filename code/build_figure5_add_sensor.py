#!/usr/bin/env python3
"""Figure 5 (author): keep the ORIGINAL Fig 5 (A-H = PERK/IRE1/ATF6/ERAD gene sets x neuron/glia) untouched,
and ADD the UPR-sensor control panel as I in the empty lower half (VECTOR-place, no raster/overlay of panels).
I = trimmed snRNA sensor donor-level violins (f4_D_sensor.pdf) — sensors unchanged (KW ns) = 'signature, not sensor'."""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
import pypdf
from pypdf import Transformation

BASE = ("/data/adbrainseq")
ANC = os.path.join(BASE, "Acta Neuropathologica Communication")
SCR = "/tmp/adbrainseq_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(HERE)
FIG5 = os.path.join(BASE, "11292025_Figure5_UPR gene set in scRNAseq_MLM_GP.ai")
SENSOR = os.path.join(HERE, "f4_D_sensor.pdf")

base = pypdf.PdfReader(FIG5).pages[0]
Wp, Hp = float(base.mediabox.width), float(base.mediabox.height)     # 612 x 792
sens = pypdf.PdfReader(SENSOR).pages[0]
ws, hs = float(sens.mediabox.width), float(sens.mediabox.height)      # 1144 x 622

# fit into the empty lower half (panels A-H end at pdf-y ~429): scale by the tighter of width/height
Htarget = 392.0; Wmax = 586.0
s = min(Wmax / ws, Htarget / hs)
wsS, hsS = ws * s, hs * s
x = (Wp - wsS) / 2                             # centred
yTop = 410.0
y = yTop - hsS

w = pypdf.PdfWriter(); pg = w.add_page(base)
pg.merge_transformed_page(sens, Transformation().scale(s).translate(x, y))

# fresh 'I' panel letter (match A-H style: Arial Bold ~18 pt), top-left of the sensor panel
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
# three sensor rows (PERK/IRE1/ATF6) get their own panel letters I / J / K, at each row's top-left
# subplot layout: top=0.965, bottom=0.05, 3 rows + hspace 0.40 -> row height 0.2408 of fig; row tops:
for frac, lab in zip([0.965, 0.628, 0.291], ["I", "J", "K"]):
    ax.text(x - 2, y + frac * hsS + 9, lab, fontsize=18, fontweight="bold", family="Arial", ha="left", va="top")
OVL = os.path.join(SCR, "f5_I_letter.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
pg.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "Figure5_REVISION_20260817.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| sensor I at x%.0f y%.0f  (%.0f x %.0f)" % (x, y, wsS, hsS))
