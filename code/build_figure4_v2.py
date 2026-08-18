#!/usr/bin/env python3
"""Figure 4 restructure (author): a · b/c · D=UPR-sensor control (unchanged) · e/f/g=TRIB3/TMED2 correlations.
The old donor-level module-score row moves to Supplementary. VECTOR-MOVE approach:
 - crop the ORIGINAL vector bands (a+b/c ; d/e/f correlations) with cropbox+pdftocairo (keeps vector, no raster),
 - the panel-letter relabel d/e/f -> e/f/g is done by EDITING the real vector text (d->e, e->f) in the content
   stream, and (since 'g' is not in the Arial-Bold subset) DELETING the old 'f' glyph and drawing a fresh 'g'
   (no cover-and-redraw), plus a fresh 'D' for the sensor panel.
 - D = trimmed snRNA sensor donor-level violins (f4_D_sensor.pdf)."""
import os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
import pypdf
from pypdf import Transformation
from pypdf.generic import RectangleObject

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
        "ADBrainSeq_V6.0/Acta Neuropathologica Communication")
SCR = "/tmp/adbrainseq_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(HERE)
ORIG = os.path.join(BASE, "main figures_05182026 4.pdf")

# ---- 1) edit the panel-letter vector text in the original content stream ----
rd = pypdf.PdfReader(ORIG); pg = rd.pages[0]
data = pg.get_contents().get_data().decode("latin1")
data = data.replace("31.9927 322.9941 Tm\n(d)Tj", "31.9927 322.9941 Tm\n(e)Tj", 1)          # panel1 d -> e
data = data.replace("186.9399 322.9941 Tm\n(e)Tj\n8.36 0 Td\n(f)Tj",
                    "186.9399 322.9941 Tm\n(f)Tj", 1)                                        # panel2 e->f ; drop old panel3 'f'
from pypdf.generic import DecodedStreamObject
newco = DecodedStreamObject(); newco.set_data(data.encode("latin1"))
w0 = pypdf.PdfWriter(); pe = w0.add_page(pg)
ref = w0._add_object(newco); pe[pypdf.generic.NameObject("/Contents")] = ref
edited = os.path.join(SCR, "fig4_edited.pdf")
with open(edited, "wb") as f: w0.write(f)

# ---- 2) crop the vector bands (cropbox + pdftocairo keeps vector) ----
def bake(src, crop, out):
    r = pypdf.PdfReader(src); p = r.pages[0]
    p.cropbox = RectangleObject(crop); p.mediabox = RectangleObject(crop)
    ww = pypdf.PdfWriter(); ww.add_page(p); raw = out + ".raw.pdf"
    with open(raw, "wb") as f: ww.write(f)
    os.system('pdftocairo -pdf "%s" "%s" 2>/dev/null' % (raw, out))
TOP = os.path.join(SCR, "f4v2_top.pdf"); CORR = os.path.join(SCR, "f4v2_corr.pdf")
bake(edited, [13, 347, 599, 794], TOP)      # a + b/c
bake(edited, [13, 103, 599, 345], CORR)     # e/f/g correlations (g letter deleted -> drawn fresh)
D = os.path.join(HERE, "f4_D_sensor.pdf")

top = pypdf.PdfReader(TOP).pages[0]; corr = pypdf.PdfReader(CORR).pages[0]; dpan = pypdf.PdfReader(D).pages[0]
def wh(p): return float(p.mediabox.width), float(p.mediabox.height)
wt, ht = wh(top); wc, hc = wh(corr); wd, hd = wh(dpan)

# ---- 3) assemble: TOP / D / CORR ----
Wt = 586.0
st = Wt / wt; sc = Wt / wc; sd = Wt / wd
htS, hcS, hdS = ht * st, hc * sc, hd * sd
M, GAP = 14.0, 14.0
Wp = 612.0; Hp = M + hcS + GAP + hdS + GAP + htS + M
xc = xd = xt = (Wp - Wt) / 2
yC = M; yD = M + hcS + GAP; yT = yD + hdS + GAP
w = pypdf.PdfWriter(); page = w.add_blank_page(width=Wp, height=Hp)
page.merge_transformed_page(corr, Transformation().scale(sc).translate(xc, yC))
page.merge_transformed_page(dpan, Transformation().scale(sd).translate(xd, yD))
page.merge_transformed_page(top,  Transformation().scale(st).translate(xt, yT))

# ---- 4) fresh vector letters: 'D' (sensor) and 'g' (correlation panel 3) ----
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
ax.text(xd + 2, yD + hdS - 2, "D", fontsize=18, fontweight="bold", family="Arial", ha="left", va="top")
# panel3 'g' : original letter at (337,323); in CORR crop [13,103] -> baked (324, 220); placed at (xc,yC) scale sc
gx = xc + (337 - 13) * sc; gy = yC + (345 - 103 - (345 - 323)) * sc   # baked-y top-ref
ax.text(gx, yC + (323 - 103) * sc + 6, "g", fontsize=18, fontweight="bold", family="Arial", ha="left", va="top")
OVL = os.path.join(SCR, "f4v2_letters.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
page.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "Figure4_REVISION_v2_20260817.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| page %.0f x %.0f" % (Wp, Hp))
