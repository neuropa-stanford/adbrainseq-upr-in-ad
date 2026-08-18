#!/usr/bin/env python3
"""Build merged Figure 1 (revision) — FULLY VECTOR (no rasterisation).
Pipeline (all vector):
  1. Figure1.pdf  --in-stream edits-->  Figure1_vedit2.pdf
       - D & G bar y-axis titles resized 11.4702 -> 12.0 pt (Tm scale) to match the 12.0 pt violin titles
       - panel letters retyped in the content stream: B->C, C->D, D->E, E->F  (glyphs present in subset)
  2. crop to the B-F region (drop original panel A) via cropbox + pdftocairo -pdf  -> fig1_BF2_vec.pdf (vector)
  3. compose: new GO consolidated panel (vector) on top + B-F (vector) below
  4. the single letter 'G' (old bars-F) is drawn as a vector overlay, because the embedded Arial-Bold
     SUBSET has no 'G' glyph (verified: an in-stream 'G' renders as .notdef box) -> one-letter vector cover.
Output: Figure1_REVISION_vector_20260813.pdf  (+ stable Figure1_MERGED_vector.pdf)
Prereqs in this dir: R1Q1_Figure_GO_consolidated.pdf. Uses pdftocairo (poppler)."""
import os, subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})
import pypdf
from pypdf import Transformation
from pypdf.generic import DecodedStreamObject, NameObject, RectangleObject

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
ASSETS = "/data/adbrainseq/Figures_ASSETS/originals_pdf/Figure1.pdf"
GO_PDF = os.path.join(HERE, "R1Q1_Figure_GO_consolidated.pdf")
SCR = "/tmp/adbrainseq_work"

# ---------- 1. in-stream vector edits (titles + letters) ----------
w = pypdf.PdfWriter(); w.append(pypdf.PdfReader(ASSETS)); pg = w.pages[0]
s = pg.get_contents().get_data().decode("latin-1")
REPL = [
    ("0 -11.4702 -11.4702 0 208.23 339.5859 Tm", "0 -12 -12 0 208.23 339.5859 Tm"),
    ("0 -6.6871 -6.6871 0 212.0493 317.6748 Tm", "0 -6.996 -6.996 0 212.2257 316.6624 Tm"),
    ("0 -11.4702 -11.4702 0 208.23 314.0156 Tm", "0 -12 -12 0 208.23 312.8344 Tm"),
    ("-21.192 -0.333 Td", "-20.355 -0.318 Td"),
    ("0 -6.6871 -6.6871 0 215.8682 535.1826 Tm", "0 -6.996 -6.996 0 216.045 534.172 Tm"),
    ("0 -11.4702 -11.4702 0 212.0488 531.5234 Tm", "0 -12 -12 0 212.049 530.341 Tm"),
    ("-3.71 -8.273 Td\n(B)Tj", "-3.71 -8.273 Td\n(C)Tj"),
    ("181.1353 179.6494 Tm\n/f-1-0 1 Tf\n(C)Tj", "181.1353 179.6494 Tm\n/f-1-0 1 Tf\n(D)Tj"),
    ("499.9238 179.6494 Tm\n/f-1-0 1 Tf\n(D)Tj", "499.9238 179.6494 Tm\n/f-1-0 1 Tf\n(E)Tj"),
    ("-0.026 -12.416 Td\n(E)Tj", "-0.026 -12.416 Td\n(F)Tj"),
]
for a, b in REPL:
    assert s.count(a) == 1, ("anchor not unique", a, s.count(a))
    s = s.replace(a, b)
co = DecodedStreamObject(); co.set_data(s.encode("latin-1")); pg[NameObject("/Contents")] = w._add_object(co)
vedit2 = os.path.join(SCR, "Figure1_vedit2.pdf")
with open(vedit2, "wb") as f: w.write(f)

# ---------- 2. crop to B-F (drop panel A), flatten to vector ----------
w2 = pypdf.PdfWriter(); w2.append(pypdf.PdfReader(vedit2))
BOX = RectangleObject([20, 198, 552, 632]); w2.pages[0].mediabox = BOX; w2.pages[0].cropbox = BOX
crop = os.path.join(SCR, "fig1_BF2_crop.pdf")
with open(crop, "wb") as f: w2.write(f)
BFV = os.path.join(SCR, "fig1_BF2_vec.pdf")
subprocess.run(["pdftocairo", "-pdf", crop, BFV], check=True)

# ---------- 3. compose GO (top) + B-F (bottom) ----------
go = pypdf.PdfReader(GO_PDF).pages[0]; gw, gh = float(go.mediabox.width), float(go.mediabox.height)
bf = pypdf.PdfReader(BFV).pages[0];    fw, fh = float(bf.mediabox.width), float(bf.mediabox.height)
M, GAP = 18, 12
Wp = max(fw, gw) + 2 * M
Hp = gh + GAP + fh + 2 * M
tx_go = M + (Wp - 2 * M - gw) / 2; ty_go = Hp - M - gh
tx_bf = M + (Wp - 2 * M - fw) / 2; ty_bf = M
w3 = pypdf.PdfWriter(); page = w3.add_blank_page(width=Wp, height=Hp)
page.merge_transformed_page(go, Transformation().translate(tx=tx_go, ty=ty_go))
page.merge_transformed_page(bf, Transformation().translate(tx=tx_bf, ty=ty_bf))

# ---------- 4. one-letter vector overlay: cover old bars-F, draw 'G' ----------
# old bars-F letter (cropped BF pts): glyph x[167,178] y[192,207], centre ~(172.5,200).
# The G y-title "...Braak 0,I,II" sits immediately to the right at x>=180, so the cover box MUST stop
# before x=180 or it erases the title's "0,I,II" (that was the bug). Tight box over the F glyph only.
gx = tx_bf + 172.5; gy = ty_bf + 200.0
cov = (tx_bf + 161, ty_bf + 189, 18, 22)  # x[161,179] y[189,211] over old F, clear of title (x>=180)
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
ax.add_patch(Rectangle((cov[0], cov[1]), cov[2], cov[3], facecolor="white", edgecolor="none"))
ax.text(gx, gy, "G", fontsize=18.26, fontweight="bold", family="Arial", color="black", ha="center", va="center")
OVL = os.path.join(SCR, "g_overlay.pdf"); fig.savefig(OVL, transparent=True); plt.close(fig)
page.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

final = os.path.join(OUTDIR, "Figure1_REVISION_vector_20260813.pdf")
with open(final, "wb") as f: w3.write(f)
with open(os.path.join(OUTDIR, "Figure1_MERGED_vector.pdf"), "wb") as f: w3.write(f)
print("wrote", final)
print("page %.1f x %.1f pt | GO %.1f x %.1f | BF %.1f x %.1f" % (Wp, Hp, gw, gh, fw, fh))
