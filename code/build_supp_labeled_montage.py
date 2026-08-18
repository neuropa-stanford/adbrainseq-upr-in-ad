#!/usr/bin/env python3
"""Per-donor IHC montage labels — COPY THE ACTUAL SOURCE VECTORS (author: copy the vector, not the image).
Channel headers: the source draws each header row as one self-contained vector text block
  (BT  <CMYK colour> (Neuron)Tj <colour> (TRIB3)Tj ... ET), font /TT0 = Arial-BoldMT, no image.
We lift that exact block from the content stream and re-emit it at every panel via `q .. cm .. Q`
(a pure translation) -> identical font/bold/colour/style, zero image pixels copied.
Braak side boxes: stamped from a white-margin vector crop of #1-1 (no image either).
Source-labelled cases (#1-1 each section, OligD #1-2) are left untouched. All vector."""
import os, re
import numpy as np
from PIL import Image
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42})
import pypdf
from pypdf import Transformation
from pypdf.generic import RectangleObject, DecodedStreamObject, NameObject, ArrayObject

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/V1 Science manuscript/"
        "ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision")
SCR = "/tmp/adbrainseq_work"
HERE = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(HERE)
MONT = os.path.join(BASE, "07222026_Neuron_MicroG_OligD.pdf")
SKIP_CH = {"Neurons": {0}, "Microglia": {0}, "OligD": {0, 1}}   # cases already carrying source headers
SKIP_BK = {0, 4}                                                # leftmost cases carry source Braak boxes

# ---------- lift the source header vector blocks from the content stream ----------
_pg = pypdf.PdfReader(MONT).pages[0]
_data = _pg.get_contents().get_data().decode("latin1")
def lift(anchor):
    m = next(mm for mm in re.finditer(r"BT.*?ET", _data, re.S) if anchor in mm.group(0))
    b = m.group(0)
    # ensure the bold header font is selected inside the block
    if "Tf" not in b: b = b.replace("BT\n", "BT\n/TT0 1 Tf\n", 1)
    # 3-column variant = drop the trailing 3.5x label (its Td + Tj)
    b3 = re.sub(r"[\d.eE+\- ]+Td\s*\n\(3\.5x[^)]*\)Tj\s*\n", "", b)
    # the block's leading Tm (first label position) = source anchor (ax, ay)
    tm = re.search(r"([\d.eE+\-]+) 0 0 ([\d.eE+\-]+) ([\d.eE+\-]+) ([\d.eE+\-]+) Tm", b)
    return b, b3, float(tm.group(3)), float(tm.group(4))
HDR = {"Neurons": lift("(Neuron)"), "Microglia": lift("(CD45)"), "OligD": lift("(MOG)")}
# lift the source Braak-box vector (two stroked rectangles) + the Braak text block (both from #1-1 Neurons)
_mbox = re.search(r"[\d.]+ 731\.\d+ [\d.]+ -34\.\d+ re", _data)
_bs = _data.rfind("\nq\n", 0, _mbox.start()); _be = _data.find("\nQ", _mbox.start())
_boxblk = _data[_bs + 3:_be]                                   # inside the q..Q
_boxblk = re.sub(r"0 792 612 -792 re\s*\nW n\s*\n", "", _boxblk)   # drop the full-page clip
_boxblk = _boxblk.replace("8 M", "0.4 w\n8 M", 1)                  # thin border to match source
BRAAK_TXT = next(m.group(0) for m in re.finditer(r"BT.*?ET", _data, re.S)
                 if "(Braak 0)" in m.group(0) and "(Braak VI)" in m.group(0))
BRAAK = _boxblk + "\n" + BRAAK_TXT                              # box rectangles + rotated text, ready to translate

# ---------- geometry ----------
img = np.array(Image.open(os.path.join(SCR, "mon200-1.png")).convert("L")); Hpx, Wpx = img.shape; ppp = Wpx / 612.0
mask = ndimage.binary_closing(img < 235, structure=np.ones((3, 7)))
lab, n = ndimage.label(mask); tiles = []
for sl in ndimage.find_objects(lab):
    if sl is None: continue
    ys, xs = sl; w = (xs.stop - xs.start) / ppp; h = (ys.stop - ys.start) / ppp
    if 85 < w < 170 and 55 < h < 110: tiles.append([xs.start / ppp, ys.start / ppp, xs.stop / ppp, ys.stop / ppp, w, h])
def section_of(y0): return "Neurons" if y0 < 199 else ("Microglia" if y0 < 398 else "OligD")
bysec = {"Neurons": [], "Microglia": [], "OligD": []}
for t in tiles: bysec[section_of(t[1])].append(t)
def rowsplit(x0, y0, x1, y1):
    a = img[int(y0 * ppp):int(y1 * ppp), int(x0 * ppp):int(x1 * ppp)]
    wy = (a > 230).mean(1); m0 = int(len(wy) * 0.35); m1 = int(len(wy) * 0.65)
    return (y0 + y1) / 2 if m1 <= m0 else y0 + (m0 + int(np.argmax(wy[m0:m1]))) / ppp
def img_left(x0, y0, y1):
    # precise image-left of a NON-leftmost panel = first dark column after the white gap near x0
    ys0 = int((y0 + 4) * ppp); ys1 = int((y1 - 4) * ppp); seen = False
    for xi in range(int((x0 - 7) * ppp), int((x0 + 6) * ppp)):
        wf = (img[ys0:ys1, xi] > 228).mean()
        if wf > 0.7: seen = True
        elif seen and wf < 0.3: return xi / ppp
    return x0
BOXW = 9.2; PH = 792.0

# source reference (Neurons #1-1) image-left ~ detected-x0 + box width
_ncases = sorted(bysec["Neurons"], key=lambda t: (round(t[1] / 45), t[0]))
sXn = _ncases[0][0] + BOXW; sY0n = _ncases[0][1]

# ---------- per-panel: build header + braak vector copies (content-stream, pure translation) ----------
hdr_ops = []
for sec, cases in bysec.items():
    cases.sort(key=lambda t: (round(t[1] / 45), t[0]))
    blk4, blk3, ax, ay = HDR[sec]
    sX = cases[0][0] + BOXW; sY0 = cases[0][1]          # source #1-1 image-left, image-top (headers)
    for idx, (x0, y0, x1, y1, w, h) in enumerate(cases):
        Xi = x0 + BOXW if idx in {0, 4} else x0
        ncol = 4 if (x1 - Xi) > 130 else 3
        ILt = img_left(x0, y0, y1)                       # PRECISE image-left of this panel
        if idx not in SKIP_CH[sec]:
            dx = ILt - sX; dy = sY0 - y0                 # header copy: align precise image-left & image-top
            hdr_ops.append("q 1 0 0 1 %.3f %.3f cm\n%s\nQ" % (dx, dy, (blk4 if ncol == 4 else blk3)))
        if idx not in SKIP_BK:                            # braak copy: box right edge sits against the image
            dx = ILt - sXn; dy = sY0n - y0
            hdr_ops.append("q 1 0 0 1 %.3f %.3f cm\n%s\nQ" % (dx, dy, BRAAK))

# ---------- compose: montage + header vector copies (FIRST, /Contents is an array) + braak stamps ----------
w0 = pypdf.PdfWriter(); pg0 = w0.add_page(pypdf.PdfReader(MONT).pages[0])
co = DecodedStreamObject(); co.set_data(("\n".join(hdr_ops)).encode("latin1")); ref = w0._add_object(co)
contents = pg0.get("/Contents")                        # fresh page -> ArrayObject; append the label copies
if isinstance(contents, ArrayObject): contents.append(ref)
else: pg0[NameObject("/Contents")] = ArrayObject([contents, ref])
with open(os.path.join(SCR, "mont_labeled_raw.pdf"), "wb") as f: w0.write(f)
r = pypdf.PdfReader(os.path.join(SCR, "mont_labeled_raw.pdf")); p = r.pages[0]
p.cropbox = RectangleObject([7, 211, 571, 787]); p.mediabox = RectangleObject([7, 211, 571, 787])
w1 = pypdf.PdfWriter(); w1.add_page(p)
with open(os.path.join(SCR, "mont_labeled_crop_raw.pdf"), "wb") as f: w1.write(f)
os.system('pdftocairo -pdf "%s" "%s" 2>/dev/null' % (os.path.join(SCR, "mont_labeled_crop_raw.pdf"), os.path.join(SCR, "mont_labeled.pdf")))
M_ = pypdf.PdfReader(os.path.join(SCR, "mont_labeled.pdf")).pages[0]
wm, hm = float(M_.mediabox.width), float(M_.mediabox.height)
LM = 22.0; Wp = wm + LM + 6; Hp = hm + 12; xm = LM; ym = 6
w = pypdf.PdfWriter(); pg = w.add_blank_page(width=Wp, height=Hp)
pg.merge_transformed_page(M_, Transformation().translate(xm, ym))
fig2 = plt.figure(figsize=(Wp / 72, Hp / 72)); ax2 = fig2.add_axes([0, 0, 1, 1]); ax2.axis("off")
ax2.set_xlim(0, Wp); ax2.set_ylim(0, Hp)
for yt, yb, new in [(5.9, 16.6, "a"), (199.4, 210.0, "b"), (398.1, 408.8, "c")]:
    ay_top = ym + (581 - yt); ay_bot = ym + (581 - yb)
    ax2.add_patch(Rectangle((xm + 2.2 - 3, ay_bot - 3), 12, (ay_top - ay_bot) + 6, facecolor="white", edgecolor="none", zorder=5))
    ax2.text(4, ay_top + 2.0, new, fontsize=18.3, fontweight="bold", family="Arial", color="black", ha="left", va="top", zorder=6)
OVL2 = os.path.join(SCR, "supp_sec_letters.pdf"); fig2.savefig(OVL2, transparent=True); plt.close(fig2)
pg.merge_transformed_page(pypdf.PdfReader(OVL2).pages[0], Transformation())
out = os.path.join(OUTDIR, "SupplementaryFigure_IHC_percase_montage_20260814.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| header vector copies:", len(hdr_ops), "", "(braak+header vector copies)")
