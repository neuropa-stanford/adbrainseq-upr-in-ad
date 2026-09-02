#!/usr/bin/env python3
"""Figure 5 (revision): keep the ORIGINAL A-H (Mathys UPR-branch gene-set violins, hand-made .ai) with its
published values, and drop the SEA-AD reproduction I-P (fig5_IJKL_seaad.pdf, from build_figure5_seaad_8panel.py)
into the empty lower half. Vector-place at scale ~1.0 so the SEA-AD fonts print identical to A-H."""
import os, shutil
import pypdf
from pypdf import Transformation

DRD = os.path.dirname(os.path.abspath(__file__))          # .../Major Revision/processed raw data
OUTDIR = os.path.dirname(DRD)                             # .../Major Revision
V6 = os.path.dirname(os.path.dirname(OUTDIR))             # .../ADBrainSeq_V6.0
FIG5 = os.path.join(V6, "11292025_Figure5_UPR gene set in scRNAseq_MLM_GP.ai")
IJKL = os.path.join(DRD, "fig5_IJKL_seaad.pdf")

base = pypdf.PdfReader(FIG5).pages[0]
Wp, Hp = float(base.mediabox.width), float(base.mediabox.height)     # 612 x 792
child = pypdf.PdfReader(IJKL).pages[0]
wc, hc = float(child.mediabox.width), float(child.mediabox.height)   # 590 x 380

Htarget, Wmax, yTop = 396.0, 590.0, 414.0
s = min(Wmax / wc, Htarget / hc)
wcS, hcS = wc * s, hc * s
x = (Wp - wcS) / 2
y = yTop - hcS

w = pypdf.PdfWriter(); pg = w.add_page(base)
pg.merge_transformed_page(child, Transformation().scale(s).translate(x, y))

out = os.path.join(OUTDIR, "Figure5_REVISION_20260827.pdf")
with open(out, "wb") as f: w.write(f)
shutil.copyfile(out, os.path.join(OUTDIR, "Figure5_REVISION_20260827.ai"))
print("wrote", out, "| IJKL at x%.0f y%.0f (%.0f x %.0f) scale %.3f" % (x, y, wcS, hcS, s))
