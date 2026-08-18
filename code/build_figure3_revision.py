#!/usr/bin/env python3
"""Assemble Figure 3 revision from the pericyte-annotated reproduced panels (R1Q3):
   top row = A (tSNE + Pericytes) | B (marker violins + PDGFRB/Per column),
   bottom  = C (Braak-split tSNE + Per in legend) over D (proportions incl. Per).
All panels are matplotlib-generated (Arial, fixes the original's Helvetica/MyriadPro mix).
Panel-letter sizes were pre-scaled in the panel scripts so they land ~18.3 pt after this layout's scaling.
Output: Figure3_REVISION_20260813.pdf (vector; tSNE scatter is high-dpi raster inside, text vector)."""
import os
import pypdf
from pypdf import Transformation

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
A = pypdf.PdfReader(os.path.join(HERE, "R1Q3_Figure3a_reproduced.pdf")).pages[0]
B = pypdf.PdfReader(os.path.join(HERE, "R1Q3_Figure3b_reproduced.pdf")).pages[0]
CD = pypdf.PdfReader(os.path.join(HERE, "R1Q3_Figure3cd_reproduced.pdf")).pages[0]
def wh(p): return float(p.mediabox.width), float(p.mediabox.height)
wa, ha = wh(A); wb, hb = wh(B); wcd, hcd = wh(CD)

M, GAP, GAP2 = 18.0, 15.0, 18.0
HT = 320.0                      # top-row height
sa, sb = HT / ha, HT / hb
wa_s, wb_s = wa * sa, wb * sb
Wt = wa_s + GAP + wb_s          # top-row width
scd = Wt / wcd                  # bottom spans the top-row width
hcd_s = hcd * scd

Wp = Wt + 2 * M
Hp = HT + GAP2 + hcd_s + 2 * M
y_top = M + hcd_s + GAP2        # bottom of the top row

w = pypdf.PdfWriter(); page = w.add_blank_page(width=Wp, height=Hp)
page.merge_transformed_page(A,  Transformation().scale(sa).translate(M, y_top))
page.merge_transformed_page(B,  Transformation().scale(sb).translate(M + wa_s + GAP, y_top))
page.merge_transformed_page(CD, Transformation().scale(scd).translate(M, M))

out = os.path.join(OUTDIR, "Figure3_REVISION_20260813.pdf")
with open(out, "wb") as fh:
    w.write(fh)
print("wrote", out)
print("page %.1f x %.1f | A %.2f B %.2f CD %.2f | Wt %.1f hcd_s %.1f" % (Wp, Hp, sa, sb, scd, Wt, hcd_s))
