#!/usr/bin/env python3
"""Assemble Figure 5 (revision) — NEW: independent replication in SEA-AD (R2.9).
  A  SEA-AD 3-group by cell type: neuronal ER-stress-associated transcripts decline; donor-level
     inhibitory-neuron trend rho=-0.24, p=0.025 (Allen MTG, 84 donors, different platform).
  B  SEA-AD internal control: UPR sensors (PERK/IRE1/ATF6) near-flat -> the R2.1 logic reproduced.
Both panels code-generated (Arial). Descriptions in the figure legend (titles trimmed). All vector."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "pdf.fonttype": 42, "ps.fonttype": 42})
import pypdf
from pypdf import Transformation

HERE = os.path.dirname(os.path.abspath(__file__)); OUTDIR = os.path.dirname(HERE)
A = pypdf.PdfReader(os.path.join(HERE, "R2Q9_SEAAD_3group_figure.pdf")).pages[0]
B = pypdf.PdfReader(os.path.join(HERE, "R2Q9_SEAAD_sensor_control.pdf")).pages[0]
def wh(p): return float(p.mediabox.width), float(p.mediabox.height)
wa, ha = wh(A); wb, hb = wh(B)
Wref = max(wa, wb)
M, GAP = 16.0, 14.0
sa = Wref / wa; sb = Wref / wb
haS, hbS = ha * sa, hb * sb
Wp = Wref + 2 * M
Hp = M + hbS + GAP + haS + M
yA = M + hbS + GAP
yB = M
w = pypdf.PdfWriter(); pg = w.add_blank_page(width=Wp, height=Hp)
pg.merge_transformed_page(A, Transformation().scale(sa).translate(M, yA))
pg.merge_transformed_page(B, Transformation().scale(sb).translate(M, yB))

# letters (proportional to this wide figure ~= 3% of width, matching Figs 1-4)
LP = round(Wp * 0.029, 1)
fig = plt.figure(figsize=(Wp / 72, Hp / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, Wp); ax.set_ylim(0, Hp)
ax.text(6, yA + haS - 2, "A", fontsize=LP, fontweight="bold", family="Arial", ha="left", va="top")
ax.text(6, yB + hbS - 2, "B", fontsize=LP, fontweight="bold", family="Arial", ha="left", va="top")
OVL = "/tmp/adbrainseq_work/f5_letters.pdf"
fig.savefig(OVL, transparent=True); plt.close(fig)
pg.merge_transformed_page(pypdf.PdfReader(OVL).pages[0], Transformation())

out = os.path.join(OUTDIR, "Figure6_REVISION_20260817.pdf")
with open(out, "wb") as f: w.write(f)
print("wrote", out, "| page %.0f x %.0f | letter %.1f pt" % (Wp, Hp, LP))
