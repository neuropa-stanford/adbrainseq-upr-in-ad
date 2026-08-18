#!/usr/bin/env python3
"""Figure 2 revision — in-stream VECTOR edits on the original Figure2.pdf (no rasterisation):
  ③ "(N gene set)" -> "(N-gene set)"  (grammar; hyphenated compound modifier; same length -> no reflow)
  ⑤ Y-axis titles 11 -> 12 pt (Tm scale) to match Figure 1's 12 pt titles; subscript "2" + anchors + the
     blk10 Td link recomputed programmatically (Log anchor kept fixed).  4 titles total.
  ⑥ panel-F "ns ns" 11.55 -> 11.00 pt to match every other annotation.
  ② verified separately: gene-set counts (30/31/73/69 vs 31/32/74/71) are dataset-specific detection
     (Nativio detects the full curated PERK31/IRE1 32/ATF6 74; Mizuno 1 fewer) -> CORRECT, left as-is.
Output: Figure2_REVISION_20260813.pdf (vector). Verify render + zero Type-3 fonts before trusting."""
import os
import pypdf
from pypdf.generic import DecodedStreamObject, NameObject

SRC = "/data/adbrainseq/Figures_ASSETS/originals_pdf/Figure2.pdf"
OUTDIR = "/data/adbrainseq/anc/Major Revision"

K = 12.0 / 11.0
def f(x): return ("%.4f" % x).rstrip("0").rstrip(".")

def title_repls(log, sub, fc, sub_old=6.413):
    """Given Log/sub/fc anchors (e,f) at size 11, return [(old_tm, new_tm), ...] resized to 12, Log fixed."""
    (le, lf), (se, sf), (fe, ff) = log, sub, fc
    out = []
    out.append((f"0 -11 -11 0 {f(le)} {f(lf)} Tm", f"0 -12 -12 0 {f(le)} {f(lf)} Tm"))
    nse, nsf, nss = le + (se - le) * K, lf + (sf - lf) * K, sub_old * K
    out.append((f"0 -6.413 -6.413 0 {f(se)} {f(sf)} Tm", f"0 -{f(nss)} -{f(nss)} 0 {f(nse)} {f(nsf)} Tm"))
    nfe, nff = le + (fe - le) * K, lf + (ff - lf) * K
    out.append((f"0 -11 -11 0 {f(fe)} {f(ff)} Tm", f"0 -12 -12 0 {f(nfe)} {f(nff)} Tm"))
    return out

w = pypdf.PdfWriter(); w.append(pypdf.PdfReader(SRC)); pg = w.pages[0]
s = pg.get_contents().get_data().decode("latin-1")

repls = []
# ③ hyphenate  " gene set" -> "-gene set"
n_gs = s.count(" gene set"); assert n_gs == 8, ("expected 8 gene-set labels, got", n_gs)
s = s.replace(" gene set", "-gene set")

# ⑥ ns size 11.55 -> 11 (keep kerning; verify position by render)
repls.append(("11.55 0 0 -11.55 212.5693 240.6479 Tm", "11 0 0 -11 212.5693 240.6479 Tm"))

# ⑤ titles -> 12 pt
# blk17 (E bars), blk21 (N bars): single titles
repls += title_repls((26.938, 391.5488), (30.6011, 370.7559), (26.938, 367.2793))
repls += title_repls((26.938, 759.7695), (30.6011, 738.9766), (26.938, 735.5))
# blk10 title1 (0,I,II,III)
repls += title_repls((51.4961, 208.4814), (55.1592, 187.4688), (51.4961, 183.959))
# blk10 title2 (0,I,II): reached by Td from title1-fc; replace Td+Log with an absolute Tm (drop fragile Td)
t1fc = (51.4961, 183.959)  # ORIGINAL title1-fc anchor
tx, ty = -35.602, 0.045
t2le = 0 * tx + (-11) * ty + t1fc[0]        # = title2 Log anchor e
t2lf = (-11) * tx + 0 * ty + t1fc[1]        # = title2 Log anchor f
repls.append(("-35.602 0.045 Td\n[( )9(L)9(o)9(g)]TJ",
              f"0 -12 -12 0 {f(t2le)} {f(t2lf)} Tm\n[( )9(L)9(o)9(g)]TJ"))
# title2 sub + fc (anchors relative to title2 Log)
def seg(anchor, orig, scale=None):
    ne, nf = t2le + (orig[0] - t2le) * K, t2lf + (orig[1] - t2lf) * K
    return ne, nf
s2 = seg(None, (54.6592, 554.5703)); ns2 = 6.413 * K
repls.append(("0 -6.413 -6.413 0 54.6592 554.5703 Tm", f"0 -{f(ns2)} -{f(ns2)} 0 {f(s2[0])} {f(s2[1])} Tm"))
f2 = seg(None, (50.9961, 551.0605))
repls.append(("0 -11 -11 0 50.9961 551.0605 Tm", f"0 -12 -12 0 {f(f2[0])} {f(f2[1])} Tm"))

for a, b in repls:
    c = s.count(a); assert c == 1, ("anchor not unique/found", repr(a), c); s = s.replace(a, b)
print("applied ③(8 hyphens) + %d Tm/Td replacements" % len(repls))

co = DecodedStreamObject(); co.set_data(s.encode("latin-1")); pg[NameObject("/Contents")] = w._add_object(co)
out = os.path.join(OUTDIR, "Figure2_REVISION_20260813.pdf")
with open(out, "wb") as fh: w.write(fh)
print("wrote", out)
