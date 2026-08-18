#!/usr/bin/env python3
"""R1.5 TABLE (English) — DOWN block first, UP block second.
Per direction: Mizuno, Nativio, Bulk-overlapped (both cohorts), snRNA total,
and the overlap between the bulk-overlapped set and snRNA (n and % of snRNA)."""
import sys, math, os, csv
sys.path.insert(0, '.')
from r1q1_gomatrix import load_mizuno, load_nativio
import openpyxl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "Arial", "font.sans-serif": ["Arial"], "axes.unicode_minus": False})

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
OUT = os.path.join(BASE, "Major Revision", "processed raw data")

GS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
      ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"),
      ("ERAD (75)", "geneset_ERAD.txt")]
sets = {n: set(l.strip() for l in open(os.path.join(OUT, f)) if l.strip()) for n, f in GS}

miz = {g: v[0] for g, v in load_mizuno().items()}
nat = {g: v[0] for g, v in load_nativio().items()}
wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD6_snRNSeqDB.xlsx"), read_only=True); ws = wb.active
CT = {"Ex": 2, "In": 4, "Mic": 8, "Oli": 10}
sn = {}
for r in ws.iter_rows(min_row=2, values_only=True):
    if not r[0]: continue
    d = {}
    for ct, ci in CT.items():
        try:
            fc = float(r[ci])
            if fc > 0 and math.isfinite(fc): d[ct] = math.log2(fc)
        except (TypeError, ValueError): pass
    if len(d) == 4: sn[str(r[0]).strip()] = d
wb.close()

def pct(a, b): return f"{a} ({100*a//max(b,1)}%)"

rows = []
for name, _ in GS:
    U = [g for g in sets[name] if g in miz and g in nat and g in sn]
    miz_d = sum(1 for g in U if miz[g] < 0); miz_u = sum(1 for g in U if miz[g] > 0)
    nat_d = sum(1 for g in U if nat[g] < 0); nat_u = sum(1 for g in U if nat[g] > 0)
    bd = set(g for g in U if miz[g] < 0 and nat[g] < 0)        # bulk overlapped down
    bu = set(g for g in U if miz[g] > 0 and nat[g] > 0)        # bulk overlapped up
    snd = set(g for g in U if sn[g]["Ex"] < 0 and sn[g]["In"] < 0)   # snRNA neuron down (total)
    snu = set(g for g in U if sn[g]["Mic"] > 0 and sn[g]["Oli"] > 0) # snRNA glia up (total)
    ov_d = bd & snd; ov_u = bu & snu
    rows.append([name,
                 miz_d, nat_d, len(bd), len(snd), pct(len(ov_d), len(bd)),
                 miz_u, nat_u, len(bu), len(snu), pct(len(ov_u), len(bu))])

headers = ["UPR\ngene set",
           "Mizuno\nDOWN", "Nativio\nDOWN", "Bulk\noverlapped\nDOWN", "snRNA\nneuron\nDOWN",
           "Reproduced\nin snRNA\n(% of bulk)",
           "Mizuno\nUP", "Nativio\nUP", "Bulk\noverlapped\nUP", "snRNA\nglia\nUP",
           "Reproduced\nin snRNA\n(% of bulk)"]

with open(os.path.join(OUT, "R1Q5_reproducibility_table.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow([h.replace("\n", " ") for h in headers]); w.writerows(rows)

colW = [0.105] + [0.0895] * 10
fig, ax = plt.subplots(figsize=(20.5, 4.3)); ax.axis("off")
fig.subplots_adjust(left=0.008, right=0.992, top=0.78, bottom=0.13)
tbl = ax.table(cellText=rows, colLabels=headers, cellLoc="center", loc="center",
               colWidths=colW, bbox=[0, 0, 1, 1])
tbl.auto_set_font_size(False); tbl.set_fontsize(9.5)
for (r, c), cell in tbl.get_celld().items():
    cell.set_edgecolor("#c8c8c8")
    if r == 0:
        cell.set_text_props(fontweight="bold", color="#1b2a3a"); cell.set_height(0.46)
        cell.set_facecolor("#dfe6ec" if c <= 5 else "#efe1df")
    else:
        cell.set_height(0.28)
        if c == 0: cell.set_text_props(fontweight="bold")
        elif c in (1, 2, 3, 4): cell.set_facecolor("#eaf1f7")
        elif c == 5: cell.set_facecolor("#d5e6f0"); cell.set_text_props(fontweight="bold", color="#12507e")
        elif c in (6, 7, 8, 9): cell.set_facecolor("#fbf0ef")
        elif c == 10: cell.set_facecolor("#f6ddd9"); cell.set_text_props(fontweight="bold", color="#9a2b1e")

xdiv = 0.105 + 5 * 0.0895
ax.plot([xdiv, xdiv], [0, 1.16], color="#333", lw=2.2, transform=ax.transAxes, clip_on=False)
ax.text(0.105 + 2.5 * 0.0895, 1.06, "DOWN-regulated", transform=ax.transAxes, ha="center",
        fontsize=12.5, fontweight="bold", color="#12507e")
ax.text(xdiv + 2.5 * 0.0895, 1.06, "UP-regulated", transform=ax.transAxes, ha="center",
        fontsize=12.5, fontweight="bold", color="#9a2b1e")

fig.suptitle("Cross-modality directional reproducibility of UPR gene sets in AD brain "
             "transcriptomics analysis", fontsize=13.5, fontweight="bold", y=0.975)
fig.text(0.5, 0.035,
         "Bulk overlapped = genes changed in the same direction in both bulk cohorts (Mizuno and Nativio). "
         "Overlap = genes both bulk-overlapped and changed in the same direction in single-nucleus RNA-seq "
         "(neurons for down, glia for up); % is of the bulk-overlapped set (i.e., reproduced in snRNA).",
         ha="center", fontsize=8.6, color="#555")
for ext in ("png", "svg"):
    fig.savefig(os.path.join(OUT, f"R1Q5_reproducibility_table.{ext}"),
                dpi=300 if ext == "png" else None, facecolor="white", bbox_inches="tight")
print("saved R1Q5_reproducibility_table.{png,svg,csv}")
for r in rows: print(r)
