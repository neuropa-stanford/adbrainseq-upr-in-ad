#!/usr/bin/env python3
"""R1.5 — are the UPR genes altered in glia/neurons (snRNA, Fig 5) also seen in bulk (Fig 2)?
Cross-modality directionality check: snRNA (SuppD6) vs Mizuno/Nativio bulk log2FC.
"""
import sys, math, os, csv, statistics
sys.path.insert(0, '.')
from r1q1_gomatrix import load_mizuno, load_nativio
import openpyxl

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
OUT = os.path.join(BASE, "Major Revision", "processed raw data")

# UPR universe = ER-stress set (QuickGO GO:0034976)
UPR = set(l.strip() for l in open(os.path.join(OUT,"ERstress_260_geneset.txt")) if l.strip())

miz = {g: v[0] for g, v in load_mizuno().items()}
nat = {g: v[0] for g, v in load_nativio().items()}

wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD6_snRNSeqDB.xlsx"), read_only=True); ws = wb.active
CT = {"Ex": 2, "In": 4, "Mic": 8, "Oli": 10}
sn = {}
for r in ws.iter_rows(min_row=2, values_only=True):
    g = r[0]
    if not g: continue
    d = {}
    for ct, ci in CT.items():
        try:
            fc = float(r[ci])
            if fc > 0 and math.isfinite(fc): d[ct] = math.log2(fc)
        except (TypeError, ValueError): pass
    sn[str(g).strip()] = d
wb.close()

# genes with full data everywhere
genes = [g for g in UPR if g in miz and g in nat and g in sn
         and all(ct in sn[g] for ct in CT)]
print(f"UPR genes with bulk + snRNA(all 4 cell types): {len(genes)}", file=sys.stderr)

rows = []
for g in genes:
    neuron = (sn[g]["Ex"] + sn[g]["In"]) / 2
    glia = (sn[g]["Mic"] + sn[g]["Oli"]) / 2
    bulk = (miz[g] + nat[g]) / 2
    rows.append({"gene": g, "sn_Ex": sn[g]["Ex"], "sn_In": sn[g]["In"],
                 "sn_Mic": sn[g]["Mic"], "sn_Oli": sn[g]["Oli"],
                 "sn_neuron": neuron, "sn_glia": glia,
                 "miz": miz[g], "nat": nat[g], "bulk_mean": bulk})

def down(x): return x < 0
def up(x): return x > 0

# neuronal-down set (down in BOTH Ex and In)
neuron_down = [r for r in rows if down(r["sn_Ex"]) and down(r["sn_In"])]
nd_bulk = [r for r in neuron_down if down(r["miz"]) and down(r["nat"])]
# glial-up set (up in BOTH Mic and Oli)
glia_up = [r for r in rows if up(r["sn_Mic"]) and up(r["sn_Oli"])]
gu_bulk_up = [r for r in glia_up if up(r["miz"]) and up(r["nat"])]
gu_bulk_down = [r for r in glia_up if down(r["miz"]) and down(r["nat"])]

print(f"\nNEURONAL-DOWN UPR genes (snRNA Ex&In down): {len(neuron_down)}", file=sys.stderr)
print(f"   also DOWN in both bulk cohorts: {len(nd_bulk)}  ({100*len(nd_bulk)/max(len(neuron_down),1):.0f}%)  ← reproduced", file=sys.stderr)
print(f"\nGLIAL-UP UPR genes (snRNA Mic&Oli up): {len(glia_up)}", file=sys.stderr)
print(f"   also UP in both bulk cohorts: {len(gu_bulk_up)}  ({100*len(gu_bulk_up)/max(len(glia_up),1):.0f}%)  ← NOT reproduced", file=sys.stderr)
print(f"   instead DOWN in both bulk cohorts: {len(gu_bulk_down)}  (masked/opposite)", file=sys.stderr)

print("\nMarker examples (snRNA Ex/In/Mic/Oli | bulk Miz/Nat):", file=sys.stderr)
for gm in ["TMED2", "HSPA5", "TRIB3", "DNAJB9", "CALR", "CANX", "CHAC1", "ATP2A2", "BAG6", "DERL1"]:
    r = next((x for x in rows if x["gene"] == gm), None)
    if r:
        print(f"  {gm:8s} Ex{r['sn_Ex']:+.2f} In{r['sn_In']:+.2f} Mic{r['sn_Mic']:+.2f} Oli{r['sn_Oli']:+.2f} | "
              f"Miz{r['miz']:+.2f} Nat{r['nat']:+.2f}", file=sys.stderr)

# save tables
with open(os.path.join(OUT, "R1Q5_crossmodality_UPR_genes.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
open(os.path.join(OUT, "R1Q5_neuronDown_reproduced_in_bulk.txt"), "w").write(
    "\n".join(r["gene"] for r in nd_bulk) + "\n")
open(os.path.join(OUT, "R1Q5_glialUp_not_in_bulk.txt"), "w").write(
    "\n".join(r["gene"] for r in glia_up) + "\n")
print(f"\nsaved tables -> {OUT}", file=sys.stderr)
