#!/usr/bin/env python3
"""R2.9 — independent-cohort replication of the neuron-down / glia-up UPR-associated direction.
Two independent human AD snRNA cohorts, from their PUBLISHED per-cell-type differential-expression tables
(no raw download; all local):
  Lau 2020 (PNAS)     : AD_vs_NC log2FC per cell type (sd03).
  Mathys 2024 (Nature): pseudobulk DE vs neurofibrillary-tangle burden (path=nft), logFC per cell subtype
                        (Suppl Table 9), aggregated to broad classes.
Metric = mean signed log2FC of the ER-stress(260) UPR set (and 4 branches) per cell type; a matched
non-UPR control set (mRNA transport, 91 genes) is included to check specificity. Sign test (genes with
consistent direction) is descriptive replication of DIRECTION only (these are DE summaries, not donor-level).
"""
import os, glob, gzip, csv, math, re, zipfile, statistics as st
from xml.etree import ElementTree as ET

# ---- minimal xlsx reader (openpyxl chokes on this file's autofilter) ----
_M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
def xlsx_sheets(path):
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(_M + "si"):
            ss.append("".join(t.text or "" for t in si.iter(_M + "t")))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    sheets = {}
    for s in ET.fromstring(z.read("xl/workbook.xml")).iter(_M + "sheet"):
        tgt = rels[s.get(_R + "id")]
        sheets[s.get("name")] = ("xl/" + tgt) if not tgt.startswith("xl/") else tgt
    def read(name):
        out = []
        for row in ET.fromstring(z.read(sheets[name])).iter(_M + "row"):
            cells = {}
            for c in row.findall(_M + "c"):
                v = c.find(_M + "v")
                if v is None or v.text is None: continue
                col = re.match(r"[A-Z]+", c.get("r")).group()
                cells[col] = ss[int(v.text)] if c.get("t") == "s" else v.text
            out.append(cells)
        return out
    return sheets, read

OUT = os.path.dirname(os.path.abspath(__file__))
def rd(f): return set(l.strip() for l in open(os.path.join(OUT, f)) if l.strip())
UPR = rd("ERstress_260_geneset.txt")
BR = {"PERK": rd("geneset_PERK.txt"), "IRE1": rd("geneset_IRE1.txt"),
      "ATF6": rd("geneset_ATF6.txt"), "ERAD": rd("geneset_ERAD.txt")}
CTRL = rd("controlset_MRNA_TRANSPORT.txt")

def summ(vals):
    vals = [v for v in vals if v is not None and math.isfinite(v)]
    if len(vals) < 5: return None
    n = len(vals); m = st.mean(vals); down = sum(v < 0 for v in vals)
    return dict(n=n, mean=m, pct_down=100 * down / n)

# ================= Lau 2020 =================
LAU = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
       "Human AD brain SEQ analysis/Single cell RNA seq/2020 Lau/pnas.2008762117.sd03_with PERK gene set.xlsx")
LAU_CT = [("Excit", "Ex (neuron)"), ("Inhit", "In (neuron)"), ("Astro", "Ast (glia)"),
          ("Mic", "Mic (glia)"), ("Oligo", "Oli (glia)"), ("Endo", "Endo (vasc)")]
_, lau_read = xlsx_sheets(LAU)
lau = {}
for sh, lab in LAU_CT:
    d = {}
    for row in lau_read(sh)[1:]:              # skip header
        g = row.get("A")
        try: d[str(g).strip()] = float(row.get("B"))
        except (TypeError, ValueError): pass
    lau[lab] = d

# ================= Mathys 2024 =================
M24 = ("/tmp/"
       "adbrainseq_work/mathys24/dereg")
PATH, REGION = "nft", "allregions"        # tangles, pooled regions
CLASS = {"Ast": "Ast (glia)", "Exc": "Ex (neuron)", "Inh": "In (neuron)", "Mic": "Mic (glia)",
         "Oli": "Oli (glia)", "Opc": "Opc (glia)"}
m24 = {v: {} for v in CLASS.values()}     # class -> gene -> [logFC across subtypes]
for fp in glob.glob(os.path.join(M24, "aggregated_fullset.*.tsv.gz")):
    pre = os.path.basename(fp).split(".")[1].split("_")[0]
    if pre not in CLASS: continue
    lab = CLASS[pre]
    with gzip.open(fp, "rt") as fh:
        rd2 = csv.DictReader(fh, delimiter="\t")
        for row in rd2:
            if row.get("path") != PATH or row.get("region") != REGION: continue
            try: fc = float(row["logFC_nb"])
            except (TypeError, ValueError, KeyError): continue
            m24[lab].setdefault(row["gene"], []).append(fc)
m24g = {lab: {g: st.mean(v) for g, v in d.items()} for lab, d in m24.items()}   # mean across subtypes

# ================= report =================
def block(title, data, order):
    print(f"\n===== {title} =====")
    print(f"{'cell type':14s}{'set':16s}{'n':>4}{'meanLog2FC':>12}{'%down':>8}")
    for lab in order:
        d = data[lab]
        for sname, s in [("ER-stress(260)", UPR)] + list(BR.items()) + [("CTRL:mRNA-transport", CTRL)]:
            r = summ([d[g] for g in s if g in d])
            if r: print(f"{lab:14s}{sname:16s}{r['n']:>4}{r['mean']:>+12.3f}{r['pct_down']:>7.0f}%")
        print()

LAU_ORDER = [l for _, l in LAU_CT]
M24_ORDER = ["Ex (neuron)", "In (neuron)", "Ast (glia)", "Mic (glia)", "Oli (glia)", "Opc (glia)"]
block("Lau 2020 (AD vs NC, log2FC)", lau, LAU_ORDER)
block(f"Mathys 2024 (vs NFT tangles, {REGION}, logFC)", m24g, M24_ORDER)

# ---- verdict: does neuron-down / glia-up reproduce (ER-stress set mean sign)? ----
def sign_pattern(data):
    def m(lab):
        r = summ([data[lab][g] for g in UPR if g in data[lab]]); return r["mean"] if r else float("nan")
    neur = [m(l) for l in data if "neuron" in l]; glia = [m(l) for l in data if "glia" in l]
    return neur, glia
print("\n===== VERDICT (ER-stress set mean) =====")
for name, data in [("Lau 2020", lau), ("Mathys 2024", m24g)]:
    neur, glia = sign_pattern(data)
    print(f"{name:14s} neurons mean={st.mean(neur):+.3f} (down={all(v<0 for v in neur)})  "
          f"glia mean={st.mean(glia):+.3f} (up={all(v>0 for v in glia)})")

# save CSV
with open(os.path.join(OUT, "R2Q9_independent_cohort_replication.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["cohort", "cell_type", "gene_set", "n_genes", "mean_log2FC", "pct_down"])
    for cohort, data, order in [("Lau2020", lau, LAU_ORDER), ("Mathys2024_nft", m24g, M24_ORDER)]:
        for lab in order:
            for sname, s in [("ER-stress(260)", UPR)] + list(BR.items()) + [("CTRL_mRNA_transport", CTRL)]:
                r = summ([data[lab][g] for g in s if g in data[lab]])
                if r: w.writerow([cohort, lab, sname, r["n"], round(r["mean"], 4), round(r["pct_down"], 1)])
print("\nsaved R2Q9_independent_cohort_replication.csv")
