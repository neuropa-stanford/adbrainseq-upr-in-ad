#!/usr/bin/env python3
"""Honest scan: do standard housekeeping/reference gene sets (ribosome, cell cycle, proteasome,
OxPhos, splicing, etc.) stay FLAT in the snRNA donor-level data, or do they move with the global
cell-type shift? Reports donor-level Cohen's d (late V-VI vs low I-II) per cell type for each."""
import csv, os, glob, math
import numpy as np

GMT = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
       "Human AD brain SEQ analysis/Alzheimer's brain disease Bulk RNA seq/"
       "2021 Mizuno_Human AD brain RNA seq_decreased PERK/Kyle_analysis.GseaPreranked.1651847739012/edb/gene_sets.gmt")
gm = {}
for line in open(GMT):
    p = line.rstrip("\n").split("\t"); gm[p[0]] = [g.strip() for g in p[2:] if g.strip()]

cats = {
 "Ribosome/translation": ["GOBP_CYTOPLASMIC_TRANSLATION", "GOBP_TRANSLATIONAL_INITIATION", "GOBP_RIBOSOME_BIOGENESIS"],
 "Cell cycle": ["GOBP_MITOTIC_CELL_CYCLE_PHASE_TRANSITION", "GOBP_DNA_REPLICATION", "GOBP_CHROMOSOME_SEGREGATION"],
 "Proteasome": ["GOBP_PROTEASOMAL_PROTEIN_CATABOLIC_PROCESS", "GOBP_PROTEASOME_ASSEMBLY"],
 "OxPhos/mito": ["GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_MITOCHONDRIAL_TRANSLATION", "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT"],
 "Splicing": ["GOBP_MRNA_SPLICING_VIA_SPLICEOSOME", "GOBP_RNA_SPLICING"],
 "Glycolysis/metab": ["GOBP_GLYCOLYTIC_PROCESS", "GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY"],
 "tRNA/aminoacyl": ["GOBP_TRNA_AMINOACYLATION_FOR_PROTEIN_TRANSLATION", "GOBP_TRNA_METABOLIC_PROCESS"],
 "DNA repair": ["GOBP_DNA_REPAIR", "GOBP_DOUBLE_STRAND_BREAK_REPAIR"],
 "Cytoskeleton": ["GOBP_ACTIN_FILAMENT_ORGANIZATION", "GOBP_MICROTUBULE_CYTOSKELETON_ORGANIZATION"],
 "Chromatin": ["GOBP_CHROMATIN_ORGANIZATION", "GOBP_HISTONE_MODIFICATION"],
}

DATA = glob.glob("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
                 "**/data_extraction", recursive=True)[0]
seen, info = [], {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    s = r["Subject"]
    if s not in info: seen.append(s); info[s] = r
braak = {s: int(float(info[s]["braaksc"])) for s in seen}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")
def cohend(a, b):
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    na, nb = len(a), len(b)
    sp = math.sqrt(((na-1)*np.var(a, ddof=1)+(nb-1)*np.var(b, ddof=1))/(na+nb-2))
    return (np.mean(a)-np.mean(b))/sp if sp > 0 else float("nan")

CT = {}
for ct in ["Ex", "In", "Mic", "Oli"]:
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dc = hdr.index("donor"); cc = hdr.index("celltype")
    gc = [j for j in range(2, len(hdr)) if j not in (dc, cc)]
    genes = [hdr[j] for j in gc]; don = [r[dc] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gc] for r in tab[1:]], float)
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0, ddof=1); ok = (sd > 0) & np.isfinite(sd)
    Z = np.full_like(X, np.nan); Z[:, ok] = (X[:, ok]-mu[ok])/sd[ok]
    CT[ct] = (Z, {g: j for j, g in enumerate(genes)}, ok, don)

def score_d(gl, ct):
    Z, gi, ok, don = CT[ct]
    cols = [gi[g] for g in gl if g in gi and ok[gi[g]]]
    if len(cols) < 8: return float("nan"), len(cols)
    sc = np.nanmean(Z[:, cols], 1)
    a = np.array([sc[i] for i, dn in enumerate(don) if grp(braak[dn]) == "late"])
    b = np.array([sc[i] for i, dn in enumerate(don) if grp(braak[dn]) == "low"])
    return cohend(a, b), len(cols)

print(f"{'category / term':54s}{'n':>4}{'Ex':>7}{'In':>7}{'Mic':>7}{'Oli':>7}")
print("--- reference: UPR sets + current control ---")
for nm, f in [("ER-stress", "ERstress_260_geneset.txt"), ("IRE1", "geneset_IRE1.txt"),
              ("ctrl: mRNA-transport", "controlset_MRNA_TRANSPORT.txt")]:
    gl = [l.strip() for l in open(f) if l.strip()]
    print(f"{nm:54s}{len(gl):>4}" + "".join(f"{score_d(gl, ct)[0]:>7.2f}" for ct in ["Ex", "In", "Mic", "Oli"]))
print("--- housekeeping candidates (sorted within category) ---")
rows = []
for cat, terms in cats.items():
    for t in terms:
        if t not in gm: continue
        gl = gm[t]
        if not (30 <= len(gl) <= 300): continue
        ds = [score_d(gl, ct) for ct in ["Ex", "In", "Mic", "Oli"]]
        rows.append((cat, t, ds))
        print(f"{(cat + ': ' + t.replace('GOBP_', ''))[:54]:54s}{ds[0][1]:>4}" + "".join(f"{x[0]:>7.2f}" for x in ds))
# flattest-in-Oli summary
print("\n--- flattest in OLIGODENDROCYTE (|d_Oli| ascending; want small while IRE1=+1.45) ---")
for cat, t, ds in sorted(rows, key=lambda r: abs(r[2][3][0]) if np.isfinite(r[2][3][0]) else 9)[:8]:
    print(f"  |Oli d|={abs(ds[3][0]):.2f}  {cat}: {t.replace('GOBP_','')}")
