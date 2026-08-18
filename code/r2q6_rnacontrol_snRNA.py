#!/usr/bin/env python3
"""R2.6 — the bulk disease-flat RNA-metabolism internal controls, run in Mathys snRNA per cell type.
Shows whether sets that are flat in bulk stay flat in snRNA (they follow the global cell-type shift).
Donor-level Cohen's d (late V-VI vs low I-II) + Wilcoxon p, all 6 cell types."""
import csv, os, glob, math
import numpy as np
from scipy.stats import mannwhitneyu

GMT = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
       "Human AD brain SEQ analysis/Alzheimer's brain disease Bulk RNA seq/"
       "2021 Mizuno_Human AD brain RNA seq_decreased PERK/Kyle_analysis.GseaPreranked.1651847739012/edb/gene_sets.gmt")
gm = {}
for line in open(GMT):
    p = line.rstrip("\n").split("\t"); gm[p[0]] = [g.strip() for g in p[2:] if g.strip()]
OUT = os.path.dirname(os.path.abspath(__file__))

# bulk disease-flat RNA-metabolism control terms
RNA_CTRL = ["GOBP_MRNA_TRANSPORT", "GOBP_RNA_DESTABILIZATION", "GOBP_RNA_EXPORT_FROM_NUCLEUS",
            "GOBP_NEGATIVE_REGULATION_OF_MRNA_METABOLIC_PROCESS", "GOBP_RNA_METHYLATION",
            "GOBP_POSITIVE_REGULATION_OF_MRNA_METABOLIC_PROCESS",
            "GOBP_NUCLEAR_TRANSCRIBED_MRNA_CATABOLIC_PROCESS"]
sets = [(t.replace("GOBP_", ""), gm[t]) for t in RNA_CTRL if t in gm]
# combined + UPR reference
sets = ([("UPR: ER-stress", [l.strip() for l in open(os.path.join(OUT, "ERstress_260_geneset.txt")) if l.strip()]),
         ("UPR: IRE1", [l.strip() for l in open(os.path.join(OUT, "geneset_IRE1.txt")) if l.strip()])]
        + sets +
        [("RNA-metab COMBINED", [l.strip() for l in open(os.path.join(OUT, "controlset_RNA_METABOLISM.txt")) if l.strip()])])

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

CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
CTZ = {}
for ct in CELLS:
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dc = hdr.index("donor"); cc = hdr.index("celltype")
    gc = [j for j in range(2, len(hdr)) if j not in (dc, cc)]
    genes = [hdr[j] for j in gc]; don = [r[dc] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gc] for r in tab[1:]], float)
    mu = np.nanmean(X, 0); sd = np.nanstd(X, 0, ddof=1); ok = (sd > 0) & np.isfinite(sd)
    Z = np.full_like(X, np.nan); Z[:, ok] = (X[:, ok]-mu[ok])/sd[ok]
    CTZ[ct] = (Z, {g: j for j, g in enumerate(genes)}, ok, don)
def score(gl, ct):
    Z, gi, ok, don = CTZ[ct]
    cols = [gi[g] for g in gl if g in gi and ok[gi[g]]]
    if len(cols) < 8: return float("nan"), float("nan"), len(cols)
    sc = np.nanmean(Z[:, cols], 1)
    a = np.array([sc[i] for i, d in enumerate(don) if grp(braak[d]) == "late"])
    b = np.array([sc[i] for i, d in enumerate(don) if grp(braak[d]) == "low"])
    a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    return cohend(a, b), mannwhitneyu(a, b, alternative="two-sided").pvalue, len(cols)

print(f"{'gene set':34s}{'n':>4}" + "".join(f"{c:>8}" for c in CELLS))
rows = []
for nm, gl in sets:
    ds = [score(gl, ct) for ct in CELLS]
    n = ds[0][2]
    print(f"{nm:34s}{n:>4}" + "".join(f"{d[0]:>8.2f}" for d in ds))
    row = {"gene_set": nm, "n_genes": n}
    for ct, d in zip(CELLS, ds):
        row[f"{ct}_d"] = round(d[0], 3); row[f"{ct}_p"] = round(d[1], 4)
    rows.append(row)
with open(os.path.join(OUT, "R2Q6_RNAcontrol_snRNA_bycell.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nsaved R2Q6_RNAcontrol_snRNA_bycell.csv  (Cohen's d, late V-VI vs low I-II)")
