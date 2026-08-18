#!/usr/bin/env python3
"""Find gene sets that move OPPOSITE to the global shift: UP in neurons (Ex,In) and DOWN in glia
(Mic,Oli) with Braak stage. Their existence shows the neuron-down/glia-up pattern is a directional
biological program, not a uniform technical artifact. Donor-level Cohen's d (late V-VI vs low I-II)."""
import csv, os, glob, math
import numpy as np

GMT = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
       "Human AD brain SEQ analysis/Alzheimer's brain disease Bulk RNA seq/"
       "2021 Mizuno_Human AD brain RNA seq_decreased PERK/Kyle_analysis.GseaPreranked.1651847739012/edb/gene_sets.gmt")
gm = {}
for line in open(GMT):
    p = line.rstrip("\n").split("\t"); gm[p[0]] = [g.strip() for g in p[2:] if g.strip()]

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
    if len(a) < 3 or len(b) < 3: return float("nan")
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
    if len(cols) < 10: return float("nan"), len(cols)
    sc = np.nanmean(Z[:, cols], 1)
    a = np.array([sc[i] for i, dn in enumerate(don) if grp(braak[dn]) == "late"])
    b = np.array([sc[i] for i, dn in enumerate(don) if grp(braak[dn]) == "low"])
    return cohend(a, b), len(cols)

rows = []
for t, gl in gm.items():
    if not (30 <= len(gl) <= 300): continue
    dex, n = score_d(gl, "Ex")
    if n < 10: continue
    din = score_d(gl, "In")[0]; dmic = score_d(gl, "Mic")[0]; doli = score_d(gl, "Oli")[0]
    if not all(np.isfinite(x) for x in (dex, din, dmic, doli)): continue
    reversal = (dex + din)/2 - (dmic + doli)/2      # neuron-up minus glia-down
    rows.append((t, n, dex, din, dmic, doli, reversal))

print(f"scanned {len(rows)} terms (30-300 genes)\n")
print("=== TOP neuron-UP / glia-DOWN gene sets (reversal score high; require Ex>0,In>0,Mic<0,Oli<0) ===")
print(f"{'term':58s}{'n':>4}{'Ex':>7}{'In':>7}{'Mic':>7}{'Oli':>7}{'rev':>7}")
strict = [r for r in rows if r[2] > 0 and r[3] > 0 and r[4] < 0 and r[5] < 0]
for r in sorted(strict, key=lambda x: -x[6])[:20]:
    print(f"{r[0].replace('GOBP_','').replace('HALLMARK_','H:')[:58]:58s}{r[1]:>4}{r[2]:>7.2f}{r[3]:>7.2f}{r[4]:>7.2f}{r[5]:>7.2f}{r[6]:>7.2f}")
print(f"\n(strict neuron-up & glia-down: {len(strict)} of {len(rows)} sets)")
print("\n=== TOP by reversal score regardless of strict signs ===")
for r in sorted(rows, key=lambda x: -x[6])[:12]:
    print(f"{r[0].replace('GOBP_','').replace('HALLMARK_','H:')[:58]:58s}{r[1]:>4}{r[2]:>7.2f}{r[3]:>7.2f}{r[4]:>7.2f}{r[5]:>7.2f}{r[6]:>7.2f}")
