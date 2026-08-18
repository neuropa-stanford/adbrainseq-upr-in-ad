#!/usr/bin/env python3
"""R1.6 — Braak-score vs gene-expression correlation, per cell type (Wenjun's analysis, recomputed).
Variable X = Braak score (braaksc, 1-6) per donor; Y = per-donor mean expression per cell type.
Pearson R + two-sided p + BH-adjusted p. Reports top +/- correlated UPR genes and TMED2/TRIB3."""
import csv, math, os, sys
import numpy as np
sys.path.insert(0, '.')
from r1q1_gomatrix import betainc

DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
        "Human AD brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

# ---- donor id (1..48, ordered unique Subject) -> braaksc ----
seen, braak = [], {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    s = r["Subject"]
    if s not in braak: seen.append(s); braak[s] = float(r["braaksc"])
id2braak = {i + 1: braak[s] for i, s in enumerate(seen)}  # fallback

# ---- UPR gene sets (for reporting universe) ----
GSDIR = OUT
upr = {}
upr["ER-stress"] = set(l.strip() for l in open(os.path.join(GSDIR, "ERstress_260_geneset.txt")) if l.strip())
for b in ["PERK", "IRE1", "ATF6", "ERAD"]:
    upr[b] = set(l.strip() for l in open(os.path.join(GSDIR, f"geneset_{b}.txt")) if l.strip())
UPR_ALL = set().union(*upr.values())

def gene_branch(g):
    return ",".join(b for b in ["PERK", "IRE1", "ATF6", "ERAD", "ER-stress"] if g in upr[b])

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(adj, 1.0); return out

def pear_p(R, n):
    if abs(R) >= 1: return 0.0
    df = n - 2; t2 = R * R * df / (1 - R * R)
    return betainc(df / 2.0, 0.5, df / (df + t2))

CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
all_rows = []
per_cell = {}
for ct in CELLS:
    rows = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = rows[0]
    # last two columns are celltype, donor
    donor_col = hdr.index("donor"); ct_col = hdr.index("celltype")
    gcols = [j for j in range(2, len(hdr)) if j not in (donor_col, ct_col)]
    genes = [hdr[j] for j in gcols]
    donors = [r[donor_col] for r in rows[1:]]
    x = np.array([braak[d] for d in donors])
    y = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gcols] for r in rows[1:]], float)
    # Pearson R vectorized over genes
    xc = x - x.mean()
    R = np.full(len(genes), np.nan)
    denomx = math.sqrt((xc ** 2).sum())
    for j in range(len(genes)):
        col = y[:, j]
        m = ~np.isnan(col)
        if m.sum() < 10 or np.nanstd(col) == 0: continue
        yc = col[m] - col[m].mean(); xx = x[m] - x[m].mean()
        d = math.sqrt((xx ** 2).sum()) * math.sqrt((yc ** 2).sum())
        if d == 0: continue
        R[j] = float((xx * yc).sum() / d)
    n = len(x)
    # restrict to UPR genes for the reported table + BH within UPR
    idxs = [j for j, g in enumerate(genes) if g in UPR_ALL and not np.isnan(R[j])]
    gg = [genes[j] for j in idxs]; RR = np.array([R[j] for j in idxs])
    pp = np.array([pear_p(r, n) for r in RR])
    qq = bh(pp)
    recs = []
    for g, r, p, q in zip(gg, RR, pp, qq):
        rec = {"cell_type": ct, "gene": g, "branch": gene_branch(g),
               "R": round(r, 3), "p": p, "adj_p_BH": round(q, 4), "n_donors": n}
        recs.append(rec); all_rows.append(rec)
    per_cell[ct] = recs
    print(f"\n=== {ct}: {len(recs)} UPR genes correlated (n={n} donors) ===", file=sys.stderr)
    recs_sorted = sorted(recs, key=lambda z: z["R"])
    print("  TOP NEGATIVE:", [(z["gene"], z["R"], f"q={z['adj_p_BH']}") for z in recs_sorted[:6]], file=sys.stderr)
    print("  TOP POSITIVE:", [(z["gene"], z["R"], f"q={z['adj_p_BH']}") for z in recs_sorted[-6:][::-1]], file=sys.stderr)
    for tgt in ["TMED2", "TRIB3"]:
        z = next((z for z in recs if z["gene"] == tgt), None)
        if z:
            rank_pos = sorted(recs, key=lambda a: -a["R"]).index(z) + 1
            rank_neg = sorted(recs, key=lambda a: a["R"]).index(z) + 1
            print(f"  >> {tgt}: R={z['R']} p={z['p']:.1e} q={z['adj_p_BH']}  "
                  f"(rank {rank_pos}/{len(recs)} most positive, {rank_neg}/{len(recs)} most negative)", file=sys.stderr)

# save full table
keys = ["cell_type", "gene", "branch", "R", "p", "adj_p_BH", "n_donors"]
with open(os.path.join(OUT, "R1Q6_Braak_correlation_UPRgenes.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(all_rows)
print(f"\nwrote R1Q6_Braak_correlation_UPRgenes.csv ({len(all_rows)} rows)", file=sys.stderr)
