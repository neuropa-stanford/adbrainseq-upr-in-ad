#!/usr/bin/env python3
"""R2.2 for BULK — donor-level pathway scores in the Nativio cohort (per-sample = per-donor).
SuppD3 has per-sample raw counts (12 AD, 10 old). CPM-normalize -> log2 -> per-donor module score
per UPR gene set (+ internal control) -> AD vs old with donor as the unit (Cohen's d, 95% CI, Wilcoxon).
Mizuno (SuppD2) holds only summary DE (log2FC/p), so per-donor Mizuno needs the GSE173955 count matrix."""
import os, sys, math, csv
import numpy as np
import openpyxl
from scipy.stats import mannwhitneyu, t as tdist

BASE = ("/data/adbrainseq/Publication/2024_scRNA seq/"
        "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication")
OUT = os.path.join(BASE, "Major Revision", "processed raw data")

# ---- Nativio per-sample counts (SuppD3): sym=col[1], AD=cols[6:18] (12), old=cols[18:28] (10) ----
wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD3_NativioSeq.xlsx"), read_only=True); ws = wb.active
genes, AD, OLD = [], [], []
for r in ws.iter_rows(min_row=3, values_only=True):
    sym = r[1]
    if not sym or str(sym).strip() in ("", "None"): continue
    try:
        ad = [float(x) for x in r[6:18]]; old = [float(x) for x in r[18:28]]
    except (TypeError, ValueError): continue
    if len(ad) != 12 or len(old) != 10: continue
    genes.append(str(sym).strip()); AD.append(ad); OLD.append(old)
wb.close()
AD = np.array(AD, float); OLD = np.array(OLD, float)                  # genes x samples
counts = np.hstack([AD, OLD])                                          # genes x 22
labels = ["AD"] * 12 + ["old"] * 10
print(f"Nativio genes {len(genes)}, samples {counts.shape[1]} (12 AD / 10 old)", file=sys.stderr)

# ---- CPM normalize per sample, log2 ----
lib = counts.sum(axis=0)
cpm = counts / lib * 1e6
logcpm = np.log2(cpm + 1.0)                                           # genes x samples
gi = {g: i for i, g in enumerate(genes)}

SETS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
        ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"),
        ("ERAD (75)", "geneset_ERAD.txt"),
        ("Internal control: mRNA transport (91)", "controlset_MRNA_TRANSPORT.txt")]
gsets = {n: [l.strip() for l in open(os.path.join(OUT, f)) if l.strip()] for n, f in SETS}

def cohens_d(a, b):
    na, nb = len(a), len(b)
    sp = math.sqrt(((na-1)*np.var(a, ddof=1)+(nb-1)*np.var(b, ddof=1))/(na+nb-2))
    return (np.mean(a)-np.mean(b))/sp if sp > 0 else float("nan")
def ci_diff(a, b, conf=0.95):
    na, nb = len(a), len(b)
    sp2 = ((na-1)*np.var(a, ddof=1)+(nb-1)*np.var(b, ddof=1))/(na+nb-2)
    se = math.sqrt(sp2*(1/na+1/nb)); tc = tdist.ppf(1-(1-conf)/2, na+nb-2)
    d = np.mean(a)-np.mean(b); return d-tc*se, d+tc*se
def bh(ps):
    ps = np.asarray(ps, float); n = len(ps); o = np.argsort(ps)
    adj = np.minimum.accumulate((ps[o]*n/(np.arange(n)+1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(adj, 1.0); return out

# z-score each gene across the 22 samples, module score = mean z of set genes per sample
mu = logcpm.mean(axis=1, keepdims=True); sd = logcpm.std(axis=1, ddof=1, keepdims=True)
Z = np.where(sd > 0, (logcpm - mu) / np.where(sd > 0, sd, 1), np.nan)   # genes x samples

rows = []
for sname, gl in gsets.items():
    idx = [gi[g] for g in gl if g in gi]
    score = np.nanmean(Z[idx, :], axis=0)                              # per-sample module score
    a = score[:12]; b = score[12:]                                     # AD vs old
    p = mannwhitneyu(a, b, alternative="two-sided").pvalue
    lo, hi = ci_diff(a, b)
    rows.append({"cohort": "Nativio (bulk)", "gene_set": sname, "n_set": len(idx),
                 "comparison": "AD vs old", "n_AD": 12, "n_old": 10,
                 "mean_AD": round(float(a.mean()), 3), "mean_old": round(float(b.mean()), 3),
                 "cohens_d": round(float(cohens_d(a, b)), 3),
                 "CI95_low": round(lo, 3), "CI95_high": round(hi, 3), "p_wilcoxon": p})
qs = bh([r["p_wilcoxon"] for r in rows])
for r, q in zip(rows, qs): r["p_adj_BH"] = round(float(q), 4)

with open(os.path.join(OUT, "R2Q2_bulk_Nativio_donorlevel.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

print("=== BULK Nativio donor-level (per-sample module score), AD (n=12) vs old (n=10) ===", file=sys.stderr)
print(f"{'gene set':40s}{'d':>7}{'95% CI':>18}{'p':>8}{'BH q':>8}", file=sys.stderr)
for r in rows:
    print(f"{r['gene_set'][:38]:40s}{r['cohens_d']:>7.2f}  [{r['CI95_low']:+.2f},{r['CI95_high']:+.2f}]"
          f"{r['p_wilcoxon']:>8.3f}{r['p_adj_BH']:>8.3f}", file=sys.stderr)
print(f"\nwrote R2Q2_bulk_Nativio_donorlevel.csv ({len(rows)} rows)", file=sys.stderr)
