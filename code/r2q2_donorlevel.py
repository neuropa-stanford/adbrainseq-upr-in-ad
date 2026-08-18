#!/usr/bin/env python3
"""R2.2 — DONOR-LEVEL statistics (the pivotal fix).
Per-donor UPR module score per cell type; donor as the unit of inference.
Module score = mean of gene-wise z-scores (z across the 48 donors within a cell type).
Groups: low Braak I-II (reference) vs III-IV vs V-VI.
Reports n donors, mean difference, Cohen's d, 95% CI, Wilcoxon rank-sum p, BH-adjusted p.
Also correlates module score with quantitative tangle score (R2.5)."""
import csv, os, sys, math
import numpy as np
from scipy.stats import mannwhitneyu, t as tdist, pearsonr

DATA = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
        "Human AD brain SEQ analysis/Single cell RNA seq/2019 Mathys/Wenjun's Braak Data Extraction/data_extraction")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

# ---- donors ----
seen, info = [], {}
for r in csv.DictReader(open(os.path.join(DATA, "donor_info.csv"))):
    s = r["Subject"]
    if s not in info: seen.append(s); info[s] = r
braak = {s: int(float(info[s]["braaksc"])) for s in seen}
tangles = {s: float(info[s]["tangles"]) for s in seen if info[s]["tangles"] not in ("", "NA")}
def grp(b): return "low" if b <= 2 else ("int" if b <= 4 else "late")

# ---- gene sets (UPR + internal control) ----
SETS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
        ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"),
        ("ERAD (75)", "geneset_ERAD.txt"),
        ("Internal control: mRNA transport (91)", "controlset_MRNA_TRANSPORT.txt")]
gsets = {n: [l.strip() for l in open(os.path.join(OUT, f)) if l.strip()] for n, f in SETS}

CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
LAB = {"Ex": "Excitatory neurons", "In": "Inhibitory neurons", "Ast": "Astrocytes",
       "Mic": "Microglia", "Oli": "Oligodendrocytes", "Opc": "OPCs"}

def cohens_d(a, b):
    na, nb = len(a), len(b)
    sp = math.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sp if sp > 0 else float("nan")

def ci_diff(a, b, conf=0.95):
    na, nb = len(a), len(b)
    sp2 = ((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2)
    se = math.sqrt(sp2 * (1 / na + 1 / nb)); df = na + nb - 2
    tc = tdist.ppf(1 - (1 - conf) / 2, df); d = np.mean(a) - np.mean(b)
    return d - tc * se, d + tc * se

def bh(ps):
    ps = np.asarray(ps, float); n = len(ps); o = np.argsort(ps)
    adj = np.minimum.accumulate((ps[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(adj, 1.0); return out

rows, scores_store = [], {}
for ct in CELLS:
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dcol = hdr.index("donor"); ccol = hdr.index("celltype")
    gcols = [j for j in range(2, len(hdr)) if j not in (dcol, ccol)]
    genes = [hdr[j] for j in gcols]
    donors = [r[dcol] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else np.nan for j in gcols] for r in tab[1:]], float)
    # z-score each gene across donors
    mu = np.nanmean(X, axis=0); sd = np.nanstd(X, axis=0, ddof=1)
    ok = (sd > 0) & np.isfinite(sd)
    Z = np.full_like(X, np.nan)
    Z[:, ok] = (X[:, ok] - mu[ok]) / sd[ok]
    gidx = {g: j for j, g in enumerate(genes)}
    for sname, glist in gsets.items():
        cols = [gidx[g] for g in glist if g in gidx and ok[gidx[g]]]
        if len(cols) < 10: continue
        score = np.nanmean(Z[:, cols], axis=1)          # per-donor module score
        scores_store[(ct, sname)] = (donors, score)
        gr = {g: np.array([score[i] for i, d in enumerate(donors) if grp(braak[d]) == g]) for g in ("low", "int", "late")}
        for comp in ("late", "int"):
            a, b = gr[comp], gr["low"]
            if len(a) < 3 or len(b) < 3: continue
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
            lo, hi = ci_diff(a, b)
            rows.append({"cell_type": ct, "cell_label": LAB[ct], "gene_set": sname,
                         "comparison": f"{comp} vs low", "n_comp": len(a), "n_low": len(b),
                         "mean_diff": round(float(np.mean(a) - np.mean(b)), 3),
                         "cohens_d": round(float(cohens_d(a, b)), 3),
                         "CI95_low": round(lo, 3), "CI95_high": round(hi, 3),
                         "p_wilcoxon": p})
        # tangles correlation (R2.5)
        tv = [(tangles[d], score[i]) for i, d in enumerate(donors) if d in tangles and np.isfinite(score[i])]
        if len(tv) > 10:
            r_, p_ = pearsonr([x[0] for x in tv], [x[1] for x in tv])
            rows.append({"cell_type": ct, "cell_label": LAB[ct], "gene_set": sname,
                         "comparison": "tangles correlation", "n_comp": len(tv), "n_low": "",
                         "mean_diff": "", "cohens_d": round(float(r_), 3), "CI95_low": "", "CI95_high": "",
                         "p_wilcoxon": float(p_)})

qs = bh([r["p_wilcoxon"] for r in rows])
for r, q in zip(rows, qs): r["p_adj_BH"] = round(float(q), 4)
with open(os.path.join(OUT, "R2Q2_donorlevel_module_scores.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
np.save(os.path.join(OUT, "R2Q2_scores.npy"), np.array(
    [(f"{ct}|{sn}", ",".join(d), ",".join(f"{x:.4f}" for x in s)) for (ct, sn), (d, s) in scores_store.items()],
    dtype=object), allow_pickle=True)

print("=== DONOR-LEVEL: late (V-VI) vs low (I-II), module score ===", file=sys.stderr)
print(f"{'cell type':22s}{'gene set':40s}{'d':>7}{'95% CI':>18}{'p':>9}{'BH q':>8}", file=sys.stderr)
for r in rows:
    if r["comparison"] != "late vs low": continue
    print(f"{r['cell_label'][:20]:22s}{r['gene_set'][:38]:40s}{r['cohens_d']:>7.2f}"
          f"  [{r['CI95_low']:+.2f},{r['CI95_high']:+.2f}]{r['p_wilcoxon']:>9.3f}{r['p_adj_BH']:>8.3f}", file=sys.stderr)
print(f"\nwrote R2Q2_donorlevel_module_scores.csv ({len(rows)} rows)", file=sys.stderr)
