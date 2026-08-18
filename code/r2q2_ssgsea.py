#!/usr/bin/env python3
"""R2.2 second method — GSVA-type single-sample enrichment (ssGSEA, Barbie et al. 2009) on donor
pseudobulk, donor = unit of inference. Convergent check against the module-score result.
Per cell type: ssGSEA score per donor per UPR gene set (+ internal control), then compare
Braak V-VI vs I-II with Cohen's d, 95% CI, Wilcoxon, BH-adjusted p; also correlate with tangles."""
import csv, os, sys, math
import numpy as np
from scipy.stats import mannwhitneyu, t as tdist, pearsonr
import glob

DATA = glob.glob("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
                 "**/data_extraction", recursive=True)[0]
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

# ---- gene sets ----
SETS = [("ER-stress (260)", "ERstress_260_geneset.txt"), ("PERK (31)", "geneset_PERK.txt"),
        ("IRE1 (32)", "geneset_IRE1.txt"), ("ATF6 (74)", "geneset_ATF6.txt"),
        ("ERAD (75)", "geneset_ERAD.txt"),
        ("Internal control: mRNA transport (91)", "controlset_MRNA_TRANSPORT.txt")]
gsets = {n: [l.strip() for l in open(os.path.join(OUT, f)) if l.strip()] for n, f in SETS}
CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "Opc"]
LAB = {"Ex": "Excitatory neurons", "In": "Inhibitory neurons", "Ast": "Astrocytes",
       "Mic": "Microglia", "Oli": "Oligodendrocytes", "Opc": "OPCs"}

def ssgsea(X, gene_names, glist, alpha=0.25):
    """X: donors x genes. Returns per-donor ssGSEA enrichment score for one gene set."""
    ndon, ngene = X.shape
    gidx = {g: j for j, g in enumerate(gene_names)}
    idx = [gidx[g] for g in glist if g in gidx]
    if len(idx) < 5: return None, len(idx)
    inset = np.zeros(ngene, bool); inset[idx] = True
    out = np.full(ndon, np.nan)
    for d in range(ndon):
        expr = X[d]
        order = np.argsort(expr, kind="mergesort")          # ascending
        ranks = np.empty(ngene); ranks[order] = np.arange(1, ngene + 1)
        desc = order[::-1]                                   # descending expression
        w = ranks[desc] ** alpha
        hit_mask = inset[desc]
        hit = np.where(hit_mask, w, 0.0)
        sh = hit.sum()
        if sh == 0: continue
        P_hit = np.cumsum(hit) / sh
        P_miss = np.cumsum(~hit_mask) / (ngene - hit_mask.sum())
        out[d] = np.sum(P_hit - P_miss)
    rng = np.nanmax(out) - np.nanmin(out)                    # range-normalize across donors
    if rng > 0: out = out / rng
    return out, len(idx)

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

rows, store = [], {}
for ct in CELLS:
    tab = list(csv.reader(open(os.path.join(DATA, f"mean_DGE_by_donor_cell_type_{ct}.csv"))))
    hdr = tab[0]; dcol = hdr.index("donor"); ccol = hdr.index("celltype")
    gcols = [j for j in range(2, len(hdr)) if j not in (dcol, ccol)]
    genes = [hdr[j] for j in gcols]; donors = [r[dcol] for r in tab[1:]]
    X = np.array([[float(r[j]) if r[j] not in ("", "NA") else 0.0 for j in gcols] for r in tab[1:]], float)
    for sname, glist in gsets.items():
        sc, nset = ssgsea(X, genes, glist)
        if sc is None: continue
        store[(ct, sname)] = (donors, sc)
        g = {k: np.array([sc[i] for i, d in enumerate(donors) if grp(braak[d]) == k and np.isfinite(sc[i])])
             for k in ("low", "int", "late")}
        for comp in ("late", "int"):
            a, b = g[comp], g["low"]
            if len(a) < 3 or len(b) < 3: continue
            p = mannwhitneyu(a, b, alternative="two-sided").pvalue
            lo, hi = ci_diff(a, b)
            rows.append({"cell_type": ct, "cell_label": LAB[ct], "gene_set": sname, "n_set": nset,
                         "comparison": f"{comp} vs low", "n_comp": len(a), "n_low": len(b),
                         "cohens_d": round(float(cohens_d(a, b)), 3),
                         "CI95_low": round(lo, 3), "CI95_high": round(hi, 3), "p_wilcoxon": p})
        tv = [(tangles[d], sc[i]) for i, d in enumerate(donors) if d in tangles and np.isfinite(sc[i])]
        if len(tv) > 10:
            r_, p_ = pearsonr([x[0] for x in tv], [x[1] for x in tv])
            rows.append({"cell_type": ct, "cell_label": LAB[ct], "gene_set": sname, "n_set": nset,
                         "comparison": "tangles correlation", "n_comp": len(tv), "n_low": "",
                         "cohens_d": round(float(r_), 3), "CI95_low": "", "CI95_high": "", "p_wilcoxon": float(p_)})

qs = bh([r["p_wilcoxon"] for r in rows])
for r, q in zip(rows, qs): r["p_adj_BH"] = round(float(q), 4)
with open(os.path.join(OUT, "R2Q2_ssGSEA_donorlevel.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---- summary: late vs low, ER-stress + IRE1 + control, vs module-score direction ----
print("=== ssGSEA donor-level: late (V-VI) vs low (I-II) ===", file=sys.stderr)
print(f"{'cell type':20s}{'gene set':40s}{'d':>7}{'p':>9}{'BH q':>8}", file=sys.stderr)
for r in rows:
    if r["comparison"] != "late vs low": continue
    print(f"{r['cell_label'][:18]:20s}{r['gene_set'][:38]:40s}{r['cohens_d']:>7.2f}"
          f"{r['p_wilcoxon']:>9.3f}{r['p_adj_BH']:>8.3f}", file=sys.stderr)
print(f"\nwrote R2Q2_ssGSEA_donorlevel.csv ({len(rows)} rows)", file=sys.stderr)
