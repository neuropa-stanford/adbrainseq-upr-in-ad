#!/usr/bin/env python3
"""Find a BULK internal-control GO gene set (40-100 genes) that does NOT change in either
Mizuno OR Nativio bulk RNA-seq — i.e., log2FC distribution centered at 0, non-significant in both."""
import sys, os, statistics
import numpy as np
sys.path.insert(0, '.')
from r1q1_gomatrix import load_mizuno, load_nativio
from scipy.stats import wilcoxon, ttest_1samp

GMT = ("/data/adbrainseq/Stanford U/PERK Human AD brain analysis/"
       "Human AD brain SEQ analysis/Alzheimer's brain disease Bulk RNA seq/"
       "2021 Mizuno_Human AD brain RNA seq_decreased PERK/Kyle_analysis.GseaPreranked.1651847739012/edb/gene_sets.gmt")
OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/Major Revision/processed raw data")

miz0 = {g: v[0] for g, v in load_mizuno().items()}
nat0 = {g: v[0] for g, v in load_nativio().items()}
import statistics as _st
_mm = _st.median(miz0.values()); _nm = _st.median(nat0.values())
miz = {g: v - _mm for g, v in miz0.items()}   # centered to genome median
nat = {g: v - _nm for g, v in nat0.items()}
gmed_m = statistics.median(miz.values()); gmed_n = statistics.median(nat.values())
print(f"global median log2FC: Mizuno {gmed_m:+.3f}  Nativio {gmed_n:+.3f}", file=sys.stderr)

# exclude UPR/proteostasis-related terms (a control must be unrelated)
EXCL = ["UNFOLDED", "ENDOPLASMIC", "ER_", "_ER", "ERAD", "PROTEIN_FOLDING", "RETICULUM",
        "TOPOLOGICALLY", "CHAPERONE", "PROTEOSTASIS", "STRESS", "HEAT"]

def flat_stats(vals):
    v = np.array(vals)
    med = float(np.median(v)); mean = float(np.mean(v))
    try:
        pw = wilcoxon(v, zero_method="wilcox").pvalue if len(v) >= 10 and np.any(v != 0) else 1.0
    except Exception:
        pw = 1.0
    pt = float(ttest_1samp(v, 0.0).pvalue) if len(v) >= 3 else 1.0
    iqr = float(np.percentile(v, 75) - np.percentile(v, 25))
    return med, mean, pw, pt, iqr

rows = []
for line in open(GMT):
    parts = line.rstrip("\n").split("\t")
    name = parts[0]
    if not name.startswith("GOBP_"):
        continue
    if any(e in name for e in EXCL):
        continue
    genes = parts[2:]
    mv = [miz[g] for g in genes if g in miz]
    nv = [nat[g] for g in genes if g in nat]
    if not (40 <= len(mv) <= 100 and 40 <= len(nv) <= 100):
        continue
    mm, mmean, mpw, mpt, miqr = flat_stats(mv)
    nm, nmean, npw, npt, niqr = flat_stats(nv)
    # flatness: small |median| relative to 0 in BOTH + non-significant one-sample in BOTH
    worst_med = max(abs(mm), abs(nm))
    worst_p = min(mpw, npw, mpt, npt)      # want this HIGH (non-sig)
    rows.append({"term": name, "n_miz": len(mv), "n_nat": len(nv),
                 "med_Miz": round(mm, 3), "p_Miz": mpw, "med_Nat": round(nm, 3), "p_Nat": npw,
                 "worst_abs_median": round(worst_med, 3), "min_p": round(worst_p, 4),
                 "iqr_Miz": round(miqr, 2), "iqr_Nat": round(niqr, 2)})

# candidates flat in BOTH: |median|<0.08 both, non-sig (p>0.2) both
# rank by flatness = small |median| in both + high min one-sample p (non-sig)
for r in rows:
    r["flatscore"] = max(abs(r["med_Miz"]), abs(r["med_Nat"])) + 0.5*(abs(r["med_Miz"])+abs(r["med_Nat"]))
flat = sorted(rows, key=lambda r: r["flatscore"])
strict = [r for r in rows if abs(r["med_Miz"])<0.05 and abs(r["med_Nat"])<0.05 and r["p_Miz"]>0.3 and r["p_Nat"]>0.3]
print(f"\n{len(rows)} terms; {len(strict)} very-flat (|median|<0.05 & p>0.3 both). Top 25 by flatness:\n", file=sys.stderr)
print(f"{'term':60s}{'nM/nN':>8}{'medMiz':>8}{'pMiz':>7}{'medNat':>8}{'pNat':>7}", file=sys.stderr)
for r in flat[:25]:
    print(f"{r['term'][:58]:60s}{r['n_miz']:>4}/{r['n_nat']:<3}{r['med_Miz']:>8.3f}{r['p_Miz']:>7.2f}"
          f"{r['med_Nat']:>8.3f}{r['p_Nat']:>7.2f}", file=sys.stderr)

import csv
with open(os.path.join(OUT, "bulk_internal_control_candidates.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    w.writerows(sorted(rows, key=lambda r: r["worst_abs_median"]))
print(f"\nsaved bulk_internal_control_candidates.csv ({len(rows)} terms)", file=sys.stderr)
