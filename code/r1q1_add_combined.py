#!/usr/bin/env python3
"""Add the Combined(Mizuno+Nativio) dataset to the 18-term GO matrix (Jonathan #2/#4).

Combined = genes concordant in DIRECTION across both cohorts.
 - no cutoff: same sign in Mizuno and Nativio (split UP/DOWN by that shared sign).
 - p<0.05  : same sign AND p<0.05 in BOTH cohorts.
Matches the definition used in JONATHAN_R1Q1_results_and_answer.md.
Appends rows to R1Q1_UPR_GO_matrix_18terms_up_down_both.csv with dataset='Combined(M+N)'.
"""
import csv, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1q1_gomatrix import load_mizuno, load_nativio, gost, GO_TERMS  # noqa

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
CSV = os.path.join(OUT, "R1Q1_UPR_GO_matrix_18terms_up_down_both.csv")

miz = load_mizuno()      # sym -> (lfc, p)
nat = load_nativio()
common = set(miz) & set(nat)
print(f"common genes: {len(common)}", file=sys.stderr)

def concordant(cut):
    up, dn = [], []
    for g in common:
        lm, pm = miz[g]; ln, pn = nat[g]
        if lm == 0 or ln == 0 or (lm > 0) != (ln > 0):
            continue
        if cut is not None and not (pm < cut and pn < cut):
            continue
        (up if lm > 0 else dn).append(g)
    return up, dn

rows = []
for thr, cut in [("no cutoff (sign only)", None), ("p<0.05", 0.05)]:
    up, dn = concordant(cut)
    both = up + dn
    print(f"Combined {thr}: UP {len(up)}  DOWN {len(dn)}  BOTH {len(both)}", file=sys.stderr)
    for direction, genes in [("UP", up), ("DOWN", dn), ("BOTH", both)]:
        tag = f"Combined__{'nocut' if cut is None else 'p05'}__{direction}"
        res = gost(genes, tag) if genes else []
        by = {r["native"]: r for r in res}
        for go, name, panel in GO_TERMS:
            r = by.get(go)
            rows.append({
                "threshold": thr, "dataset": "Combined(M∩N)", "direction": direction,
                "n_query_genes": len(genes), "GO": go, "term": name,
                "figure_panel": panel,
                "adjP": r["p_value"] if r else "",
                "neglog10_adjP": (-math.log10(r["p_value"]) if r and r["p_value"] > 0 else
                                  (330 if r else "")),
                "genes_in_term": r["intersection_size"] if r else "",
                "term_size": r["term_size"] if r else "",
                "significant": (r["p_value"] < 0.05) if r else False,
            })

existing = list(csv.DictReader(open(CSV)))
existing = [r for r in existing if r["dataset"] != "Combined(M∩N)"]
fields = list(existing[0].keys())
with open(CSV, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields)
    w.writeheader(); w.writerows(existing)
    for r in rows:
        w.writerow({k: r.get(k, "") for k in fields})
print(f"appended {len(rows)} Combined rows -> {CSV}", file=sys.stderr)
