#!/usr/bin/env python3
"""R2.6 — gene-set overlap matrix + per-dataset detection table (documents the text-vs-figure
gene-set-size differences the reviewer cited). Saves two CSVs for the Supplementary."""
import os, sys, csv
sys.path.insert(0, ".")
from r1q1_gomatrix import load_thap, load_mizuno, load_nativio
OUT = os.path.dirname(os.path.abspath(__file__))
S = [("PERK", "geneset_PERK.txt"), ("IRE1", "geneset_IRE1.txt"), ("ATF6", "geneset_ATF6.txt"),
     ("ERAD", "geneset_ERAD.txt"), ("ER-stress (260)", "ERstress_260_geneset.txt")]
sets = {n: set(l.strip() for l in open(os.path.join(OUT, f)) if l.strip()) for n, f in S}
names = [n for n, _ in S]

with open(os.path.join(OUT, "R2Q6_geneset_overlap_matrix.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["set(rows) ∩ set(cols)"] + names + [f"Jaccard%_{n}" for n in names])
    for a in names:
        inter = [len(sets[a] & sets[b]) for b in names]
        jac = [round(100 * len(sets[a] & sets[b]) / len(sets[a] | sets[b]), 1) for b in names]
        w.writerow([a] + inter + jac)

thap, miz, nat = load_thap(), load_mizuno(), load_nativio()
with open(os.path.join(OUT, "R2Q6_geneset_detection_by_dataset.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["gene_set", "curated_size", "detected_Thapsigargin", "detected_Mizuno", "detected_Nativio",
                "note"])
    for n in names:
        s = sets[n]
        w.writerow([n, len(s), sum(g in thap for g in s), sum(g in miz for g in s), sum(g in nat for g in s),
                    "smaller per-cohort numbers = genes not detected in that cohort's expressed universe"])
# shared marker genes the reviewer named
with open(os.path.join(OUT, "R2Q6_shared_ISR_genes.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["gene", "in_sets", "note"])
    for g, note in [("ASNS", "ATF4/ISR target"), ("TRIB3", "ATF4/ISR target"), ("DDIT3", "CHOP; ISR/UPR-shared"),
                    ("XBP1", "IRE1 output; ER-chaperone-shared"), ("ATF4", "ISR master TF"), ("HSPA5", "BiP; ER chaperone")]:
        w.writerow([g, ";".join(n for n in names if g in sets[n]), note])
print("saved R2Q6_geneset_overlap_matrix.csv, R2Q6_geneset_detection_by_dataset.csv, R2Q6_shared_ISR_genes.csv")
