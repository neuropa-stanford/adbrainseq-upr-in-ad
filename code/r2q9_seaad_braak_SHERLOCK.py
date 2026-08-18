#!/usr/bin/env python3
"""SEA-AD (MTG) TAU-specific reanalysis -- RUN ON SHERLOCK (Lin lab allocation; scanpy env).

Goal: Braak V,VI  vs  Braak 0,I,II  per-cell-type log2FC of the UPR gene set -- the SAME tau contrast
as our manuscript Fig 4 -- in the fully independent SEA-AD cohort (Allen, MTG, 84 donors).

Why Sherlock: combined h5ad is ~30 GB; local python has a numpy2/h5py clash. Sherlock conda scanpy +
in-network S3 solve both. Data is PUBLIC (AWS open data, --no-sign-request; no DUC/approval needed).

Works for BOTH cohorts (Braak V,VI vs 0,I,II per cell type):
  SEA-AD  (MTG-only, public):   python thisscript.py SEAAD_MTG_...h5ad  ERstress_260_geneset.txt  out  SEAAD
  Mathys24 (all-regions, ROSMAP): python thisscript.py all_brain_regions_...h5ad ERstress_260_geneset.txt out MATHYS24 PFC
Args: <h5ad> <geneset.txt> <outdir> [tag] [region_filter]
  region_filter (optional): keep only cells whose region/brain-region obs column == this value (e.g. PFC).
Output (small CSVs to bring back):
  {tag}_Braak56vs012_log2FC_perCellType.csv   (gene x cell type x log2FC)
  {tag}_Braak56vs012_donor_pseudobulk.csv     (donor-level pseudobulk, for donor-level stats/plot)
"""
import sys, os
import numpy as np
import anndata as ad          # lighter than scanpy; read_h5ad is all we need

H5AD, GENESET, OUTDIR = sys.argv[1], sys.argv[2], sys.argv[3]
TAG = sys.argv[4] if len(sys.argv) > 4 else "COHORT"
REGION = sys.argv[5] if len(sys.argv) > 5 else None
os.makedirs(OUTDIR, exist_ok=True)
WANT = set(l.strip() for l in open(GENESET) if l.strip()) | {"EIF2AK3", "ERN1", "ATF6"}
SENS = {"EIF2AK3", "ERN1", "ATF6"}

print(">>> read h5ad (full load; sparse) ...", flush=True)
A = ad.read_h5ad(H5AD)          # bigmem node; sparse load. (backed .raw subsetting caused OOM)
print("shape:", A.shape)
print("obs columns:", list(A.obs.columns))

def pick(cols, *cands):
    for c in cands:
        for x in cols:
            if str(x).lower() == c.lower(): return x
    for c in cands:
        for x in cols:
            if c.lower() in str(x).lower(): return x
    return None
oc = A.obs.columns
COL_SUB = pick(oc, "Subclass", "subclass", "cell_type", "broad.cell.type", "cell_type_high_resolution", "major.celltype")
COL_DONOR = pick(oc, "Donor ID", "donor_id", "donor", "external_donor_name", "sample_id", "subject", "projid", "individualID")
COL_BRAAK = pick(oc, "Braak stage", "Braak", "braak", "braaksc")
COL_REGION = pick(oc, "region", "brain_region", "Region", "roi")
print(f">>> Subclass={COL_SUB} | Donor={COL_DONOR} | Braak={COL_BRAAK} | Region={COL_REGION}")
print(">>> unique Subclass:", sorted(map(str, A.obs[COL_SUB].unique()))[:60])
print(">>> unique Braak:", sorted(map(str, A.obs[COL_BRAAK].unique())))
if COL_REGION is not None:
    print(">>> unique Region:", sorted(map(str, A.obs[COL_REGION].unique())))
# optional region subset (e.g. Mathys 2024 all-regions -> keep PFC only)
if REGION is not None and COL_REGION is not None:
    keepr = (A.obs[COL_REGION].astype(str) == REGION).values
    print(f">>> region filter '{REGION}': keeping {keepr.sum()} / {len(keepr)} cells")
    A = A[keepr].copy() if not A.isbacked else A[keepr]

# ---- gene subset (memory-safe): subset the ~260 UPR genes from the in-memory matrix ----
vn = np.array(A.var_names)
gmask = np.isin(vn, list(WANT))
Xall = A.X
if gmask.sum() < 5 and A.raw is not None:      # genes may live in .raw
    print(">>> few genes in .X var; using .raw var")
    vn = np.array(A.raw.var_names); gmask = np.isin(vn, list(WANT)); Xall = A.raw.X
gidx = np.where(gmask)[0]; gnames = vn[gidx]
print(f">>> {len(gidx)} UPR/sensor genes found")
Xsub = Xall[:, gidx]                            # small: n_cells x ~260
X = np.asarray(Xsub.todense()) if hasattr(Xsub, "todense") else np.asarray(Xsub)
del Xall                                         # free the big matrix
# detect counts vs normalized
is_counts = np.allclose(X[:50], np.round(X[:50])) and X.max() > 30
print(f">>> matrix looks like {'RAW COUNTS' if is_counts else 'normalized/log'} (max={X.max():.2f})")

# ---- broad cell class ----
def broad(s):
    s = str(s)
    # Mathys-style short labels
    if s in ("Exc", "Ex", "Excitatory", "Excitatory neuron", "Excitatory neurons"): return "Ex"
    if s in ("Inh", "In", "Inhibitory", "Inhibitory neuron", "Inhibitory neurons"): return "In"
    if s in ("Ast", "Astro", "Astrocyte", "Astrocytes"): return "Ast"
    if s in ("Oli", "Oligo", "Oligodendrocyte", "Oligodendrocytes", "ODC", "OligD"): return "Oli"
    if s in ("Opc", "OPC"): return "OPC"
    if s in ("Mic", "MG", "Microglia", "Micro-PVM", "Immune"): return "Mic"
    # SEA-AD subclass labels
    if s in ("Lamp5", "Lamp5 Lhx6", "Sncg", "Vip", "Sst", "Sst Chodl", "Pvalb", "Chandelier", "Pax6"): return "In"
    if s.startswith("L") and any(k in s for k in ("IT", "ET", "CT", "NP", "6b")): return "Ex"
    if "Astro" in s: return "Ast"
    if "Oligo" in s: return "Oli"
    if "OPC" in s: return "OPC"
    if "Micro" in s or "PVM" in s: return "Mic"
    return "other"
bclass = np.array([broad(x) for x in A.obs[COL_SUB].values])
donor = A.obs[COL_DONOR].astype(str).values
def bnum(s):
    s = str(s).upper().replace("BRAAK", "").strip()
    return {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6,
            "1": 1, "2": 2, "3": 3, "4": 4, "5": 5, "6": 6}.get(s, np.nan)
bk = np.array([bnum(x) for x in A.obs[COL_BRAAK].values])
grp = np.where(np.isin(bk, [5, 6]), "high", np.where(np.isin(bk, [0, 1, 2]), "low", "mid"))

# ---- per donor x cell type pseudobulk (log-CPM if counts; else mean) ----
CELLS = ["Ex", "In", "Ast", "Mic", "Oli", "OPC"]
import csv, collections
donor_pb = {}   # (cell, donor) -> (grp, vector over gnames)
for c in CELLS:
    for d in np.unique(donor):
        m = (bclass == c) & (donor == d)
        if m.sum() < 10: continue
        g = grp[m][0]
        # keep ALL Braak groups (low=0,I,II ; mid=III,IV ; high=V,VI) for the 3-group Fig-4 contrasts
        sx = X[m]
        if is_counts:
            v = sx.sum(0); vec = np.log2(v / max(v.sum(), 1) * 1e4 + 1)
        else:
            vec = sx.mean(0)
        bn = bk[m][0]
        donor_pb[(c, d)] = (g, vec, int(bn) if bn == bn else -1)   # bn==bn filters NaN

with open(os.path.join(OUTDIR, f"{TAG}_Braak56vs012_donor_pseudobulk.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["cell_type", "donor", "braak_group", "braak_num"] + list(gnames))
    for (c, d), (g, vec, bn) in donor_pb.items(): w.writerow([c, d, g, bn] + [round(float(x), 4) for x in vec])

rows = []
for c in CELLS:
    hi = [vec for (cc, d), (g, vec, bn) in donor_pb.items() if cc == c and g == "high"]
    lo = [vec for (cc, d), (g, vec, bn) in donor_pb.items() if cc == c and g == "low"]
    if len(hi) < 2 or len(lo) < 2:
        print(f"  {c}: hi={len(hi)} lo={len(lo)} (skip)"); continue
    fc = np.vstack(hi).mean(0) - np.vstack(lo).mean(0)
    for gi, gn in enumerate(gnames):
        rows.append([c, gn, "UPR-sensor" if gn in SENS else "UPR-target", round(float(fc[gi]), 4), len(hi), len(lo)])
with open(os.path.join(OUTDIR, f"{TAG}_Braak56vs012_log2FC_perCellType.csv"), "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["cell_type", "gene", "gene_set", "log2FC_Braak56_vs_012", "n_high", "n_low"]); w.writerows(rows)

print("\n=== direction (UPR-target, log2FC Braak 5,6 / 0,I,II) ===")
agg = collections.defaultdict(list)
for c, gn, gs, fc, nh, nl in rows:
    if gs == "UPR-target": agg[c].append(fc)
for c in CELLS:
    v = agg.get(c, [])
    if v: print(f"  {c:4s} n={len(v)} mean={np.mean(v):+.3f} {100*sum(x<0 for x in v)//len(v)}% down  (donors hi/lo shown in CSV)")
print("\nDONE -> bring back", OUTDIR, "*.csv")
