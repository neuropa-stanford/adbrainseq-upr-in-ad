#!/usr/bin/env python3
"""R1.1 / R2.3 — preranked GSEA on the SAME 18 GO terms used for the ORA figures.

Rank metric: sign(log2FC) x -log10(P)  (Nativio: Welch t statistic, same sign convention).
Gene sets:   QuickGO human protein annotations, GO term + descendants (is_a/part_of/occurs_in).
Statistic:   classic weighted (p=1) Kolmogorov-Smirnov enrichment score, NES and P by
             gene-set permutation (size-matched random sets), BH-adjusted across terms.

Two rank-list variants, mirroring the ORA "cutoff" axis:
  full     — every detected gene (this is the analysis GSEA is designed for)
  p<0.05   — ranking restricted to genes passing p<0.05 (shown to demonstrate that the
             AD signal survives in GSEA where ORA lost it entirely)
Direction is read from the sign of NES (positive = up in AD / treated), so GSEA needs no
UP/DOWN split.
"""
import csv, json, math, os, sys, urllib.request
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r1q1_gomatrix import load_thap, load_mizuno, GO_TERMS, BASE, welch, clean  # noqa
import openpyxl

OUT = ("/data/adbrainseq/Publication/2024_scRNA seq/"
       "V1 Science manuscript/ADBrainSeq_V6.0/Acta Neuropathologica Communication/"
       "Major Revision/processed raw data")
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "quickgo")
os.makedirs(CACHE, exist_ok=True)
NPERM = 10000
rng = np.random.default_rng(1)

# ------------------------------------------------------------------ gene sets
def gene_set(go):
    fn = os.path.join(CACHE, go.replace(":", "_") + ".txt")
    if not os.path.exists(fn):
        url = ("https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?"
               f"goId={go}&goUsage=descendants&goUsageRelationships=is_a,part_of,occurs_in"
               "&taxonId=9606&taxonUsage=exact&geneProductType=protein&downloadLimit=50000")
        req = urllib.request.Request(url, headers={"Accept": "text/tsv"})
        try:
            txt = urllib.request.urlopen(req, timeout=600).read().decode()
        except urllib.error.HTTPError:
            # obsolete / merged terms (GO:0030433, GO:0042886) return 500 — same two
            # terms that g:Profiler could not test in the ORA analysis
            txt = ""
        syms = sorted({l.split("\t")[2] for l in txt.splitlines()[1:]
                       if len(l.split("\t")) > 2 and l.split("\t")[2]})
        open(fn, "w").write("\n".join(syms))
    return [s for s in open(fn).read().split("\n") if s]

# ------------------------------------------------------------------ Nativio with t stat
def load_nativio_t():
    wb = openpyxl.load_workbook(os.path.join(BASE, "SuppD3_NativioSeq.xlsx"), read_only=True)
    ws = wb.active
    out = {}
    for r in ws.iter_rows(min_row=3, values_only=True):
        sym = clean(r[1])
        if not sym or sym == "None":
            continue
        try:
            ad = [float(x) for x in r[6:18] if x is not None]
            old = [float(x) for x in r[18:28] if x is not None]
        except (TypeError, ValueError):
            continue
        w = welch(ad, old)
        if w is None:
            continue
        ma, mb, t, p = w
        if ma <= 0 or mb <= 0:
            continue
        if sym not in out or p < out[sym][1]:
            out[sym] = (math.log2(ma / mb), p, t)
    wb.close()
    return out

# ------------------------------------------------------------------ GSEA core
def es_from_positions(pos, w, N, K):
    """Weighted KS enrichment score. pos: sorted 0-based hit positions (…,K), w: |rank metric|."""
    NR = w.sum()
    if NR <= 0:
        return 0.0
    phit = np.cumsum(w) / NR
    i = np.arange(1, K + 1)
    pmiss = (pos + 1 - i) / (N - K)
    top = np.max(phit - pmiss)
    bot = np.min(np.concatenate(([0.0], phit[:-1])) - (pos - (i - 1)) / (N - K))
    return top if top >= -bot else bot

def gsea(ranks, members, nperm=NPERM):
    """ranks: dict gene -> metric. members: list of genes. Returns dict."""
    genes = np.array(list(ranks.keys()))
    vals = np.array([ranks[g] for g in genes], float)
    order = np.argsort(-vals)
    genes, vals = genes[order], vals[order]
    N = len(genes)
    gi = {g: i for i, g in enumerate(genes)}
    pos = np.array(sorted(gi[m] for m in members if m in gi))
    K = len(pos)
    if K < 10:
        return {"K": K, "ES": np.nan, "NES": np.nan, "p": np.nan}
    absv = np.abs(vals)
    es = es_from_positions(pos, absv[pos], N, K)

    # null: size-matched random gene sets
    rp = np.sort(rng.integers(0, N, size=(nperm, K)), axis=1)
    NRr = absv[rp].sum(axis=1)
    phit = np.cumsum(absv[rp], axis=1) / NRr[:, None]
    i = np.arange(1, K + 1)
    pmiss = (rp + 1 - i) / (N - K)
    dev = phit - pmiss
    dev_pre = np.concatenate([np.zeros((nperm, 1)), phit[:, :-1]], axis=1) - (rp - (i - 1)) / (N - K)
    null = np.where(dev.max(axis=1) >= -dev_pre.min(axis=1), dev.max(axis=1), dev_pre.min(axis=1))

    if es >= 0:
        pn = null[null >= 0]
        p = (np.sum(pn >= es) + 1) / (len(pn) + 1)
        nes = es / np.mean(pn) if len(pn) and np.mean(pn) > 0 else np.nan
    else:
        pn = null[null < 0]
        p = (np.sum(pn <= es) + 1) / (len(pn) + 1)
        nes = -es / np.mean(pn) if len(pn) and np.mean(pn) < 0 else np.nan
    # leading edge
    le = int(np.sum(pos <= np.argmax(np.abs(
        np.cumsum(np.where(np.isin(np.arange(N), pos), absv, 0)) / absv[pos].sum()
        - np.cumsum(~np.isin(np.arange(N), pos)) / (N - K))))) if K else 0
    return {"K": K, "ES": float(es), "NES": float(nes), "p": float(p), "leading_edge": le}

def bh(ps):
    ps = np.array(ps, float)
    ok = ~np.isnan(ps)
    out = np.full(len(ps), np.nan)
    v = ps[ok]; n = len(v)
    o = np.argsort(v)
    adj = np.minimum.accumulate((v[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    tmp = np.empty(n); tmp[o] = np.minimum(adj, 1.0)
    out[ok] = tmp
    return out

# ------------------------------------------------------------------ main
def main():
    sets = {go: gene_set(go) for go, _, _ in GO_TERMS}
    for go, name, _ in GO_TERMS:
        print(f"{go} {name[:42]:44s} {len(sets[go]):5d} genes", file=sys.stderr)

    thap = {g: (v[0], v[1]) for g, v in load_thap().items()}
    miz = {g: (v[0], v[1]) for g, v in load_mizuno().items()}
    nat = load_nativio_t()

    def metric(d, use_t=False):
        out = {}
        for g, v in d.items():
            p = max(v[1], 1e-300)
            if use_t:
                out[g] = float(v[2])
            else:
                out[g] = math.copysign(-math.log10(p), v[0])
        return out

    data = {
        "Thapsigargin": (metric(thap), thap),
        "Mizuno": (metric(miz), miz),
        "Nativio": (metric(nat, use_t=True), nat),
    }

    rows = []
    for ds, (ranks, raw) in data.items():
        for variant in ["full ranked list", "p<0.05 subset"]:
            r = ranks if variant.startswith("full") else {
                g: v for g, v in ranks.items() if raw[g][1] < 0.05}
            print(f"GSEA {ds} / {variant}: {len(r)} genes", file=sys.stderr)
            res = []
            for go, name, panel in GO_TERMS:
                out = gsea(r, sets[go])
                out.update({"dataset": ds, "variant": variant, "GO": go, "term": name,
                            "n_ranked_genes": len(r), "set_size_GO": len(sets[go])})
                res.append(out)
            for rr, q in zip(res, bh([x["p"] for x in res])):
                rr["padj"] = float(q) if not np.isnan(q) else ""
                rr["direction"] = ("" if np.isnan(rr["NES"]) else
                                   ("UP" if rr["NES"] > 0 else "DOWN"))
                rr["significant"] = bool(rr["padj"] != "" and rr["padj"] < 0.05)
            rows += res

    fn = os.path.join(OUT, "R1Q1_UPR_GSEA_18terms.csv")
    keys = ["dataset", "variant", "GO", "term", "set_size_GO", "K", "n_ranked_genes",
            "ES", "NES", "direction", "p", "padj", "leading_edge", "significant"]
    with open(fn, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print("wrote", fn, len(rows), "rows", file=sys.stderr)

if __name__ == "__main__":
    main()
