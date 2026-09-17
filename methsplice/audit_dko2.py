#!/usr/bin/env python3
"""Settle the DKO concentration null.

The manuscript computes the 62% by assuming exons in novel-junction genes
carry the population-average relabeling null (6.97%). audit_dko.py showed the
relabelings' OWN novel-junction classes run at ~11%, but those classes are
smaller and differently selected, so they are not a clean comparison.

The clean test: FIX the gene partition at the observed one, then ask what
fraction of exons in each class exceeds |dPSI|>0.05 under each of the four
exact relabelings. That is the null rate for the exact sets the 62% uses.
2026-09-11
"""
import csv, itertools, os
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"
MS   = os.path.join(ROOT, "methsplice")
JUNC = os.path.join(ROOT, "rna", "junc")
PSI  = os.path.join(MS, "psi_hs")

def log(*a): print(*a, flush=True)

def load_label(l):
    out = {}
    for r in csv.DictReader(open(os.path.join(PSI, l + ".psi.tsv")), delimiter="\t"):
        tot = float(r["inc"]) + float(r["skip"])
        if tot >= 10:
            out[(r["chrom"], r["start"], r["end"])] = (float(r["inc"])/tot, r["gene"])
    return out

def juncs(l, m=5):
    s = set()
    for r in csv.DictReader(open(os.path.join(JUNC, l + ".junctions.tsv")), delimiter="\t"):
        try:
            if float(r["n_reads"]) >= m:
                s.add((r["chrom"], int(r["intron_start"]), int(r["intron_end"])))
        except (ValueError, KeyError): continue
    return s

row = next(r for r in csv.DictReader(open(os.path.join(ROOT, "contrasts.tsv")),
                                     delimiter="\t") if r["contrast"] == "HCT116_DKO_vs_WT")
g1 = row["group1"].split(","); g2 = row["group2"].split(","); labs = g1 + g2; n1 = len(g1)
data = {l: load_label(l) for l in labs}
common = sorted(set.intersection(*(set(d) for d in data.values())))
J = {l: juncs(l) for l in labs}
log(f"exons tested in all four: {len(common)}")

by_chrom = defaultdict(list)
for k in common:
    by_chrom[k[0]].append((int(k[1]), int(k[2]), data[g1[0]][k][1]))

novel = set.intersection(*(J[l] for l in g1)) - set.union(*(J[l] for l in g2))
log(f"observed DKO-novel junctions: {len(novel)}")

def genes_for(rule):
    out = set()
    for c, js, je in novel:
        for s, e, gene in by_chrom.get(c, []):
            if rule == "published":
                if (js < e and je > s) or (s - 10000 < js < e + 10000): out.add(gene)
            else:
                if js < e and je > s: out.add(gene)
    return out

idx = list(range(len(labs)))
splits = []
for comb in itertools.combinations(idx, n1):
    rest = [i for i in idx if i not in comb]
    if list(comb) == idx[:n1] or rest == idx[:n1]: continue
    splits.append((list(comb), rest))

def rate(ks, a, b):
    n = 0
    for k in ks:
        p1 = sum(data[labs[i]][k][0] for i in a)/len(a)
        p2 = sum(data[labs[i]][k][0] for i in b)/len(b)
        if abs(p1-p2) > .05: n += 1
    return n/len(ks)

for rule in ("published", "overlap-only"):
    ng = genes_for(rule)
    ing = [k for k in common if data[g1[0]][k][1] in ng]
    outg = [k for k in common if data[g1[0]][k][1] not in ng]
    log("\n" + "="*70); log(f"RULE: {rule}   genes={len(ng)}  in={len(ing)} "
                            f"({100*len(ing)/len(common):.1f}%)  out={len(outg)}")
    ri = rate(ing, idx[:n1], idx[n1:]); ro = rate(outg, idx[:n1], idx[n1:])
    log(f"  OBSERVED   in {100*ri:.2f}%   out {100*ro:.2f}%   "
        f"enrichment {ri/ro:.2f}x")
    nis, nos = [], []
    for a, b in splits:
        x = rate(ing, a, b); y = rate(outg, a, b); nis.append(x); nos.append(y)
        log(f"    relabeling {a}v{b}:  in {100*x:.2f}%   out {100*y:.2f}%")
    mi, mo = float(np.mean(nis)), float(np.mean(nos))
    log(f"  NULL, SAME FIXED PARTITION:  in {100*mi:.2f}%  out {100*mo:.2f}%")
    log(f"  manuscript assumes both classes have the pooled null 6.97%")
    for tag, ni_, no_ in (("manuscript's pooled null", .0697, .0697),
                          ("measured, fixed partition", mi, mo)):
        ein = (ri-ni_)*len(ing); eout = (ro-no_)*len(outg)
        if ein+eout <= 0: log(f"  {tag}: no net excess"); continue
        log(f"  {tag}: share in novel-junction genes = "
            f"{100*ein/(ein+eout):.1f}%   (excess {ein:.0f} vs {eout:.0f} exons)")
log("\nDONE")
