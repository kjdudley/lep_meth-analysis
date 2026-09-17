#!/usr/bin/env python3
"""Audit of the DKO novel-junction-gene concentration claim.

The manuscript states: exons in genes that gained novel junctions carry the
dPSI excess at 13.77% against 8.53% elsewhere, and "measured against the
contrast's own relabeling null, approximately 62% of the genuine excess
concentrates in the 27% of exons whose genes gained novel junctions."

ms1_dko_overlap.py computes 13.77/8.53 but no null and no 62%. This
recomputes both, reports the null across all 4 exact 2v2 relabelings, and
repeats everything under a containment-only gene-assignment rule (the
published rule also accepts any junction starting within 10 kb of a tested
exon).  2026-09-11
"""
import csv, itertools, os, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"
MS   = os.path.join(ROOT, "methsplice")
JUNC = os.path.join(ROOT, "rna", "junc")
PSI  = os.path.join(MS, "psi_hs")

def log(*a): print(*a, flush=True)

def load_label(label):
    out = {}
    for r in csv.DictReader(open(os.path.join(PSI, label + ".psi.tsv")),
                            delimiter="\t"):
        tot = float(r["inc"]) + float(r["skip"])
        if tot >= 10:
            out[(r["chrom"], r["start"], r["end"])] = (float(r["inc"])/tot, r["gene"])
    return out

def juncs(label, min_reads=5):
    s = set()
    for r in csv.DictReader(open(os.path.join(JUNC, label + ".junctions.tsv")),
                            delimiter="\t"):
        try:
            if float(r["n_reads"]) >= min_reads:
                s.add((r["chrom"], int(r["intron_start"]), int(r["intron_end"])))
        except (ValueError, KeyError):
            continue
    return s

def main():
    row = next(r for r in csv.DictReader(
        open(os.path.join(ROOT, "contrasts.tsv")), delimiter="\t")
        if r["contrast"] == "HCT116_DKO_vs_WT")
    g1 = row["group1"].split(","); g2 = row["group2"].split(",")
    labs = g1 + g2; n1 = len(g1)
    log(f"HCT116_DKO_vs_WT: g1(DKO)={g1}  g2(WT)={g2}")
    data = {l: load_label(l) for l in labs}
    common = set.intersection(*(set(d) for d in data.values()))
    log(f"exons tested in all four libraries: {len(common)}")
    J = {l: juncs(l) for l in labs}
    for l in labs: log(f"  {l}: {len(J[l])} junctions at >=5 reads")

    by_chrom = defaultdict(list)
    for k in common:
        by_chrom[k[0]].append((int(k[1]), int(k[2]), data[g1[0]][k][1]))

    def genes_for(novel, rule):
        out = set()
        for c, js, je in novel:
            for s, e, gene in by_chrom.get(c, []):
                if rule == "published":
                    if (js < e and je > s) or (s - 10000 < js < e + 10000):
                        out.add(gene)
                else:                      # containment / overlap only
                    if js < e and je > s:
                        out.add(gene)
        return out

    def split_stats(order, rule):
        a = [labs[i] for i in order[:n1]]; b = [labs[i] for i in order[n1:]]
        novel = set.intersection(*(J[l] for l in a)) - set.union(*(J[l] for l in b))
        ng = genes_for(novel, rule)
        dp = {}
        for k in common:
            dp[k] = abs(sum(data[l][k][0] for l in a)/len(a)
                        - sum(data[l][k][0] for l in b)/len(b))
        ing = [k for k in common if data[g1[0]][k][1] in ng]
        outg = [k for k in common if data[g1[0]][k][1] not in ng]
        r_in = sum(1 for k in ing if dp[k] > .05)/max(len(ing), 1)
        r_out = sum(1 for k in outg if dp[k] > .05)/max(len(outg), 1)
        r_all = sum(1 for k in common if dp[k] > .05)/len(common)
        return dict(novel=len(novel), genes=len(ng), n_in=len(ing), n_out=len(outg),
                    r_in=r_in, r_out=r_out, r_all=r_all)

    idx = list(range(len(labs)))
    for rule in ("published", "overlap-only"):
        log("\n" + "="*72); log(f"GENE-ASSIGNMENT RULE: {rule}"); log("="*72)
        o = split_stats(idx, rule)
        log(f"observed: {o['novel']} novel junctions -> {o['genes']} genes")
        log(f"  exons in novel-junction genes: {o['n_in']} "
            f"({100*o['n_in']/(o['n_in']+o['n_out']):.1f}% of tested)"
            f"   manuscript says 27%")
        log(f"  |dPSI|>0.05 in-gene {100*o['r_in']:.2f}%  out-gene "
            f"{100*o['r_out']:.2f}%   manuscript says 13.77% / 8.53%")
        nulls = []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if list(comb) == idx[:n1] or rest == idx[:n1]: continue
            nulls.append(split_stats(list(comb)+rest, rule))
        log(f"  {len(nulls)} exact non-identity relabelings:")
        for i, s in enumerate(nulls):
            log(f"    split {i+1}: overall rate {100*s['r_all']:.2f}%  "
                f"in {100*s['r_in']:.2f}%  out {100*s['r_out']:.2f}%  "
                f"({s['novel']} novel -> {s['genes']} genes)")
        r0s = [s["r_all"] for s in nulls]
        log(f"  null overall rate: mean {100*np.mean(r0s):.2f}%  "
            f"range {100*min(r0s):.2f}-{100*max(r0s):.2f}%")

        def share(r0):
            ein = (o["r_in"] - r0) * o["n_in"]
            eout = (o["r_out"] - r0) * o["n_out"]
            tot = ein + eout
            return float("nan") if tot <= 0 else ein / tot
        log(f"\n  share of excess in novel-junction genes, as a function of the null:")
        for r0 in sorted(set([np.mean(r0s), min(r0s), max(r0s)] +
                             [x/100 for x in (5, 6, 7, 8, 9)])):
            s = share(r0)
            tag = ""
            if abs(r0 - np.mean(r0s)) < 1e-12: tag = "  <-- null mean"
            if abs(r0 - min(r0s)) < 1e-12: tag += "  <-- null min"
            if abs(r0 - max(r0s)) < 1e-12: tag += "  <-- null max"
            log(f"    null={100*r0:5.2f}%  ->  share = "
                + (f"{100*s:5.1f}%" if s == s else "  n/a (no excess)") + tag)
        log(f"  manuscript asserts approximately 62%")

if __name__ == "__main__":
    main()
