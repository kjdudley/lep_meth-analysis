#!/usr/bin/env python3
"""Where does a genome's methylation actually sit: gene bodies, or everything else?

    meth_partition.py --bed X.combined.bed.gz --gff X.gff3 --label X

Molanna reads 10.63% of CpGs at >=70% in the 30-79x band, against 0.29% for
H. armigera and 0.14% for Xylocopa. That is an order of magnitude past typical
insect gene-body methylation, and the obvious question is whether it IS
gene-body methylation. Xu et al. 2026 (Nat Ecol Evol, doi:10.1038/s41559-026-03090-6)
found that cnidarian gene-body methylation is a genome-defence mechanism
concentrated at transposable elements inside gene bodies, so "abundant" and
"gene-body" are not the same claim.

Classes follow exon_psi.py's convention: GFF3 is 1-based inclusive, the bed
start is 0-based, so a bed row at start s is 1-based position s+1. Exonic
wins over intronic; intronic is inside a gene but outside any exon;
intergenic is outside every gene.

Counted at the project's reporting standard, not a convenient one: >=70%
methylation within a 30x-79x coverage band, with site counts carried.
"""
import argparse, gzip, sys
from bisect import bisect_right
from collections import defaultdict

def merge(iv):
    iv.sort(); out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

def index(d):
    """chrom -> (starts[], ends[]) from merged intervals, for bisect lookup."""
    return {c: ([s for s, _ in m], [e for _, e in m]) for c, m in d.items()}

def hit(idx, c, p):
    t = idx.get(c)
    if not t: return False
    starts, ends = t
    i = bisect_right(starts, p) - 1
    return i >= 0 and p <= ends[i]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--threshold", type=float, default=70.0)
    ap.add_argument("--cov-lo", type=int, default=30)
    ap.add_argument("--cov-hi", type=int, default=79)
    a = ap.parse_args()

    genes, exons = defaultdict(list), defaultdict(list)
    for line in open(a.gff):
        if line[0] == "#": continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8: continue
        if f[2] == "gene":   genes[f[0]].append((int(f[3]), int(f[4])))
        elif f[2] == "exon": exons[f[0]].append((int(f[3]), int(f[4])))
    gi = index({c: merge(v) for c, v in genes.items()})
    ei = index({c: merge(v) for c, v in exons.items()})
    print(f"  annotation: {sum(len(v) for v in genes.values())} gene rows, "
          f"{sum(len(v) for v in exons.values())} exon rows", flush=True)

    n = defaultdict(int); m = defaultdict(int); allsites = defaultdict(int)
    rows = 0
    with gzip.open(a.bed, "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            f = line.split("\t", 7)
            if len(f) < 7: continue
            try:
                p = int(f[1]) + 1; pct = float(f[3]); cov = int(f[5])
            except ValueError:
                continue
            rows += 1
            c = f[0]
            if   hit(ei, c, p): k = "exonic"
            elif hit(gi, c, p): k = "intronic"
            else:               k = "intergenic"
            allsites[k] += 1
            if a.cov_lo <= cov <= a.cov_hi:
                n[k] += 1
                if pct >= a.threshold: m[k] += 1

    tot_n = sum(n.values()); tot_m = sum(m.values())
    print(f"\n=== {a.label}: methylation by class, >={a.threshold:.0f}% in "
          f"{a.cov_lo}-{a.cov_hi}x band")
    print(f"  {rows:,} CpGs total; {tot_n:,} in band; {tot_m:,} methylated in band")
    if not tot_n: sys.exit("FATAL: no sites in band")
    print()
    print(f"  {'class':11} {'sites in band':>14} {'% of band':>10} "
          f"{'methylated':>11} {'% of class':>11} {'% of meth':>10} {'enrich':>7}")
    for k in ("exonic", "intronic", "intergenic"):
        share_sites = 100 * n[k] / tot_n
        share_meth = 100 * m[k] / tot_m if tot_m else 0
        pct_class = 100 * m[k] / n[k] if n[k] else 0
        enr = (share_meth / share_sites) if share_sites else 0
        print(f"  {k:11} {n[k]:>14,} {share_sites:>9.2f}% {m[k]:>11,} "
              f"{pct_class:>10.3f}% {share_meth:>9.2f}% {enr:>6.2f}x")
    print()
    print("  enrich > 1 means methylation is concentrated there relative to")
    print("  how many CpGs that class contributes. Gene-body methylation should")
    print("  put exonic and intronic above 1 and intergenic below it.")

if __name__ == "__main__":
    main()
