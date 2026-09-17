#!/usr/bin/env python3
"""Bound the off-by-one in the bed-streaming exon-methylation route.

exon_psi.py converts the bed's 0-based start to 1-based and tests
`s <= pos <= e`, which covers the exon exactly. ms1_coral_contrast.py and
ms1_degron_dose.py tested `s <= pos < e` on the raw 0-based start, which
covers the exon minus its first base. Both are fixed as of 2026-09-11, but
the published coral and degron numbers were produced by the old window.

This measures the difference rather than arguing about it: exon methylation
computed both ways over the same exon set and the same bed, reporting how
many exons gain a CpG, how many cross the MIN_CPG floor as a result, and how
far the per-exon mean moves.  2026-09-11
"""
import csv, gzip, os, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"; MS = os.path.join(ROOT, "methsplice")


def log(*a): print(*a, flush=True)


def sweep(bed, exons, off):
    """off=0 reproduces the OLD window (pos 0-based, s <= pos < e);
       off=1 reproduces the NEW/correct one (pos 1-based, s <= pos <= e)."""
    by_c = defaultdict(list)
    for c, s, e in exons:
        by_c[c].append((s, e))
    for c in by_c:
        by_c[c].sort()
    acc = defaultdict(lambda: [0.0, 0]); idx, active = {}, {}
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.split("\t", 5)
            c = p[0]; lst = by_c.get(c)
            if lst is None:
                continue
            pos = int(p[1]) + off; score = float(p[3])
            if c not in idx:
                idx[c], active[c] = 0, []
            i = idx[c]; act = active[c]
            while i < len(lst) and lst[i][0] <= pos:
                act.append(lst[i]); i += 1
            idx[c] = i
            hi = (lambda s, e: e) if off else (lambda s, e: e - 1)
            if act and hi(*act[0]) < pos:
                active[c] = act = [r for r in act if hi(*r) >= pos]
            for s, e in act:
                if s <= pos <= (e if off else e - 1):
                    a = acc[(c, s, e)]; a[0] += score; a[1] += 1
    return acc


def main():
    runs = {}
    for r in csv.DictReader(open(os.path.join(ROOT, "pipeline",
                                              "methsplice_runs_coral.tsv")),
                            delimiter="\t"):
        if r["assay"] == "wgbs":
            runs[r["label"].removeprefix("apb_")] = r["run"]
    pdir = os.path.join(MS, "psi_apc_rs")
    code = sorted(c for c in
                  (f.removeprefix("apc_").removesuffix(".psi.tsv")
                   for f in os.listdir(pdir) if f.endswith(".psi.tsv"))
                  if c in runs)[0]
    log(f"using coral fragment {code} ({runs[code]}), RefSeq exons")
    exons = set()
    for r in csv.DictReader(open(os.path.join(pdir, f"apc_{code}.psi.tsv")),
                            delimiter="\t"):
        if float(r["inc"]) + float(r["skip"]) >= 10:
            exons.add((r["chrom"], int(r["start"]), int(r["end"])))
    log(f"exons at MIN_TOT=10: {len(exons)}")
    bed = os.path.join(MS, "bs_ap", runs[code] + ".combined.bed.gz")

    log("sweeping with the OLD window ...")
    old = sweep(bed, exons, 0)
    log("sweeping with the CORRECTED window ...")
    new = sweep(bed, exons, 1)

    gained = 0; tot_old = 0; tot_new = 0; moved = []
    crossed = 0
    for k in exons:
        o = old.get(k, [0.0, 0]); n = new.get(k, [0.0, 0])
        tot_old += o[1]; tot_new += n[1]
        if n[1] > o[1]:
            gained += 1
        if o[1] < 3 <= n[1]:
            crossed += 1
        if o[1] >= 3 and n[1] >= 3:
            moved.append(abs(n[0]/n[1] - o[0]/o[1]))
    moved = np.array(moved)
    log("")
    log(f"CpGs counted   old {tot_old:,}   corrected {tot_new:,}   "
        f"(+{tot_new - tot_old:,}, {100*(tot_new-tot_old)/max(tot_old,1):.3f}%)")
    log(f"exons gaining a CpG: {gained} of {len(exons)} "
        f"({100*gained/len(exons):.2f}%)")
    log(f"exons crossing the 3-CpG floor because of it: {crossed}")
    log(f"per-exon mean methylation shift, exons above the floor in both: "
        f"n={len(moved)}  mean {moved.mean():.4f} pp  max {moved.max():.4f} pp  "
        f"exons moving >1 pp: {int((moved > 1).sum())}")


if __name__ == "__main__":
    main()
