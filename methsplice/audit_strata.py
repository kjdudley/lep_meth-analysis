#!/usr/bin/env python3
"""Gene-clustered CIs on the coverage-stratified degron dose slopes.

audit_hs_degron.py part 2 found the measured-derasure slope negative in the
sparsely covered strata and POSITIVE in the best-covered stratum, in all six
degron contrasts. That was printed without CIs. This supplies them, for both
coupling definitions, and prints per-stratum exon characteristics so the
measurement-quality/exon-class confound is visible.  2026-09-11
"""
import csv, gzip, os, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"; MS = os.path.join(ROOT, "methsplice")
def log(*a): print(*a, flush=True)

def load_psi(p, mt=10):
    out = {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        t = float(r["inc"]) + float(r["skip"])
        if t >= mt: out[(r["chrom"], int(r["start"]), int(r["end"]))] = (float(r["inc"])/t, r["gene"])
    return out

def exon_meth(bed, exons, min_cpg=3):
    by_c = defaultdict(list)
    for c, s, e in exons: by_c[c].append((s, e))
    for c in by_c: by_c[c].sort()
    acc = defaultdict(lambda: [0.0, 0]); idx, active = {}, {}
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            p = line.split("\t", 5); c = p[0]; lst = by_c.get(c)
            if lst is None: continue
            pos = int(p[1]); sc = float(p[3])
            if c not in idx: idx[c], active[c] = 0, []
            i = idx[c]; act = active[c]
            while i < len(lst) and lst[i][0] <= pos: act.append(lst[i]); i += 1
            idx[c] = i
            if act and act[0][1] <= pos: active[c] = act = [r for r in act if r[1] > pos]
            for s, e in act:
                if s <= pos < e: a = acc[(c, s, e)]; a[0] += sc; a[1] += 1
    return {k: (v[0]/v[1], v[1]) for k, v in acc.items() if v[1] >= min_cpg}

def slope_ci(x, y, g, nboot=2000, seed=7):
    x = np.asarray(x, float); y = np.asarray(y, float); g = np.asarray(g)
    o = np.argsort(g, kind="stable"); x, y, g = x[o], y[o], g[o]
    _, st = np.unique(g, return_index=True); b = np.append(np.sort(st), len(g))
    def seg(v):
        cs = np.concatenate(([0.0], np.cumsum(v))); return cs[b[1:]] - cs[b[:-1]]
    n_g = np.diff(b).astype(float); Sx, Sy, Sxx, Sxy = seg(x), seg(y), seg(x*x), seg(x*y)
    def sl(n, sx, sy, sxx, sxy):
        den = sxx - sx*sx/n; num = sxy - sx*sy/n
        return np.where(den > 0, num/np.where(den > 0, den, 1.0), 0.0)
    obs = float(sl(n_g.sum(), Sx.sum(), Sy.sum(), Sxx.sum(), Sxy.sum()))
    rng = np.random.default_rng(seed); G = len(n_g)
    d = rng.integers(0, G, size=(nboot, G))
    bs = np.sort(sl(n_g[d].sum(1), Sx[d].sum(1), Sy[d].sum(1), Sxx[d].sum(1), Sxy[d].sum(1)))
    return obs, bs[int(.025*nboot)], bs[int(.975*nboot)]

contrasts = {r["contrast"]: (r["group1"].split(","), r["group2"].split(","))
             for r in csv.DictReader(open(os.path.join(ROOT, "contrasts.tsv")), delimiter="\t")
             if r["contrast"].startswith("degron_")}
def cond_of(title):
    """Copied verbatim from ms1_degron_dose.py: dual must be matched FIRST."""
    day = "Day6" if "Day6" in title else ("Day8" if "Day8" in title else None)
    if day is None: return None
    if "UHRF1-AID/DNMT1-AID" in title: return ("dual", day)
    if "DNMT1-AID" in title: return ("d1aid", day)
    if "UHRF1-AID" in title: return ("u1aid", day)
    return None

beds = defaultdict(list)
for r in csv.DictReader(open(os.path.join(MS, "runs_bs2.tsv")), delimiter="\t"):
    if r["outdir"] != "bs_degron": continue
    c = cond_of(r["title"])
    if c:
        pth = os.path.join(MS, "bs_degron", r["run"] + ".combined.bed.gz")
        if os.path.exists(pth): beds[c].append(pth)
for c in sorted(beds): log(f"  {c}: {len(beds[c])} methylomes")
if not beds: sys.exit("REFUSING: no degron beds found")

labels = sorted({l for g in contrasts.values() for gg in g for l in gg})
psi, gene_of = {}, {}
for l in labels:
    d = load_psi(os.path.join(MS, "psi_hs", l + ".psi.tsv"))
    psi[l] = {k: v[0] for k, v in d.items()}
    for k, v in d.items(): gene_of[k] = v[1]
common = set.intersection(*(set(psi[l]) for l in labels))
log(f"exons at depth in all {len(labels)} labels: {len(common)}")
base = exon_meth(os.path.join(MS, "bs_hs", "hct116_wt.hcg.combined.bed.gz"), common)
log(f"WT NOMe baseline exons (cov>=1, >=3 CpG): {len(base)}")
cond = {}
for key, ps in sorted(beds.items()):
    acc = defaultdict(list)
    for p in ps:
        for k, v in exon_meth(p, common).items(): acc[k].append(v[0])
    cond[key] = {k: float(np.mean(v)) for k, v in acc.items() if len(v) == len(ps)}
    log(f"  {key}: {len(cond[key])} exons over {len(ps)} methylomes")

log("\n" + "="*100)
log("COVERAGE-STRATIFIED DOSE SLOPES WITH GENE-CLUSTERED CIs")
log("="*100)
for name in sorted(contrasts):
    g1, g2 = contrasts[name]
    line, day = name.split("_")[1], name.split("_")[2]
    cm = cond.get((line, day))
    if not cm: continue
    usable = sorted(k for k in common if k in base and k in cm)
    ncpg = np.array([base[k][1] for k in usable])
    dpsi = np.array([abs(np.mean([psi[l][k] for l in g1]) - np.mean([psi[l][k] for l in g2]))
                     for k in usable])
    genes = np.array([gene_of[k] for k in usable])
    xb = np.array([base[k][0] for k in usable])
    xd = np.array([abs(cm[k] - base[k][0]) for k in usable])
    log(f"\n{name}   n={len(usable)}")
    for tag, xv in (("baseline_meth", xb), ("measured_derasure", xd)):
        row = []
        for lo_, hi_, nm in ((3, 5, "3-5"), (6, 15, "6-15"), (16, 10**9, "16+")):
            m = (ncpg >= lo_) & (ncpg <= hi_)
            if m.sum() < 300: row.append(f"{nm}: n/a"); continue
            s, lo, hi = slope_ci(xv[m], dpsi[m], genes[m])
            star = "*" if (lo > 0 or hi < 0) else " "
            row.append(f"{nm}: {s:+.6f}{star} [{lo:+.6f},{hi:+.6f}] n={int(m.sum())}")
        log(f"  {tag:18s} " + "   ".join(row))
    # exon-class confound: what ARE the high-CpG exons?
    log("  exon character by stratum: " + "   ".join(
        f"{nm}: len={int(np.mean([k[2]-k[1] for k,mm in zip(usable,(ncpg>=lo_)&(ncpg<=hi_)) if mm])) if ((ncpg>=lo_)&(ncpg<=hi_)).sum() else 0}bp "
        f"meth={np.mean(xb[(ncpg>=lo_)&(ncpg<=hi_)]):.1f}% "
        f"|dPSI|={np.mean(dpsi[(ncpg>=lo_)&(ncpg<=hi_)]):.4f}"
        for lo_, hi_, nm in ((3,5,"3-5"),(6,15,"6-15"),(16,10**9,"16+"))
        if ((ncpg>=lo_)&(ncpg<=hi_)).sum() > 0))
log("\n* = CI excludes zero")
log("DONE")
