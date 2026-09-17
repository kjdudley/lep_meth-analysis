#!/usr/bin/env python3
"""Independent source-level audit of the METHSPLICE human + degron arms.

Written 2026-09-11 as a check on ms1_contrasts.py / ms1_degron_dose.py /
ms1_dko_dose.py. Reimplements the statistics from the PSI tables and the
methylation beds rather than importing or re-reading their outputs, and adds
three sensitivity tests the originals do not run:

  S1  per-CpG coverage floor  (originals apply NONE; Methods claims one)
  S2  baseline-coverage strata (attenuation from the borrowed WT baseline)
  S3  exonic mean methylation of every bed (characterises that baseline)

Prints a report; writes nothing into the repo.
"""
import csv, gzip, itertools, os, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"
MS   = os.path.join(ROOT, "methsplice")
MIN_JUNC = 10
MIN_CPG  = 3
NBOOT    = 2000
SEED     = 13

def log(*a):
    print(*a, flush=True)

def find(name):
    """contrasts.tsv lives at the project root on Aqua but in methsplice/ in
    the repo; the result tables are in methsplice/ in both. Resolve either."""
    for d in (MS, ROOT):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    sys.exit("FATAL: cannot find " + name)

# ---------------------------------------------------------------- PSI loading
def load_psi(label):
    for d in ("psi_hs", "psi_hs_meth", "psi_mm"):
        p = os.path.join(MS, d, label + ".psi.tsv")
        if os.path.exists(p):
            out = {}
            with open(p) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    tot = float(r["inc"]) + float(r["skip"])
                    if tot >= MIN_JUNC:
                        out[(r["chrom"], int(r["start"]), int(r["end"]))] = (
                            float(r["inc"]) / tot, r["gene"])
            return out
    sys.exit("no psi table for " + label)

# ------------------------------------------------------- bed -> per-exon meth
def exon_meth(bed, exons, min_cov=1):
    """Sweep one bed. Returns {(c,s,e): (mean_score, n_cpg)} for n_cpg>=MIN_CPG."""
    by_c = defaultdict(list)
    for c, s, e in exons:
        by_c[c].append((s, e))
    for c in by_c:
        by_c[c].sort()
    acc = defaultdict(lambda: [0.0, 0])
    idx, active = {}, {}
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.split("\t", 6)
            c = p[0]
            lst = by_c.get(c)
            if lst is None:
                continue
            if min_cov > 1 and float(p[5]) < min_cov:
                continue
            pos = int(p[1]); score = float(p[3])
            if c not in idx:
                idx[c], active[c] = 0, []
            i = idx[c]; act = active[c]
            while i < len(lst) and lst[i][0] <= pos:
                act.append(lst[i]); i += 1
            idx[c] = i
            if act and act[0][1] <= pos:
                active[c] = act = [r for r in act if r[1] > pos]
            for s, e in act:
                if s <= pos < e:
                    a = acc[(c, s, e)]
                    a[0] += score; a[1] += 1
    return {k: (v[0] / v[1], v[1]) for k, v in acc.items() if v[1] >= MIN_CPG}

# ------------------------------------------------- gene-clustered bootstrap CI
def slope_ci(x, y, gene_ids, nboot=NBOOT, seed=SEED):
    """Exact gene-clustered bootstrap via per-gene sufficient statistics."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    g = np.asarray(gene_ids)
    order = np.argsort(g, kind="stable")
    x, y, g = x[order], y[order], g[order]
    _, start = np.unique(g, return_index=True)
    bounds = np.append(start, len(g))
    def seg(v):
        cs = np.concatenate(([0.0], np.cumsum(v)))
        return cs[bounds[1:]] - cs[bounds[:-1]]
    n_g  = np.diff(bounds).astype(float)
    Sx, Sy = seg(x), seg(y)
    Sxx, Sxy = seg(x * x), seg(x * y)
    def sl(n, sx, sy, sxx, sxy):
        den = sxx - sx * sx / n
        num = sxy - sx * sy / n
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
    obs = float(sl(n_g.sum(), Sx.sum(), Sy.sum(), Sxx.sum(), Sxy.sum()))
    rng = np.random.default_rng(seed)
    G = len(n_g)
    draws = rng.integers(0, G, size=(nboot, G))
    bs = sl(n_g[draws].sum(1), Sx[draws].sum(1), Sy[draws].sum(1),
            Sxx[draws].sum(1), Sxy[draws].sum(1))
    bs = np.sort(bs)
    return obs, bs[int(0.025 * nboot)], bs[int(0.975 * nboot)]

# ------------------------------------------------------------- part 1: contrasts
def part1(contrasts):
    log("\n" + "=" * 74)
    log("PART 1  contrasts_results.tsv recomputed from the PSI tables")
    log("=" * 74)
    stored = {r["contrast"]: r for r in
              csv.DictReader(open(find("contrasts_results.tsv")),
                             delimiter="\t")}
    log(f"{'contrast':32s} {'n_exons':>8s} {'>0.05 obs':>10s} {'null':>9s} "
        f"{'excess':>7s} {'stored':>7s} {'relab':>7s}")
    for name, (g1, g2) in contrasts.items():
        labs = g1 + g2; n1 = len(g1)
        data = {l: load_psi(l) for l in labs}
        common = set.intersection(*(set(d) for d in data.values()))
        keys = sorted(common)
        V = np.array([[data[l][k][0] for l in labs] for k in keys])
        idx = list(range(len(labs)))
        def dp(order):
            o = list(order)
            return np.abs(V[:, o[:n1]].mean(1) - V[:, o[n1:]].mean(1))
        obs = dp(idx)
        n05, n10, nmean = [], [], []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if list(comb) == idx[:n1] or rest == idx[:n1]:
                continue
            d = dp(list(comb) + rest)
            n05.append(int((d > 0.05).sum())); n10.append(int((d > 0.10).sum()))
            nmean.append(float(d.mean()))
        o05 = int((obs > 0.05).sum()); o10 = int((obs > 0.10).sum())
        e05 = float(np.mean(n05)); e10 = float(np.mean(n10))
        ge = sum(1 for c in n05 if c >= o05)
        st = stored[name]
        flag = "" if (int(st["n_exons"]) == len(keys)
                      and int(st["n_gt05_obs"]) == o05
                      and abs(float(st["excess05"]) - o05 / e05) < 0.015) else "  <-- MISMATCH"
        log(f"{name:32s} {len(keys):8d} {o05:10d} {e05:9.1f} "
            f"{o05/e05:7.2f} {st['excess05']:>7s} {ge:3d}/{len(n05):<3d}{flag}")
        log(f"{'':32s} mean|dPSI| {obs.mean():.5f} (stored {st['mean_dpsi_obs']})"
            f"  null {np.mean(nmean):.5f} (stored {st['mean_dpsi_null']})"
            f"  >0.10 {o10} vs {e10:.1f} (stored {st['n_gt10_obs']}/{st['n_gt10_null']})")

# ------------------------------------------------------------ part 2: degron dose
def part2():
    log("\n" + "=" * 74)
    log("PART 2  degron_dose.tsv recomputed, plus coverage/attenuation tests")
    log("=" * 74)
    conds = {}
    for r in csv.DictReader(open(find("runs_bs2.tsv")), delimiter="\t"):
        if r["outdir"] != "bs_degron":
            continue
        t = r["title"]
        day = "Day6" if "Day6" in t else ("Day8" if "Day8" in t else None)
        if day is None:
            continue
        line = ("dual" if "UHRF1-AID/DNMT1-AID" in t else
                "d1aid" if "DNMT1-AID" in t else
                "u1aid" if "UHRF1-AID" in t else None)
        if line:
            conds.setdefault((line, day), []).append(
                os.path.join(MS, "bs_degron", r["run"] + ".combined.bed.gz"))
    for c in sorted(conds):
        log(f"  {c}: {len(conds[c])} methylomes")

    contrasts = {}
    for r in csv.DictReader(open(find("contrasts.tsv")), delimiter="\t"):
        n = r["contrast"]
        if n.startswith("degron_") and ("Day6" in n or "Day8" in n):
            contrasts[n] = (r["group1"].split(","), r["group2"].split(","))
    labels = sorted({l for g in contrasts.values() for l in g[0] + g[1]})
    psi, gene_of = {}, {}
    for l in labels:
        d = load_psi(l)
        psi[l] = {k: v[0] for k, v in d.items()}
        for k, v in d.items():
            gene_of[k] = v[1]
    common = set.intersection(*(set(psi[l]) for l in labels))
    log(f"exons at depth>={MIN_JUNC} in all {len(labels)} labels: {len(common)}")

    wt_bed = os.path.join(MS, "bs_hs", "hct116_wt.hcg.combined.bed.gz")
    log("streaming WT NOMe baseline (cov>=1) ...")
    base1 = exon_meth(wt_bed, common, 1)
    log("streaming WT NOMe baseline (cov>=5) ...")
    base5 = exon_meth(wt_bed, common, 5)
    log(f"  baseline exons: cov>=1 {len(base1)}   cov>=5 {len(base5)}")
    log(f"  WT NOMe exonic mean methylation: "
        f"cov>=1 {np.mean([v[0] for v in base1.values()]):.2f}%  "
        f"cov>=5 {np.mean([v[0] for v in base5.values()]):.2f}%")

    cond1, cond5 = {}, {}
    for c in sorted(conds):
        a1, a5 = defaultdict(list), defaultdict(list)
        for p in conds[c]:
            log(f"streaming {os.path.basename(p)} ...")
            for k, v in exon_meth(p, common, 1).items():
                a1[k].append(v[0])
            for k, v in exon_meth(p, common, 5).items():
                a5[k].append(v[0])
        cond1[c] = {k: float(np.mean(v)) for k, v in a1.items() if len(v) == len(conds[c])}
        cond5[c] = {k: float(np.mean(v)) for k, v in a5.items() if len(v) == len(conds[c])}
        log(f"  {c}: exonic mean {np.mean(list(cond1[c].values())):.2f}% "
            f"({len(cond1[c])} exons)   cov>=5 {np.mean(list(cond5[c].values())):.2f}% "
            f"({len(cond5[c])} exons)")

    stored = defaultdict(dict)
    for r in csv.DictReader(open(find("degron_dose.tsv")), delimiter="\t"):
        stored[r["contrast"]][r["test"]] = r

    log("\n-- recomputed against degron_dose.tsv (cov>=1, as implemented) --")
    rows = []
    for name in sorted(contrasts):
        g1, g2 = contrasts[name]
        line, day = name.split("_")[1], name.split("_")[2]
        cm = cond1.get((line, day))
        if not cm:
            continue
        usable = sorted(k for k in common if k in base1 and k in cm)
        er = float(np.mean([cm[k] - base1[k][0] for k in usable]))
        dpsi = np.array([abs(np.mean([psi[l][k] for l in g1])
                             - np.mean([psi[l][k] for l in g2])) for k in usable])
        genes = [gene_of[k] for k in usable]
        for test, xv in (("baseline_meth", np.array([base1[k][0] for k in usable])),
                         ("measured_derasure",
                          np.array([abs(cm[k] - base1[k][0]) for k in usable]))):
            s, lo, hi = slope_ci(xv, dpsi, genes)
            st = stored[name][test]
            ok = (abs(s - float(st["slope"])) < 3e-6
                  and abs(er - float(st["mean_erasure_pct"])) < 0.05
                  and len(usable) == int(st["n_exon"]))
            log(f"{name:20s} {test:18s} n={len(usable):6d} (stored {st['n_exon']:>6s}) "
                f"erasure {er:+7.2f} (stored {st['mean_erasure_pct']:>7s})  "
                f"slope {s:+.6f} [{lo:+.6f},{hi:+.6f}]  stored {st['slope']}"
                f" [{st['ci_lo']},{st['ci_hi']}] {'OK' if ok else '<-- MISMATCH'}")
        rows.append((name, er))

    log("\n-- S1: same, with a per-CpG coverage floor of 5 --")
    for name in sorted(contrasts):
        g1, g2 = contrasts[name]
        line, day = name.split("_")[1], name.split("_")[2]
        cm = cond5.get((line, day))
        if not cm:
            continue
        usable = sorted(k for k in common if k in base5 and k in cm)
        if len(usable) < 200:
            log(f"{name:20s} only {len(usable):5d} exons survive cov>=5 -- FLOOR NOT APPLICABLE")
            continue
        er = float(np.mean([cm[k] - base5[k][0] for k in usable]))
        dpsi = np.array([abs(np.mean([psi[l][k] for l in g1])
                             - np.mean([psi[l][k] for l in g2])) for k in usable])
        genes = [gene_of[k] for k in usable]
        xv = np.array([abs(cm[k] - base5[k][0]) for k in usable])
        s, lo, hi = slope_ci(xv, dpsi, genes)
        xb = np.array([base5[k][0] for k in usable])
        sb, lob, hib = slope_ci(xb, dpsi, genes)
        log(f"{name:20s} n={len(usable):6d} erasure {er:+7.2f}  "
            f"derasure slope {s:+.6f} [{lo:+.6f},{hi:+.6f}]  "
            f"baseline slope {sb:+.6f} [{lob:+.6f},{hib:+.6f}]")

    log("\n-- S2: attenuation, measured-derasure slope by WT baseline CpG count --")
    for name in sorted(contrasts):
        g1, g2 = contrasts[name]
        line, day = name.split("_")[1], name.split("_")[2]
        cm = cond1.get((line, day))
        if not cm:
            continue
        usable = sorted(k for k in common if k in base1 and k in cm)
        ncpg = np.array([base1[k][1] for k in usable])
        dpsi = np.array([abs(np.mean([psi[l][k] for l in g1])
                             - np.mean([psi[l][k] for l in g2])) for k in usable])
        xv = np.array([abs(cm[k] - base1[k][0]) for k in usable])
        genes = np.array([gene_of[k] for k in usable])
        out = []
        for lo_, hi_ in ((3, 5), (6, 15), (16, 10**9)):
            m = (ncpg >= lo_) & (ncpg <= hi_)
            if m.sum() < 500:
                out.append(f"{lo_}-{hi_}:n/a")
                continue
            s, a, b = slope_ci(xv[m], dpsi[m], genes[m])
            out.append(f"{lo_}-{hi_ if hi_ < 10**8 else '+'}: {s:+.6f} (n={int(m.sum())})")
        log(f"{name:20s} " + "   ".join(out))

    log("\n-- slope-of-slopes on the recomputed erasures --")
    for test in ("baseline_meth", "measured_derasure"):
        xs = np.array([float(stored[n][test]["mean_erasure_pct"]) for n, _ in rows])
        ys = np.array([float(stored[n][test]["slope"]) for n, _ in rows])
        b = np.polyfit(-xs, ys, 1)[0]
        log(f"  {test:18s} d(slope)/d(erasure depth) = {b:+.2e} per point")

# ------------------------------------------------------------ part 3: DKO dose
def part3():
    log("\n" + "=" * 74)
    log("PART 3  HCT116 DKO dose-response recomputed")
    log("=" * 74)
    g1 = ["hsk_hct116dko_a", "hsk_hct116dko_b"]
    g2 = ["hsk_hct116wt_a", "hsk_hct116wt_b"]
    data = {l: load_psi(l) for l in g1 + g2}
    common = set.intersection(*(set(d) for d in data.values()))
    log(f"exons at depth in all four labels: {len(common)}")
    wt = os.path.join(MS, "bs_hs", "hct116_wt.hcg.combined.bed.gz")
    dko = os.path.join(MS, "bs_hs", "hct116_dko1.hcg.combined.bed.gz")
    log("streaming WT NOMe ...");  bw = exon_meth(wt, common, 1)
    log("streaming DKO NOMe ..."); bd = exon_meth(dko, common, 1)
    usable = sorted(k for k in common if k in bw and k in bd)
    mw = np.mean([bw[k][0] for k in usable]); md = np.mean([bd[k][0] for k in usable])
    log(f"exonic mean methylation  WT {mw:.2f}%   DKO {md:.2f}%   erasure {md-mw:+.2f} pp")
    dpsi = np.array([abs(np.mean([data[l][k][0] for l in g1])
                         - np.mean([data[l][k][0] for l in g2])) for k in usable])
    genes = [data[g1[0]][k][1] for k in usable]
    xb = np.array([bw[k][0] for k in usable])
    s, lo, hi = slope_ci(xb, dpsi, genes)
    log(f"slope(|dPSI| ~ WT baseline meth) = {s:+.6f} [{lo:+.6f}, {hi:+.6f}]  n={len(usable)}")
    log("  DESIGN records ALL -0.000100 [-0.000138, -0.000063]")

def main():
    contrasts = {}
    for r in csv.DictReader(open(find("contrasts.tsv")), delimiter="\t"):
        n = r["contrast"]
        if n.startswith("degron_") or n == "HCT116_DKO_vs_WT":
            contrasts[n] = (r["group1"].split(","), r["group2"].split(","))
    part1(contrasts)
    part2()
    part3()
    log("\nAUDIT COMPLETE")

if __name__ == "__main__":
    main()
