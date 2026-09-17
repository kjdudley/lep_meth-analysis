#!/usr/bin/env python3
"""Source-level audit of the remaining METHSPLICE arms: mouse, coral (both
annotation routes) and the positive control (assulta + armigera).
Recomputes each statistic from the PSI tables and methylation beds rather
than re-reading stored outputs. Prints a report; writes nothing to the repo.

    audit_rest.py mouse|coral_rs|coral_denovo|pos|all
2026-09-11
"""
import csv, gzip, itertools, math, os, random, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"
MS   = os.path.join(ROOT, "methsplice")
HAS  = os.path.join(ROOT, "hassulta")
PIPE = os.path.join(ROOT, "pipeline")

def log(*a): print(*a, flush=True)

def find(name):
    for d in (MS, ROOT, HAS, PIPE):
        p = os.path.join(d, name)
        if os.path.exists(p): return p
    sys.exit("FATAL: cannot find " + name)

def load_psi(path, min_tot):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tot = float(r["inc"]) + float(r["skip"])
            if tot >= min_tot:
                out[(r["chrom"], int(r["start"]), int(r["end"]))] = (
                    float(r["inc"]) / tot, r["gene"])
    return out

def exon_meth(bed, exons, min_cpg=3):
    """Independent reimplementation of stream_exon_meth."""
    by_c = defaultdict(list)
    for c, s, e in exons: by_c[c].append((s, e))
    for c in by_c: by_c[c].sort()
    acc = defaultdict(lambda: [0.0, 0]); idx, active = {}, {}
    with gzip.open(bed, "rt") as fh:
        for line in fh:
            if line[0] == "#": continue
            p = line.split("\t", 5)
            c = p[0]; lst = by_c.get(c)
            if lst is None: continue
            pos = int(p[1]); score = float(p[3])
            if c not in idx: idx[c], active[c] = 0, []
            i = idx[c]; act = active[c]
            while i < len(lst) and lst[i][0] <= pos:
                act.append(lst[i]); i += 1
            idx[c] = i
            if act and act[0][1] <= pos:
                active[c] = act = [r for r in act if r[1] > pos]
            for s, e in act:
                if s <= pos < e:
                    a = acc[(c, s, e)]; a[0] += score; a[1] += 1
    return {k: v[0]/v[1] for k, v in acc.items() if v[1] >= min_cpg}

def slope_ci(x, y, gene_ids, nboot=2000, seed=7):
    """Gene-clustered bootstrap via per-gene sufficient statistics."""
    x = np.asarray(x, float); y = np.asarray(y, float); g = np.asarray(gene_ids)
    o = np.argsort(g, kind="stable"); x, y, g = x[o], y[o], g[o]
    _, st = np.unique(g, return_index=True); b = np.append(np.sort(st), len(g))
    def seg(v):
        cs = np.concatenate(([0.0], np.cumsum(v))); return cs[b[1:]] - cs[b[:-1]]
    n_g = np.diff(b).astype(float)
    Sx, Sy, Sxx, Sxy = seg(x), seg(y), seg(x*x), seg(x*y)
    def sl(n, sx, sy, sxx, sxy):
        den = sxx - sx*sx/n; num = sxy - sx*sy/n
        return np.where(den > 0, num/np.where(den > 0, den, 1.0), 0.0)
    obs = float(sl(n_g.sum(), Sx.sum(), Sy.sum(), Sxx.sum(), Sxy.sum()))
    rng = np.random.default_rng(seed); G = len(n_g)
    d = rng.integers(0, G, size=(nboot, G))
    bs = np.sort(sl(n_g[d].sum(1), Sx[d].sum(1), Sy[d].sum(1),
                    Sxx[d].sum(1), Sxy[d].sum(1)))
    return obs, bs[int(.025*nboot)], bs[int(.975*nboot)]

# ------------------------------------------------------------------ A: mouse
def part_mouse():
    log("\n" + "="*74); log("PART A  mouse contrasts recomputed from psi_mm"); log("="*74)
    stored = {r["contrast"]: r for r in
              csv.DictReader(open(find("contrasts_results.tsv")), delimiter="\t")}
    cfg = [(r["contrast"], r["group1"].split(","), r["group2"].split(","))
           for r in csv.DictReader(open(find("contrasts.tsv")), delimiter="\t")
           if r["contrast"].startswith("mouse_")]
    for name, g1, g2 in cfg:
        labs = g1 + g2; n1 = len(g1)
        data = {l: load_psi(os.path.join(MS, "psi_mm", l + ".psi.tsv"), 10) for l in labs}
        keys = sorted(set.intersection(*(set(d) for d in data.values())))
        V = np.array([[data[l][k][0] for l in labs] for k in keys])
        idx = list(range(len(labs)))
        def dp(o): o = list(o); return np.abs(V[:, o[:n1]].mean(1) - V[:, o[n1:]].mean(1))
        obs = dp(idx); n05, n10, nm = [], [], []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if list(comb) == idx[:n1] or rest == idx[:n1]: continue
            d = dp(list(comb) + rest)
            n05.append(int((d > .05).sum())); n10.append(int((d > .10).sum()))
            nm.append(float(d.mean()))
        o05 = int((obs > .05).sum()); o10 = int((obs > .10).sum())
        e05, e10 = float(np.mean(n05)), float(np.mean(n10))
        ge = sum(1 for v in nm if v >= obs.mean())
        st = stored[name]
        ok = (int(st["n_exons"]) == len(keys) and int(st["n_gt05_obs"]) == o05
              and abs(float(st["excess05"]) - o05/e05) < .015
              and int(st["n_gt10_obs"]) == o10
              and abs(float(st["mean_dpsi_obs"]) - obs.mean()) < 5e-5
              and abs(float(st["n_gt05_null"]) - e05) < .5)
        log(f"{name}")
        log(f"   n_exon      {len(keys):7d}   stored {st['n_exons']:>7s}")
        log(f"   mean|dPSI|  {obs.mean():.5f}   stored {st['mean_dpsi_obs']}"
            f"   null {np.mean(nm):.5f} stored {st['mean_dpsi_null']}")
        log(f"   >0.05       {o05:7d} vs {e05:8.1f}  x{o05/e05:.2f}"
            f"   stored {st['n_gt05_obs']} vs {st['n_gt05_null']} x{st['excess05']}")
        log(f"   >0.10       {o10:7d} vs {e10:8.1f}  x{o10/e10:.2f}"
            f"   stored {st['n_gt10_obs']} vs {st['n_gt10_null']}")
        log(f"   relabelings {len(nm)}, >= obs mean: {ge}   stored {st['relabelings_ge_obs']}")
        log(f"   ==> {'OK' if ok else '<-- MISMATCH'}")

# ------------------------------------------------------------------ B: coral
def part_coral(psi_dirname, stored_name, label):
    log("\n" + "="*74); log(f"PART B  coral recomputed: {label} ({psi_dirname})"); log("="*74)
    stored = {(r["test"], r["cell"]): r for r in
              csv.DictReader(open(find(stored_name)), delimiter="\t")}
    wgbs = {}
    for r in csv.DictReader(open(find("methsplice_runs_coral.tsv")), delimiter="\t"):
        if r["assay"] == "wgbs":
            wgbs[r["label"].removeprefix("apb_")] = r["run"]
    pdir = os.path.join(MS, psi_dirname)
    psi, gene_of, unpaired = {}, {}, []
    for f in sorted(os.listdir(pdir)):
        if not f.endswith(".psi.tsv"): continue
        code = f.removeprefix("apc_").removesuffix(".psi.tsv")
        if code not in wgbs: unpaired.append(code); continue
        d = load_psi(os.path.join(pdir, f), 10)
        psi[code] = {k: v[0] for k, v in d.items()}
        for k, v in d.items(): gene_of[k] = v[1]
    codes = sorted(psi)
    log(f"RNA codes with a WGBS mate: {len(codes)}  -> {' '.join(codes)}")
    log(f"RNA codes with no WGBS mate: {unpaired}")
    log(f"WGBS codes with no RNA mate: {sorted(set(wgbs) - set(codes))}")
    common = set.intersection(*(set(psi[c]) for c in codes))
    log(f"exons passing MIN_TOT=10 in all {len(codes)}: {len(common)}")
    meth = {}
    for c in codes:
        meth[c] = exon_meth(os.path.join(MS, "bs_ap", wgbs[c] + ".combined.bed.gz"), common)
        log(f"  {c} ({wgbs[c]}): {len(meth[c])}/{len(common)}")
    usable = sorted(k for k in common if all(k in meth[c] for c in codes))
    log(f"exons with >=3 covered CpGs in every fragment: {len(usable)}")

    cells = defaultdict(lambda: {"c": [], "t": []})
    for c in codes: cells[(c[:2], c[3:6])][c[2]].append(c)
    log("cell membership:")
    for k, v in sorted(cells.items()):
        log(f"   {k[0]}{k[1]}  control={v['c']}  treat={v['t']}"
            + ("   <-- DROPPED, one arm empty" if not v["c"] or not v["t"] else ""))
    genes = [gene_of[k] for k in usable]
    for (site, genet), grp in sorted(cells.items()):
        if not grp["c"] or not grp["t"]: continue
        dm = np.array([np.mean([meth[c][k] for c in grp["t"]]) -
                       np.mean([meth[c][k] for c in grp["c"]]) for k in usable])
        dp = np.array([np.mean([psi[c][k] for c in grp["t"]]) -
                       np.mean([psi[c][k] for c in grp["c"]]) for k in usable])
        s, lo, hi = slope_ci(dm, dp, genes)
        st = stored.get(("dose_response", f"{site}{genet}"))
        if st is None:
            log(f"  dose_response {site}{genet}: slope {s:+.9f}  <-- NOT IN STORED TABLE")
            continue
        ok = abs(s - float(st["value"])) < 1e-6 and len(usable) == int(st["n_exon"])
        log(f"  dose_response {site}{genet}: slope {s:+.9f} [{lo:+.6f},{hi:+.6f}] "
            f"n={len(usable)}")
        log(f"                stored  {st['value']} [{st['ci_lo']},{st['ci_hi']}] "
            f"n={st['n_exon']}   {'OK' if ok else '<-- MISMATCH'}")

    M = np.array([[meth[c][k] for c in codes] for k in usable])
    P = np.array([[psi[c][k] for c in codes] for k in usable])
    def rowr(A, B):
        a = A - A.mean(1, keepdims=True); b = B - B.mean(1, keepdims=True)
        den = np.sqrt((a*a).sum(1) * (b*b).sum(1))
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (a*b).sum(1) / den
        return r[np.isfinite(r)]
    r = rowr(M, P); st = stored[("within_exon_r", "ALL")]
    ok = abs(r.mean() - float(st["value"])) < 1e-4 and len(r) == int(st["n_exon"])
    log(f"  within_exon_r: mean r {r.mean():+.4f} over {len(r)} exons "
        f"  stored {st['value']} n={st['n_exon']}  {'OK' if ok else '<-- MISMATCH'}")
    rng = np.random.default_rng(7); n_us = len(usable)
    for tag, nsamp in (("as-published (2,000-exon subsample)", 2000),
                       ("correctly matched (all exons)", n_us)):
        nulls = []
        for _ in range(200):
            perm = rng.permutation(len(codes))
            sel = rng.choice(n_us, size=min(nsamp, n_us), replace=False) \
                  if nsamp < n_us else slice(None)
            nulls.append(rowr(M[sel], P[sel][:, perm]).mean())
        nulls = np.sort(np.array(nulls))
        log(f"    null, {tag}: mean {nulls.mean():+.4f} "
            f"band [{nulls[4]:+.4f}, {nulls[194]:+.4f}] width {nulls[194]-nulls[4]:.4f}")
    log(f"    stored row reports band [{st['ci_lo']}, {st['ci_hi']}] in the "
        f"ci_lo/ci_hi COLUMNS and null mean {st['null']}")

    mm = M.mean(1); pm = P.mean(1)
    for cut, sel in (("meth>=50", mm >= 50), ("meth<10", mm < 10)):
        frac = float((pm[sel] > .98).mean()); st = stored[("ceiling", cut)]
        ok = abs(frac - float(st["value"])) < 2e-4 and int(sel.sum()) == int(st["n_exon"])
        log(f"  ceiling {cut}: {frac:.4f} n={int(sel.sum())}  stored {st['value']} "
            f"n={st['n_exon']}  {'OK' if ok else '<-- MISMATCH'}")
    log(f"  exons in neither band (10<=meth<50): "
        f"{int(((mm >= 10) & (mm < 50)).sum())}")

# ----------------------------------------------------------- C: pos control
def part_pos():
    log("\n" + "="*74); log("PART C  positive control recomputed"); log("="*74)
    MIN_TOT, DPSI = 20, 0.10

    def contrast(pdir, g1, g2, name, stored, key_obs, key_med):
        labs = g1 + g2; n1 = len(g1)
        psi = {l: {k: v[0] for k, v in
                   load_psi(os.path.join(pdir, l + ".psi.tsv"), MIN_TOT).items()}
               for l in labs}
        ks = sorted(set.intersection(*(set(psi[l]) for l in labs)))
        V = np.array([[psi[l][k] for l in labs] for k in ks])
        idx = list(range(len(labs)))
        def cnt(o):
            o = list(o)
            return int((np.abs(V[:, o[:n1]].mean(1) - V[:, o[n1:]].mean(1)) >= DPSI).sum())
        obs = cnt(idx); nulls = []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if set(comb) == set(idx[:n1]) or set(rest) == set(idx[:n1]): continue
            nulls.append(cnt(list(comb) + rest))
        nulls.sort(); med = nulls[len(nulls)//2]
        p = (1 + sum(1 for v in nulls if v >= obs)) / (1 + len(nulls))
        mean_null = float(np.mean(nulls))
        ok = (obs == int(stored[key_obs]) and med == int(stored[key_med]))
        log(f"{name}: n_exon={len(ks)} obs={obs} null_med={med} null_mean={mean_null:.1f} "
            f"null_max={nulls[-1]} splits={len(nulls)} p={p:.4f}")
        log(f"    ratio vs MEDIAN {obs/med:.3f}   vs MEAN {obs/mean_null:.3f}"
            f"   stored {stored}")
        log(f"    ==> {'OK' if ok else '<-- MISMATCH'}")

    dev = list(csv.DictReader(open(find("runs_dev.tsv")), delimiter="\t"))
    by = defaultdict(list)
    for r in dev: by[(r["stage"], r["sex"])].append(r["label"])
    L4 = by[("4_instar_larva", "NA")]; L5 = by[("5_instar_larva", "NA")]
    PM, PF = by[("pupa", "M")], by[("pupa", "F")]
    AM, AF = by[("adult", "M")], by[("adult", "F")]
    sd = {r["contrast"]: r for r in
          csv.DictReader(open(find("dev_contrast.tsv")), delimiter="\t")}
    pdir = os.path.join(HAS, "psi")
    for nm, g1, g2 in (("larva4_vs_larva5", L4, L5),
                       ("larva5_vs_pupa", L5, PM + PF),
                       ("pupa_vs_adult", PM + PF, AM + AF),
                       ("pupaM_vs_pupaF", PM, PF),
                       ("adultM_vs_adultF", AM, AF)):
        s = sd[nm]
        contrast(pdir, g1, g2, nm,
                 {"observed": s["observed"], "null_median": s["null_median"],
                  "n_exon": s["n_exon"], "n_splits": s["n_splits"],
                  "p": s["p_exact"], "enr": s["enrichment"]},
                 "observed", "null_median")

    hm = list(csv.DictReader(open(find("runs_harm.tsv")), delimiter="\t"))
    hl = [r["label"] for r in hm if r["stage"] == "larva"]
    ha = [r["label"] for r in hm if r["stage"] == "adult"]
    hs = {r["test"]: r for r in
          csv.DictReader(open(find("harm_contrast.tsv")), delimiter="\t")}["metamorphosis_splicing"]
    contrast(os.path.join(HAS, "psi_harm"), hl, ha, "harm metamorphosis",
             {"observed": hs["value"], "null_median": hs["null_median"],
              "ratio": hs["ratio"], "p": hs["p"]}, "observed", "null_median")

def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "mouse"): part_mouse()
    if what in ("all", "coral_rs"):
        part_coral("psi_apc_rs", "coral_contrast_rs.tsv", "RefSeq (PRIMARY)")
    if what in ("all", "coral_denovo"):
        part_coral("psi_apc", "coral_contrast_denovo.tsv", "StringTie denovo")
    if what in ("all", "pos"): part_pos()
    log("\nAUDIT COMPLETE")

if __name__ == "__main__":
    main()
