#!/usr/bin/env python3
"""Audit of the two claims the METHSPLICE manuscript still makes positively:
the novel-junction asymmetry (S3 of ms1_cryptic.py) and the coral pairing.

  audit_cryptic.py s3      novel-junction asymmetry: depth, null, all contrasts
  audit_cryptic.py coral   re-pair the coral arm with the .1 suffix normalised

S3 as published is a raw count ratio: junctions at >=5 reads in EVERY
replicate of one group and NO replicate of the other. It carries no
relabeling null and no depth normalisation, unlike S1 and S2 in the same
script. Both are supplied here.  2026-09-11
"""
import csv, gzip, itertools, os, sys
import numpy as np
from collections import defaultdict

ROOT = "/work/cyberomics/lep_meth"
MS   = os.path.join(ROOT, "methsplice")
JUNC = os.path.join(ROOT, "rna", "junc")
PIPE = os.path.join(ROOT, "pipeline")

def log(*a): print(*a, flush=True)

def find(name):
    for d in (MS, ROOT, PIPE):
        p = os.path.join(d, name)
        if os.path.exists(p): return p
    sys.exit("FATAL: cannot find " + name)

# --------------------------------------------------------------- S3 audit
def junc_table(label):
    """(dict junction->reads, total reads over all junctions)"""
    d = {}; tot = 0.0
    with open(os.path.join(JUNC, label + ".junctions.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try: n = float(r["n_reads"])
            except (ValueError, KeyError): continue
            d[(r["chrom"], r["intron_start"], r["intron_end"])] = n
            tot += n
    return d, tot

def gain_loss(sets_g1, sets_g2):
    g = len(set.intersection(*sets_g1) - set.union(*sets_g2))
    l = len(set.intersection(*sets_g2) - set.union(*sets_g1))
    return g, l

def part_s3():
    log("="*78); log("S3 NOVEL-JUNCTION ASYMMETRY: depth, calibration, direction")
    log("="*78)
    cfg = {r["contrast"]: (r["group1"].split(","), r["group2"].split(","))
           for r in csv.DictReader(open(find("contrasts.tsv")), delimiter="\t")}
    stored = {}
    for f in ("cryptic.tsv", "cryptic_degron.tsv"):
        for r in csv.DictReader(open(find(f)), delimiter="\t"):
            stored[r["contrast"]] = r
    rows = []
    for name in [c for c in cfg if c in stored]:
        g1, g2 = cfg[name]; labs = g1 + g2; n1 = len(g1)
        try:
            tabs = {l: junc_table(l) for l in labs}
        except FileNotFoundError as e:
            log(f"{name}: missing junction file ({e.filename}); skipped"); continue
        depth = {l: tabs[l][1] for l in labs}
        nj5 = {l: sum(1 for v in tabs[l][0].values() if v >= 5) for l in labs}
        sets5 = {l: {k for k, v in tabs[l][0].items() if v >= 5} for l in labs}
        g, l = gain_loss([sets5[x] for x in g1], [sets5[x] for x in g2])
        st = stored[name]
        ok = (g == int(st["novel_g1"]) and l == int(st["novel_g2"]))
        d1 = np.mean([depth[x] for x in g1]); d2 = np.mean([depth[x] for x in g2])
        j1 = np.mean([nj5[x] for x in g1]);   j2 = np.mean([nj5[x] for x in g2])

        # exact balanced relabeling null for the ratio
        idx = list(range(len(labs))); nulls = []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if list(comb) == idx[:n1] or rest == idx[:n1]: continue
            a = [sets5[labs[i]] for i in comb]; b = [sets5[labs[i]] for i in rest]
            gg, ll = gain_loss(a, b)
            nulls.append(gg / max(ll, 1))
        nulls = np.sort(np.array(nulls))
        obs_ratio = g / max(l, 1)
        p = (1 + int((nulls >= obs_ratio).sum())) / (1 + len(nulls))

        # depth-normalised: rescale each library's read counts to the shallowest
        ref = min(depth[x] for x in labs)
        setsN = {x: {k for k, v in tabs[x][0].items() if v * ref / depth[x] >= 5}
                 for x in labs}
        gN, lN = gain_loss([setsN[x] for x in g1], [setsN[x] for x in g2])
        nullsN = []
        for comb in itertools.combinations(idx, n1):
            rest = [i for i in idx if i not in comb]
            if list(comb) == idx[:n1] or rest == idx[:n1]: continue
            gg, ll = gain_loss([setsN[labs[i]] for i in comb],
                               [setsN[labs[i]] for i in rest])
            nullsN.append(gg / max(ll, 1))
        nullsN = np.sort(np.array(nullsN))
        obsN = gN / max(lN, 1)
        pN = (1 + int((nullsN >= obsN).sum())) / (1 + len(nullsN))

        log(f"\n{name}   ({n1}v{len(g2)})")
        log(f"  reproduce stored : {g}/{l} = x{obs_ratio:.2f}   "
            f"stored {st['novel_g1']}/{st['novel_g2']} x{st['novel_ratio']}  "
            f"{'OK' if ok else '<-- MISMATCH'}")
        log(f"  library depth    : g1 {d1:,.0f} junction-reads, g2 {d2:,.0f}  "
            f"ratio g1/g2 = {d1/d2:.2f}")
        log(f"  junctions >=5rd  : g1 {j1:,.0f}, g2 {j2:,.0f}  ratio {j1/j2:.2f}")
        log(f"  relabeling null  : median x{np.median(nulls):.2f}  "
            f"range x{nulls[0]:.2f}-x{nulls[-1]:.2f}  ({len(nulls)} splits)  "
            f"p={p:.3f}")
        log(f"  depth-normalised : {gN}/{lN} = x{obsN:.2f}   null median "
            f"x{np.median(nullsN):.2f} range x{nullsN[0]:.2f}-x{nullsN[-1]:.2f}  p={pN:.3f}")
        rows.append((name, obs_ratio, d1/d2, j1/j2, obsN, p, pN))

    log("\n" + "="*78); log("SUMMARY: does the asymmetry track demethylation, or depth?")
    log("="*78)
    log(f"{'contrast':34s}{'x pub':>8s}{'depth g1/g2':>13s}{'j>=5 g1/g2':>12s}"
        f"{'x norm':>9s}{'p raw':>8s}{'p norm':>8s}")
    for r in rows:
        log(f"{r[0]:34s}{r[1]:8.2f}{r[2]:13.2f}{r[3]:12.2f}{r[4]:9.2f}{r[5]:8.3f}{r[6]:8.3f}")
    if len(rows) > 2:
        a = np.log([r[1] for r in rows]); b = np.log([r[3] for r in rows])
        cc = float(np.corrcoef(a, b)[0, 1])
        log(f"\ncorrelation of log(published ratio) with log(junction-count ratio): "
            f"r = {cc:+.3f} over {len(rows)} contrasts")
        n_fav = sum(1 for r in rows if r[1] > 1)
        log(f"published ratio favours the demethylated arm in {n_fav}/{len(rows)} contrasts")
        n_favN = sum(1 for r in rows if r[4] > 1)
        log(f"depth-normalised ratio favours it in {n_favN}/{len(rows)}")

# ------------------------------------------------------------- coral repair
def load_psi(path, min_tot=10):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tot = float(r["inc"]) + float(r["skip"])
            if tot >= min_tot:
                out[(r["chrom"], int(r["start"]), int(r["end"]))] = (
                    float(r["inc"]) / tot, r["gene"])
    return out

def exon_meth(bed, exons, min_cpg=3):
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

def part_coral():
    log("="*78)
    log("CORAL RE-PAIRED: technical-replicate suffix normalised")
    log("="*78)
    wgbs = {}
    for r in csv.DictReader(open(find("methsplice_runs_coral.tsv")), delimiter="\t"):
        if r["assay"] == "wgbs":
            wgbs[r["label"].removeprefix("apb_")] = r["run"]
    for psi_dirname, tag in (("psi_apc_rs", "RefSeq (PRIMARY)"),
                             ("psi_apc", "StringTie denovo")):
        log(f"\n---- {tag} ----")
        pdir = os.path.join(MS, psi_dirname)
        # code -> list of psi files (a fragment may have >1 RNA run)
        frag = defaultdict(list)
        for f in sorted(os.listdir(pdir)):
            if not f.endswith(".psi.tsv"): continue
            code = f.removeprefix("apc_").removesuffix(".psi.tsv")
            base = code[:8] if len(code) > 8 else code
            frag[base].append(f)
        paired = {c: fs for c, fs in frag.items() if c in wgbs}
        log(f"fragments after suffix normalisation: {len(paired)} "
            f"(published pairing used {sum(1 for c in frag if len(c)==8 and c in wgbs)})")
        for c, fs in sorted(paired.items()):
            if len(fs) > 1: log(f"   {c}: {len(fs)} RNA runs {fs}")
        psi, gene_of = {}, {}
        for c, fs in paired.items():
            # a fragment with two runs: average PSI over runs that pass MIN_TOT
            acc = defaultdict(list)
            for f in fs:
                for k, v in load_psi(os.path.join(pdir, f)).items():
                    acc[k].append(v[0]); gene_of[k] = v[1]
            psi[c] = {k: float(np.mean(v)) for k, v in acc.items()}
        codes = sorted(psi)
        common = set.intersection(*(set(psi[c]) for c in codes))
        log(f"exons passing MIN_TOT=10 in all {len(codes)}: {len(common)}")
        meth = {}
        for c in codes:
            meth[c] = exon_meth(os.path.join(MS, "bs_ap",
                                             wgbs[c] + ".combined.bed.gz"), common)
        usable = sorted(k for k in common if all(k in meth[c] for c in codes))
        log(f"exons with >=3 covered CpGs in every fragment: {len(usable)}")
        genes = [gene_of[k] for k in usable]
        cells = defaultdict(lambda: {"c": [], "t": []})
        for c in codes: cells[(c[:2], c[3:6])][c[2]].append(c)
        for (site, genet), grp in sorted(cells.items()):
            if not grp["c"] or not grp["t"]:
                log(f"  {site}{genet}: STILL DROPPED  c={grp['c']} t={grp['t']}")
                continue
            dm = np.array([np.mean([meth[c][k] for c in grp["t"]]) -
                           np.mean([meth[c][k] for c in grp["c"]]) for k in usable])
            dp = np.array([np.mean([psi[c][k] for c in grp["t"]]) -
                           np.mean([psi[c][k] for c in grp["c"]]) for k in usable])
            s, lo, hi = slope_ci(dm, dp, genes)
            log(f"  dose_response {site}{genet} ({len(grp['t'])}t v {len(grp['c'])}c): "
                f"slope {s:+.9f} [{lo:+.6f}, {hi:+.6f}]  n={len(usable)}")

def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "s3"
    if what == "s3": part_s3()
    elif what == "coral": part_coral()
    else: sys.exit("usage: audit_cryptic.py s3|coral")
    log("\nDONE")

if __name__ == "__main__":
    main()
