#!/usr/bin/env python3
"""Size-matched top-K variant of the cross-species overlap test, on the real
Helicoverpa data (methods paper, completes the simulation-grid arc).

    harm_topk.py  (same arguments as harm_overlap_null.py)

The simulation grid showed that under within-species label permutation,
threshold-defined changing sets shrink (the permutation destroys real
within-species regulation), and any statistic sensitive to set size fails:
the enrichment ratio degenerates low (87% FPR in the collapse regime) and
the permuted hypergeometric p pins at ~1 (100% FPR). The two controls that
hold grid-wide are rank correlation (harm_spearman.py, already run:
p=0.69) and this one: rank genes by their maximum per-exon |dPSI|, take the
top K in each species with K FIXED at the observed threshold-set size, and
score the overlap of those size-matched sets against the same overlap
recomputed under label permutation. Permuted sets have exactly the same
sizes by construction, so the marginal-shrink pathology cannot occur.

Imports all data handling from harm_overlap_null.py so both tests read the
data identically.
"""
import argparse, csv, itertools, os, random, re, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import harm_overlap_null as base


def gene_mags(psi, ks, gene_of, g1, g2, gid_of):
    mags = {}
    for k in ks:
        a = sum(psi[l][k] for l in g1) / len(g1)
        b = sum(psi[l][k] for l in g2) / len(g2)
        g = gid_of(gene_of[k])
        if g is not None:
            d = abs(a - b)
            if d > mags.get(g, -1.0):
                mags[g] = d
    return mags


def topk(mags, shared, K):
    ranked = sorted(((m, g) for g, m in mags.items() if g in shared),
                    reverse=True)
    return {g for _, g in ranked[:K]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi-assulta", required=True)
    ap.add_argument("--runs-assulta", required=True)
    ap.add_argument("--psi-harm", required=True)
    ap.add_argument("--runs-harm", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--k-factors", default="1",
                    help="comma list; each scales the observed K (e.g. 0.5,1,2)")
    a = ap.parse_args()

    to_gid = base.gff_maps(a.gff)
    best = {}
    with open(a.paf) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 12 or int(p[11]) < 10:
                continue
            q, t, nm = p[0], p[5], int(p[9])
            g = q.rsplit(".", 1)[0] if q.startswith("MSTRG") else q
            if g not in best or nm > best[g][1]:
                best[g] = (t, nm)
    a2gid = {}
    for g, (t, _) in best.items():
        gid = to_gid.get(t) or to_gid.get(re.sub(r"^rna-", "", t))
        if gid:
            a2gid[g] = gid

    am = list(csv.DictReader(open(a.runs_assulta), delimiter="\t"))
    aby = defaultdict(list)
    for r in am:
        aby[(r["stage"], r["sex"])].append(r["label"])
    a_g1 = aby[("pupa", "M")] + aby[("pupa", "F")]
    a_g2 = aby[("adult", "M")] + aby[("adult", "F")]
    apsi, agene = base.load_psi(a.psi_assulta, a_g1 + a_g2)
    aks = sorted(set.intersection(*(set(apsi[l]) for l in a_g1 + a_g2)))

    hm = list(csv.DictReader(open(a.runs_harm), delimiter="\t"))
    h_g1 = [r["label"] for r in hm if r["stage"] == "larva"]
    h_g2 = [r["label"] for r in hm if r["stage"] == "adult"]
    hpsi, hgene = base.load_psi(a.psi_harm, h_g1 + h_g2)
    hks = sorted(set.intersection(*(set(hpsi[l]) for l in h_g1 + h_g2)))

    a_gid = lambda g: a2gid.get(g)
    h_gid = lambda g: to_gid.get(g) or to_gid.get(re.sub(r"^(gene|rna)-", "", g))
    a_bg = {a_gid(agene[k]) for k in aks} - {None}
    h_bg = {h_gid(hgene[k]) for k in hks} - {None}
    shared = a_bg & h_bg
    print(f"shared testable background: {len(shared)} genes")

    # K fixed at the observed threshold-set sizes (same definition as the
    # published test), then everything below is size-matched.
    obs_a_thr = base.changing_genes(apsi, aks, agene, a_g1, a_g2, a_gid) & shared
    obs_h_thr = base.changing_genes(hpsi, hks, hgene, h_g1, h_g2, h_gid) & shared
    Ka, Kh = len(obs_a_thr), len(obs_h_thr)
    print(f"K fixed at observed sizes: assulta {Ka}, armigera {Kh}")

    factors = [float(x) for x in a.k_factors.split(",")]
    obs_mag_a = gene_mags(apsi, aks, agene, a_g1, a_g2, a_gid)
    obs_mag_h = gene_mags(hpsi, hks, hgene, h_g1, h_g2, h_gid)

    a_pool, h_pool = a_g1 + a_g2, h_g1 + h_g2
    a_splits = list(itertools.combinations(range(len(a_pool)), len(a_g1)))
    h_splits = list(itertools.combinations(range(len(h_pool)), len(h_g1)))
    rng = random.Random(23)
    Ks = [(max(1, int(round(Ka * f))), max(1, int(round(Kh * f)))) for f in factors]
    obs_by_f = {}
    for f, (ka, kh) in zip(factors, Ks):
        oa = topk(obs_mag_a, shared, ka); oh = topk(obs_mag_h, shared, kh)
        obs_by_f[f] = len(oa & oh)
        print(f"K-factor {f:g} (K={ka}/{kh}): observed overlap {obs_by_f[f]} "
              f"(independence expectation {ka*kh/len(shared):.1f})")
    nulls_by_f = {f: [] for f in factors}
    for _ in range(a.n_perm):
        while True:
            sa = rng.choice(a_splits)
            g1 = [a_pool[i] for i in sa]
            if set(g1) != set(a_g1) and set(g1) != set(a_g2):
                break
        g2 = [a_pool[i] for i in range(len(a_pool)) if i not in sa]
        pm_a = gene_mags(apsi, aks, agene, g1, g2, a_gid)
        while True:
            sh = rng.choice(h_splits)
            k1 = [h_pool[i] for i in sh]
            if set(k1) != set(h_g1) and set(k1) != set(h_g2):
                break
        k2 = [h_pool[i] for i in range(len(h_pool)) if i not in sh]
        pm_h = gene_mags(hpsi, hks, hgene, k1, k2, h_gid)
        for f, (ka, kh) in zip(factors, Ks):
            ca = topk(pm_a, shared, ka); ch = topk(pm_h, shared, kh)
            nulls_by_f[f].append(len(ca & ch))
    with open(a.out, "w") as o:
        o.write("k_factor\tK_assulta\tK_armigera\tobserved\texpectation\t"
                "null_median\tnull_max\tp\n")
        for f, (ka, kh) in zip(factors, Ks):
            nl = sorted(nulls_by_f[f])
            obs = obs_by_f[f]
            med = nl[len(nl) // 2]
            p = (1 + sum(1 for v in nl if v >= obs)) / (1 + len(nl))
            o.write(f"{f:g}\t{ka}\t{kh}\t{obs}\t{ka*kh/len(shared):.2f}\t"
                    f"{med}\t{nl[-1]}\t{p:.5f}\n")
            with open(f"{a.out}.nulls.f{f:g}", "w") as nf:
                nf.write("null_topk_overlap\n")
                for v in nulls_by_f[f]:
                    nf.write(f"{v}\n")
            print(f"K-factor {f:g}: obs {obs} vs null median {med} "
                  f"(max {nl[-1]}) -> p = {p:.5f}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
