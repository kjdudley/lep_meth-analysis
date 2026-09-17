#!/usr/bin/env python3
"""Is the cross-species gene overlap real, or manufactured by shared exon
properties?

    harm_overlap_null.py --psi-assulta hassulta/psi --runs-assulta ... \
        --psi-harm hassulta/psi_harm --runs-harm ... --gff ... --paf ... \
        --out hassulta/harm_overlap_null.tsv

The concern the hypergeometric cannot address: exons sitting at intermediate
PSI have more room to move, so they exceed |dPSI| >= 0.1 more readily in
BOTH species. If that propensity is a conserved property of a gene, the
overlap inflates with no conserved REGULATION behind it.

The fix is a null that keeps every exon property and destroys only the
developmental structure: **relabel the samples within each species**
(6v6 -> a random non-identity 6/6 split, from the same 924), recompute each
species' changing set, and re-measure the overlap. Under this null the PSI
distributions, depths, gene sets, annotation quality and detectability are
all exactly as observed; only which samples count as "pupa/larva" versus
"adult" changes. If the observed overlap survives, it is developmental
regulation, not exon geometry.

Both species are permuted independently, so the null also destroys any
shared batch structure.
"""
import argparse, csv, itertools, math, os, random, re, sys
from collections import defaultdict

MIN_TOT = 20
DPSI = 0.10


def load_psi(psi_dir, labels):
    psi, gene_of = {}, {}
    for l in labels:
        d = {}
        p = os.path.join(psi_dir, f"{l}.psi.tsv")
        if not os.path.exists(p):
            sys.exit(f"REFUSING: missing {p}")
        for r in csv.DictReader(open(p), delimiter="\t"):
            tot = float(r["inc"]) + float(r["skip"])
            if tot < MIN_TOT:
                continue
            k = (r["chrom"], int(r["start"]), int(r["end"]))
            d[k] = float(r["inc"]) / tot
            gene_of[k] = r["gene"]
        psi[l] = d
    return psi, gene_of


def changing_genes(psi, ks, gene_of, g1, g2, gid_of):
    out = set()
    for k in ks:
        a = sum(psi[l][k] for l in g1) / len(g1)
        b = sum(psi[l][k] for l in g2) / len(g2)
        if abs(a - b) >= DPSI:
            g = gid_of(gene_of[k])
            if g:
                out.add(g)
    return out


def gff_maps(gff):
    to_gid = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] not in ("gene", "mRNA", "transcript"):
                continue
            gid = re.search(r"GeneID:(\d+)", p[8])
            if not gid:
                continue
            i = re.search(r"(?:^|;)ID=([^;]+)", p[8])
            nm = re.search(r"(?:^|;)gene=([^;]+)", p[8])
            keys = []
            if i:
                keys += [i.group(1), re.sub(r"^(gene|rna)-", "", i.group(1))]
            if nm:
                keys.append(nm.group(1))
            for key in keys:
                to_gid[key] = gid.group(1)
    return to_gid


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
    a = ap.parse_args()

    to_gid = gff_maps(a.gff)
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
    apsi, agene = load_psi(a.psi_assulta, a_g1 + a_g2)
    aks = sorted(set.intersection(*(set(apsi[l]) for l in a_g1 + a_g2)))

    hm = list(csv.DictReader(open(a.runs_harm), delimiter="\t"))
    h_g1 = [r["label"] for r in hm if r["stage"] == "larva"]
    h_g2 = [r["label"] for r in hm if r["stage"] == "adult"]
    hpsi, hgene = load_psi(a.psi_harm, h_g1 + h_g2)
    hks = sorted(set.intersection(*(set(hpsi[l]) for l in h_g1 + h_g2)))

    a_gid = lambda g: a2gid.get(g)
    h_gid = lambda g: to_gid.get(g) or to_gid.get(re.sub(r"^(gene|rna)-", "", g))
    a_bg = {a_gid(agene[k]) for k in aks} - {None}
    h_bg = {h_gid(hgene[k]) for k in hks} - {None}
    shared = a_bg & h_bg
    print(f"shared testable background: {len(shared)} genes")

    obs_a = changing_genes(apsi, aks, agene, a_g1, a_g2, a_gid) & shared
    obs_h = changing_genes(hpsi, hks, hgene, h_g1, h_g2, h_gid) & shared
    obs = len(obs_a & obs_h)
    print(f"observed: assulta {len(obs_a)}, armigera {len(obs_h)}, "
          f"overlap {obs}")

    a_pool, h_pool = a_g1 + a_g2, h_g1 + h_g2
    a_splits = [s for s in itertools.combinations(range(len(a_pool)), len(a_g1))]
    h_splits = [s for s in itertools.combinations(range(len(h_pool)), len(h_g1))]
    rng = random.Random(23)
    nulls_raw, sizes = [], []
    for _ in range(a.n_perm):
        while True:
            sa = rng.choice(a_splits)
            g1 = [a_pool[i] for i in sa]
            if set(g1) != set(a_g1) and set(g1) != set(a_g2):
                break
        g2 = [a_pool[i] for i in range(len(a_pool)) if i not in sa]
        ca = changing_genes(apsi, aks, agene, g1, g2, a_gid) & shared
        while True:
            sh = rng.choice(h_splits)
            k1 = [h_pool[i] for i in sh]
            if set(k1) != set(h_g1) and set(k1) != set(h_g2):
                break
        k2 = [h_pool[i] for i in range(len(h_pool)) if i not in sh]
        ch = changing_genes(hpsi, hks, hgene, k1, k2, h_gid) & shared
        nulls_raw.append(len(ca & ch))
        sizes.append((len(ca), len(ch)))
    # RAW overlap is confounded by set size: the true contrast yields far
    # bigger changing sets than a relabeling, and bigger sets overlap more
    # mechanically. The statistic must be ENRICHMENT OVER EACH SPLIT'S OWN
    # independent-draw expectation. (Found 2026-08-27: the raw comparison
    # gave x4.82 / p=0.001 while the enrichment comparison does not.)
    Nsh = len(shared)
    obs_exp = len(obs_a) * len(obs_h) / Nsh
    obs_enr = obs / obs_exp if obs_exp else float("nan")
    enrs = []
    for (sa_n, sh_n), ovl in zip(sizes, nulls_raw):
        e = sa_n * sh_n / Nsh
        if e > 0:
            enrs.append(ovl / e)
    enrs.sort()
    p_enr = (1 + sum(1 for v in enrs if v >= obs_enr)) / (1 + len(enrs))
    print(f"ENRICHMENT over own expectation: observed x{obs_enr:.2f} "
          f"(={obs}/{obs_exp:.1f}); null median x{enrs[len(enrs)//2]:.2f}, "
          f"max x{enrs[-1]:.2f}; p = {p_enr:.4f}")
    nulls = sorted(nulls_raw)
    med = nulls[len(nulls) // 2]
    p = (1 + sum(1 for v in nulls if v >= obs)) / (1 + len(nulls))
    ma = sum(s[0] for s in sizes) / len(sizes)
    mh = sum(s[1] for s in sizes) / len(sizes)
    print(f"null overlap: median {med}, max {nulls[-1]}, "
          f"mean set sizes assulta {ma:.0f} / armigera {mh:.0f}")
    print(f"observed {obs} vs null median {med}  ->  x{obs/med if med else float('inf'):.2f}"
          f"  p = {p:.4f}  ({len(nulls)} permutations)")
    with open(a.out, "w") as o:
        o.write("statistic\tvalue\n")
        o.write(f"shared_background\t{len(shared)}\n")
        o.write(f"observed_assulta_changing\t{len(obs_a)}\n")
        o.write(f"observed_armigera_changing\t{len(obs_h)}\n")
        o.write(f"observed_overlap\t{obs}\n")
        o.write(f"null_median_overlap\t{med}\n")
        o.write(f"null_max_overlap\t{nulls[-1]}\n")
        o.write(f"null_mean_assulta_set\t{ma:.1f}\n")
        o.write(f"null_mean_armigera_set\t{mh:.1f}\n")
        o.write(f"ratio_vs_null\t{obs/med if med else ''}\n")
        o.write(f"p_permutation_RAW_confounded\t{p:.4f}\n")
        o.write(f"observed_enrichment_over_expectation\t{obs_enr:.4f}\n")
        o.write(f"null_median_enrichment\t{enrs[len(enrs)//2]:.4f}\n")
        o.write(f"null_max_enrichment\t{enrs[-1]:.4f}\n")
        o.write(f"p_PRIMARY_enrichment\t{p_enr:.4f}\n")
        o.write(f"n_permutations\t{len(nulls)}\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
