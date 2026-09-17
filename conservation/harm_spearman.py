#!/usr/bin/env python3
"""The threshold-free conservation test: do orthologous genes change by the
same AMOUNT across metamorphosis in the two Helicoverpa species?

    harm_spearman.py --psi-assulta ... --psi-harm ... --gff ... --paf ... \
        --out hassulta/harm_spearman.tsv

Why this replaces the overlap test. Binarising at |dPSI| >= 0.1 created two
artefacts at once: set-size confounding (a real contrast yields bigger
changing sets than a relabeling, and bigger sets overlap more) and shared
exon propensity (exons at intermediate PSI have room to move in BOTH
species). Together they produced x2.09 at p = 6.5e-28 for what the proper
control showed to be a complete null (p = 0.81). See DESIGN.md retraction.

A rank correlation of the CONTINUOUS per-gene change magnitude has neither
problem: no threshold, no sets, no size to confound.

Gene-level statistic: mean |dPSI| over the gene's tested exons. Direction
is deliberately NOT compared -- without exon-level orthology between the
species, the sign of a given exon's shift is not comparable across them,
and pretending otherwise would invent a correspondence that the data does
not establish.

Null: within-species relabeling (each species permuted independently from
its own 924 splits), recomputing both magnitude vectors and the
correlation. This absorbs everything static about the genes -- exon count,
expression, depth, noise, annotation quality -- because under relabeling
those are unchanged and only the developmental grouping is destroyed.
"""
import argparse, csv, itertools, math, os, random, re, sys
from collections import defaultdict

MIN_TOT = 20


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


def gene_magnitude(psi, ks, gene_of, g1, g2, gid_of):
    acc = defaultdict(lambda: [0.0, 0])
    for k in ks:
        g = gid_of(gene_of[k])
        if not g:
            continue
        a = sum(psi[l][k] for l in g1) / len(g1)
        b = sum(psi[l][k] for l in g2) / len(g2)
        acc[g][0] += abs(a - b); acc[g][1] += 1
    return {g: v[0] / v[1] for g, v in acc.items()}


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = float(pos)
        # average ties
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            if j > i:
                avg = sum(r[order[t]] for t in range(i, j + 1)) / (j - i + 1)
                for t in range(i, j + 1):
                    r[order[t]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(rx)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    return num / den if den else 0.0


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

    ma = gene_magnitude(apsi, aks, agene, a_g1, a_g2, a_gid)
    mh = gene_magnitude(hpsi, hks, hgene, h_g1, h_g2, h_gid)
    shared = sorted(set(ma) & set(mh))
    print(f"ortholog pairs with a magnitude in both species: {len(shared)}")
    if len(shared) < 200:
        sys.exit("REFUSING: too few ortholog pairs")
    obs = spearman([ma[g] for g in shared], [mh[g] for g in shared])
    print(f"observed Spearman(|dPSI| assulta, |dPSI| armigera) = {obs:+.4f}")

    a_pool, h_pool = a_g1 + a_g2, h_g1 + h_g2
    a_splits = list(itertools.combinations(range(len(a_pool)), len(a_g1)))
    h_splits = list(itertools.combinations(range(len(h_pool)), len(h_g1)))
    rng = random.Random(31)
    nulls = []
    for _ in range(a.n_perm):
        while True:
            sa = rng.choice(a_splits)
            g1 = [a_pool[i] for i in sa]
            if set(g1) != set(a_g1) and set(g1) != set(a_g2):
                break
        g2 = [a_pool[i] for i in range(len(a_pool)) if i not in sa]
        while True:
            sh = rng.choice(h_splits)
            k1 = [h_pool[i] for i in sh]
            if set(k1) != set(h_g1) and set(k1) != set(h_g2):
                break
        k2 = [h_pool[i] for i in range(len(h_pool)) if i not in sh]
        na = gene_magnitude(apsi, aks, agene, g1, g2, a_gid)
        nh = gene_magnitude(hpsi, hks, hgene, k1, k2, h_gid)
        sh_g = [g for g in shared if g in na and g in nh]
        nulls.append(spearman([na[g] for g in sh_g], [nh[g] for g in sh_g]))
    nulls.sort()
    med = nulls[len(nulls) // 2]
    p = (1 + sum(1 for v in nulls if v >= obs)) / (1 + len(nulls))
    print(f"null Spearman: median {med:+.4f}, 95th pct "
          f"{nulls[int(0.95*len(nulls))]:+.4f}, max {nulls[-1]:+.4f}")
    print(f"p = {p:.4f} ({len(nulls)} within-species relabelings)")
    with open(a.out, "w") as o:
        o.write("statistic\tvalue\n")
        o.write(f"ortholog_pairs\t{len(shared)}\n")
        o.write(f"observed_spearman\t{obs:.4f}\n")
        o.write(f"null_median\t{med:.4f}\n")
        o.write(f"null_p95\t{nulls[int(0.95*len(nulls))]:.4f}\n")
        o.write(f"null_max\t{nulls[-1]:.4f}\n")
        o.write(f"p_permutation\t{p:.4f}\n")
        o.write(f"n_permutations\t{len(nulls)}\n")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
