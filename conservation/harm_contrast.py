#!/usr/bin/env python3
"""The cross-species test: does metamorphosis restructure splicing in
H. armigera, and does the adhesion hypothesis from H. assulta hold up?

    harm_contrast.py --psi-harm hassulta/psi_harm --runs hassulta/runs_harm.tsv \
        --gff hassulta/ref/Harm.gff --gene2go hassulta/ref/harm_gene2go.tsv \
        --psi-assulta hassulta/psi --runs-assulta hassulta/runs_dev.tsv \
        --paf hassulta/ref/has_vs_harm.paf --out-prefix hassulta/harm

Three tests, ALL pre-registered in hassulta/DESIGN.md before the armigera
data was touched:

 1. Does metamorphosis restructure splicing here too? Same statistic as
    assulta (exons at |dPSI| >= 0.1 vs exact relabeling calibration).
    Prediction on record: yes, ratio > 1.5.
 2. THE KEY TEST: are the eight adhesion/cytoskeleton/muscle-attachment GO
    terms enriched? Those terms topped the assulta ranking but none
    survived FDR there, so they are a HYPOTHESIS. Testing them in a second
    species is legitimate; re-testing in assulta would be circular. The
    term list is hard-coded below and was fixed in advance.
 3. Do the same genes change? Overlap of the armigera changing set with the
    assulta changing set mapped through the ortholog map, hypergeometric.

Treatment (HAC) is balanced across stages, so it is variance, not
confounding, and is not modelled.
"""
import argparse, csv, itertools, math, os, re, sys
from collections import defaultdict

MIN_TOT = 20
DPSI = 0.10
# fixed in advance -- see DESIGN.md "Pre-registered, before the armigera
# data has been touched"
ADHESION_GO = ["GO:0005925", "GO:0098609", "GO:0031941", "GO:0005178",
               "GO:0030018", "GO:0008305", "GO:0007229", "GO:0007160"]


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    if k <= 0:
        return 1.0
    tot = float("-inf")
    den = log_choose(N, n)
    for i in range(k, min(K, n) + 1):
        t = log_choose(K, i) + log_choose(N - K, n - i) - den
        tot = t if tot == float("-inf") else (max(tot, t) + math.log1p(
            math.exp(min(tot, t) - max(tot, t))))
    return min(1.0, math.exp(tot)) if tot > float("-inf") else 0.0


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


def count_changed(psi, ks, g1, g2):
    n = 0
    for k in ks:
        a = sum(psi[l][k] for l in g1) / len(g1)
        b = sum(psi[l][k] for l in g2) / len(g2)
        if abs(a - b) >= DPSI:
            n += 1
    return n


def gff_maps(gff):
    """gene-feature id / transcript id / gene name -> NCBI GeneID."""
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
            g = gid.group(1)
            i = re.search(r"(?:^|;)ID=([^;]+)", p[8])
            nm = re.search(r"(?:^|;)gene=([^;]+)", p[8])
            for key in filter(None, [i.group(1) if i else None,
                                     re.sub(r"^(gene|rna)-", "", i.group(1)) if i else None,
                                     nm.group(1) if nm else None]):
                to_gid[key] = g
    return to_gid


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--psi-harm", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--gene2go", required=True)
    ap.add_argument("--psi-assulta", required=True)
    ap.add_argument("--runs-assulta", required=True)
    ap.add_argument("--paf", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    meta = list(csv.DictReader(open(a.runs), delimiter="\t"))
    larva = [r["label"] for r in meta if r["stage"] == "larva"]
    adult = [r["label"] for r in meta if r["stage"] == "adult"]
    psi, gene_of = load_psi(a.psi_harm, larva + adult)
    ks = sorted(set.intersection(*(set(psi[l]) for l in larva + adult)))
    print(f"armigera: {len(larva)}v{len(adult)}, exons tested in all: {len(ks)}")

    out = open(a.out_prefix + "_contrast.tsv", "w")
    out.write("test\tvalue\tnull_median\tratio\tp\tdetail\n")

    # --- TEST 1: does metamorphosis restructure splicing here too?
    obs = count_changed(psi, ks, larva, adult)
    pool = larva + adult
    nulls = []
    for s in itertools.combinations(range(len(pool)), len(larva)):
        h1 = [pool[i] for i in s]
        if set(h1) == set(larva) or set(h1) == set(adult):
            continue
        nulls.append(count_changed(psi, ks, h1,
                                   [pool[i] for i in range(len(pool)) if i not in s]))
    nulls.sort()
    med = nulls[len(nulls) // 2]
    p1 = (1 + sum(1 for v in nulls if v >= obs)) / (1 + len(nulls))
    print(f"TEST1 metamorphosis: obs={obs} null_med={med} x{obs/med:.2f} "
          f"p={p1:.4f} ({len(nulls)} relabelings)")
    out.write(f"metamorphosis_splicing\t{obs}\t{med}\t{obs/med:.3f}\t{p1:.4f}\t"
              f"{len(nulls)} exact relabelings\n")

    # --- TEST 2: the pre-registered adhesion hypothesis
    to_gid = gff_maps(a.gff)
    gid2go = defaultdict(set)
    go_name = {}
    with open(a.gene2go) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 6:
                gid2go[p[1]].add(p[2])
                go_name[p[2]] = p[5]

    def gid_of(k):
        g = gene_of[k]
        return to_gid.get(g) or to_gid.get(re.sub(r"^(gene|rna)-", "", g))
    bg_genes = {gid_of(k) for k in ks}
    bg_genes = {g for g in bg_genes if g and g in gid2go}
    ch_genes = {gid_of(k) for k in ks
                if abs(sum(psi[l][k] for l in larva) / len(larva)
                       - sum(psi[l][k] for l in adult) / len(adult)) >= DPSI}
    ch_genes = {g for g in ch_genes if g and g in gid2go}
    N, n = len(bg_genes), len(ch_genes)
    print(f"armigera genes with GO: background {N}, changing {n}")
    hits_bg = {g for g in bg_genes if gid2go[g] & set(ADHESION_GO)}
    hits_ch = {g for g in ch_genes if gid2go[g] & set(ADHESION_GO)}
    K, k = len(hits_bg), len(hits_ch)
    if K >= 3:
        p2 = hypergeom_sf(k, N, K, n)
        ratio = (k / n) / (K / N) if K and n else float("nan")
        print(f"TEST2 adhesion set (PRE-REGISTERED): {k}/{n} changing vs "
              f"{K}/{N} background  x{ratio:.2f}  p={p2:.4f}")
        out.write(f"adhesion_PREREGISTERED\t{k}\t{K}\t{ratio:.3f}\t{p2:.4f}\t"
                  f"union of 8 GO terms; N={N} n={n}\n")
        for t in ADHESION_GO:
            Kt = sum(1 for g in bg_genes if t in gid2go[g])
            kt = sum(1 for g in ch_genes if t in gid2go[g])
            if Kt >= 3:
                pt = hypergeom_sf(kt, N, Kt, n)
                out.write(f"  {t}\t{kt}\t{Kt}\t"
                          f"{((kt/n)/(Kt/N)) if Kt and n else 0:.3f}\t{pt:.4f}\t"
                          f"{go_name.get(t,'')}\n")
    else:
        print(f"TEST2: only {K} background genes carry the adhesion terms — "
              f"UNTESTABLE, reported as such")
        out.write(f"adhesion_PREREGISTERED\tNA\t{K}\tNA\tNA\tuntestable: "
                  f"only {K} background genes\n")

    # --- TEST 3: do the same genes change?
    ameta = list(csv.DictReader(open(a.runs_assulta), delimiter="\t"))
    aby = defaultdict(list)
    for r in ameta:
        aby[(r["stage"], r["sex"])].append(r["label"])
    ag1 = aby[("pupa", "M")] + aby[("pupa", "F")]
    ag2 = aby[("adult", "M")] + aby[("adult", "F")]
    apsi, agene = load_psi(a.psi_assulta, ag1 + ag2)
    aks = sorted(set.intersection(*(set(apsi[l]) for l in ag1 + ag2)))
    a_changed = {agene[k] for k in aks
                 if abs(sum(apsi[l][k] for l in ag1) / len(ag1)
                        - sum(apsi[l][k] for l in ag2) / len(ag2)) >= DPSI}
    a_tested = {agene[k] for k in aks}
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
    shared_bg = {a2gid[g] for g in a_tested if g in a2gid} & bg_genes
    a_ch_gid = {a2gid[g] for g in a_changed if g in a2gid} & shared_bg
    h_ch_gid = ch_genes & shared_bg
    ov = a_ch_gid & h_ch_gid
    N3, K3, n3, k3 = len(shared_bg), len(a_ch_gid), len(h_ch_gid), len(ov)
    if N3 >= 100 and K3 and n3:
        p3 = hypergeom_sf(k3, N3, K3, n3)
        exp = K3 * n3 / N3
        print(f"TEST3 shared changing genes: {k3} observed vs {exp:.1f} "
              f"expected (x{k3/exp:.2f}), p={p3:.4g}  "
              f"[assulta {K3}, armigera {n3}, shared background {N3}]")
        out.write(f"gene_overlap\t{k3}\t{exp:.2f}\t{k3/exp:.3f}\t{p3:.4g}\t"
                  f"assulta={K3} armigera={n3} background={N3}\n")
    else:
        print(f"TEST3: shared background too small ({N3})")
        out.write(f"gene_overlap\tNA\tNA\tNA\tNA\tshared background {N3}\n")
    out.close()
    print(f"wrote {a.out_prefix}_contrast.tsv")


if __name__ == "__main__":
    main()
