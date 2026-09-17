#!/usr/bin/env python3
"""Simulation: why cross-species splicing-conservation tests return
spectacular false positives, and what fixes them.

    sim_conservation_artefact.py --reps 50 --out sim_conservation.tsv

Ground truth in every simulated dataset: **there is NO conserved regulation
whatsoever.** Each species' developmental effect is applied to an
INDEPENDENTLY drawn subset of genes. Any cross-species signal a statistic
reports is therefore, by construction, false.

The one thing shared between species is what is shared in real data:
STATIC gene properties. Orthologues have correlated expression, so they
have correlated read depth, so their PSI estimates have correlated noise;
and they have similar exon counts. Genes that are noisy in one species are
noisy in the other, so they cross a |dPSI| threshold in both and they carry
large |dPSI| in both -- with no regulation involved.

Statistics compared, on identical data:
  1. hypergeometric overlap of threshold-defined "changing gene" sets
  2. Spearman correlation of per-gene mean |dPSI|
  3. the same two, scored against a WITHIN-SPECIES LABEL PERMUTATION that
     recomputes the whole statistic (the fix)

The depth-heterogeneity sweep shows the artefact is driven by exactly the
property real transcriptomes have most of.
"""
import argparse, math, sys
import numpy as np


def log_choose(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    if k <= 0:
        return 1.0
    tot, den = float("-inf"), log_choose(N, n)
    for i in range(k, min(K, n) + 1):
        t = log_choose(K, i) + log_choose(N - K, n - i) - den
        tot = t if tot == float("-inf") else (max(tot, t) + math.log1p(
            math.exp(min(tot, t) - max(tot, t))))
    return min(1.0, math.exp(tot)) if tot > float("-inf") else 0.0


def rankdata(a):
    order = np.argsort(a, kind="mergesort")
    r = np.empty(len(a), float)
    r[order] = np.arange(len(a), dtype=float)
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt)); np.add.at(sums, inv, r)
    return (sums / cnt)[inv]


def spear(x, y):
    rx, ry = rankdata(np.asarray(x)), rankdata(np.asarray(y))
    rx = rx - rx.mean(); ry = ry - ry.mean()
    d = math.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else 0.0


def simulate(rng, G, n_per_group, sigma_shared, frac_regulated):
    """Per species: a padded (G, 2n, Kmax) sample array + exon mask.
    NO conserved regulation -- each species regulates its own random subset.
    The ONLY thing shared is static: orthologue expression is correlated,
    so log depth shares a component, so PSI noise is correlated."""
    n_ex = 1 + rng.poisson(2.5, size=G)
    Kmax = int(n_ex.max())
    shared_depth = rng.normal(0, sigma_shared, size=G)
    species = []
    for sp in range(2):
        depth = np.exp(3.0 + shared_depth + rng.normal(0, 0.3, size=G))
        regulated = rng.random(G) < frac_regulated   # INDEPENDENT per species
        X = np.zeros((G, 2 * n_per_group, Kmax))
        mask = np.zeros((G, Kmax), bool)
        for g in range(G):
            k = int(n_ex[g])
            mask[g, :k] = True
            base = rng.beta(1.2, 1.2, size=k)
            se = 1.0 / np.sqrt(depth[g] * 4.0)
            eff = rng.normal(0, 0.25, size=k) if regulated[g] else np.zeros(k)
            mu0 = np.clip(base, 0.01, 0.99)
            mu1 = np.clip(base + eff, 0.01, 0.99)
            X[g, :n_per_group, :k] = np.clip(
                mu0[None, :] + rng.normal(0, se, size=(n_per_group, k)), 0, 1)
            X[g, n_per_group:, :k] = np.clip(
                mu1[None, :] + rng.normal(0, se, size=(n_per_group, k)), 0, 1)
        species.append((X, mask))
    return species


def magnitudes(X, mask, idx, n):
    """Per-gene mean |dPSI| for a sample split. idx is ONE permutation of
    the 2n samples applied to EVERY gene -- which is what a real label
    permutation does. (Permuting independently per gene, as a first draft
    of this script did, destroys more structure than relabeling and makes
    the null too lenient.)"""
    a = X[:, idx[:n], :].mean(axis=1)
    b = X[:, idx[n:], :].mean(axis=1)
    d = np.abs(b - a)
    return (d * mask).sum(axis=1) / mask.sum(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=50)
    ap.add_argument("--genes", type=int, default=3000)
    ap.add_argument("--n-per-group", type=int, default=6)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    thr, n = 0.10, a.n_per_group
    rows = []
    for sigma in (0.0, 0.4, 0.8, 1.2):
        stats = {"hyper_fp": 0, "rho_fp": 0, "perm_ov_fp": 0, "perm_rho_fp": 0}
        enrs, rhos, hyper_ps = [], [], []
        for rep in range(a.reps):
            rng = np.random.default_rng(1000 * rep + int(sigma * 10))
            (X0, m0), (X1, m1) = simulate(rng, a.genes, n, sigma, 0.15)
            true_idx = np.arange(2 * n)
            g0 = magnitudes(X0, m0, true_idx, n)
            g1 = magnitudes(X1, m1, true_idx, n)
            c0 = set(np.where(g0 >= thr)[0]); c1 = set(np.where(g1 >= thr)[0])
            K, nn, k, N = len(c0), len(c1), len(c0 & c1), a.genes
            p_hyper = hypergeom_sf(k, N, K, nn) if K and nn else 1.0
            exp = K * nn / N if N else 0
            enr = k / exp if exp else float("nan")
            rho = spear(g0, g1)
            enrs.append(enr); rhos.append(rho); hyper_ps.append(p_hyper)
            ov_null, rho_null = [], []
            for _ in range(a.n_perm):
                i0 = rng.permutation(2 * n); i1 = rng.permutation(2 * n)
                p0 = magnitudes(X0, m0, i0, n); p1 = magnitudes(X1, m1, i1, n)
                a0 = set(np.where(p0 >= thr)[0]); a1 = set(np.where(p1 >= thr)[0])
                e = len(a0) * len(a1) / N if N else 0
                if e:
                    ov_null.append(len(a0 & a1) / e)
                rho_null.append(spear(p0, p1))
            ov_null = np.array(ov_null); rho_null = np.array(rho_null)
            p_ov = (1 + np.sum(ov_null >= enr)) / (1 + len(ov_null))
            p_rho = (1 + np.sum(rho_null >= rho)) / (1 + len(rho_null))
            stats["hyper_fp"] += p_hyper < 0.05
            stats["rho_fp"] += rho > 0.10
            stats["perm_ov_fp"] += p_ov < 0.05
            stats["perm_rho_fp"] += p_rho < 0.05
        r = a.reps
        med_hyper = float(np.median(hyper_ps))
        rows.append((sigma, float(np.mean(enrs)), float(np.mean(rhos)),
                     med_hyper, stats["hyper_fp"] / r, stats["rho_fp"] / r,
                     stats["perm_ov_fp"] / r, stats["perm_rho_fp"] / r))
        print(f"sigma={sigma:.1f} | overlap enr x{np.mean(enrs):.2f}  "
              f"rho {np.mean(rhos):+.3f}  median hypergeom p {med_hyper:.2e}"
              f" || NAIVE FPR: hypergeom {stats['hyper_fp']/r:.0%}, "
              f"rho>0.10 {stats['rho_fp']/r:.0%}"
              f" || CORRECTED FPR: {stats['perm_ov_fp']/r:.0%} / "
              f"{stats['perm_rho_fp']/r:.0%}")
    with open(a.out, "w") as o:
        o.write("sigma_shared_depth\tmean_overlap_enrichment\tmean_spearman\t"
                "median_hypergeom_p\tnaive_hypergeom_FPR\tnaive_rho_FPR\t"
                "permutation_overlap_FPR\tpermutation_rho_FPR\n")
        for t in rows:
            o.write("\t".join(f"{v:.6g}" for v in t) + "\n")
    print(f"wrote {a.out}  (ground truth in EVERY dataset: no conservation)")


if __name__ == "__main__":
    main()
