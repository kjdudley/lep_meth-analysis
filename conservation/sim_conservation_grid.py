#!/usr/bin/env python3
"""Grid extension of the conservation-artefact simulation (methods paper).

    sim_conservation_grid.py --reps 30 --genes 2859 --n-perm 200 \
        --out sim_conservation_grid.tsv

Reuses simulate()/magnitudes()/spear()/hypergeom_sf() from
sim_conservation_artefact.py unchanged (imported, not copied), so the grid
is guaranteed to run the same generative model as the published single-cell
of the grid (n=6, thr=0.10). Sweeps:

  n_per_group in {3, 6, 12}   (typical RNA-seq designs)
  dpsi_thr    in {0.05, 0.10, 0.20}
  sigma       in {0.0, 0.4, 0.8, 1.2}  (shared-depth heterogeneity)

Ground truth everywhere: no conserved regulation. Output adds n_per_group
and dpsi_thr columns ahead of the original schema.

v2: the first grid run exposed a set-size pathology in the
enrichment-ratio permutation itself (at thr=0.20, n=6, sigma=0 the
corrected-overlap FPR hit 87% at observed enrichment ~1.0): destroying the
real within-species regulation shrinks the permuted threshold sets toward
empty, the ratio statistic degenerates, and the observed ratio beats a
null centred below 1. This version therefore also scores the
HYPERGEOMETRIC P-VALUE against its own permutation distribution
(perm_pval_FPR), and records mean observed/permuted set sizes, so the
paper can show all three regimes: naive (anticonservative with shared
depth), ratio-permutation (fails when permuted sets collapse), and
p-value-permutation (fails harder: near-empty permuted sets pin the
permuted hypergeometric p at ~1, so any unremarkable observed p scores as
extreme; 100% FPR at n=12/thr=0.20 in v2).

v3 diagnosis: the label permutation destroys real within-species
regulation along with any conservation, so the permuted threshold sets are
systematically smaller than the observed ones, and every statistic whose
marginals depend on set size inherits the mismatch. The fix is a statistic
with permutation-invariant marginals: rank correlation (validated in v2,
worst FPR 13%), or TOP-K OVERLAP, where each species' set is its top K
genes by magnitude with K fixed at the observed threshold-set size, so
permuted sets are size-matched by construction (permutation_topk_FPR).
"""
import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_conservation_artefact as sim


def one_cell(a, n, thr, sigma):
    st = {"hyper_fp": 0, "rho_fp": 0, "perm_ov_fp": 0, "perm_rho_fp": 0,
          "perm_pp_fp": 0, "perm_tk_fp": 0}
    enrs, rhos, hyper_ps, obs_sizes, perm_sizes = [], [], [], [], []
    for rep in range(a.reps):
        rng = np.random.default_rng(
            7_000_000 + 100_000 * n + int(thr * 1000) * 100 + 10 * rep
            + int(sigma * 10))
        (X0, m0), (X1, m1) = sim.simulate(rng, a.genes, n, sigma, 0.15)
        true_idx = np.arange(2 * n)
        g0 = sim.magnitudes(X0, m0, true_idx, n)
        g1 = sim.magnitudes(X1, m1, true_idx, n)
        c0 = set(np.where(g0 >= thr)[0]); c1 = set(np.where(g1 >= thr)[0])
        K, nn, k, N = len(c0), len(c1), len(c0 & c1), a.genes
        K0, K1 = max(K, 1), max(nn, 1)
        t0 = set(np.argsort(g0)[-K0:]); t1 = set(np.argsort(g1)[-K1:])
        obs_topk = len(t0 & t1)
        p_hyper = sim.hypergeom_sf(k, N, K, nn) if K and nn else 1.0
        exp = K * nn / N if N else 0
        enr = k / exp if exp else float("nan")
        rho = sim.spear(g0, g1)
        enrs.append(enr); rhos.append(rho); hyper_ps.append(p_hyper)
        obs_sizes.append(0.5 * (K + nn))
        ov_null, rho_null, pp_null, tk_null, psz = [], [], [], [], []
        for _ in range(a.n_perm):
            i0 = rng.permutation(2 * n); i1 = rng.permutation(2 * n)
            p0 = sim.magnitudes(X0, m0, i0, n)
            p1 = sim.magnitudes(X1, m1, i1, n)
            a0 = set(np.where(p0 >= thr)[0]); a1 = set(np.where(p1 >= thr)[0])
            e = len(a0) * len(a1) / N if N else 0
            if e:
                ov_null.append(len(a0 & a1) / e)
            pp_null.append(sim.hypergeom_sf(len(a0 & a1), N, len(a0), len(a1))
                           if a0 and a1 else 1.0)
            q0 = set(np.argsort(p0)[-K0:]); q1 = set(np.argsort(p1)[-K1:])
            tk_null.append(len(q0 & q1))
            psz.append(0.5 * (len(a0) + len(a1)))
            rho_null.append(sim.spear(p0, p1))
        ov_null = np.array(ov_null); rho_null = np.array(rho_null)
        pp_null = np.array(pp_null); tk_null = np.array(tk_null)
        perm_sizes.append(float(np.mean(psz)))
        p_ov = (1 + np.sum(ov_null >= enr)) / (1 + len(ov_null))
        p_rho = (1 + np.sum(rho_null >= rho)) / (1 + len(rho_null))
        p_pp = (1 + np.sum(pp_null <= p_hyper)) / (1 + len(pp_null))
        p_tk = (1 + np.sum(tk_null >= obs_topk)) / (1 + len(tk_null))
        st["hyper_fp"] += p_hyper < 0.05
        st["rho_fp"] += rho > 0.10
        st["perm_ov_fp"] += p_ov < 0.05
        st["perm_rho_fp"] += p_rho < 0.05
        st["perm_pp_fp"] += p_pp < 0.05
        st["perm_tk_fp"] += p_tk < 0.05
    r = a.reps
    return (n, thr, sigma, float(np.mean(enrs)), float(np.mean(rhos)),
            float(np.median(hyper_ps)), st["hyper_fp"] / r, st["rho_fp"] / r,
            st["perm_ov_fp"] / r, st["perm_rho_fp"] / r,
            st["perm_pp_fp"] / r, st["perm_tk_fp"] / r,
            float(np.mean(obs_sizes)), float(np.mean(perm_sizes)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--genes", type=int, default=2859)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    hdr = ("n_per_group\tdpsi_thr\tsigma_shared_depth\t"
           "mean_overlap_enrichment\tmean_spearman\tmedian_hypergeom_p\t"
           "naive_hypergeom_FPR\tnaive_rho_FPR\t"
           "permutation_overlap_FPR\tpermutation_rho_FPR\t"
           "permutation_pval_FPR\tpermutation_topk_FPR\t"
           "mean_obs_set\tmean_perm_set\n")
    # resumable: skip cells already in the output
    done = set()
    if os.path.exists(a.out):
        with open(a.out) as fh:
            first = next(fh, "")
            if first != hdr:
                sys.exit(f"FATAL: {a.out} exists with a different schema; "
                         "refusing to resume into it")
            for line in fh:
                f = line.split("\t")
                done.add((int(f[0]), float(f[1]), float(f[2])))
    else:
        with open(a.out, "w") as o:
            o.write(hdr)
    for n in (3, 6, 12):
        for thr in (0.05, 0.10, 0.20):
            for sigma in (0.0, 0.4, 0.8, 1.2):
                if (n, thr, sigma) in done:
                    print(f"have n={n} thr={thr} sigma={sigma}")
                    continue
                row = one_cell(a, n, thr, sigma)
                with open(a.out, "a") as o:
                    o.write("\t".join(f"{v:.6g}" for v in row) + "\n")
                print(f"n={n} thr={thr:.2f} sigma={sigma:.1f} | "
                      f"enr x{row[3]:.2f} rho {row[4]:+.3f} | "
                      f"naive FPR {row[6]:.0%}/{row[7]:.0%} | "
                      f"ratio {row[8]:.0%} rho {row[9]:.0%} "
                      f"pval {row[10]:.0%} topk {row[11]:.0%} | sets "
                      f"{row[12]:.0f}/{row[13]:.0f}", flush=True)
    # completeness assertion on the artifact itself
    with open(a.out) as fh:
        nrows = sum(1 for _ in fh) - 1
    assert nrows == 36, f"grid has {nrows} rows, expected 36"
    print(f"grid complete: 36 cells -> {a.out}")


if __name__ == "__main__":
    main()
