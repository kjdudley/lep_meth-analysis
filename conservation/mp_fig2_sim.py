import os
#!/usr/bin/env python3
"""Methods paper, Figure 2: the simulation.

(a) Naive false-positive rates vs shared depth heterogeneity at the real
    comparison's dimensions (n=6, thr=0.10, 2,859 genes) — ground truth is
    zero conservation in every replicate.
(b) The 36-cell grid as four small heatmaps (naive hypergeom, ratio-perm,
    pval-perm, topk-perm), FPR by (n, thr) at the worst sigma per cell —
    the two failed fixes and the control that holds, side by side.

Reads hassulta/sim_conservation.tsv (the original sigma sweep) and
hassulta/sim_conservation_grid.tsv (v3). Grayscale-survivable.
"""
import csv, os, sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.size"] = 8
import matplotlib.pyplot as plt
import numpy as np

# Repo-relative: tables ship in conservation/tables/ in this repository.
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tables")
sweep = list(csv.DictReader(open(os.path.join(ROOT, "sim_conservation.tsv")), delimiter="\t"))
grid = list(csv.DictReader(open(os.path.join(ROOT, "sim_conservation_grid.tsv")), delimiter="\t"))
assert len(grid) == 36, len(grid)

INK = "#1a1a2e"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#2e8b57"

fig = plt.figure(figsize=(7.0, 4.6))
gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1.0], hspace=0.62, wspace=0.35,
                      left=0.09, right=0.98, top=0.90, bottom=0.13)

# ---- (a) the sigma sweep at real dimensions
ax = fig.add_subplot(gs[0, :2])
sig = [float(r["sigma_shared_depth"]) for r in sweep]
ax.plot(sig, [100*float(r["naive_hypergeom_FPR"]) for r in sweep], "o-",
        color=BLUE, label="naive hypergeometric overlap", lw=1.6)
ax.plot(sig, [100*float(r["naive_rho_FPR"]) for r in sweep], "s-",
        color=ORANGE, label="naive rank correlation", lw=1.6)
ax.plot(sig, [100*float(r["permutation_rho_FPR"]) for r in sweep], "^-",
        color=GREEN, label="rank corr. + label-permutation recompute", lw=1.6)
ax.axhline(5, color=INK, lw=0.7, ls=":")
ax.text(0.02, 7, "nominal 5%", fontsize=6.5, color=INK)
ax.set_xlabel("shared depth heterogeneity (σ)")
ax.set_ylabel("false-positive rate (%)")
ax.set_ylim(-4, 104)
ax.set_title("a", loc="left", fontweight="bold")
ax.text(0.5, 1.02, "zero true conservation in every dataset",
        transform=ax.transAxes, ha="center", fontsize=7, style="italic")
ax.legend(fontsize=6.2, frameon=False, loc="center right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- (a2) observed vs simulated effect sizes
ax2 = fig.add_subplot(gs[0, 2:])
ax2.plot(sig, [float(r["mean_overlap_enrichment"]) for r in sweep], "o-",
         color=BLUE, lw=1.6, label="simulated overlap enrichment")
ax2.plot(sig, [float(r["mean_spearman"]) * 4 for r in sweep], "s-",
         color=ORANGE, lw=1.6, label="simulated Spearman ρ (×4 scale)")
ax2.axhline(2.02, color=BLUE, lw=1.0, ls="--")
ax2.text(0.02, 2.06, "observed enrichment ×2.02", fontsize=6.2, color=BLUE)
ax2.axhline(0.3417 * 4, color=ORANGE, lw=1.0, ls="--")
ax2.text(0.02, 0.3417 * 4 + 0.04, "observed ρ +0.34", fontsize=6.2, color=ORANGE)
ax2.set_xlabel("shared depth heterogeneity (σ)")
ax2.set_ylabel("effect size")
ax2.set_title("b", loc="left", fontweight="bold")
ax2.text(0.5, 1.02, "static properties alone approach the observed values",
         transform=ax2.transAxes, ha="center", fontsize=7, style="italic")
ax2.legend(fontsize=6.2, frameon=False, loc="lower right")
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)

# ---- (c) grid heatmaps: worst-sigma FPR per (n, thr), four statistics
cells = {}
for r in grid:
    key = (int(r["n_per_group"]), float(r["dpsi_thr"]))
    for col in ("naive_hypergeom_FPR", "permutation_overlap_FPR",
                "permutation_pval_FPR", "permutation_topk_FPR"):
        cells.setdefault(col, {}).setdefault(key, []).append(float(r[col]))
ns, thrs = [3, 6, 12], [0.05, 0.10, 0.20]
panels = [("naive_hypergeom_FPR", "naive\nhypergeometric"),
          ("permutation_overlap_FPR", "permuted\nenrichment ratio"),
          ("permutation_pval_FPR", "permuted\nhypergeometric p"),
          ("permutation_topk_FPR", "permuted top-K\n(size-matched)")]
for i, (col, label) in enumerate(panels):
    ax = fig.add_subplot(gs[1, i])
    M = np.array([[max(cells[col][(n, t)]) for t in thrs] for n in ns])
    im = ax.imshow(M * 100, cmap="Greys", vmin=0, vmax=100, aspect="auto")
    for y in range(3):
        for x in range(3):
            v = M[y, x] * 100
            ax.text(x, y, f"{v:.0f}", ha="center", va="center", fontsize=6.6,
                    color="white" if v > 55 else INK)
    ax.set_xticks(range(3), [f"{t:g}" for t in thrs], fontsize=6.2)
    ax.set_yticks(range(3), [str(n) for n in ns], fontsize=6.2)
    if i == 0:
        ax.set_ylabel("n per group", fontsize=7)
    ax.set_xlabel("|dPSI| threshold", fontsize=7)
    ax.set_title(label, fontsize=7, pad=3)
    if i == 0:
        ax.text(-0.32, 1.32, "c", transform=ax.transAxes, fontsize=11,
                fontweight="bold", color=INK)
fig.text(0.5, 0.012, "c: worst-case false-positive rate (%) across σ, per design cell",
         fontsize=7, style="italic", ha="center", color=INK)
out = os.path.join(ROOT, "figures", "mp_fig2_sim")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out + ".{pdf,png}")
