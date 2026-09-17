import os
#!/usr/bin/env python3
"""Methods paper, Figure 3: the mechanism.

(a) Schematic: the causal chain from shared static properties to false
    conservation signal, and where each null cuts the chain.
(b) The marginal-shrink pathology: mean observed vs mean permuted
    threshold-set size per grid cell (log-log). Cells collapse toward the
    x-axis at the strict threshold, which is where the ratio and p-value
    permutation statistics fail.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.size"] = 8
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# Repo-relative: tables ship in conservation/tables/ in this repository.
_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(_HERE, "tables")      # inputs ship here
OUTDIR = _HERE                             # figures are written beside the script
os.makedirs(os.path.join(OUTDIR, "figures"), exist_ok=True)
grid = list(csv.DictReader(open(os.path.join(ROOT, "sim_conservation_grid.tsv")), delimiter="\t"))
INK = "#1a1a2e"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#2e8b57"

fig = plt.figure(figsize=(7.0, 3.1))
gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1.0], wspace=0.30,
                      left=0.03, right=0.97, top=0.90, bottom=0.15)

# ---- (a) schematic
ax = fig.add_subplot(gs[0, 0]); ax.set_axis_off()
ax.set_xlim(0, 10); ax.set_ylim(0, 10)
def box(x, y, w, h, text, fc="#eef2f8"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12",
                                fc=fc, ec=INK, lw=0.8))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=6.6,
            color=INK)
def arrow(x0, y0, x1, y1, color=INK, ls="-"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                 mutation_scale=9, color=color, lw=1.1,
                                 linestyle=ls))
box(0.2, 7.6, 2.9, 1.6, "orthologs share\nexpression level,\nexon count, baseline PSI")
box(3.9, 7.6, 2.6, 1.6, "correlated read depth,\ncorrelated noise\nin per-gene statistics")
box(7.2, 7.6, 2.6, 1.6, "same genes cross\nthresholds in\nBOTH species")
arrow(3.1, 8.4, 3.85, 8.4); arrow(6.5, 8.4, 7.15, 8.4)
box(7.2, 4.6, 2.6, 1.5, "\"conserved\nregulation\"\np = 6.5e-28", fc="#fbe9e2")
arrow(8.5, 7.5, 8.5, 6.2)
box(0.2, 4.6, 2.9, 1.5, "independence null:\nblind to the whole chain", fc="#f5f5f5")
box(3.9, 4.6, 2.6, 1.5, "random-gene null:\nbreaks the pairing,\nnot the chain", fc="#f5f5f5")
box(0.3, 1.0, 5.2, 2.0,
    "label permutation: keeps the chain intact,\ndestroys only treatment structure\n"
    "(statistic must be size-invariant:\nrank correlation, top-K)", fc="#e8f3ec")
arrow(5.7, 2.6, 7.5, 4.5, color=GREEN)
ax.text(7.0, 2.6, "measures the chain's\ncontribution", fontsize=6.2,
        color=GREEN, ha="left")
ax.set_title("a", loc="left", fontweight="bold")

# ---- (b) set-size collapse
ax2 = fig.add_subplot(gs[0, 1])
mk = {0.05: "o", 0.10: "s", 0.20: "^"}
col = {0.05: BLUE, 0.10: GREEN, 0.20: ORANGE}
for r in grid:
    thr = float(r["dpsi_thr"])
    obs = float(r["mean_obs_set"]); perm = max(float(r["mean_perm_set"]), 0.5)
    ax2.plot(obs, perm, mk[thr], color=col[thr], ms=4.5, alpha=0.85)
lo, hi = 0.4, 3000
ax2.plot([lo, hi], [lo, hi], color=INK, lw=0.7, ls=":")
ax2.text(500, 700, "equal size", fontsize=6, rotation=38, color=INK)
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlim(50, 3000); ax2.set_ylim(lo, 3000)
ax2.set_xlabel("mean observed set size (genes)")
ax2.set_ylabel("mean permuted set size")
for thr in (0.05, 0.10, 0.20):
    ax2.plot([], [], mk[thr], color=col[thr], label=f"|dPSI| ≥ {thr:g}")
ax2.legend(fontsize=6.2, frameon=False, loc="upper left")
ax2.set_title("b", loc="left", fontweight="bold")
ax2.text(0.97, 0.06, "permuted sets collapse at\nstrict thresholds; size-dependent\nstatistics fail here",
         transform=ax2.transAxes, ha="right", fontsize=6.2, color=ORANGE)
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)
out = os.path.join(ROOT, "figures", "mp_fig3_mechanism")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
