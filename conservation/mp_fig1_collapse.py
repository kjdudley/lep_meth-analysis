#!/usr/bin/env python3
"""Methods paper, Figure 1: the real-data false positive, decomposed.

(a) Cross-species overlap of the top-K splicing genes (K size-matched at
    507/511): the independence expectation (what the hypergeometric tests
    against), the label-permutation null (the chain's contribution), and
    the observed value. The gap between expectation and null median is the
    artefact; the gap between null and observed is the residual signal.
(b) The threshold-free check: Spearman correlation of per-gene |dPSI|,
    observed vs the label-permutation null (summary intervals; the naive
    p-value against zero is meaningless).
"""
import os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.size"] = 8
import matplotlib.pyplot as plt
import numpy as np

# Repo-relative: tables ship in conservation/tables/ in this repository.
_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(_HERE, "tables")      # inputs ship here
OUTDIR = _HERE                             # figures are written beside the script
os.makedirs(os.path.join(OUTDIR, "figures"), exist_ok=True)
nulls = np.array([int(l) for l in
                  open(os.path.join(ROOT, "harm_topk.nulls.tsv")).read().split()[1:]])
OBS, EXP = 183, 90.6
INK = "#1a1a2e"; BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#2e8b57"

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0),
                              gridspec_kw=dict(width_ratios=[1.5, 1.0],
                                               wspace=0.32, left=0.08,
                                               right=0.97, top=0.86, bottom=0.17))

# (a) decomposition
ax.hist(nulls, bins=np.arange(min(nulls)-0.5, max(nulls)+1.5, 2),
        color="#b9cbe4", edgecolor="none",
        label="label-permutation null (n = 1,000)")
ax.axvline(EXP, color=INK, lw=1.2, ls=":")
ax.text(EXP, ax.get_ylim()[1]*0.97, " independence\n expectation (90.6)",
        fontsize=6.4, va="top", color=INK)
med = float(np.median(nulls))
ax.axvline(med, color=BLUE, lw=1.2)
ax.text(med-2, ax.get_ylim()[1]*0.97, "null median (163) ", fontsize=6.4,
        va="top", ha="right", color=BLUE)
ax.axvline(OBS, color=ORANGE, lw=1.6)
ax.text(OBS-1.5, ax.get_ylim()[1]*0.72, "observed (183)\np = 0.011",
        fontsize=6.6, color=ORANGE, ha="right")
ymid = ax.get_ylim()[1]*0.44
ax.annotate("", xy=(med, ymid), xytext=(EXP, ymid),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
ax.text((EXP+med)/2, ymid*1.06, "the artefact:\nshared gene properties\n(hypergeometric p = 6.5e-28)",
        ha="center", fontsize=6.2, color=INK)
ax.annotate("", xy=(OBS, ymid*0.55), xytext=(med, ymid*0.55),
            arrowprops=dict(arrowstyle="<->", color=GREEN, lw=0.9))
ax.text((med+OBS)/2, ymid*0.62, "residual\nsignal", ha="center",
        fontsize=6.2, color=GREEN)
ax.set_xlabel("cross-species overlap of size-matched top-K splicing genes")
ax.set_ylabel("permutations")
ax.set_title("a", loc="left", fontweight="bold")
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)

# (b) spearman summary
obs_r, nmed, np95, nmax = 0.3417, 0.3484, 0.3709, 0.3894
ax2.errorbar([0], [nmed], yerr=[[nmed-0.30], [nmax-nmed]], fmt="none",
             ecolor="#b9cbe4", elinewidth=10, capsize=0)
ax2.plot([0], [nmed], "_", color=BLUE, ms=22, mew=2, label="null median")
ax2.plot([0], [np95], "_", color=BLUE, ms=14, mew=1, alpha=0.6)
ax2.plot([0], [obs_r], "o", color=ORANGE, ms=7, label="observed")
ax2.text(0.09, obs_r - 0.012, "observed +0.342\np = 0.69", fontsize=6.6,
         color=ORANGE, va="center")
ax2.text(0.09, nmax, "null max +0.389", fontsize=6.2, color=BLUE, va="center")
ax2.text(0.09, nmed + 0.008, "null median +0.348", fontsize=6.2, color=BLUE,
         va="center")
ax2.set_xlim(-0.35, 0.75); ax2.set_ylim(0.29, 0.40)
ax2.set_xticks([])
ax2.set_ylabel("Spearman ρ of per-gene |dPSI|")
ax2.set_title("b", loc="left", fontweight="bold")
ax2.text(0.5, 1.02, "naive p < 1e-70 against zero;\nrelabelled data reproduce it",
         transform=ax2.transAxes, ha="center", fontsize=6.4, style="italic")
for sp in ("top", "right"):
    ax2.spines[sp].set_visible(False)
out = os.path.join(OUTDIR, "figures", "mp_fig1_collapse")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
