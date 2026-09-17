#!/usr/bin/env python3
"""methsplice Figure 2 (RNA style: capital panel tags, no in-figure
titles, Helvetica/Arial 8-10pt (RNA minimum enforced), lines >=0.5pt, RGB).

(A) Measured global erasure per degron line and day.
(B) Coupling slopes (dPSI on baseline exon methylation, and on measured
    per-exon erasure) with gene-clustered CIs, per contrast.
(C) Coupling slope against erasure depth across contrasts: the
    slope-of-slopes (+2.4e-07) showing no dose scaling.
Data: methsplice/degron_dose.tsv.
"""
import csv, os
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.6, "lines.linewidth": 1.2})
import matplotlib.pyplot as plt
import numpy as np

# Repo-relative: inputs ship in methsplice/tables/, figures are written beside
# this script. The original walked three directories up from a different layout.
_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = _HERE
os.makedirs(os.path.join(_HERE, "figures"), exist_ok=True)
rows = list(csv.DictReader(open(os.path.join(ROOT, "tables/degron_dose.tsv")), delimiter="\t"))
INK = "#1a1a2e"; BLUE, ORANGE = "#2a78d6", "#eb6834"
LINE = {"d1aid": "DNMT1-AID", "u1aid": "UHRF1-AID", "dual": "dual"}

def parse(c):
    p = c.split("_")           # degron_d1aid_Day6
    return p[1], int(p[2][3:])

contrasts = {}
for r in rows:
    line, day = parse(r["contrast"])
    d = contrasts.setdefault((line, day), {"erasure": -float(r["mean_erasure_pct"])})
    d[r["test"]] = (float(r["slope"]), float(r["ci_lo"]), float(r["ci_hi"]))

order = [("d1aid", 6), ("d1aid", 8), ("u1aid", 6), ("u1aid", 8),
         ("dual", 6), ("dual", 8)]
labels = [f"{LINE[l]}\nd{d}" for l, d in order]

fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.4),
                         gridspec_kw=dict(wspace=0.45, left=0.09,
                                          right=0.98, top=0.90, bottom=0.34))
A, B, C = axes
for ax, tag in zip(axes, "ABC"):
    ax.text(-0.18, 1.06, tag, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

# A: erasure
x = np.arange(6)
A.bar(x, [contrasts[k]["erasure"] for k in order], color=BLUE, width=0.6)
A.set_xticks(x, labels, fontsize=8, rotation=90, ha="center")
A.set_ylabel("global methylation lost\n(percentage points)", fontsize=8)
for i, k in enumerate(order):
    A.text(i, contrasts[k]["erasure"] + 0.8, f"{contrasts[k]['erasure']:.0f}",
           ha="center", fontsize=8, color=INK)

# B: slopes with CIs
for off, key, col, lab in ((-0.17, "baseline_meth", BLUE, "on baseline methylation"),
                           (0.17, "measured_derasure", ORANGE, "on measured erasure")):
    y = [contrasts[k][key][0] * 1e6 for k in order]
    lo = [contrasts[k][key][0]*1e6 - contrasts[k][key][1]*1e6 for k in order]
    hi = [contrasts[k][key][2]*1e6 - contrasts[k][key][0]*1e6 for k in order]
    B.errorbar(x + off, y, yerr=[lo, hi], fmt="o", ms=3.6, color=col,
               elinewidth=0.9, capsize=1.6, label=lab)
B.axhline(0, color=INK, lw=0.6, ls=":")
B.set_xticks(x, labels, fontsize=8, rotation=90, ha="center")
B.set_ylabel("coupling slope\n(dPSI per methylation point, x1e-6)", fontsize=8)
B.legend(fontsize=8, frameon=False, loc="lower left")

# C: slope vs erasure depth, both coupling definitions
er = np.array([contrasts[k]["erasure"] for k in order])
xs = np.array([min(er)-2, max(er)+2])
worst = 0.0
for key, col in (("baseline_meth", BLUE), ("measured_derasure", ORANGE)):
    sl = np.array([contrasts[k][key][0] * 1e6 for k in order])
    C.plot(er, sl, "o", color=col, ms=4.2)
    b, a = np.polyfit(er, sl, 1)
    worst = max(worst, abs(b * 1e-6))
    C.plot(xs, a + b*xs, color=col, lw=0.8, ls="--")
C.axhline(0, color=INK, lw=0.6, ls=":")
C.set_xlabel("erasure depth (percentage points)", fontsize=8)
C.set_ylabel("coupling slope (x1e-6)", fontsize=8)
C.text(0.04, 0.06, f"|slope of slopes| < {worst:.1e}\nper point, both definitions",
       transform=C.transAxes, fontsize=8, color=INK)

out = os.path.join(_HERE, "figures", "msf_fig2_dose")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
