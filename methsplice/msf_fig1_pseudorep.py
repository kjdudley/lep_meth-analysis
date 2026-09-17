#!/usr/bin/env python3
"""methsplice Figure 1: the founding result as a pseudo-replication
artefact. RNA style.

(A) Individual-level permutation outcomes (6 vs 6 bees): observed counts
    at the permutation floor and at p <= 0.05 against uniform expectation;
    significant exons after FDR: zero.
(B) The same processed data under the original read-pooled Fisher
    procedure regenerates the published scale (exon mode: 469 genes vs
    their 524; intron retention: 387 introns vs their 27 genes; units
    per legend).
(C) Within-group standard deviation of individual exon PSIs, computed
    from methsplice/contrast.tsv: the variance the pooled test ignores.
"""
import csv, os, re, statistics
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.linewidth": 0.6, "lines.linewidth": 1.2})
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
INK = "#1a1a2e"; BLUE, ORANGE = "#2a78d6", "#eb6834"

# (C) data: per-exon within-group SD over the six control individuals
ctrl_cols = ["psi_112c", "psi_114c", "psi_115c", "psi_131c", "psi_136c",
             "psi_138c"]
sds = []
with open(os.path.join(ROOT, "methsplice/contrast.tsv")) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        try:
            v = [float(r[c]) for c in ctrl_cols]
        except (ValueError, KeyError):
            continue
        if 0.05 < statistics.mean(v) < 0.95:      # informative exons
            sds.append(statistics.stdev(v))
print(f"informative exons: {len(sds)}, median within-group SD "
      f"{statistics.median(sds):.3f}")

fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.0),
                         gridspec_kw=dict(wspace=0.50, left=0.08,
                                          right=0.98, top=0.88, bottom=0.26))
A, B, C = axes
for ax, tag in zip(axes, "ABC"):
    ax.text(-0.24, 1.07, tag, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

# A: permutation outcomes vs expectation
x = np.arange(3)
obs = [87, 888, 0]
exp = [81, 1871, None]
A.bar(x - 0.19, [87, 888, 0], width=0.36, color=BLUE, label="observed")
A.bar(x[:2] + 0.19, [81, 1871], width=0.36, color="#b9cbe4",
      label="uniform expectation")
A.set_yscale("symlog", linthresh=10)
A.set_xticks(x, ["at permutation\nfloor", "p <= 0.05", "FDR < 0.05"],
             fontsize=8)
A.set_ylabel("exons (n = 37,416 tested)", fontsize=8)
A.text(2 - 0.19, 1.5, "0", ha="center", fontsize=8, color=INK)
A.legend(fontsize=8, frameon=False, loc="upper right")

# B: scale reconstruction. Our two counts are read from their artefacts:
# the exon-mode gene count is the number of distinct genes in
# pooled_fisher_sig.tsv, and the intron count is the "2013 SHAPE" line of
# ir_summary.txt. The two published values are literals because they come
# from Li-Byarlay et al. 2013 Table 1, not from any artefact of ours.
PUBLISHED_EXON_GENES = 524      # Li-Byarlay et al. 2013
PUBLISHED_IR_GENES = 27         # Li-Byarlay et al. 2013

_pf = os.path.join(ROOT, "methsplice/pooled_fisher_sig.tsv")
if not os.path.exists(_pf):
    raise SystemExit("FATAL: pooled_fisher_sig.tsv missing")
with open(_pf) as fh:
    _rows = list(csv.DictReader(fh, delimiter="\t"))
ours_exon_genes = len({r["gene"] for r in _rows})

_ir = open(os.path.join(ROOT, "methsplice/ir_summary.txt")).read()
_m = re.search(r"FDR<0\.1\):\s*(\d+)\s+of\s+\d+\s+introns", _ir)
if not _m:
    raise SystemExit("FATAL: cannot parse the 2013-shape count from ir_summary.txt")
ours_ir = int(_m.group(1))
print(f"panel B: {len(_rows)} exons in {ours_exon_genes} genes; "
      f"{ours_ir} introns")

x = np.arange(2)
ours = [ours_exon_genes, ours_ir]
pub = [PUBLISHED_EXON_GENES, PUBLISHED_IR_GENES]
B.bar(x - 0.19, ours, width=0.36, color=ORANGE,
      label="read-pooled Fisher,\nour reprocessing")
B.bar(x + 0.19, pub, width=0.36, color=BLUE,
      label="published claim")
B.set_xticks(x, ["exon mode\n(genes)", "intron retention\n(introns / genes)"],
             fontsize=8)
B.set_ylabel("reported significant units", fontsize=8)
for xi, v in ((x[0]-0.19, ours[0]), (x[0]+0.19, pub[0]),
              (x[1]-0.19, ours[1]), (x[1]+0.19, pub[1])):
    B.text(xi, v + 12, str(v), ha="center", fontsize=8, color=INK)
B.legend(fontsize=8, frameon=False, loc="upper right",
         bbox_to_anchor=(1.02, 1.13))
B.set_ylim(0, 640)

# C: within-group SD histogram
C.hist(sds, bins=np.arange(0, 0.31, 0.01), color=BLUE, edgecolor="none")
C.axvline(statistics.median(sds), color=ORANGE, lw=1.2)
C.text(statistics.median(sds) + 0.008, C.get_ylim()[1]*0.9,
       f"median {statistics.median(sds):.3f}", fontsize=8, color=ORANGE)
C.set_xlabel("within-group SD of individual PSI", fontsize=8)
C.set_ylabel("exons (0.05 < mean PSI < 0.95)", fontsize=8)

out = os.path.join(ROOT, "figures", "msf_fig1_pseudorep")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
