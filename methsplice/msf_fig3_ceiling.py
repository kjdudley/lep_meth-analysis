#!/usr/bin/env python3
"""methsplice Figure 3: the constitutive-exon ceiling. RNA style.

(A) Honeybee (computed from methsplice/contrast.tsv): fraction of exons
    at control PSI > 0.98 by baseline methylation stratum; sparse-gbM
    genomes show the mark labelling constitutive exons.
(B) Coral (Acropora palmata), read from coral_contrast_rs.tsv and
    coral_contrast_denovo.tsv:
    the same fraction for heavily (>=50%) vs sparsely (<10%) methylated
    exons, under both annotations; the separation collapses where
    methylation blankets gene bodies.
Human evidence for the ceiling is the negative dose-response slopes
(Figure 2), cited in text.
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
INK = "#1a1a2e"; BLUE, ORANGE = "#2a78d6", "#eb6834"

# bee strata from the committed per-exon table
strata = [(0, 0.001, "0"), (0.001, 5, "0-5"), (5, 20, "5-20"),
          (20, 50, "20-50"), (50, 101, ">=50")]
counts = [[0, 0] for _ in strata]     # [n, n_ceiling]
with open(os.path.join(ROOT, "tables/contrast.tsv")) as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    for r in rdr:
        try:
            m = float(r["meth_c"]); psi = float(r["psi_c"])
        except ValueError:
            continue
        for i, (lo, hi, _) in enumerate(strata):
            if lo <= m < hi:
                counts[i][0] += 1
                counts[i][1] += psi > 0.98
                break
fr = [c[1]/c[0] if c[0] else float("nan") for c in counts]
ns = [c[0] for c in counts]
print("bee strata:", [f"{l}:{f*100:.1f}% (n={n})" for (_,_,l), f, n in zip(strata, fr, ns)])

fig, (A, B) = plt.subplots(1, 2, figsize=(7.0, 3.0),
                           gridspec_kw=dict(wspace=0.42, left=0.11,
                                            right=0.98, top=0.88, bottom=0.26))
for ax, tag in zip((A, B), "AB"):
    ax.text(-0.22, 1.07, tag, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

x = np.arange(len(strata))
A.bar(x, [f*100 for f in fr], color=BLUE, width=0.62)
A.set_xticks(x, [l for _, _, l in strata], fontsize=8)
A.set_xlabel("baseline exon methylation (%)", fontsize=8)
A.set_ylabel("exons at PSI > 0.98 (%)", fontsize=8)
A.set_ylim(0, 122)
for i, (f, n) in enumerate(zip(fr, ns)):
    A.text(i, f*100 + 2, f"{f*100:.0f}", ha="center", fontsize=8, color=INK)

# read from the coral contrast tables rather than carrying the values: the
# ceiling rows are `ceiling <meth>=50|meth<10>` with the fraction in `value`.
def _ceiling(fn):
    rows = list(csv.DictReader(
        open(os.path.join(ROOT, "tables", fn)), delimiter="\t"))
    g = {r["cell"]: float(r["value"]) for r in rows if r["test"] == "ceiling"}
    if "meth>=50" not in g or "meth<10" not in g:
        raise SystemExit(f"FATAL: {fn} has no ceiling rows")
    return (100 * g["meth>=50"], 100 * g["meth<10"])

vals = {"RefSeq": _ceiling("coral_contrast_rs.tsv"),
        "de novo": _ceiling("coral_contrast_denovo.tsv")}
print(f"coral ceiling: RefSeq {vals['RefSeq']}, de novo {vals['de novo']}")
x = np.arange(2)
B.bar(x - 0.19, [vals["RefSeq"][0], vals["de novo"][0]], width=0.36,
      color=ORANGE, label="methylated (>=50%)")
B.bar(x + 0.19, [vals["RefSeq"][1], vals["de novo"][1]], width=0.36,
      color=BLUE, label="sparse (<10%)")
B.set_xticks(x, ["RefSeq\nannotation", "de novo\nannotation"], fontsize=8)
B.set_ylabel("exons at PSI > 0.98 (%)", fontsize=8)
B.set_ylim(0, 122)
# Legend moved above the bars: at 8 pt (the RNA minimum) it no longer fits
# inside a panel whose bars reach ~85%, and the color key must stay in the
# artwork rather than move to the caption.
B.legend(fontsize=8, frameon=False, loc="upper center", ncol=2,
         handlelength=1.1, columnspacing=0.9, borderaxespad=0.1)

out = os.path.join(_HERE, "figures", "msf_fig3_ceiling")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
