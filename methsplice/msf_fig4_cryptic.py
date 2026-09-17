#!/usr/bin/env python3
"""methsplice Figure 4: chronic demethylation produces transcripts, not
exon choices. RNA style (capital tags, no in-figure titles).

EVERY NUMBER IN THIS FIGURE IS READ FROM AN ARTEFACT. Rewritten 2026-09-11
after the source audit found panels A, B and D hard-coding their values and
panel B sourcing them from a DESIGN.md note rather than a table, with no way
for ms_consistency.py to see a drift. Panel C read cryptic_degron.tsv, which
holds 5 of the 9 degron contrasts. A missing artefact is now a hard failure,
not a silent fallback.

(A) Condition-exclusive novel junction counts, DKO vs WT, with the acute
    DNMT1i control            <- cryptic_all.tsv, else cryptic.tsv
(B) dPSI excess concentrating in novel-junction genes, each class against
    its OWN relabeling null   <- dko_overlap.tsv (containment rule, primary)
(C) The within-gene 3'-coverage tilt across all nine degron contrasts
                              <- cryptic_all.tsv, else cryptic_degron.tsv
                                 + cryptic_missing.tsv
(D) Novel-junction ends vs wild-type PRO-cap peaks   <- latent_tss.txt
"""
import csv, os, re, sys
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
MS = os.path.join(ROOT, "methsplice")
INK = "#1a1a2e"; INK2 = "#52514e"
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#2e8b57"


def need(name):
    p = os.path.join(MS, name)
    if not os.path.exists(p):
        sys.exit(f"FATAL: {name} missing. This figure reads every number from "
                 f"an artefact and will not substitute a literal.")
    return p


def tsv(name):
    return list(csv.DictReader(open(need(name)), delimiter="\t"))


def cryptic_rows():
    """Prefer the complete calibrated table; else assemble the partial ones.
    ms1_cryptic2.py writes cryptic_all.tsv over every contrast in the config;
    before it existed the arm lived in cryptic.tsv (5 rows), cryptic_degron.tsv
    (5 of 9 degron rows) and cryptic_missing.tsv (the 4 that were never run)."""
    if os.path.exists(os.path.join(MS, "cryptic_all.tsv")):
        return {r["contrast"]: r for r in tsv("cryptic_all.tsv")}, True
    out = {}
    for f in ("cryptic.tsv", "cryptic_degron.tsv", "cryptic_missing.tsv"):
        for r in tsv(f):
            out[r["contrast"]] = r
    return out, False


CR, COMPLETE = cryptic_rows()
OV = {(r["rule"], r["class"]): r for r in tsv("dko_overlap.tsv")}
print(f"cryptic table: {len(CR)} contrasts"
      f"{' (cryptic_all.tsv)' if COMPLETE else ' (assembled from partials)'}")

fig, axgrid = plt.subplots(2, 2, figsize=(7.0, 5.4),
                           gridspec_kw=dict(wspace=0.34, hspace=0.55,
                                            left=0.10, right=0.98,
                                            top=0.93, bottom=0.09))
axes = axgrid.ravel()
A, B, C, D = axes
for ax, tag in zip(axes, "ABCD"):
    ax.text(-0.16, 1.06, tag, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color=INK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

# ---- A: junction asymmetry, from the cryptic table
pick = ["HCT116_DKO_vs_WT", "HCT116_DNMT1i_3d"]
g1 = [int(CR[c]["novel_g1"]) for c in pick]
g2 = [int(CR[c]["novel_g2"]) for c in pick]
rat = [float(CR[c]["novel_ratio"]) for c in pick]
x = np.arange(len(pick))
A.bar(x - 0.19, g1, width=0.36, color=ORANGE, label="perturbed-only")
A.bar(x + 0.19, g2, width=0.36, color=BLUE, label="WT-only")
A.set_yscale("log")
A.set_ylim(200, max(g1 + g2) * 6)
A.set_xticks(x, ["DKO vs WT\n(chronic)", "DNMT1i 3 d\n(acute)"], fontsize=8)
A.set_ylabel("condition-exclusive novel junctions", fontsize=8)
for i, (v, r) in enumerate(zip(g1, rat)):
    A.text(i - 0.19, v * 1.35, f"x{r:.2f}", ha="center", fontsize=8, color=INK)
# The relabeling null, when the complete table carries it. Placed on the
# panel floor: the headroom above the bars belongs to the ratio labels and
# the legend, and a log axis makes anything in between collide.
if COMPLETE and "novel_null_med" in CR[pick[0]]:
    for i, c in enumerate(pick):
        # The DKO's interval is sequencing variance only -- its four libraries
        # are two runs of one library per genotype -- so it is labelled a
        # relabeling range, not a null, and carries no p. See DESIGN
        # "Finding 9".
        if c == "HCT116_DKO_vs_WT":
            A.text(i, 245, f"relabeling range x{float(CR[c]['novel_null_lo']):.2f}"
                   f"-{float(CR[c]['novel_null_hi']):.2f}\n(runs, not replicates)",
                   ha="center", fontsize=6, color=INK2)
        else:
            A.text(i, 245,
                   f"null x{float(CR[c]['novel_null_med']):.2f}  p={CR[c]['novel_p']}",
                   ha="center", fontsize=6.5, color=INK2)
A.legend(fontsize=8, frameon=False, loc="upper right", handlelength=1.2,
         borderaxespad=0.1)

# ---- B: concentration, containment rule, each class against its own null
oin = OV[("overlap", "novel_junction_genes")]
oout = OV[("overlap", "other_genes")]
vals = [100 * float(oin["obs_rate"]), 100 * float(oout["obs_rate"])]
nulls = [100 * float(oin["null_rate"]), 100 * float(oout["null_rate"])]
frac_in = int(oin["n_exon"]) / (int(oin["n_exon"]) + int(oout["n_exon"]))
B.bar([0, 1], vals, color=[ORANGE, BLUE], width=0.55, zorder=2)
for i, nv in enumerate(nulls):
    B.plot([i - 0.34, i + 0.34], [nv, nv], color=INK, lw=1.1, ls=":", zorder=3)
    B.text(i + 0.37, nv, f"null {nv:.2f}", fontsize=6.5, color=INK,
           va="center", ha="left")
B.set_xticks([0, 1], [f"novel-junction\ngenes ({100*frac_in:.1f}%)",
                      "other\ngenes"], fontsize=8)
B.set_ylabel("exons with |dPSI| > 0.05 (%)", fontsize=8)
for i, v in enumerate(vals):
    B.text(i, v + 0.4, f"{v:.2f}", ha="center", fontsize=8, color=INK)
B.text(0.5, max(vals) * 1.13,
       f"{vals[0]/vals[1]:.2f}x enrichment; "
       f"{100*float(oin['share']):.0f}% of excess",
       ha="center", fontsize=7, color=INK)
B.set_xlim(-0.7, 1.9)
B.set_ylim(0, max(vals) * 1.25)

# ---- C: 3'-tilt across ALL nine degron contrasts
series = {}
for c, r in CR.items():
    p = c.split("_")
    if p[0] != "degron":
        continue
    series.setdefault(p[1], []).append((int(p[2][3:]), float(r["s1_dslope"])))
n_deg = sum(len(v) for v in series.values())
if n_deg != 9:
    sys.exit(f"FATAL: {n_deg} degron contrasts available, expected 9. "
             f"cryptic_missing.tsv or cryptic_all.tsv is absent.")
lab = {"d1aid": "DNMT1-AID", "u1aid": "UHRF1-AID", "dual": "dual degron"}
col = {"d1aid": BLUE, "u1aid": GREEN, "dual": ORANGE}
for k in [x for x in ("d1aid", "u1aid", "dual") if x in series]:
    pts = sorted(series[k])
    C.plot([d for d, _ in pts], [v for _, v in pts], "o-", ms=4,
           color=col.get(k, INK), label=lab.get(k, k))
C.axhline(0, color=INK, lw=0.6, ls=":")
C.set_xlabel("days of depletion", fontsize=8)
C.set_ylabel("3'-coverage tilt\n(within-gene downstream gain)", fontsize=8)
C.set_xticks([2, 6, 8])
C.legend(fontsize=7.5, frameon=False, loc="upper left")

# ---- D: latent TSS test, parsed from latent_tss.txt
txt = open(need("latent_tss.txt")).read()
obs, nul = [], []
for line in txt.splitlines():
    m = re.search(r"\(([\d.]+)%\).*?null \(\+50kb\): ([\d.]+)%", line)
    if m:
        obs.append(float(m.group(1))); nul.append(float(m.group(2)))
if len(obs) != 2:
    sys.exit(f"FATAL: parsed {len(obs)} rows from latent_tss.txt, expected 2")
x = np.arange(2)
D.bar(x - 0.19, obs, width=0.36, color=ORANGE, label="observed")
D.bar(x + 0.19, nul, width=0.36, color="#b9cbe4", label="positional null")
D.set_xticks(x, ["DKO-novel\njunction ends", "expressed\ncontrols"],
             fontsize=8)
D.set_ylabel("ends within 300 bp of\nWT PRO-cap peak (%)", fontsize=8)
for i, (o, n) in enumerate(zip(obs, nul)):
    D.text(i, o + 0.25, f"x{o/n:.1f}", ha="center", fontsize=8, color=INK)
D.legend(fontsize=8, frameon=False, loc="upper left")
D.set_ylim(0, max(obs) * 1.35)

out = os.path.join(ROOT, "figures", "msf_fig4_cryptic")
fig.savefig(out + ".pdf"); fig.savefig(out + ".png", dpi=170)
print("wrote", out)
print(f"  A: {g1[0]}/{g2[0]} x{rat[0]}   B: {vals[0]:.2f}/{vals[1]:.2f} "
      f"nulls {nulls[0]:.2f}/{nulls[1]:.2f}   C: {n_deg} degron contrasts   "
      f"D: {obs} vs {nul}")
