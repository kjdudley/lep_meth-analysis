#!/usr/bin/env python3
"""Manuscript figures for the COMPARATIVE draft.

Conventions are Molecular Ecology / comparative genomics standard: sans-serif,
7-8 pt, despined axes, no gridlines heavier than hairline, vector PDF plus a
PNG proof, single column 85 mm and double column 170 mm.

Palette is Okabe-Ito, the field standard for colourblind-safe scientific
figures. The three caddisfly hues were validated (all six checks pass, worst
adjacent CVD dE 11.0 deutan). H. armigera is deliberately NEUTRAL grey: it is
the reference species, not a fourth category, so it must not read as a peer.
That is why the four-colour set fails a categorical chroma check and the
three-colour set passes.

Every figure carries a legend and direct species labels, so identity is never
carried by colour alone.

NO DUAL AXES. The first plan for Figure 3 put enrichment and absolute rate on
twinned y-axes; that is the most common serious chart error and it is replaced
here by two panels sharing a category axis.
"""
import csv, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

S = "/tmp/claude-1000/-home-eldudy-lep-meth/14f2dac2-87c2-49f7-bfc8-867bf4195ecc/scratchpad"
R = "/home/eldudy/lep_meth/results/lowmeth/hi70"
REG = "/home/eldudy/lep_meth/results/qc/panel_registry.tsv"
MM = 1/25.4

plt.rcParams.update({
    "hatch.linewidth": 0.6,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    # Heredity house style: "Do not make rules thinner than 1pt (0.36mm)."
    "axes.linewidth": 1.0, "xtick.major.width": 1.0, "ytick.major.width": 1.0,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "legend.frameon": False, "savefig.dpi": 600, "figure.dpi": 150,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
C = {"molanna": "#0072B2", "ceraclea": "#009E73", "philopotamus": "#D55E00", "harmi": "#595959"}
DASH = {"molanna": (None, None), "philopotamus": (4.5, 1.6),
        "ceraclea": (1.4, 1.4), "harmi": (6, 1.6, 1.4, 1.6)}
MARK = {"molanna": "o", "philopotamus": "s", "ceraclea": "^", "harmi": "D"}
HATCH = {"molanna": "", "philopotamus": "///", "ceraclea": "...", "harmi": "\\\\"}
def blend(hexcol, frac, bg=(1.0, 1.0, 1.0)):
    """Pre-blend towards the page so EPS, which has no transparency, matches the PDF.

    matplotlib's PostScript backend renders partially transparent artists OPAQUE,
    which turned Figure 4's interquartile band into a solid block over its own
    median line. Blending here keeps all three output formats identical.
    """
    r, g, b = (int(hexcol[i:i+2], 16)/255 for i in (1, 3, 5))
    return tuple(c*frac + s*(1-frac) for c, s in zip((r, g, b), bg))

def dashed(sp):
    d = DASH[sp]
    return {} if d[0] is None else {"dashes": list(d)}
NAME = {"molanna": "M. angustata", "ceraclea": "C. dissimilis",
        "philopotamus": "P. montanus", "harmi": "H. armigera"}
ORDER = ["molanna", "philopotamus", "ceraclea", "harmi"]   # fixed, never cycled
INK, MUTED = "#1a1a1a", "#666666"

def despine(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelcolor=INK)

def panel(ax, letter, x=-0.16):
    ax.text(x, 1.04, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color=INK)

def tsv(p):
    """Read a headed TSV, tolerating a stray header row anywhere in the file.

    figdata.pbs sorted its own output with `sort -k1,1`, which moved the header
    line to line 65 of fig_covbins.tsv. The data is fine; the header just is not
    first. Rather than trust position, take the first line as the field names and
    drop any later row that repeats them.
    """
    rows = [l.rstrip("\n").split("\t") for l in open(p) if l.strip()]
    hdr = max(rows, key=lambda r: sum(1 for c in r if c and not c.replace(".", "").isdigit()))
    hdr = rows[0] if len(set(rows[0])) == len(rows[0]) and not rows[0][1].replace(".", "").isdigit() else hdr
    return [dict(zip(hdr, r)) for r in rows if r != hdr and len(r) == len(hdr)]

def save(fig, name):
    for ext in ("pdf", "eps", "png"):
        fig.savefig(f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {name}.pdf / .eps / .png")

# ---------------------------------------------------------------- Figure 1
dec, cov = tsv(f"{S}/fig_deciles.tsv"), tsv(f"{S}/fig_covbins.tsv")
MODAL = {"molanna": 20, "philopotamus": 50, "ceraclea": 50, "harmi": 80}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(169*MM, 62*MM))
for sp in ORDER:
    d = sorted([r for r in dec if r["species"] == sp], key=lambda r: int(r["decile"]))
    x = [int(r["decile"]) + 5 for r in d]
    y = [float(r["pct_of_sites"]) for r in d]
    a1.plot(x, y, "-", marker=MARK[sp], color=C[sp], lw=1.4, ms=3.4, mew=0,
            label=NAME[sp], zorder=3 if sp != "harmi" else 2, **dashed(sp))
a1.set_yscale("log")
a1.set_xlabel("Site methylation (%)"); a1.set_ylabel("Sites in bin (% of all CpGs)")
a1.set_xticks(range(0, 101, 20)); a1.set_xlim(-3, 103)
a1.set_ylim(0.03, 400)
a1.legend(loc="upper right", ncol=1, handlelength=1.4, columnspacing=1.0,
          labelcolor=INK, borderpad=0.2, bbox_to_anchor=(1.02, 1.03))
despine(a1); panel(a1, "(a)")

for sp in ORDER:
    d = sorted([r for r in cov if r["species"] == sp], key=lambda r: int(r["cov_bin"]))
    d = [r for r in d if int(r["n_sites"]) >= 50_000 and int(r["cov_bin"]) <= 130
         and float(r["pct_ge70"]) > 0]        # log axis cannot show an exact zero
    x = [int(r["cov_bin"]) + 5 for r in d]
    y = [float(r["pct_ge70"]) for r in d]
    a2.plot(x, y, "-", color=C[sp], lw=1.4, zorder=3, **dashed(sp))
    m = MODAL[sp]
    mr = [r for r in d if int(r["cov_bin"]) == m]
    if mr:
        a2.plot(m+5, float(mr[0]["pct_ge70"]), MARK[sp], color=C[sp], ms=6, mew=1.1,
                mec="white", zorder=4)
# floor panel: native Diptera, ULI and C. loewi excluded
reg = {r["panel_label"]: r for r in tsv(REG)}
nat = [l for l, r in reg.items() if r["clade"] == "Diptera"
       and r["amplified_uli"] == "no" and l != "Calliphora_loewi"]
bands, med = [], []
for b in range(20, 130, 10):
    v = []
    for l in nat:
        rows = {r[1]: r[2] for r in csv.reader(open(f"{R}/{l}.floor.tsv"), delimiter="\t") if len(r) > 2}
        n, p = rows.get(f"cov{b}_n"), rows.get(f"cov{b}_pct")
        if n and p and float(n) >= 200_000: v.append(float(p))
    if len(v) >= 4: bands.append(b+5); med.append(st.median(v))
a2.plot(bands, med, "--", color=MUTED, lw=1.2, zorder=2)
a2.set_ylim(3e-3, 40)
a2.annotate("floor panel\n(11 native Diptera, median)", xy=(bands[2], med[2]),
            xytext=(58, 0.030), fontsize=6.5, color=MUTED, ha="left", va="center",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=1.0,
                            connectionstyle="arc3,rad=-0.2"))
for sp in ORDER:
    m = MODAL[sp]
    rr = [r for r in cov if r["species"] == sp and int(r["cov_bin"]) == m]
    if rr:
        off = {"molanna": (4, -13), "philopotamus": (7, 5),
               "ceraclea": (7, -9), "harmi": (7, 4)}[sp]
        a2.annotate(NAME[sp], xy=(m+5, float(rr[0]["pct_ge70"])), xytext=off,
                    textcoords="offset points", fontsize=6.5, color=INK,
                    style="italic", ha="left")
a2.set_yscale("log"); a2.set_xlabel("Coverage (x)"); a2.set_ylabel("CpGs methylated $\\geq$70% (%)")
a2.set_xlim(15, 132)
despine(a2); panel(a2, "(b)")
fig.text(0.5, -0.07, "Filled markers mark each genome's modal coverage band, where its rate is quoted.",
         ha="center", fontsize=6.5, color=MUTED)
fig.tight_layout(w_pad=2.4)
save(fig, "Figure1")

# ---------------------------------------------------------------- Figure 2
reg = {r["panel_label"]: r for r in tsv(REG)}

def modal_rate(lbl):
    n, p = {}, {}
    for r in csv.reader(open(f"{R}/{lbl}.floor.tsv"), delimiter="\t"):
        if len(r) < 3: continue
        m = r[1]
        if m.startswith("cov") and m.endswith("_n"):   n[int(m[3:-2])] = int(float(r[2]))
        if m.startswith("cov") and m.endswith("_pct"): p[int(m[3:-4])] = float(r[2])
    if not n: return None
    b = max(n, key=lambda k: n[k])
    return (b, n[b], p.get(b))

CLADE_COL = {"Mollusca": "#CC79A7", "Cnidaria": "#E69F00", "Trichoptera": "#0072B2"}
pts = []
for lbl, r in reg.items():
    if r["amplified_uli"] != "no" or lbl == "Calliphora_loewi": continue
    m = modal_rate(lbl)
    if not m or m[2] is None or m[2] <= 0: continue
    pts.append((r["clade"], lbl, m[2]))
CADDIS = [("Trichoptera", "M. angustata", 9.55), ("Trichoptera", "P. montanus", 8.68),
          ("Trichoptera", "C. dissimilis", 5.98)]
pts += CADDIS
order = [c for c in ["Mollusca", "Cnidaria", "Trichoptera", "Lepidoptera",
                     "Hymenoptera", "Diptera"] if any(x[0] == c for x in pts)]

fig, ax = plt.subplots(figsize=(112*MM, 80*MM))
rng = np.random.default_rng(0)
for i, cl in enumerate(order):
    v = [p for c, l, p in pts if c == cl]
    isc = cl == "Trichoptera"
    col = CLADE_COL.get(cl, "#8c8c8c")
    jit = rng.uniform(-0.15, 0.15, len(v)) if len(v) > 3 else np.zeros(len(v))
    ax.scatter(v, [i]*len(v) + jit, s=40 if isc else 17,
               marker="o" if not isc else "D",
               edgecolor="white" if isc else "none",
               linewidth=1.0 if isc else 0,
               facecolors=col if isc else blend(col, 0.7),
               zorder=4 if isc else 3)
    if len(v) >= 4:
        ax.plot([st.median(v)]*2, [i-0.28, i+0.28], color=INK, lw=1.3, zorder=5)

# caddisfly labels: leader lines at staggered heights so three points in one
# decade do not collide
ti = order.index("Trichoptera")
for (c, l, p), dy in zip(sorted(CADDIS, key=lambda x: -x[2]), (0.46, 0.72, 0.98)):
    ax.annotate(l, xy=(p, ti), xytext=(p*1.25, ti - dy), fontsize=6.4, style="italic",
                color=INK, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=1.0, shrinkA=0, shrinkB=2))
hm = [p for c, l, p in pts if l == "ilHelArmi9"][0]
ax.annotate("H. armigera", xy=(hm, order.index("Lepidoptera")), xytext=(hm*0.35, order.index("Lepidoptera") + 0.55),
            fontsize=6.4, style="italic", color=INK, ha="center",
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=1.0, shrinkA=0, shrinkB=2))

ax.set_xscale("log")
ax.set_yticks(range(len(order)))
ax.set_yticklabels([f"{c} ({sum(1 for x in pts if x[0]==c)})" for c in order])
ax.invert_yaxis()
ax.set_ylim(len(order)-0.45, -0.7)
ax.set_xlim(1.2e-3, 60)
ax.set_xlabel("CpGs methylated $\\geq$70% in modal coverage band (%)")
lo, hi = min(p for c, l, p in pts if c == "Diptera"), max(p for c, l, p in pts if c == "Diptera")
ax.axvspan(lo, hi, color="#f2f2f2", zorder=0)
ax.text((lo*hi)**0.5, -0.55, "floor panel", fontsize=6.2, color=MUTED, ha="center")
despine(ax); ax.tick_params(axis="y", length=0)
fig.text(0.5, -0.045, "Each point is one genome at its own modal coverage band. Vertical rules are clade medians\n"
         "where n $\\geq$ 4; clade sizes in parentheses. Amplified and ULI libraries excluded throughout.",
         ha="center", fontsize=6.2, color=MUTED)
fig.tight_layout()
save(fig, "Figure2")

# ---------------------------------------------------------------- Figure 3
# Compartment partition, from meth_partition.py (jobs 25381429 molanna/xylo,
# 25401425 ceraclea, 25401472 harmi, 25401976 philopotamus).
# TWO PANELS, NOT TWINNED AXES. Absolute rate and enrichment are different
# measures on different scales; putting them on one pair of axes would be a
# dual-axis chart, which is the most common serious chart error.
PART = {   # species: (exonic%, intronic%, intergenic%), (exonic x, intronic x, intergenic x)
    "philopotamus": ((51.18, 10.05, 4.34), (6.24, 1.23, 0.53)),
    "molanna":      ((46.01, 14.80, 4.26), (4.33, 1.39, 0.40)),
    "ceraclea":     ((31.93, 6.67, 3.62),  (5.58, 1.17, 0.63)),
    "harmi":        ((1.142, 0.370, 0.128),(3.86, 1.25, 0.43)),
}
CLASSES = ["exonic", "intronic", "intergenic"]
fig, (b1, b2) = plt.subplots(1, 2, figsize=(169*MM, 62*MM))
w, xs = 0.19, np.arange(3)
for k, sp in enumerate(ORDER):
    rates, enr = PART[sp]
    off = (k - 1.5) * w
    b1.bar(xs + off, rates, w*0.88, color=C[sp], label=NAME[sp], zorder=3,
           hatch=HATCH[sp], edgecolor="white", linewidth=1.0)
    b2.bar(xs + off, enr,   w*0.88, color=C[sp], zorder=3,
           hatch=HATCH[sp], edgecolor="white", linewidth=1.0)
b1.set_yscale("log")
b1.set_ylabel("CpGs methylated $\\geq$70% in class (%)")
b1.set_ylim(0.05, 200)
b1.legend(loc="upper right", ncol=2, handlelength=1.0, columnspacing=0.8,
          labelcolor=INK, borderpad=0.2, handletextpad=0.5)
b2.axhline(1.0, color=MUTED, lw=1.0, ls="--", zorder=2)
b2.set_ylabel("Enrichment (class share of methylation /\nclass share of sites)")
b2.set_ylim(0, 7.2)
b2.text(2.42, 1.12, "no enrichment", fontsize=6.2, color=MUTED, ha="right")
for ax, letter in ((b1, "(a)"), (b2, "(b)")):
    ax.set_xticks(xs); ax.set_xticklabels(CLASSES); despine(ax)
    panel(ax, letter, x=-0.22)
    ax.set_axisbelow(True)
fig.text(0.5, -0.06, "(a) absolute rate per compartment, log scale. (b) enrichment relative to each class's share of "
         "covered sites.\nThe profiles coincide; the levels differ roughly thirty to forty fold.",
         ha="center", fontsize=6.3, color=MUTED)
fig.tight_layout(w_pad=2.6)
save(fig, "Figure3")

# ---------------------------------------------------------------- Figure 4
fig, axs = plt.subplots(1, 4, figsize=(169*MM, 52*MM), sharey=True)
for ax, sp, letter in zip(axs, ORDER, ("(a)", "(b)", "(c)", "(d)")):
    rows = []
    for r in csv.DictReader(open(f"{S}/{sp}.gene_oe_meth.tsv"), delimiter="\t"):
        try:
            oe, pc = float(r["cpg_oe"]), float(r["pct_conf"])
            if oe == oe and oe < 5: rows.append((pc, oe))
        except (ValueError, KeyError): pass
    rows.sort()
    n, nb = len(rows), 8
    xs_, med, q1, q3 = [], [], [], []
    for i in range(nb):
        chunk = rows[i*n//nb:(i+1)*n//nb]
        if len(chunk) < 30: continue
        v = sorted(o for _, o in chunk)
        xs_.append(st.median([p for p, _ in chunk]))
        med.append(st.median(v)); q1.append(v[len(v)//4]); q3.append(v[3*len(v)//4])
    ax.fill_between(range(len(xs_)), q1, q3, facecolor="none", edgecolor=C[sp],
                    hatch="////", linewidth=0.0, alpha=0.55, zorder=2)
    ax.fill_between(range(len(xs_)), q1, q3, facecolor=blend(C[sp], 0.10), lw=0, zorder=1)
    ax.plot(range(len(xs_)), med, "-", marker=MARK[sp], color=C[sp], lw=1.5,
            ms=3.6, mew=0, zorder=3, **dashed(sp))
    ax.axhline(1.0, color=MUTED, lw=1.0, ls="--", zorder=1)
    ax.set_title(NAME[sp], style="italic", fontsize=7.5, color=INK, pad=3)
    ax.set_xticks([0, len(xs_)-1]); ax.set_xticklabels(["least", "most"])
    ax.set_xlabel("gene methylation octile")
    despine(ax); panel(ax, letter)
axs[0].set_ylabel("CpG observed / expected")
axs[0].set_ylim(0.4, 1.6)
fig.text(0.5, -0.10, "Median CpG o/e with interquartile range, genes binned by their own measured methylation. "
         "A real historical-methylation\nsignal slopes DOWN: more methylation, more CpG loss. "
         "Dashed line marks the neutral expectation.", ha="center", fontsize=6.3, color=MUTED)
fig.tight_layout(w_pad=1.4)
save(fig, "Figure4")
