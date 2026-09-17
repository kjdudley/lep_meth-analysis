#!/usr/bin/env python3
"""Does CpG o/e track MEASURED methylation, within one genome?

WHY THIS EXISTS. `cpgoe_gc_control.py` showed that per-gene o/e over BUSCO
orthologues fails GC control. But BUSCO genes are a poor substrate for this
metric and that result cannot distinguish two explanations:

  (a) CpG o/e does not track methylation, or
  (b) BUSCO genes cannot show it, because they are under strong purifying
      selection (which removes the CpG->TpG mutations o/e counts) and are, in
      insects, the methylated class (so the unmethylated comparison class is
      absent by construction).

This test separates them. It stays inside ONE genome, uses ALL annotated genes
rather than universal orthologues, and classifies them by **our own measured
methylation** rather than by species identity. Cross-species composition never
enters, and both methylation classes are present across a range of selective
constraint.

  IF o/e separates measured-methylated from measured-unmethylated genes at
  matched GC -> the metric works, and BUSCO genes were the wrong substrate.
  IF it does not -> the metric does not track methylation, and the BUSCO result
  was not an artefact of its gene set.

Either answer is worth having and neither is currently known.

TWO SUBCOMMANDS:

  bed   --gff X.gff3.gz --fai ref.fa.fai --out genes.bed
        Gene spans in the FASTA/BED namespace, ready for cpg_oe.py --bed.

  join  --regions X.regions.tsv --cpg X.combined.bed.gz --out table.tsv
        Joins per-gene o/e to per-gene measured methylation, then analyses.

NAMESPACE. Ensembl Rapid Release renames placed chromosomes to bare 1..30/W/Z
while the CpG bed and FASTA use INSDC accessions (OX424490.1). Ignoring this
once made an entire exon analysis read 507 CpGs out of 28.8 M and still exit 0.
Sequences are matched on unambiguous LENGTH, exactly as exon_enrichment.py does.
"""
import argparse, gzip, sys
from collections import defaultdict
import numpy as np


def op(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def gff_lengths(gff):
    lens = {}
    for line in op(gff):
        if line.startswith("##sequence-region"):
            f = line.split()
            if len(f) >= 4:
                try: lens[f[1]] = int(f[3])
                except ValueError: pass
        elif line and line[0] != "#":
            break
    return lens


def alias_gff_to_fasta(fai, gff_lens):
    """GFF seqid -> FASTA name, matched on unambiguous length (both directions)."""
    fa = {}
    for line in open(fai):
        f = line.split("\t")
        if len(f) >= 2:
            fa.setdefault(int(f[1]), []).append(f[0])
    bylen = defaultdict(list)
    for sid, L in gff_lens.items():
        bylen[L].append(sid)
    out = {}
    for L, sids in bylen.items():
        names = fa.get(L, [])
        if len(sids) == 1 and len(names) == 1:      # unambiguous only
            out[sids[0]] = names[0]
    return out


def cmd_bed(a):
    lens = gff_lengths(a.gff)
    alias = alias_gff_to_fasta(a.fai, lens) if a.fai else {}
    n_out = n_skip = 0
    with open(a.out, "w") as w:
        for line in op(a.gff):
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            seq = alias.get(f[0], f[0])
            gid = "NA"
            for kv in f[8].split(";"):
                if kv.startswith("ID="):
                    gid = kv[3:]; break
            try:
                s, e = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            if e - s < a.min_len:
                n_skip += 1; continue
            w.write(f"{seq}\t{s}\t{e}\t{gid}\n")
            n_out += 1
    sys.stderr.write(f"genes written {n_out}, skipped short {n_skip}, "
                     f"aliased seqids {len(alias)}\n")
    if alias and n_out == 0:
        sys.exit("no genes written -- check the GFF")


def cmd_join(a):
    # ---- per-gene o/e and GC from cpg_oe.py --bed
    oe, gc, spans = {}, {}, {}
    for line in open(a.regions):
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[0] == "label":
            continue
        try:
            gc[f[1]] = float(f[6]); oe[f[1]] = float(f[7])
            spans[f[1]] = (f[2], int(f[3]), int(f[4]))
        except ValueError:
            continue
    sys.stderr.write(f"genes with o/e: {len(oe)}\n")

    # ---- per-gene measured methylation from the CpG bed
    # Interval index per sequence, so the bed is streamed once.
    idx = defaultdict(list)
    for g, (s, st, en) in spans.items():
        idx[s].append((st, en, g))
    for s in idx:
        idx[s].sort()
    starts = {s: np.array([x[0] for x in v]) for s, v in idx.items()}

    cov = defaultdict(int); hi = defaultdict(int)
    for line in op(a.cpg):
        if line[0] == "#":
            continue
        f = line.split("\t")
        if len(f) < 6:
            continue
        s = f[0]
        if s not in idx:
            continue
        try:
            pos = int(f[1]); pct = float(f[3]); depth = float(f[5])
        except ValueError:
            continue
        # OUR OWN §2.2 RULE: only sites in the validated coverage band count.
        if not (a.cov_lo <= depth < a.cov_hi):
            continue
        j = int(np.searchsorted(starts[s], pos, side="right")) - 1
        for k in range(max(0, j - 3), min(len(idx[s]), j + 1)):
            st, en, g = idx[s][k]
            if st <= pos < en:
                cov[g] += 1
                if pct >= 50: hi[g] += 1
                break

    rows = []
    for g in oe:
        if cov[g] >= a.min_cpg:
            rows.append((g, gc[g], oe[g], cov[g], hi[g], 100.0 * hi[g] / cov[g]))
    sys.stderr.write(f"genes with >={a.min_cpg} CpGs in {a.cov_lo}-{a.cov_hi}x: {len(rows)}\n")
    if len(rows) < 200:
        sys.exit("too few genes pass the CpG/coverage filter to analyse")

    with open(a.out, "w") as w:
        w.write("gene\tgc\tcpg_oe\tn_cpg\tn_conf\tpct_conf\n")
        for r in rows:
            w.write(f"{r[0]}\t{r[1]:.5f}\t{r[2]:.5f}\t{r[3]}\t{r[4]}\t{r[5]:.4f}\n")

    G = np.array([r[1] for r in rows]); O = np.array([r[2] for r in rows])
    P = np.array([r[5] for r in rows])

    print("=" * 74)
    print(f"WITHIN-GENOME: does o/e track MEASURED methylation?   n={len(rows)} genes")
    print("=" * 74)
    print(f"  corr(o/e, GC)            = {np.corrcoef(O, G)[0,1]:+.3f}   <- the known confound")
    print(f"  corr(o/e, %conf-meth)    = {np.corrcoef(O, P)[0,1]:+.3f}   <- the signal")
    print("  A real methylation signal is NEGATIVE: more methylation, more")
    print("  historical CpG loss, lower o/e.")

    # dose-response across methylation, GC held in bins
    print(f"\n  median o/e by methylation quartile, WITHIN GC bins")
    qs = np.percentile(P, [25, 50, 75])
    qlab = ["Q1 (least)", "Q2", "Q3", "Q4 (most)"]
    gbins = [(0.30, 0.38), (0.38, 0.44), (0.44, 0.52)]
    print(f"  {'GC bin':<14}" + "".join(f"{q:>13}" for q in qlab) + f"{'Q4-Q1':>10}")
    for lo, hi_ in gbins:
        gm = (G >= lo) & (G < hi_)
        if gm.sum() < 100: continue
        cells, meds = [], []
        for i in range(4):
            if i == 0:   pm = P <= qs[0]
            elif i == 3: pm = P > qs[2]
            else:        pm = (P > qs[i-1]) & (P <= qs[i])
            m = gm & pm
            if m.sum() >= 25:
                v = float(np.median(O[m])); meds.append(v)
                cells.append(f"{v:.3f}({m.sum()})")
            else:
                meds.append(np.nan); cells.append("   -   ")
        d = (meds[3] - meds[0]) if np.isfinite(meds[0]) and np.isfinite(meds[3]) else np.nan
        print(f"  {f'{lo:.2f}-{hi_:.2f}':<14}" + "".join(f"{c:>13}" for c in cells) + f"{d:>+10.3f}")
    print("\n  Q4-Q1 NEGATIVE in every GC bin => o/e tracks methylation with GC held.")
    print("  Q4-Q1 near zero or inconsistent => it does not, and the BUSCO result")
    print("  was not an artefact of that gene set.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("bed"); b.add_argument("--gff", required=True)
    b.add_argument("--fai"); b.add_argument("--out", required=True)
    b.add_argument("--min-len", type=int, default=200); b.set_defaults(f=cmd_bed)
    j = sub.add_parser("join"); j.add_argument("--regions", required=True)
    j.add_argument("--cpg", required=True); j.add_argument("--out", required=True)
    j.add_argument("--min-cpg", type=int, default=20)
    j.add_argument("--cov-lo", type=float, default=30); j.add_argument("--cov-hi", type=float, default=80)
    j.set_defaults(f=cmd_join)
    a = ap.parse_args(); a.f(a)
