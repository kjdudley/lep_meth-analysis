#!/usr/bin/env python3
"""Paper 1 section 2.3: the false-positive floor as a function of THRESHOLD.

Section 2.3 fits the floor against coverage at a fixed >=50% cut. That cut has
never been justified by measurement -- it was chosen because it has an observed
zero false-positive rate, which is an argument for specificity and says nothing
about what it costs in sensitivity. This sweeps the threshold so the choice
becomes a derived optimum rather than an assertion, and so the paper can answer a
question a reader will certainly ask: CAN this chemistry see intermediate
methylation at all?

WHY THE ANSWER IS NOT OBVIOUS. At coverage c a site clears threshold t when
t*c/100 of its reads are called modified. At >=50% and 30x that needs 15 of 30
miscalls, which at the vendor's 0.11 per-read FPR is ~6.8 SD out and effectively
never happens. At >=20% and 30x it needs 6 of 30 -- ~1.6 SD -- which is order
several percent. So the floor rises steeply as the threshold falls, while the
real signal is a fixed set of sites that does not. Where those two curves cross
is an empirical question, it depends on coverage, and this measures it.

THE DISCRIMINATOR THAT MATTERS. Lowering the threshold does two different things
depending on the biology, and they look completely different here:

  * if methylation is genuinely bimodal (heavily methylated sites plus
    unmethylated ones), a lower cut adds floor and no signal -- the moth rate and
    the unmethylated rate rise together
  * if there is a real INTERMEDIATE population being missed, the moth rate rises
    FASTER than the floor

That distinction cannot be made from the >=50% number alone, and it is the
substantive reason to run this rather than a QC exercise.

TWO SCORE COLUMNS, DELIBERATELY BOTH. pb-CpG-tools emits a modelled score
(column 4) and the raw counts (columns 7 and 8). They disagree, and they disagree
most where a low threshold lives. A real row from ilNocJanh1:

    OY755051.1  12182  12183  13.3  Total  13  1  12  7.7
                                    ^model            ^raw = 1/13

One modified read of thirteen reads as 7.7% raw and 13.3% modelled. At >=50% that
gap is irrelevant because either way the evidence is overwhelming. At >=20% the
model's prior may be doing most of the work, so a conclusion that holds for the
model and not the counts is a conclusion about pb-CpG-tools, not about the
genome. `lowmeth_test.py` reads the model column only; this reads both.

RAW IS RECOMPUTED FROM COUNTS, not taken from column 9, so the two are
cross-checked rather than assumed equal. Disagreements are counted and reported.

MIN COVERAGE IS A TRAP AT LOW THRESHOLDS. At >=20% a 4x site needs 1 read, so a
single miscall makes it "methylated". Output is coverage-binned so those bins are
separable downstream, but do not pool across bins at low thresholds.

No numpy: this runs under any python3, which avoids the env2/env_phase trap that
cost three jobs on 2026-08-12.

Usage:
  floor_sweep.py --bed X.combined.bed.gz --label X --mt CONTIG --out X.sweep.tsv
"""
import argparse
import gzip
from collections import defaultdict

THRESHOLDS = [10, 20, 30, 40, 50, 60, 70]
COV_CAP = 200          # top bin is "200x and above"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mt", default="", help="mitochondrial contig, reported apart")
    ap.add_argument("--min-cov", type=int, default=4)
    ap.add_argument("--thresholds", default=",".join(str(t) for t in THRESHOLDS))
    a = ap.parse_args()
    ths = [float(t) for t in a.thresholds.split(",")]

    # (score_type, threshold, cov_bin) -> hits ;  cov_bin -> n
    hits = defaultdict(int)
    n_by_bin = defaultdict(int)
    mt_hits = defaultdict(int)
    mt_n = 0
    n_rows = n_skipped = n_disagree = n_noraw = 0

    op = gzip.open if a.bed.endswith(".gz") else open
    with op(a.bed, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                n_skipped += 1
                continue
            try:
                model = float(f[3])
                cov = float(f[5])
            except ValueError:
                n_skipped += 1
                continue
            if cov < a.min_cov:
                continue
            n_rows += 1

            # Recompute raw from the counts rather than trusting column 9, and
            # cross-check. A silent format change would otherwise be invisible.
            raw = None
            if len(f) >= 9:
                try:
                    mod, unmod = float(f[6]), float(f[7])
                    tot = mod + unmod
                    if tot > 0:
                        raw = 100.0 * mod / tot
                        if abs(raw - float(f[8])) > 0.15:
                            n_disagree += 1
                except ValueError:
                    raw = None
            if raw is None:
                n_noraw += 1

            b = min(int(cov // 10) * 10, COV_CAP)

            if a.mt and f[0] == a.mt:
                mt_n += 1
                for t in ths:
                    if model >= t:
                        mt_hits[("model", t)] += 1
                    if raw is not None and raw >= t:
                        mt_hits[("raw", t)] += 1
                continue

            n_by_bin[b] += 1
            for t in ths:
                if model >= t:
                    hits[("model", t, b)] += 1
                if raw is not None and raw >= t:
                    hits[("raw", t, b)] += 1

    with open(a.out, "w") as w:
        w.write("label\tscore_type\tthreshold\tcov_bin\tn\thits\tpct\n")
        for st in ("model", "raw"):
            for t in ths:
                for b in sorted(n_by_bin):
                    n = n_by_bin[b]
                    h = hits[(st, t, b)]
                    w.write(f"{a.label}\t{st}\t{t:g}\t{b}\t{n}\t{h}"
                            f"\t{100.0 * h / n if n else 0:.6f}\n")
        # mitochondrion apart: it is a sanity check, never the floor (2.5)
        if a.mt and mt_n:
            for st in ("model", "raw"):
                for t in ths:
                    h = mt_hits[(st, t)]
                    w.write(f"{a.label}\t{st}_mt\t{t:g}\tMT\t{mt_n}\t{h}"
                            f"\t{100.0 * h / mt_n if mt_n else 0:.6f}\n")

    print(f"{a.label}: {n_rows} sites at >={a.min_cov}x, "
          f"{len(n_by_bin)} coverage bins, mito {mt_n}")
    if n_noraw:
        print(f"  !! {n_noraw} rows had no usable counts -- raw rows are NOT "
              f"comparable to model rows for this species")
    if n_disagree:
        print(f"  !! {n_disagree} rows where recomputed raw disagrees with "
              f"column 9 by >0.15pp -- check the bed format before using raw")
    if n_skipped:
        print(f"  ({n_skipped} malformed rows skipped)")
    print(f"  wrote {a.out}")


if __name__ == "__main__":
    main()
