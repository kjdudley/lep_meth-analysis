#!/usr/bin/env python3
"""Read the per-species sweeps and show where signal and floor cross.

The sweep files are per-species and long-format; this is the table that answers
the question. For each threshold and coverage band it prints the mean floor
across the unmethylated genomes beside the rate in a known-methylated moth, and
their ratio.

HOW TO READ IT. The ratio is signal-to-floor. It should be large at >=50%, where
the floor is essentially zero, and collapse as the threshold falls. The threshold
where it stops being large is the sensitivity limit of this chemistry, and it is
coverage-dependent -- which is the whole point of reporting it per band.

AND THE DISCRIMINATOR. Watch the moth column on its own as the threshold falls:

  * rises no faster than the floor  -> methylation here is bimodal, and the >=50%
    cut is not costing real signal. The intermediate sites are not there.
  * rises FASTER than the floor     -> there IS an intermediate population that
    the >=50% cut discards, and section 2.2's estimates are lower bounds by more
    than the paper currently says.

Both are publishable; they are different claims and the sweep separates them.

No numpy, for the same reason as floor_sweep.py.

Usage:
  floor_sweep_summary.py --dir lowmeth/sweep --unmeth A,B,C --meth ilNocJanh1
"""
import argparse
import glob
import os
from collections import defaultdict

BANDS = [30, 50, 70]


def load(path):
    """-> {(score_type, threshold, cov_bin): (n, hits, pct)}"""
    out = {}
    with open(path) as fh:
        hdr = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if hdr is None:
                hdr = f
                continue
            if len(f) < 7:
                continue
            try:
                out[(f[1], float(f[2]), f[3])] = (int(f[4]), int(f[5]), float(f[6]))
            except ValueError:
                continue
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--unmeth", required=True, help="comma-separated labels")
    ap.add_argument("--meth", required=True, help="comma-separated labels")
    ap.add_argument("--score-type", default="model", choices=["model", "raw"])
    a = ap.parse_args()

    data = {}
    for f in sorted(glob.glob(os.path.join(a.dir, "*.sweep.tsv"))):
        data[os.path.basename(f).replace(".sweep.tsv", "")] = load(f)

    unmeth = [s for s in a.unmeth.split(",") if s in data]
    meth = [s for s in a.meth.split(",") if s in data]
    missing = [s for s in a.unmeth.split(",") + a.meth.split(",") if s not in data]
    print(f"sweeps loaded: {len(data)}   unmethylated: {len(unmeth)}   "
          f"methylated: {len(meth)}   score column: {a.score_type}")
    if missing:
        print(f"  WARNING absent: {','.join(missing)}")
    if not unmeth or not meth:
        raise SystemExit("need at least one of each to compute a ratio")

    ths = sorted({k[1] for d in data.values() for k in d if k[0] == a.score_type})

    for band in BANDS:
        b = str(band)
        print(f"\n=== coverage {band}-{band + 9}x "
              f"({a.score_type} score) ===")
        print(f"{'thresh':>7}{'floor %':>11}{'floor SD':>11}"
              f"{'moth %':>11}{'ratio':>10}{'genomes':>9}")
        for t in ths:
            fl = [data[s][(a.score_type, t, b)][2]
                  for s in unmeth if (a.score_type, t, b) in data[s]]
            mo = [data[s][(a.score_type, t, b)][2]
                  for s in meth if (a.score_type, t, b) in data[s]]
            if not fl or not mo:
                continue
            fm = sum(fl) / len(fl)
            mm = sum(mo) / len(mo)
            sd = (sum((x - fm) ** 2 for x in fl) / len(fl)) ** 0.5
            ratio = (mm / fm) if fm > 0 else float("inf")
            rs = f"{ratio:>10.1f}" if ratio != float("inf") else f"{'inf':>10}"
            print(f"{t:>7.0f}{fm:>11.4f}{sd:>11.4f}{mm:>11.4f}{rs}{len(fl):>9}")

    # Between-genome spread is what sets the limit of detection (see 2.3), so
    # report the LoD at each threshold rather than only the mean floor.
    print(f"\n=== limit of detection = mean floor + 3 SD, {a.score_type} score ===")
    print(f"{'thresh':>7}" + "".join(f"{str(b) + '-' + str(b + 9) + 'x':>14}"
                                     for b in BANDS))
    for t in ths:
        row = f"{t:>7.0f}"
        for band in BANDS:
            b = str(band)
            fl = [data[s][(a.score_type, t, b)][2]
                  for s in unmeth if (a.score_type, t, b) in data[s]]
            if not fl:
                row += f"{'-':>14}"
                continue
            fm = sum(fl) / len(fl)
            sd = (sum((x - fm) ** 2 for x in fl) / len(fl)) ** 0.5
            row += f"{fm + 3 * sd:>14.4f}"
        print(row)


if __name__ == "__main__":
    main()
