#!/usr/bin/env python3
"""Call DNMT1 / DNMT3 presence per genome from an existing miniprot PAF.

WHY THIS EXISTS. dnmt_queries_v2.faa answers "is there any DNA
methyltransferase?", which is what the calibration paper needs. It cannot
answer "DNMT1 or DNMT3?", because its headers carry no family: every
methyltransferase is `DNMT|<acc>|<species>|<len>`. Comparing methylation motifs
against methyltransferase complement needs the family.

Two things make the family call harder than reading a best hit.

  1. LENGTH DOES NOT WORK. UniProt assigns Tribolium A0A1Z2TI83 (729 aa) to
     DNMT1 and Apis A0A7M7MWZ4 (1185 aa) to Dnmt3. Splitting the same query set
     at 1000 aa puts two genuine DNMT1 sequences in the DNMT3 bucket, and then
     every DNMT1-bearing genome scores as DNMT1+DNMT3. That is what happened on
     the first attempt: all 37 Lepidoptera came back with DNMT3, which
     contradicts the literature and the premise of the analysis being scaled.

  2. ONE BAD QUERY IS ENOUGH. Even with correct labels, Lepidoptera still
     scored DNMT3-positive on a single query, A0ACB9SVW9, a TrEMBL entry
     auto-annotated "methyltransferase 3-related" with no real gene symbol. It
     hits a locus no other DNMT query touches. A real DNMT3 is hit by the Apis,
     Blattella, Locusta, Bemisia and Reticulitermes DNMT3 queries together.

So presence requires TWO INDEPENDENT QUERIES OF THE SAME FAMILY AGREEING ON A
LOCUS. Two queries hitting two different places is not corroboration:
Achroia grisella has two DNMT3 hits at three loci and no DNMT3.

Validated against known biology on the 65 screened panel genomes: Lepidoptera
DNMT1 only, Hymenoptera DNMT1+DNMT3, Diptera and Nematoda neither.

    dnmt_family_call.py --paf-dir results/qc/dnmt --queries hpc/queries/dnmt_queries_v3.faa
"""
import argparse, collections, glob, os


def families(faa):
    fam = {}
    for line in open(faa):
        if line.startswith(">"):
            p = line[1:].strip().split("|")
            if p[0] in ("DNMT1", "DNMT3"):
                fam[p[1]] = p[0]
    return fam


def call(paf, fam, min_qcov=0.5, window=50000, min_queries=2):
    """Distinct queries per family clustered by locus; a family is present if
    any locus is supported by >= min_queries distinct queries."""
    hits = collections.defaultdict(list)
    for line in open(paf):
        p = line.rstrip("\n").split("\t")
        if len(p) < 12 or not p[0].startswith("DNMT"):
            continue
        acc = p[0].split("|")[1]
        f = fam.get(acc)
        if f is None:
            continue
        if (int(p[3]) - int(p[2])) / int(p[1]) < min_qcov:
            continue
        hits[f].append((p[5], int(p[7]), acc))
    out = {}
    for f, hs in hits.items():
        best = 0
        for ctg, pos, _ in hs:
            near = {a for c, q, a in hs if c == ctg and abs(q - pos) <= window}
            best = max(best, len(near))
        out[f] = best
    return out.get("DNMT1", 0), out.get("DNMT3", 0), min_queries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf-dir", required=True)
    ap.add_argument("--queries", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    fam = families(a.queries)
    if not fam:
        raise SystemExit("FATAL: no family-labelled queries; is this v3?")
    rows = []
    for f in sorted(glob.glob(os.path.join(a.paf_dir, "*.dnmt.paf"))):
        lab = os.path.basename(f)[:-len(".dnmt.paf")]
        d1, d3, m = call(f, fam)
        g = ("DNMT1+DNMT3" if d1 >= m and d3 >= m else
             "DNMT1 only" if d1 >= m else
             "DNMT3 only" if d3 >= m else "neither")
        rows.append((lab, d1, d3, g))
    fh = open(a.out, "w") if a.out else None
    hdr = "label\tdnmt1_queries_at_locus\tdnmt3_queries_at_locus\tgenotype"
    print(hdr) if not fh else fh.write(hdr + "\n")
    for lab, d1, d3, g in rows:
        line = f"{lab}\t{d1}\t{d3}\t{g}"
        print(line) if not fh else fh.write(line + "\n")
    if fh:
        fh.close()
        print(f"wrote {a.out}  ({len(rows)} genomes)")


if __name__ == "__main__":
    main()
