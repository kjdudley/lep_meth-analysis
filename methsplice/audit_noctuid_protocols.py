#!/usr/bin/env python3
"""Calibration audit of a deposited-methylation species set.

    audit_noctuid_protocols.py --runs runs.tsv --out audit.tsv

For each run: pull ENA library_construction_protocol and read counts, and
classify the library as NATIVE (5mC preserved) or AMPLIFIED/ULI (5mC erased
by PCR, but MM/ML tags deposited anyway). The amplified case is the
calibration paper's §2.1b trap: nothing in the structured metadata
distinguishes it, and an erased genome reads as biologically unmethylated.

Written for the noctuid splice-site methylation set, where one species in
the comparison turned out to be amplified.
"""
import argparse, csv, re, sys, urllib.request

FIELDS = ("run_accession,scientific_name,library_construction_protocol,"
          "instrument_model,read_count,base_count,study_accession")
AMPLIFIED = re.compile(
    r"\b(ULI|ultra[- ]?low|low[- ]?input|amplif|PCR[- ]amplif|"
    r"whole[- ]genome amplif|WGA|Ultra-Low Input)\b", re.I)
NATIVE_HINT = re.compile(r"\bnative|non[- ]amplif|no PCR|PCR[- ]free\b", re.I)


def fetch(run):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={run}"
           f"&result=read_run&fields={FIELDS}&format=tsv")
    for _ in range(4):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                rows = r.read().decode().rstrip("\n").split("\n")
            if len(rows) >= 2:
                return dict(zip(rows[0].split("\t"), rows[1].split("\t")))
        except Exception:
            continue
    return None


def classify(protocol):
    if not protocol.strip():
        return "UNKNOWN", "no protocol text deposited"
    if AMPLIFIED.search(protocol):
        return "AMPLIFIED", "PCR/ULI wording present — 5mC erased"
    if NATIVE_HINT.search(protocol):
        return "NATIVE", "explicitly native/PCR-free"
    return "NATIVE(assumed)", "no amplification wording found"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True,
                    help="TSV with a 'run' column (or run in col 1)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs = []
    with open(a.runs) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        head = next(rdr)
        idx = head.index("run") if "run" in head else 0
        if not head[idx].lower().startswith("run"):
            runs.append(head[idx])
        for row in rdr:
            if row and row[idx].strip():
                runs.append(row[idx].strip())
    print(f"auditing {len(runs)} runs")
    out = open(a.out, "w")
    out.write("run\tspecies\tverdict\treason\tinstrument\tgbase\tprotocol\n")
    counts = {}
    for r in runs:
        rec = fetch(r)
        if rec is None:
            print(f"  {r}: ENA lookup FAILED")
            out.write(f"{r}\tNA\tLOOKUP_FAILED\tENA did not respond\tNA\tNA\t\n")
            continue
        proto = rec.get("library_construction_protocol", "") or ""
        verdict, why = classify(proto)
        counts[verdict] = counts.get(verdict, 0) + 1
        try:
            gb = f"{int(rec.get('base_count') or 0)/1e9:.1f}"
        except ValueError:
            gb = "NA"
        sp = rec.get("scientific_name", "NA")
        flag = "  <<< EXCLUDE" if verdict == "AMPLIFIED" else ""
        print(f"  {r}  {sp:32s} {verdict:16s} {gb:>6s} Gb{flag}")
        out.write(f"{r}\t{sp}\t{verdict}\t{why}\t"
                  f"{rec.get('instrument_model','NA')}\t{gb}\t"
                  f"{proto[:300].replace(chr(9),' ')}\n")
    out.close()
    print("\nsummary: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
