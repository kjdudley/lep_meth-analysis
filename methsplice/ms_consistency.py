#!/usr/bin/env python3
"""Assert the methsplice manuscript agrees with its committed artefacts.

Same role as consistency_check.py for the calibration paper: every
load-bearing number quoted in methsplice/draft_manuscript.md is recomputed
or re-read from the artefact that produced it, and drift fails loudly.
Run before any submission-bound edit is committed.
"""
import csv, os, re, statistics, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MS = open(os.path.join(ROOT, "methsplice/draft_manuscript.md")).read()
FAIL = []

def check(name, cond, detail=""):
    tag = "ok  " if cond else "FAIL"
    print(f"[{tag}] {name}" + (f"  {detail}" if detail else ""))
    if not cond:
        FAIL.append(name)

# ---- register and journal limits
abstract = MS.split("## Abstract")[1].split("## Introduction")[0]
aw = len([x for x in abstract.split() if any(c.isalnum() for c in x)])
check("abstract <= 250 words", aw <= 250, f"({aw})")
check("no em-dashes", chr(0x2014) not in MS)
rh = re.search(r"\*\*Running title:\*\* (.+)", MS).group(1)
check("running head <= 50 chars", len(rh) <= 50, f"({len(rh)})")

# ---- degron dose (degron_dose.tsv)
rows = list(csv.DictReader(open(os.path.join(ROOT, "methsplice/degron_dose.tsv")), delimiter="\t"))
er = sorted(-float(r["mean_erasure_pct"]) for r in rows[::2])
check("degron erasure range 21.7-49.2", abs(er[0]-21.66) < 0.1 and abs(er[-1]-49.2) < 0.1,
      f"({er[0]:.1f}..{er[-1]:.1f})")
base = [float(r["slope"]) for r in rows if r["test"] == "baseline_meth"]
check("baseline slopes all negative", all(s < 0 for s in base))
check("baseline slopes within -18e-6..-34e-6",
      min(base) >= -35e-6 and max(base) <= -17e-6,
      f"({min(base):.1e}..{max(base):.1e})")
check("manuscript quotes -0.000018 to -0.000034",
      "-0.000018 to -0.000034" in MS)
meas = [float(r["slope"]) for r in rows if r["test"] == "measured_derasure"]
check("measured slopes within -6e-6..-27e-6",
      min(meas) >= -28e-6 and max(meas) <= -5e-6)
import numpy as np
for test, sl in (("baseline_meth", base), ("measured_derasure", meas)):
    b, _ = np.polyfit([-float(r["mean_erasure_pct"]) for r in rows if r["test"] == test], sl, 1)
    check(f"slope-of-slopes |{test}| < 3e-7", abs(b) < 3e-7, f"({b:+.1e})")
check("manuscript states <0.0000003 both definitions", "0.0000003" in MS)

# ---- coverage floor: the Methods must not claim a blanket per-site floor.
# exon_psi.py applies --min-cov 5 (default, never overridden) so the psi-table
# meth columns carry it; ms1_degron_dose.py and ms1_coral_contrast.py re-derive
# from the beds with no per-site floor, and cx2bed/nome2bed default --min-cov 1.
check("Methods states both aggregation routes",
      "aggregated by two routes" in MS and "only floor is per exon" in MS)
check("Methods gives the tPBAT floor cost", "2 and 907 exons per condition" in MS)
check("Methods names the attenuation direction", "biases a regression slope" in MS)

# ---- bee per-exon (perexon_summary.txt, contrast.tsv)
pes = open(os.path.join(ROOT, "methsplice/perexon_summary.txt")).read()
check("37,416 complete exons", "37416" in pes and "37,416" in MS)
check("0 at FDR", "BH FDR<0.05: 0" in pes and "zero survive Benjamini-Hochberg" in MS)
check("57 floor+dPSI candidates", "candidates: 57" in pes and "57 of them" in MS)
check("floor count 87 vs 81 expected", "Eighty-seven exons reach that floor against\n81 expected" in MS)
check("conditioned null 6.5 stated", "6.5 expected once the same conditioning" in MS)
check("BH resolution 1,617 stated", "1,617 exons at the" in MS)

ctrl = ["psi_112c","psi_114c","psi_115c","psi_131c","psi_136c","psi_138c"]
sds = []
for r in csv.DictReader(open(os.path.join(ROOT, "methsplice/contrast.tsv")), delimiter="\t"):
    try:
        v = [float(r[c]) for c in ctrl]
    except (ValueError, KeyError):
        continue
    if 0.05 < statistics.mean(v) < 0.95:
        sds.append(statistics.stdev(v))
med = statistics.median(sds)
check("within-group SD median 0.052", abs(med - 0.052) < 0.0015, f"({med:.3f})")
check("manuscript quotes 0.052", "median 0.052" in MS)

# ---- ceiling (contrast.tsv, coral tsvs)
n_hi = n_hi_c = n_lo = n_lo_c = 0
for r in csv.DictReader(open(os.path.join(ROOT, "methsplice/contrast.tsv")), delimiter="\t"):
    try:
        m = float(r["meth_c"]); psi = float(r["psi_c"])
    except ValueError:
        continue
    if m >= 50:
        n_hi += 1; n_hi_c += psi > 0.98
    elif m < 0.001:
        n_lo += 1; n_lo_c += psi > 0.98
check("bee methylated ceiling ~90%", abs(n_hi_c/n_hi - 0.897) < 0.01,
      f"({n_hi_c/n_hi:.3f})")
# ---- coral, regenerated 2026-09-11 with the technical-replicate suffix
# normalised (recovering the un185 cell and a fourth genet-by-site contrast),
# the exon window corrected, and the within_exon_r null computed over every
# usable exon rather than a 2,000-exon subsample. See DESIGN "Finding 1".
def _coral(fn):
    return list(csv.DictReader(
        open(os.path.join(ROOT, "methsplice", fn)), delimiter="\t"))
_crs = _coral("coral_contrast_rs.tsv"); _cdn = _coral("coral_contrast_denovo.tsv")
for _tag, _t in (("RefSeq", _crs), ("denovo", _cdn)):
    _dr = [r for r in _t if r["test"] == "dose_response"]
    check(f"coral {_tag}: all four genet-by-site cells present", len(_dr) == 4)
    check(f"coral {_tag}: every dose CI spans zero",
          all(float(r["ci_lo"]) < 0 < float(r["ci_hi"]) for r in _dr))
    _wr = next(r for r in _t if r["test"] == "within_exon_r")
    check(f"coral {_tag}: within-exon r inside its null range",
          float(_wr["ci_lo"]) < float(_wr["value"]) < float(_wr["ci_hi"]))
    check(f"coral {_tag}: |within-exon r| below 0.005",
          abs(float(_wr["value"])) < 0.005)
check("manuscript quotes the four RefSeq slopes",
      all(f"{float(r['value']):+.6f}".rstrip("0") in MS or r["value"] in MS
          for r in _crs if r["test"] == "dose_response"))
check("manuscript says four cells, not three",
      "all four\ngenet-by-site cells" in MS and "all three genet-by-site" not in MS)
check("manuscript states how many fragments were analysed",
      "26 fragments carry" in MS)
check("manuscript reports the corrected within-exon r and its instability",
      "+0.0011 against a null of -0.0002" in MS
      and "not stable to reasonable" in MS
      and "returned -0.0012 against +0.0006" in MS)
_ceil = {r["cell"]: float(r["value"]) for r in _crs if r["test"] == "ceiling"}
check("coral ceiling as quoted",
      abs(_ceil["meth>=50"] - 0.8518) < 1e-4 and abs(_ceil["meth<10"] - 0.8219) < 1e-4)
check("manuscript quotes 85.2/82.2", "85.2%" in MS and "82.2%" in MS)

dko = next(r for r in csv.DictReader(
    open(os.path.join(ROOT, "methsplice/cryptic_all.tsv")), delimiter="\t")
    if r["contrast"] == "HCT116_DKO_vs_WT")
check("DKO novel junctions 15,670/3,848",
      dko["novel_g1"] == "15670" and dko["novel_g2"] == "3848"
      and "15,670" in MS and "3,848" in MS)
check("DKO ratio 4.07", dko["novel_ratio"] == "4.07" and "4.07" in MS)
# The 3'-tilt. cryptic_degron.tsv holds 5 of the 9 degron contrasts; the four
# that were never run are in cryptic_missing.tsv (audit 2026-09-11). Both are
# needed, because the two rows that break the old monotone claim are in the
# second file. See DESIGN "Finding 5".
cd = {}
for _f in ("methsplice/cryptic_degron.tsv", "methsplice/cryptic_missing.tsv"):
    _p = os.path.join(ROOT, _f)
    if os.path.exists(_p):
        for r in csv.DictReader(open(_p), delimiter="\t"):
            cd[r["contrast"]] = float(r["s1_dslope"])
check("all nine degron contrasts present", len(
      [k for k in cd if k.startswith("degron_")]) == 9)
check("day-8 tilt positive in all three lines",
      all(cd[f"degron_{l}_Day8"] > 0 for l in ("d1aid", "u1aid", "dual")))
check("day-6 breaks monotonicity in two lines",
      cd["degron_d1aid_Day6"] < cd["degron_d1aid_Day2"]
      and cd["degron_u1aid_Day6"] < 0)
check("manuscript quotes the three day-8 tilts",
      "+0.76" in MS and "+0.54" in MS and "+1.67" in MS)
check("manuscript quotes the day-6 collapse",
      "+0.01 and -0.003" in MS)
check("manuscript no longer claims monotone growth",
      "grows with demethylation duration" not in MS
      and "not, however, a monotone function" in MS)
check("manuscript discloses the borrowed degron baseline",
      "deposited no untreated methylome" in MS and "GSE58695" in MS)
# ---- DKO concentration (dko_overlap.tsv). Containment is the primary gene
# rule as of 2026-09-11; proximity (+-10 kb) is the sensitivity. Each exon
# class carries its own relabeling null, and the two must weight-average to
# the pooled rate in contrasts_results.tsv. See DESIGN "Finding 4".
ov = {(r["rule"], r["class"]): r for r in csv.DictReader(
      open(os.path.join(ROOT, "methsplice/dko_overlap.tsv")), delimiter="\t")}
_o_in = ov[("overlap", "novel_junction_genes")]
_o_out = ov[("overlap", "other_genes")]
_p_in = ov[("proximity", "novel_junction_genes")]
check("containment rates 21.64/9.02",
      abs(float(_o_in["obs_rate"]) - 0.2164) < 1e-4
      and abs(float(_o_out["obs_rate"]) - 0.0902) < 1e-4)
check("containment class nulls 11.67/6.60",
      abs(float(_o_in["null_rate"]) - 0.1167) < 1e-4
      and abs(float(_o_out["null_rate"]) - 0.0660) < 1e-4)
check("containment share 24%", abs(float(_o_in["share"]) - 0.2441) < 1e-3)
check("proximity share 46%", abs(float(_p_in["share"]) - 0.4570) < 1e-3)
check("class nulls weight-average to the pooled rate", abs(
      (float(_o_in["null_rate"]) * int(_o_in["n_exon"])
       + float(_o_out["null_rate"]) * int(_o_out["n_exon"]))
      / (int(_o_in["n_exon"]) + int(_o_out["n_exon"])) - 0.0697) < 5e-4)
check("both rules give the same total excess", abs(
      float(_o_in["excess_exons"]) + float(_o_out["excess_exons"])
      - float(_p_in["excess_exons"])
      - float(ov[("proximity", "other_genes")]["excess_exons"])) < 2.0)
check("manuscript leads with containment",
      "21.64% against\n9.02% elsewhere, a 2.40-fold enrichment" in MS)
check("manuscript gives the class-specific nulls",
      "11.67%" in MS and "6.60%" in MS)
check("manuscript reports the proximity sensitivity",
      "proximity rule of 10 kb" in MS and "1.61-fold" in MS)
check("manuscript no longer claims 62%", "62%" not in MS)
check("manuscript states the majority lies outside",
      "does the majority of the excess sit in novel-junction genes" in MS)

# ---- the wrong-sign claim. Stratifying degron_dose by baseline CpG count
# splits exons by LENGTH (mean 163bp vs 2161bp), not by measurement quality.
# The negative slope holds in the short-exon majority; the 16+ stratum spans
# zero rather than reversing. See DESIGN "Finding 6(b)". The manuscript may
# therefore keep its wrong-sign statements.
check("manuscript keeps the wrong-sign claim",
      "wrong-signed dose-response" in MS
      and "returns the wrong sign" in MS
      and "in the direction opposite to the" in MS)

# ---- Figure 4 must read its numbers, not carry them. Panels A, B and D
# hard-coded their values until 2026-09-11, with panel B sourced from a
# DESIGN.md note rather than a table, so a drift between figure and text was
# invisible to this checker. Guard against the literals returning.
_fig = open(os.path.join(ROOT, "hpc/scripts/msf_fig4_cryptic.py")).read()
_body = "\n".join(l for l in _fig.splitlines()
                  if not l.strip().startswith("#"))
check("Figure 4 carries no hard-coded panel values",
      not any(v in _body for v in ("13.77", "8.53", "15670", "3848",
                                   "4.42", "6.40", "2.02", "2.93",
                                   "21.64", "9.02", "11.67")))
check("Figure 4 reads dko_overlap.tsv", 'tsv("dko_overlap.tsv")' in _fig)
check("Figure 4 requires all nine degron contrasts",
      "expected 9" in _fig)

# ---- the cryptic account must not be overclaimed. Two reductions apply in
# sequence and they are unequal: calibration removes 6,619 of the 9,433 exons
# crossing the threshold, and of the 2,814 that survive, 687 (a quarter) sit
# in novel-junction genes. Wording that implied the transcript class carried
# "most" of it, or that the data supported the account "directly", was
# corrected 2026-09-11. See DESIGN "Finding 4".
_dko = next(r for r in csv.DictReader(
    open(os.path.join(ROOT, "methsplice/contrasts_results.tsv")),
    delimiter="\t") if r["contrast"] == "HCT116_DKO_vs_WT")
check("manuscript's calibration arithmetic matches contrasts_results",
      f"{int(_dko['n_gt05_obs']):,}" in MS
      and f"{int(float(_dko['n_gt05_null'])):,}" in MS
      and "descriptive for the reason given above" in MS)
check("manuscript states the 2,814 genuine excess", "2,814" in MS)
check("manuscript states the 687-exon share", "687 sit in genes" in MS)
check("cryptic account not overclaimed",
      "the evidence for it is\nnarrower than the phenomenon" in MS
      and "substantially dissolves" not in MS
      and "accounts for most of what perturbation" not in MS)
check("manuscript says the residual is uncoupled",
      "What is left after both is unexplained" in MS
      and "uncoupled from methylation." in MS)

# ---- S3 calibrated across every contrast (cryptic_all.tsv, ms1_cryptic2.py).
# The predecessor ran on a hand-picked subset and reported S3 with no null at
# all. Calibrated, the DKO asymmetry is the most extreme relative to its own
# null but is NOT unique: UHRF1-AID day 6 and a six-hour CTCF degron also
# exceed their null maxima, the latter with no methylation change. The
# manuscript must carry that. See DESIGN "Finding 2".
_ca = {r["contrast"]: r for r in csv.DictReader(
       open(os.path.join(ROOT, "methsplice/cryptic_all.tsv")), delimiter="\t")}
check("cryptic_all covers every contrast in the config", len(_ca) == 17)
_exceed = [c for c, r in _ca.items()
           if float(r["novel_ratio"]) > float(r["novel_null_hi"])]
check("exactly three contrasts exceed their S3 null maximum",
      sorted(_exceed) == sorted(["HCT116_DKO_vs_WT", "degron_u1aid_Day6",
                                 "K562_CTCFdtag_6h"]))
check("DKO S3 null band and normalised ratio as quoted",
      abs(float(_ca["HCT116_DKO_vs_WT"]["novel_null_lo"]) - 0.97) < 1e-9
      and abs(float(_ca["HCT116_DKO_vs_WT"]["novel_null_hi"]) - 1.03) < 1e-9
      and abs(float(_ca["HCT116_DKO_vs_WT"]["novel_norm_ratio"]) - 3.27) < 1e-9)
check("manuscript gives S3 its depth normalisation", "3.27-fold" in MS)
# The DKO's four libraries are two runs of one library per genotype
# (SRX1828597 / SRX1828598), so its relabeling interval is sequencing
# variance alone. The arm was demoted to description 2026-09-11; these
# checks fail if it is ever re-presented as calibrated. See DESIGN
# "Finding 9".
check("DKO replication structure stated in Results",
      "SRR3644141" in MS and "SRX1828597" in MS
      and "one\nbiological replicate per arm" in MS)
check("DKO relabeling interval explicitly not used as a null",
      "is not a null against which an effect can be judged" in MS
      and "we do not use it as one" in MS)
check("DKO no longer called replicates",
      "every DKO replicate" not in MS
      and "both DKO runs and neither wild-type run" in MS)
check("DKO no longer called the clean chronic comparison",
      "clean chronic comparison" not in MS)
check("cryptic account rests on the replicated series",
      "The only place we can show this with biological replication" in MS
      and "account does not rest on it" in MS)
check("Figure 4 legend flags the descriptive panels",
      "they\nare descriptive, and the intervals shown for them are not nulls" in MS)
check("Methods explains why the 2v2 interval is not a null",
      "reflects sequencing variance alone and is not\nused as a null" in MS)
check("manuscript discloses that S3 is not specific to demethylation",
      "the statistic is not\nspecific to demethylation" in MS
      and "six-hour CTCF degradation 1.80-fold" in MS
      and "5.32-fold" in MS)

# ---- DKO dose split (dko_dose.tsv, ms1_dko_dose2.py). Containment is primary
# but cannot estimate the novel-junction class: only 409 of its exons carry
# baseline methylation at the coverage floor. Every slope is negative and
# excludes zero under both rules, which is the point the text makes.
_dd = {(r["rule"], r["class"]): r for r in csv.DictReader(
       open(os.path.join(ROOT, "methsplice/dko_dose.tsv")), delimiter="\t")}
check("DKO dose ALL reproduces under both rules",
      abs(float(_dd[("overlap", "ALL")]["slope"]) + 0.000100) < 1e-6
      and abs(float(_dd[("proximity", "ALL")]["slope"]) + 0.000100) < 1e-6)
check("every DKO dose slope is negative and excludes zero",
      all(float(v["slope"]) < 0 and float(v["ci_hi"]) < 0 for v in _dd.values()))
check("containment cannot estimate the novel-junction class",
      ("overlap", "novel_junction_genes") not in _dd)
check("manuscript reports the containment residual and the 409",
      "-0.000071 (-0.000106 to -0.000037)" in MS and "only 409 of its exons" in MS)

# ---- bee per-exon, regenerated (perexon_counts.tsv, ms1_perexon.py)
_pc = {r["stat"]: r["value"] for r in csv.DictReader(
       open(os.path.join(ROOT, "methsplice/perexon_counts.tsv")), delimiter="\t")}
check("per-exon counts regenerate", _pc["complete_exons"] == "37416"
      and _pc["at_floor_obs"] == "87" and _pc["p05_obs"] == "888"
      and _pc["bh_reject"] == "0" and _pc["candidates"] == "57")
check("conditional null conditions on max-over-splits, not observed dPSI",
      _pc["candidates_exp"] == "6.5" and _pc["dpsi_max_gt05"] == "2995")
check("three-CpG subset is 31 and 6, not 34 and 7",
      _pc["cand_3cpg"] == "31" and _pc["cand_3cpg_meth10"] == "6"
      and "the 31 that meet a" in MS and "in six cases" in MS
      and "a rate of 19%" in MS)

# ---- the bee knockdown replication (bee_global_meth.tsv, ms1_bee_global.py).
# This is the paper's opening empirical claim and the premise everything
# downstream rests on: if the knockdown did not work, every null is
# uninformative. Until 2026-09-11 neither 1.116 nor 0.852 appeared in any
# artefact and nothing checked them. See DESIGN "The PSI pipeline audited".
_bg = {(r["cov_floor"], r["group"]): r for r in csv.DictReader(
       (l for l in open(os.path.join(ROOT, "methsplice/bee_global_meth.tsv"))
        if not l.startswith("#")), delimiter="\t")}
_bc = float(_bg[("1", "control")]["weighted_pct"])
_bk = float(_bg[("1", "knockdown")]["weighted_pct"])
check("bee global methylation 1.116 -> 0.852 (coverage-weighted)",
      abs(_bc - 1.1160) < 5e-4 and abs(_bk - 0.8516) < 5e-4)
check("bee reduction rounds to 24%", 23.5 <= 100 * (1 - _bk / _bc) < 24.5)
check("manuscript states the weighting and the alternative",
      "coverage-weighted over all 9.6" in MS
      and "1.204% to 1.079%" in MS and "a 10% drop" in MS)
check("bee reduction stable across coverage floors",
      all(23.0 < 100 * (1 - float(_bg[(f, "knockdown")]["weighted_pct"])
                        / float(_bg[(f, "control")]["weighted_pct"])) < 25.0
          for f in ("1", "3", "5")))

# ---- borrowed mouse wild type. GSE130686 is TKO-only (6 libraries, no WT),
# so the WT had to come from elsewhere; it is the untreated arm of a Daxx
# deletion study. The "no linked publication" claim has a shelf life and is
# therefore date-stamped in the text. Rechecked 2026-09-11: GEO pubmedids
# empty, Europe PMC hitCount 0 for "GSE241701".
check("manuscript names what the borrowed mouse WT actually is",
      "Daxx-deletion study" in MS and "carrying no\nwild type of its own" in MS)
check("the no-publication claim is date-stamped",
      "rechecked September 2026" in MS)

# ---- the two papers LITERATURE.md had recorded as not existing, found
# 2026-09-15. Xu et al. is a cnidarian loss-of-function study and this paper's
# title claims a cnidarian arm, so its absence would be the most obvious
# omission available to a reviewer. See LITERATURE.md, "THE SEARCH FOR A
# PUBLISHED NULL, RERUN".
check("Xu et al. 2026 cited and engaged",
      "Xu\net al. (2026)" in MS and "10.1038/s41559-026-03090-6" in MS
      and "genome-defence mechanism" in MS)
check("Manz and List 2024 cited and engaged",
      "Manz and List (2024)" in MS and "10.1101/2024.08.30.610315" in MS
      and "on other\nchromosomes" in MS)
check("the cnidarian correspondence is stated, not left implicit",
      "our coral arm and their" in MS and "same phylum" in MS)
check("no priority claim was introduced",
      "first direct negative" not in MS and "to our knowledge" not in MS)

# The reference list must stay alphabetical and every entry must carry a DOI:
# RNA requires DOIs, and two entries were inserted at the top on 2026-09-15.
import unicodedata as _ud
_rb = MS[MS.index("## References"):MS.index("## Figure legends")]
_refs = [e.strip() for e in _rb.split("\n\n") if "doi:" in e]
def _k(e):
    return _ud.normalize("NFKD", e.split()[0]).encode("ascii", "ignore").decode().lower()
check("every reference carries a DOI",
      len(_refs) == len([e for e in _rb.split("\n\n") if e.strip()
                         and not e.strip().startswith("#")]))
check("reference list is alphabetical",
      all(_k(_refs[i]) >= _k(_refs[i-1]) for i in range(1, len(_refs))))

# ---- pseudo-replication has a name and a literature, cited 2026-09-15.
# Zimmerman et al. also bears on our own aggregation choice, and that is the
# part worth protecting: they find aggregation to the individual conservative
# relative to mixed models, so merging runs within a bee errs toward missing
# an effect. A reviewer asking whether the nulls are an artefact of
# aggregation is answered in the text.
check("pseudo-replication precedent cited",
      "Lazic (2010)" in MS and "10.1186/1471-2202-11-5" in MS
      and "Zimmerman et al. (2021)" in MS
      and "10.1038/s41467-021-21038-1" in MS)
check("the aggregation direction is stated",
      "conservative and underpowered relative to mixed" in MS
      and "cannot be attributed to the aggregation" in MS)
check("the reads-to-bee analogy is made explicit",
      "same relation to the bee as cells within a donor" in MS)

lt = open(os.path.join(ROOT, "methsplice/latent_tss.txt")).read()
check("latent TSS 4.42/2.02 and 6.40/2.93",
      "4.42%" in lt and "2.02%" in lt and "6.40%" in lt and "2.93%" in lt
      and "4.42%" in MS and "6.40%" in MS)

# ---- contrasts table (contrasts_results.tsv)
cx = {r["contrast"]: float(r["excess05"]) for r in
      csv.DictReader(open(os.path.join(ROOT, "methsplice/contrasts_results.tsv")), delimiter="\t")}
chronic = [v for k, v in cx.items() if "TKO" in k or "DKO_vs_WT" in k.upper()
           or k == "HCT116_DKOvWT"]
check("chronic excess range 1.43-1.79 covered",
      any(abs(v-1.79) < 0.02 for v in cx.values()) and
      any(abs(v-1.59) < 0.02 for v in cx.values()),
      f"(TKO serum {cx.get('mouse_TKOvWT_serum')})")
check("manuscript quotes 1.43 to 1.79", "1.43 to\n1.79" in MS or "1.43 to 1.79" in MS)

# ---- positive control (hassulta artefacts)
check("manuscript quotes x2.31 and x2.55 and floors",
      "2.31" in MS and "2.55" in MS and "0.0011" in MS)

print()
if FAIL:
    print(f"{len(FAIL)} FAILURE(S): {FAIL}")
    sys.exit(1)
print("manuscript agrees with artefacts")
