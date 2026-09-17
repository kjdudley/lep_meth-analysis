<!-- Extracted from the project parameter ledger for public release.
     Sections 9 and 11, which cover the allele-specific phased arm, are omitted:
     that work belongs to a separate multi-author manuscript and is not this
     project's to publish. Everything the methylation-and-splicing paper relies
     on is retained. -->

# Methods

Every parameter actually used, with its justification — or with an explicit note
that it is conventional and unjustified. Written 2026-08-12 to be
publication-tight; a reviewer should be able to reproduce any number from this
file plus `PROVENANCE.tsv` (1,467 job records) and `SCRIPT_HISTORY.md`.

**Thresholds are separated into three kinds throughout: derived from evidence in
this project, taken from a tool default, or chosen by convention.** Only the
first kind should be defended as a finding.

---

## 1. Data sources and access dates

All three sources are live services. Dates below are the mtimes of the resulting
tables, which bound the access date from above; **exact query timestamps were not
recorded and cannot be recovered.** State the table dates in any manuscript, not
"accessed 2026".

| source | query | table | date |
|---|---|---|---|
| NCBI Datasets v2alpha | all GenBank assemblies, Noctuidae | `assemblies_screened.tsv` (351 ToLIDs) | 2026-08-07 |
| ENA portal | `tax_tree(7100) AND instrument_platform="PACBIO_SMRT"` | `runs_screened.tsv` | 2026-08-07 |
| ENA portal | `tax_tree(7100) AND library_source="TRANSCRIPTOMIC"` | `rna_raw.tsv` (5,272 runs) | 2026-08-07 |
| ENA + BioSamples | per-run sample, per-sample attributes | `rna_match.tsv`, `hpc/samples.tsv` (73 spp.) | 2026-08-07 |
| Ensembl Rapid Release | per-species geneset URLs | `ensembl_manifest.tsv` (121 rows) | 2026-08-10 |

Annotation **provider and version vary by species** (`braker` vs `ensembl`;
2022_03 through 2024_04) and are recorded per species in `ensembl_manifest.tsv`.
This is a real heterogeneity, not a detail: it was implicated in the CpG o/e
work and must be carried as a covariate in any cross-species comparison.

## 2. Software

Captured 2026-08-12 (`SOFTWARE_VERSIONS.txt`):

| tool | version | used for |
|---|---|---|
| `jasmine` | 26.1.3 | re-calling only |
| `pbmm2` | 26.2.99 | all alignment |
| `aligned_bam_to_cpg_scores` (pb-CpG-tools) | 3.0.0 | all CpG aggregation |
| `hisat2` | 2.2.3 | RNA alignment |
| `stringtie` | 3.0.3 | *Xylocopa* annotation |
| `minimap2` | 2.31-r1302 | — |
| `ccs` | 6.4.0 | *Xylocopa* only (`--hifi-kinetics`) |
| Python | 3.11.15 | all analysis |
| `samtools` | **1.24 in `env`, 1.23.1 in `env2`** | mixed |

**The samtools discrepancy is real** — analyses run from `env` and `env2` used
different versions. Nothing here is known to depend on it, but it must be
disclosed rather than reported as a single version.

Deposited 5mC calls were produced by the depositors, not by us: primrose 1.2.0
and 1.4.0, jasmine 2.2.0 and 2.3.0. Caller version per dataset is recorded.

## 3. Methylation aggregation

`aligned_bam_to_cpg_scores --modsites-mode reference --pileup-mode model
--threads 16`, i.e. **tool defaults except `--modsites-mode reference`**.

* `--min-coverage 4` — **tool default**, not chosen. Every bed inherits it, and
  it is visible in each file's `##min-coverage=4` header. It sets the floor of
  the coverage-stratification analysis.
* Alignment: `pbmm2 --preset HIFI --sort -j 16`. `-j 8` for *P. flavicincta*
  only, after `pbmm2` aborted twice at `-j 16` with a thread-teardown crash on
  its 120 GB input. **Recorded because it is a per-species deviation.**

## 4. Confident methylation, and the floor

* **≥50% methylation at a site = "confident".** *Derived*: rule 2 of this
  project — the genome-wide mean is dominated by the ~3% noise floor, whereas
  the ≥50% fraction separates species. Not a tool default.
* **Coverage stratification in 10× bins to 200×+.** *Derived*: the confident-call
  rate is strongly coverage-dependent (job `24926958`) — 5.4% at 0–9× falling to
  0.000% at 100×+ in unmethylated lineages, while a methylated moth is flat at
  0.33–0.48% from 30× to 99×. **Any methylation figure must state its coverage.**
* **The mitochondrial control is reported but is NOT the floor.** Mitochondrial
  coverage is ~1,900× where the false-positive rate is 0.000% regardless, so it
  systematically understates the nuclear rate. Retained as a per-dataset sanity
  check only. Mitochondrion identified by `grep -i -m1` of the FASTA headers for
  `mitochondrion|mitochondrial|organelle: mitochondrion` — **a heuristic, not an
  annotation lookup**; first match wins.
* Genome mean is **unweighted** across sites (not coverage-weighted).

## 5. Sequence-namespace reconciliation

CpG beds use INSDC accessions (`OX424490.1`); Ensembl Rapid Release renames
placed chromosomes to bare `1`..`30`/`W`/`Z`. Reconciled by matching **sequence
length** between the `.fa.fai` and the GFF's `##sequence-region` pragmas,
accepting **only lengths that are unique on both sides**.

Guard: analyses exit non-zero if **>5% of CpGs** fall outside the annotation.
*Conventional* threshold, but its purpose is derived — job `24906658` exited 0
having binned all 28.8 M chromosomal CpGs as intergenic through exactly this
mismatch.

## 6. RNA and PSI

* HISAT2 index built **without `--ss/--exon`**. *Derived*: supplying annotated
  splice sites biases junction discovery toward the annotation, and the question
  is whether an annotated exon is skipped.
* `hisat2 -p 14 --dta --no-unal`; `samtools sort -@ 4 -m 2G`.
* FASTQ verified against **ENA's published md5**, not byte count.
* Junctions from CIGAR `N` operations. Filters: `-F 3844` (drop unmapped,
  secondary, supplementary, QC-fail, duplicate), **MAPQ ≥ 10** (*conventional*),
  **anchor ≥ 8 aligned bases on both sides** (*conventional*; short anchors are
  the standard source of spurious junctions).
* PSI = inclusion / (inclusion + skip), inclusion averaged over the two flanks.
* **Classes: constitutive PSI ≥ 0.95; alternative 0.05–0.95; rarely-included
  < 0.05.** *Conventional.* The 0.95 boundary is not derived and the result's
  sensitivity to it has **not** been tested — a known gap.
* **Minimum 10 inclusion+skip reads** to classify an exon. *Conventional.*
* **Internal exons only** — internal in ≥1 transcript. *Derived*: first and last
  exons have no upstream/downstream junction, so their PSI is undefined by this
  method and would be spuriously low.
* CpGs require **coverage ≥ 5** (`--min-cov 5`), on top of the bed's own ≥4.

## 7. Statistics

* **Gene clustering throughout.** Exons within a gene are not independent. Genes
  are assigned **positionally**, never by the PSI table's gene column, which
  falls back to a transcript ID when the GFF carries no `gene_id` on exon lines.
* **Cluster bootstrap**, 2,000 resamples of whole genes, percentile CI.
* **≥20 CpGs per class** for the within-gene test. *Derived, and this one is a
  result*: the sensitivity curve (`skew_sensitivity.py`) shows the mean
  within-gene difference flipping sign as the threshold rises, stabilising
  negative at ≥20 in all three species. Below it, per-class rates from a handful
  of CpGs create a spurious positive tail.
* **Sign tests exclude ties.** *Derived*: `pct_ge50` is heavily zero-inflated —
  479 of 699 *A. litura* genes tie at exactly zero — and counting ties as
  non-negative produced p ≈ 1e-48 against an effect that was present.
* **Paired analyses use the continuous statistic (`mean_meth`), not
  `pct_ge50`.** *Derived*, same reason; the zero-inflated version left only 6–9
  untied pairs of 181–266.
* No p-values are reported for pooled comparisons — effect sizes with
  cluster-bootstrap intervals instead.
* **No phylogenetic correction is applied.** Species are not independent and
  PLAN warns the naive pooled regression is anticonservative. **Outstanding gap**
  — required before any cross-species claim.

## 8. Sequence-composition control

* GC fraction and CpG per kb computed from exon sequence, **N-aware**, requiring
  ≥20 informative bases; sequence retrieved in a single `samtools faidx -r` call.
* Stratified within **4×4 GC × CpG-density quantile cells**; cells require ≥30
  constitutive and ≥15 alternative exons. *Conventional.*
* Within-gene matching pairs each alternative exon with the constitutive exon in
  the same gene minimising `(ΔGC/0.05)² + (ΔCpG_kb/10)²`. **The scaling constants
  are arbitrary** and the sensitivity of the result to them has not been tested —
  a known gap, and material, because this control returned a negative.

## 10. Known methodological gaps — status

Listed so they are disclosed rather than discovered. Four were addressed
2026-08-12; three more (9–11) were **added** the same day, two of them found by
questioning analyses that had already been written up. The list growing is the
list working.

| # | gap | status |
|---|---|---|
| 1 | No phylogenetic correction | **Partially closed.** `genus_cluster.py` |
| 2 | PSI class boundaries untested | **CLOSED** — the null is boundary-independent |
| 3 | Matching constants arbitrary | **CLOSED, with a caveat that matters** |
| 6 | Mitochondrion heuristic | **CLOSED** — 23/24 verified, 1 flagged |
| 4 | Access timestamps unrecorded | **Addressed** — access dates are recorded in section 1; re-accessed and date before submission |
| 5 | samtools differs between envs | Disclose; no dependence known |
| 7 | NUMTs not excluded | Disclose — no longer threatening |
| 8 | No cross-platform anchor | **Closing**, job `24928274` |
| 9 | Coverage curve confounds depth with site difficulty | **OPEN, under test** — jobs `24932728`–`29` |
| 10 | Reference divergence affects the floor, unquantified | **OPEN** — measurable in *C. elegans* |
| 11 | No positive control; false-negative rate unverified | Disclose — estimates are lower bounds |

### Gap 9 — the calibration curve may be measuring the wrong thing

Raised 2026-08-12. `FPR(coverage)` in paper 1 §2.3 is built from **natural**
coverage variation: sites within one genome binned by their own read depth. Depth
within a genome is **not assigned at random** — shallow sites are
disproportionately repetitive, low-complexity or hard to map, which is to say
disproportionately the sites most likely to be miscalled for reasons unrelated to
depth. The curve therefore confounds *shallow* with *difficult*. *A. cognata*
reads 5.4% at 0–9× against 0.05% at 40–49×, and that 100-fold difference is
currently attributed entirely to depth.

`subsample_floor.pbs` separates them by holding the site set fixed — sites
ordinary at native depth, 40–59× — and subsampling the alignment, so difficulty
is constant and only depth varies. If the subsampled and natural curves agree,
§2.3 stands; if natural is much higher, the curve overstates the depth effect at
the low end and the recommended band may move. **Do not publish the curve before
this reports.**

### Gap 10 — reference divergence

Every reuse study maps reads to a reference of imperfect match, and the effect on
the 5mC floor has never been quantified. It becomes measurable in *C. elegans*,
where ToL reads from wild isolate 73214 are mapped to WBcel235 (N2): stratifying
the floor by local divergence turns the confound into a calibration of reference
mismatch itself. Until then, note that the panel's floors are measured against
assemblies of varying relatedness to the sequenced specimen.

### Gap 11 — no positive control

There is **no fully-methylated positive control anywhere for the HiFi calls**.
The only positive control in the project is the SssI in the *H. armigera*
bisulfite dataset, which failed its own assay and belongs to a different
technology. The vendor's 11% false-negative rate is cited but never
independently verified. **Every methylation estimate in paper 1 is therefore a
lower bound of unknown tightness**, and the cross-platform anchor (§2.9) is the
only external check that bears on it.

### Gap 1 — phylogenetic non-independence (partially closed)

**CORRECTED 2026-08-12. The BUSCO run does not produce a tree, and was never
going to.** This was stated repeatedly as "a PGLS awaits the orthologue tree from
the BUSCO run". Checked directly: `Snakefile.busco` has exactly two rules,
`busco_oe` and `summarise_pergene`. Its output `regions.tsv` carries *coordinates
and per-gene statistics* — `label, region, seq, start, end, n_acgt, gc, cpg_oe` —
not sequences. And **no alignment or phylogenetics software is installed in any
environment**: mafft, muscle, iqtree, raxml, FastTree, astral and trimal are all
absent from `env`, `env2`, `env_busco` and `env_phase`.

So the BUSCO run supplies a *substrate* — orthologous gene coordinates keyed by
BUSCO id across species — from which a tree could be built by extracting
sequences, aligning, concatenating and inferring. **None of those steps exist.**
Gap 1 is therefore **not on track**, and waiting for BUSCO to finish will not
close it.

**Two routes, and the cheap one is probably also the better one:**

1. **Use a published phylogeny.** For 24 measured species, insect relationships
   are well established and a dated tree (TimeTree, or a published Lepidoptera
   phylogeny) gives the branch lengths PGLS needs. Standard practice, hours not
   weeks, and not reviewer-contentious.
2. **Build one from the BUSCO orthologues.** Defensible and self-contained, but
   requires installing a toolchain, writing an extract/align/infer pipeline, and
   running it — against a benefit that is mostly aesthetic.

**Also ask whether paper 1 needs PGLS at all.** Its central claims are
*within-genome* (coverage curves, subsampling); species pooling enters only in
the floor summary, where genus-clustered bootstrap is already in hand and widens
the CI by 1.12×. PGLS matters most for **paper 2's** comparative claims. Do not
let gap 1 block paper 1 by default.

A full PGLS needs a tree, and the orthologue substrate (`24927989`) is not
finished. Genus is available now, and the panel deliberately contains
congeners. `genus_cluster.py` resamples **genera** rather than species — the
same clustering logic already applied to exons within genes, one level up.

Applied to confident methylation across the 24 measured species:

* 24 species in **20 genera**, 3 genera with >1 species
* mean 0.3704%; naive species bootstrap 95% CI **[0.3125, 0.4388]**
* genus-clustered 95% CI **[0.3075, 0.4492]** — **1.12× wider**
* mean |deviation| **within** genus **0.0295** vs **between** genera **0.1106**

**Congeners are ~3.7× more similar to each other than genera are to the panel
mean.** Species are therefore demonstrably non-independent, and the naive
interval is demonstrably too narrow — at this n the inflation is 12%.

**Genus must be taken from the species binomial, not the ToLID prefix.** The
first version used the ToLID, which encodes only three letters of the genus:
*Polymixis lichenea* (ilPolLich1) and *Polia nebulosa* (ilPolNebu2) both give
`ilPol` and were wrongly treated as congeners. Correcting it removed one of four
apparent multi-species genera — and **strengthened** the result, since a false
pair had been inflating within-genus deviation (0.0352 → 0.0295). Report genus-clustered intervals for any cross-species statistic. The
honest manuscript sentence: *genus-clustered intervals are reported; a
phylogenetically explicit model awaits the orthologue tree.*

### Gap 2 — PSI class boundary (closed)

Swept at 0.90 / 0.95 / 0.98 / 0.99 (job `24928146`). The **pooled** ratio does
move — *A. litura* 0.73 → 0.44 as the threshold rises — so a pooled figure
should never be quoted without its boundary. But the **within-gene** statistic
sits at **46.9–58.2% of genes with alt < con across every species and every
boundary**, i.e. at or near chance throughout. **The within-gene null is not a
threshold artefact.**

### Gap 3 — composition-matching constants (closed, with a caveat that matters)

Swept three scale settings × three callipers. The negative is stable: 44.9–48.7%
of pairs negative wherever enough pairs exist.

**But the calliper exposed a design limitation.** Requiring matches to be within
a stated tolerance collapses the sample: *A. litura* falls from 181 pairs to
**7** at the tightest setting, *N. janthina* from 266 to **10**. Only the
loosest combinations retain a testable number.

**Within a gene, constitutive and alternative exons are rarely close in
composition.** The original test appeared to work only because it accepted the
nearest constitutive exon however poor the match. So the composition-matched
design may be **unable to answer the question at all** — the matched pairs
barely exist. That is a limitation of the design, not evidence either way, and
it is a further argument for phasing, where matching is exact by construction
because it is the same exon.

### Gap 6 — mitochondrion identification (closed)

All 24 species checked against expected mitogenome size (job `24928147`).
**23 OK.** Every Lepidopteran falls in 15.3–15.7 kb; the Dipteran and Coleopteran
at 18.4 and 18.6 kb.

**One flagged: *Aethecerus discolor*, 30.0 kb — SUSPECT.** Ichneumonid wasps are
known for expanded mitogenomes so this may be genuine, but it is twice the
Lepidopteran size and must be confirmed manually before that species' mito
control is quoted. It is also one of the three non-Lepidopterans whose apparent
methylation is now known to be floor.

### Gap 8 — cross-platform anchor (closing)

An external anchor **does** exist, separate from the failed in-house WGBS:
**Jones et al. 2018, *G3*** — a single-nucleotide-resolution WGBS methylome for
*Helicoverpa armigera*, reporting **~0.9% of CpG sites methylated**, data at
ArrayExpress **E-MTAB-4779**.

*H. armigera* has deposited native HiFi with MM tags (`ERR12102456`, Sequel IIe,
30.4 Gb) and was in PLAN's gate test 3 for exactly this purpose. Queued as
`ilHelArmi9` (job `24928274`) against GCA_963930815.1.

The comparison: does our confident-methylation estimate agree with 0.9%, and
does the genomic distribution agree? Our Lepidopteran range is 0.23–0.99%, so
0.9% is squarely inside it. **Caveat: different individuals and populations, so
per-site concordance is bounded by biological variation — genome-wide level and
distribution are the comparable quantities, not individual sites.**

### Gaps 4, 5, 7 — disclosed, not closed

**4. Access timestamps.** The queries were not logged and cannot be recovered.
Report table mtimes as an upper bound.

**5. samtools version.** `view`, `mpileup` and `faidx` behave identically
between 1.23.1 and 1.24 for the operations used. Re-running to unify would
consume compute to change nothing.

**7. NUMTs.** This has stopped being a threat. Nuclear mitochondrial insertions
mismapping to the mitochondrion would carry nuclear methylation and **raise**
the apparent mitochondrial signal. We observe 0.00%. Possible NUMT contamination
therefore makes that zero *more* robust, and reinforces rather than undermines
the conclusion that the mitochondrial control is too optimistic.

### Gap 8 — open

The only external validation the programme could have. Needs published WGBS or
EM-seq for a species we have measured; the *H. armigera* WGBS cannot serve, as
it failed its own SssI positive control. A literature and data search, not a
compute task.
