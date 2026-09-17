# Analysis code for insect DNA methylation and comparative transcriptomics

Kevin J. Dudley, Queensland University of Technology

This repository holds the analysis code, derived tables and software provenance
for three manuscripts. It is a curated extract, not a working tree: it contains
what the papers cite and nothing else.

All three studies are reanalyses of public archive data. **No new sequence data
were generated**, so the primary data are obtained from the accessions given in
each paper and reproduced below.

---

## What is here, by paper

### `conservation/` — Shared gene properties manufacture cross-species conservation signal

The claim is that intersecting thresholded gene sets across orthologs produces
extreme false positives, because orthologs share static properties that make
the same genes cross thresholds in both species without conserved regulation.

| script | what it does |
|---|---|
| `harm_overlap_null.py` | threshold-set overlap, with the within-species label permutation |
| `harm_spearman.py` | the threshold-free variant: rank correlation of per-gene effect magnitude |
| `harm_topk.py` | top-K overlap at K fixed to the observed set size |
| `harm_contrast.py` | per-exon PSI and the group contrasts the above consume |
| `sim_conservation_artefact.py` | simulation with zero conserved regulation |
| `sim_conservation_grid.py` | the 36-cell design grid across group size, threshold and heterogeneity |

**Seeding.** `harm_overlap_null.py` and `harm_topk.py` use `random.Random(23)`,
`harm_spearman.py` uses `random.Random(31)`, and the simulations derive their
seed from the replicate index and the heterogeneity parameter
(`np.random.default_rng(1000 * rep + int(sigma * 10))`). Every reported p-value
and false-positive rate is therefore exactly reproducible.

**Data.** *H. assulta* RNA-seq PRJEB70911, 18 runs ERR12364891 to ERR12364908;
*H. armigera* RNA-seq PRJNA1061802, 12 runs SRR27472011 to SRR27472022.
Assemblies GCA_050947285.1 and GCF_030705265.1. Run manifests are in
`conservation/tables/`.

Derived tables in `conservation/tables/` are sufficient to reproduce every
figure and statistic **without re-running alignment**, which is the expensive
stage.

### `comparative/` — Caddisfly gene body methylation

Direct methylome measurement from deposited PacBio HiFi base-modification
calls, for three Trichoptera spanning both suborders, against an eleven-genome
native dipteran floor panel.

| script | what it does |
|---|---|
| `caddis_methyl.pbs` | fetch, align with pbmm2, call CpG scores, summarise |
| `meth_summary.pbs` | derives the methylome summary from an existing bed |
| `caddis_floor.pbs`, `floor_sweep.py` | the banded false-positive floor |
| `miniprot_annot.pbs` | protein-homology annotation for genomes without RNA |
| `meth_partition.py` / `.pbs` | assigns CpGs to exonic, intronic, intergenic |
| `gene_oe_meth.py` / `.pbs` | CpG observed/expected against measured methylation, per gene |
| `philo_oe_checks.py` / `.pbs` | the exon-only and repeat-masked o/e controls |
| `dnmt_prescreen.pbs`, `dnmt_family_call.py` | methyltransferase family calls |
| `figdata.pbs`, `make_figures.py` | figure data and the figures |

`comparative/tables/floor/` holds the per-genome floor tables for the whole
77-genome panel; `panel_registry.tsv` carries the clade and library class for
each, which determines panel membership.

### `methsplice/` — No evidence that DNA methylation instructs exon inclusion

Audit and figure code for the perturbation, dose-response and dynamics
reanalysis. `ms_consistency.py` is the manuscript consistency checker: it tests
that the paper matches its own artefacts, and it carries guards that fail if
withdrawn claims are reintroduced. It proves transcription, not correctness.

---

## `provenance/`

`tool_versions.txt` records the software versions the analyses ran on,
interrogated from the binaries directly rather than inferred from package
metadata:

    samtools 1.24      bcftools 1.24     hisat2 2.2.3      stringtie 3.0.3
    minimap2 2.31-r1302    miniprot 0.18-r281
    pbmm2 26.2.99      pb-CpG-tools 3.0.0    Python 3.10.13

The `*.explicit.txt` files are conda explicit-spec exports of the five
environments, suitable for `conda create --file`.

---

## Running this elsewhere

The PBS scripts were written for a specific cluster and hardcode paths under
`/work/cyberomics/lep_meth`, along with PBS directives for that scheduler. They
are included as the exact provenance of what was run, not as a portable
pipeline. The Python analysis scripts take their inputs as arguments and will
run anywhere.

## Licence and citation

Code is MIT licensed. If you use it, please cite the corresponding paper; each
manuscript names the scripts it relies on in its Software availability section.
