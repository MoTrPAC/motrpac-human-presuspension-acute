# Landscape manuscript figures

Figure and supplementary-table code for the MoTrPAC Human Pre-Suspension **landscape**
manuscript, the cross-tissue integrative analysis. The other manuscripts live beside this
directory under `figures/blood/`, `figures/muscle/`, `figures/adipose/` and
`figures/splicing/`.

**One figure, one script.** Every manuscript figure is a single runnable R script, named for the
figure and living in the folder of the figure group it belongs to: `figure_2/FIG2.R`,
`figure_4/ED6.R`, `supplementary_figure_1/SF1.R`, and so on. Each script builds every panel of
that figure, one vector PDF per panel, and nothing else. Supplementary tables that are not a
panel's own numbers are built the same way by an `ST*.R` script in the same folder.

A main figure and its Extended Data supplement share a folder, because they are usually one
computation split across two figure numbers. The **figure index below gives the script path for
every figure**; you should not have to guess it.

- **94 panels across 17 figures.** Five further panels are made outside R and are not built
  here.
- **31 of the 36 numbered supplementary sub-tables.** The other five belong to other analysts.
- `config/panel_map.json` and `config/table_map.json` are authoritative for the numbering. A
  script reads its own output path, page size and RNG seed from the manifest, so renumbering a
  panel is an edit to one JSON file and nowhere else.

Outputs land at declared paths:

```
outputs/panels/<figure_dir>/<PANEL>_<slug>.pdf
outputs/tables/<TABLE>_<slug>.tsv
```

Note that `<figure_dir>` is the **manuscript** figure, not the folder the script lives in:
`figure_4/ED6.R` writes to `outputs/panels/extended_data_6/`. The mapping is in
`config/panel_map.json`.

Every script finds this directory by walking up from its own path, so no environment variable is
needed. `LANDSCAPE_ROOT` overrides that if you need it to.

## Figure index

| Figure | Script | Panels | Title |
|---|---|---|
| FIG1 | `figure_1/FIG1.R` | FIG1C | Study design and cohort overview |
| ED1 | `figure_1/ED1.R` | ED1A, ED1B, ED1C, ED1D | Extended Data 1 - cohort and data overview |
| FIG2 | `figure_2/FIG2.R` | FIG2A, FIG2B, FIG2Bii, FIG2C, FIG2D, FIG2E, FIG2F | Differential abundance landscape |
| ED2 | `figure_2/ED2.R` | ED2A, ED2B, ED2C, ED2D, ED2E, ED2F | Extended Data 2 - differential abundance supplement |
| FIG3 | `figure_3/FIG3.R` | FIG3A, FIG3B, FIG3C, FIG3D, FIG3EFG | Resistance vs endurance exercise |
| ED3 | `extended_data_3/ED3.R` | ED3A, ED3B, ED3C, ED3D, ED3E | Extended Data 3 - sex sensitivity and cellular deconvolution |
| FIG4 | `figure_4/FIG4.R` | FIG4A, FIG4B, FIG4C, FIG4D, FIG4E, FIG4F, FIG4G, FIG4H, FIG4I | Clustering and PLIER latent variables |
| ED4 | `figure_3/ED4.R` | ED4A, ED4B, ED4C, ED4D, ED4E, ED4F, ED4G | Extended Data 4 - RE vs EE supplement |
| ED5 | `figure_4/ED5.R` | ED5A, ED5B, ED5C, ED5D, ED5E | Extended Data 5 - fuzzy c-means supplement |
| ED6 | `figure_4/ED6.R` | ED6A, ED6Aii, ED6B, ED6C, ED6D, ED6E, ED6F, ED6G, ED6H, ED6I | Extended Data 6 - cross-tissue PLIER supplement |
| FIG5 | `figure_5/FIG5.R` | FIG5A, FIG5B, FIG5C, FIG5D, FIG5E, FIG5F, FIG5G, FIG5H | Transcription-factor motif enrichment and TF exercise responses |
| ED7 | `figure_5/ED7.R` | ED7A, ED7B | Extended Data 7 - transcription-factor supplement |
| FIG6 | `figure_6/FIG6.R` | FIG6C, FIG6D, FIG6E, FIG6F | Multi-omic regulatory network analysis (SC-ION) |
| ED8 | `figure_6/ED8.R` | ED8A, ED8B, ED8C, ED8D, ED8E, ED8F | Extended Data 8 - TFEB pathway regulation in response to EE and RE |
| FIG7 | `figure_7/FIG7.R` | FIG7A, FIG7B, FIG7C, FIG7D | Secretome |
| ED9 | `figure_7/ED9.R` | ED9A, ED9B, ED9C, ED9D | Extended Data 9 - secretome supplement |
| SF1 | `supplementary_figure_1/SF1.R` | SF1A, SF1B, SF1C, SF1D, SF1E, SF1F, SF1G | Supplementary Figure 1 - skeletal muscle ATAC-seq QC |

Ninety-four panels. Two panel ids are not manuscript letters: `FIG2Bii` and `ED6Aii` are second
graphics from a computation whose figure has no letter left, and `FIG3EFG` is one panel spanning
three letters, with one legend and one colour scale.

### Five panels are made outside R

They are recorded in each figure's `not_built` block in `config/panel_map.json` and are not
built here.

| Panel | What it is | Made with |
|---|---|---|
| FIG1A | Study design: screening, randomization and the three intervention arms | BioRender |
| FIG1B | Sampling timeline: the seven collection points across the two days | BioRender |
| FIG6A | SC-ION workflow schematic | Illustrator |
| FIG6B | Merged EE+RE network and the TFEB first-neighbour subnetwork | Cytoscape |
| FIG6G | TFEB signalling module schematic | BioRender |

Every other panel letter in the manuscript is built here, so a reader who finds no `FIG4A` on
disk is reading a build that failed rather than a panel this directory leaves to someone else.

### Single-feature plots

Sixteen panels are the same kind of plot: one feature's abundance over the acute timepoints,
exercise-group mean with a 95% confidence interval. FIG2D, ED2B, ED2D, ED3B, ED3C, FIG3D,
FIG4I, ED4D, ED4E, FIG6F, ED8B, ED8D, ED8F, FIG7B, FIG7C and ED9C are all defined once in
`single_feature_plots.R`, and their figure scripts call into it.

## Supplementary tables

The manuscript numbers its supplementary tables ST1 to ST7. Each number is one table made of
several sub-tables, and the sub-table is the unit of work: one row in `config/table_map.json`,
one TSV. The letter suffix is a stable id, not an ordinal recomputed when a neighbour moves.
Extended Data Table 1 is in the manuscript rather than the supplement, has no sub-tables, and is
carried in the same manifest under the id `EDT1`.

The manifest holds all 36 numbered entries, the 35 supplementary sub-tables and `EDT1`,
including the five this directory does not build. The numbering only makes sense whole: a reader
who finds `ST2d` on disk and no `ST2e` needs to see that `ST2e` exists and is someone else's
rather than conclude the build dropped it.

Sub-tables are written as one tab-delimited file each, not as assembled workbooks. Assembling
ST1 to ST7 into the shape a journal wants is a submission step. The submitted workbooks, and the
submitted figures, are kept under `assembled/`; see `assembled/README.md` for how they relate to
what the scripts emit.

**A sub-table whose numbers are a panel's numbers is written by that panel's figure script**, in
the same run, from the same in-memory object. Nineteen of the 31 built here are produced that
way. The remaining twelve come from the standalone `ST*.R` scripts, each in its figure's folder.

| Table | Sub-table | Built here | Written by |
|---|---|---|---|
| ST1 | `ST1a` Participants per randomization group | yes | `figure_1/ST1.R` |
| | `ST1b` Participants per group and sex | yes | `figure_1/ST1.R` |
| | `ST1c` Participants per tissue and omics assay | yes | `figure_1/FIG1.R` (FIG1C) |
| | `ST1d` Outliers | yes | `figure_1/ST1.R` |
| | `ST1e` Named features | yes | `figure_1/ED1.R` (ED1C) |
| ST2 | `ST2a` Covariates included | yes | `figure_2/ST2.R` |
| | `ST2b` Total number of features detected | yes | `figure_1/FIG1.R` (FIG1C) |
| | `ST2c` Number of DA features | yes | `figure_2/FIG2.R` (FIG2A) |
| | `ST2d` DA percent overlap | yes | `figure_2/ST2.R` |
| | `ST2e` Cross-tissue ORA | yes | `figure_2/FIG2.R` (FIG2B) |
| | `ST2f` Phospho ORA (PTMsigDB) | yes | `figure_2/ED2.R` (ED2E) |
| | `ST2g` Sex-specific DA estimate correlation | yes | `extended_data_3/ED3.R` (ED3A) |
| | `ST2h` Sex-specific differences enrichment | yes | `extended_data_3/ED3.R` (ED3A) |
| ST3 | `ST3a` EE muscle-only ORA features | yes | `figure_3/FIG3.R` (FIG3B) |
| | `ST3b` RE muscle-only ORA features | yes | `figure_3/FIG3.R` (FIG3B) |
| | `ST3c` All-muscle RE/EE-only ORA | yes | `figure_3/FIG3.R` (FIG3B) |
| | `ST3d` Triangle selected pathways | yes | `figure_3/ED4.R` (ED4B) |
| | `ST3e` DA features with opposite EE/RE direction | no | not in this repository |
| | `ST3f` Clinical x omics all significant associations | yes | `figure_3/ED4.R` (ED4G) |
| | `ST3g` Over-representation analysis of muscle transcripts associated with baseline VO2peak | yes | `figure_3/ST3.R` |
| | `ST3h` Over-representation analysis of adipose transcripts associated with baseline HOMA-IR | yes | `figure_3/ST3.R` |
| ST4 | `ST4a` ORA of fuzzy c-means clusters, filtered | yes | `figure_4/ST4.R` |
| | `ST4b` Cross-tissue RNA-seq PLIER | no | not in this repository |
| | `ST4c` Transcript affinity | no | not in this repository |
| | `ST4d` Cross-tissue metabolomics PLIER | no | not in this repository |
| | `ST4e` Metabolite affinity | no | not in this repository |
| ST5 | `ST5a` Transcription factor enrichment and expression | yes | `figure_5/ST5.R` |
| ST6 | `ST6a` Combined EE+RE edge table | yes | `figure_6/ST6.R` |
| | `ST6b` Combined EE+RE node table | yes | `figure_6/ST6.R` |
| | `ST6c` TFEB ChIP targets | yes | `figure_6/FIG6.R` (FIG6D) |
| | `ST6d` EE 3.5-4 h up | yes | `figure_6/FIG6.R` (FIG6D) |
| | `ST6e` EE 3.5-4 h down | yes | `figure_6/FIG6.R` (FIG6D) |
| | `ST6f` RE 3.5-4 h up | yes | `figure_6/FIG6.R` (FIG6D) |
| | `ST6g` RE 3.5-4 h down | yes | `figure_6/FIG6.R` (FIG6D) |
| ST7 | `ST7a` Exerkine candidate list | yes | `figure_7/FIG7.R` (FIG7A) |
| EDT1 | `EDT1` Selected molecular features associated with baseline clinical traits and acute exercise | yes | `figure_3/EDT1.R` |

`ST5a` is standalone even though its numbers are Figure 5's: no panel computes its grid, and it
reads the same helpers the FIG5 panels do. `ST2d` is standalone because the cross-tissue
comparison it tabulates is not a panel.

One threshold is worth stating here rather than leaving in a script. **`ST2h` selects its ORA
input at two different FDR cuts**, the primary analysis at 0.05 and the male-vs-female
sex-difference contrast at 0.10, because the sensitivity analysis is fitted on half the samples
per cell and 0.05 on both leaves 18 input genes against 95 at 0.10. This departs from the
published legend, which states 0.05 for both.

## What cannot be run or released

Three input files carry values attributable to individual MoTrPAC participants and may not be
published. **They are not in this repository, and no copy of them belongs in it.** A reader
with consortium access points an environment variable at a copy held outside the checkout;
there is no ignore rule for them, so a copy placed inside could be committed.

| File | What it is | Read by |
|---|---|---|
| `1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS` | Genotype PCA over 2,742 samples x 94,839 LD-pruned SNPs. `sample.id` carries 238 MoTrPAC participant IDs alongside the 1000 Genomes samples, and `eigenvect` holds 32 genotype principal components per individual. | **ED1A** |
| `CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv` | CIBERSORTx deconvolution of MoTrPAC blood RNA-seq: 951 rows keyed by vial label, 22 immune cell-type proportions per sample. | **ED3D**, and `extended_data_3/celltype_da_fit.R`, whose output **ED3E** reads |

Running **ED1A**, **ED3D** or **ED3E** requires an approved MoTrPAC consortium data-access
request and a local copy of the file, supplied through `ED1A_PCA_RDS`, `ED3D_CIBERSORTX_CSV` or
`CELLTYPE_DA_CIBERSORTX_CSV`. Without it, those panel functions stop with that message rather
than failing on a missing path, so `Rscript figures/landscape/figure_1/ED1.R` still produces ED1B, ED1C
and ED1D, and `Rscript figures/landscape/extended_data_3/ED3.R` still produces ED3A, ED3B and ED3C. **Every
other panel of ED1 and ED3 runs without these files.**

ED1A needs the genotype PCA even when the participants' own coordinates come from the released
data-hub resource: that release publishes no coordinate for any of the 2,504 reference samples,
so the 1000 Genomes cloud the participants are projected onto comes from the gated object alone.

`1kg_samples.tsv` is individual-level and ships here. It is the public 1000 Genomes sample
table, one row per 1000 Genomes individual with sex and population, freely redistributable, and
it contains no MoTrPAC participants. Individual-level and publishable are independent
properties, and the three files above are withheld because they are both.

### This is not the same question as data access

Two separate things:

| Question | Answer |
|---|---|
| Can this code be **run**? | Most of it needs `MotrpacHumanPreSuspensionData`, which is subject-level and requires a formal MoTrPAC consortium data-access request. That is true of this manuscript as a whole, not just the three files above. Several panels also read epigenomics QC and DA tables from consortium GCS prefixes that neither package ships. |
| Can a file in this repository be **published**? | Yes for everything present. A script that reads restricted data holds none of it, which is what makes the code publishable while its input is not. The three files above are the exception, and they are absent. |

A panel that declares only `MotrpacHumanPreSuspensionAnalysis` runs from the public aggregate
layer. `data_packages` in `config/panel_map.json` records which per panel.

## How this repository fits with the others

The Pre-Suspension human work is split across several repositories. Individual-level molecular
and phenotypic data cannot be distributed publicly, so they live in a separate, access-gated
package, and everything that can be released publicly, meaning aggregate results and all
analysis code, is kept clear of them.

```
                MoTrPAC BIC, consortium GCS buckets
                     (raw assay data, gated)
                                │
            motrpac-human-presuspension-repro
        normalizes omics data, applies statistical models,
        builds every data object, versions it, uploads it,
                and carries it into both packages
                                │
              ┌─────────────────┴─────────────────┐
              ▼                                   ▼
 MotrpacHumanPreSuspensionData      MotrpacHumanPreSuspensionAnalysis
 subject-level data, access-gated   aggregate results, public
              └─────────────────┬─────────────────┘
                                ▼
               motrpac-human-presuspension-acute
               manuscript figure code + QC vignettes
                       (this repository)
                                ▼
                          manuscripts
```

| Repository | What it holds | Access |
|---|---|---|
| [`motrpac-human-presuspension-repro`](https://github.com/MoTrPAC/motrpac-human-presuspension-repro) | the end-to-end rebuild pipeline and its pinned software environment | public; a full run needs consortium bucket access |
| [`MotrpacHumanPreSuspensionData`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionData) | subject-level molecular and phenotypic data objects | formal data-access request to the consortium |
| [`MotrpacHumanPreSuspensionAnalysis`](https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis) | differential analysis, group summary statistics, enrichment, clustering, feature-to-gene map, plotting functions | public |
| [`MotrpacBicQC`](https://github.com/MoTrPAC/MotrpacBicQC) | consortium QC conventions and the GCS read helpers the packages call | public |
| [`motrpac-human-presuspension-acute`](https://github.com/MoTrPAC/motrpac-human-presuspension-acute) | per-manuscript figure code and QC vignettes; this repository | code public; three panels need Data access, see above |

Every palette, timepoint order, ome name, tissue name, exercise-group label and sex label comes
from `MotrpacHumanPreSuspensionAnalysis` rather than from a figure script.

## Running things

### Before the first run

You need all five of these. Nothing here checks them for you.

| | |
|---|---|
| **R 4.4.1** | the version every panel in `assembled/` was drawn under |
| **The two data packages** | `MotrpacHumanPreSuspensionData` 2.0.3 and `MotrpacHumanPreSuspensionAnalysis` 2.0.7, installed from source. `Data` needs an approved consortium data-access request; `Analysis` is public |
| **The packages in `config/panel_packages.txt`** | 37 CRAN and Bioconductor packages. `config/required_packages.tsv` pins the five whose version changes the numbers: the two data packages, `MotrpacBicQC`, and `ggpubr >= 1.0.0` with the `rstatix` it requires |
| **`MotrpacBicQC` >= 1.8.1** | from the `develop` branch, not the `v1.8.0` tag. Analysis 2.0.7 reads it live rather than from a vendored snapshot |
| **`gsutil` on `PATH`**, authenticated | several figures download epigenomics QC and differential-analysis tables at run time and cache them under `staging/` |

Scripts must be run with **`Rscript`**, not sourced from an R console or RStudio. Each one locates
this directory from its own file path, which a console session does not have.

### Build everything

```bash
Rscript figures/landscape/run_all.R              # every figure and table
Rscript figures/landscape/run_all.R --fits       # run the four fits first
Rscript figures/landscape/run_all.R FIG2 ED6     # only those
```

Twenty-four scripts: 17 figures and 7 standalone table scripts. Each runs in its own `Rscript`
process, so a figure that fails or that grows to several GB loading its differential analysis
cannot affect the next one. The script list comes from the two manifests, so a figure added there
is picked up without editing the runner.

**The four fits are not run unless you ask.** They are hours of compute and their output is
cached, so the usual pattern is to run them once by hand and then rebuild figures freely. Without
them the sixteen fit-dependent panels report SKIPPED and name the fit to run, which is not a
failure and does not change the exit status.

It prints a per-script roll-up, lists any failures, writes `outputs/run_report.tsv` with one row
per panel and sub-table, and exits non-zero if anything failed.

`single_feature_plots.R` is not part of a run-all: all sixteen of its panels are drawn by the
figure scripts that own them.

### Build one figure

Run from the repository root:

```bash
Rscript figures/landscape/figure_2/FIG2.R            # every panel of Figure 2
Rscript figures/landscape/figure_2/FIG2.R FIG2B      # just that panel
Rscript figures/landscape/figure_4/ED6.R             # every panel of Extended Data 6
```

Each panel is a function inside its figure's script, so a whole-figure run loads the shared data
once instead of once per panel, and a single panel can still be rebuilt on its own. A run prints
one line per panel and carries on past a failure:

```
[ OK ] FIG2A    outputs/panels/figure_2/FIG2A_da_feature_matrix.pdf
[SKIP] ED1A     input not in this repository: 1kg.motrpac...pca.RDS
[FAIL] FIG2B    <the error>
        7 panel(s): 6 ok, 1 failed, 0 skipped
```

An unknown panel id is an error that lists the ids the script does build.

### Build a supplementary table

```bash
Rscript figures/landscape/figure_3/ST3.R        # every ST3 sub-table built here
Rscript figures/landscape/figure_3/ST3.R ST3g   # one sub-table
Rscript figures/landscape/figure_3/EDT1.R       # Extended Data Table 1
```

**Nineteen of the 31 sub-tables are written by a figure script, not a table script**, because
their numbers are a panel's numbers. `ST2e` arrives when `figure_2/FIG2.R` runs; there is no
`ST2e` in `figure_2/ST2.R`. `config/table_map.json` names the owning script for every sub-table.

### If a panel fails

The `[FAIL]` line carries the id and the error. Where the cause is one of the
recognisable ones, an indented `why:` follows it saying what to do:

```
[FAIL] FIG4A    there is no package called 'ComplexHeatmap'
        why: the package 'ComplexHeatmap' is not installed.
             config/panel_packages.txt lists every package the figures need.
```

It recognises a missing package, a function not in scope, a data-package object not on the
search path, a bucket read that failed for want of `gsutil` or credentials, a vendored input that
is not where it should be, an `.env` that was not found, a path that does not exist, an input
whose shape does not match the pinned package version, and data the manifest does not declare.
Anything it does not recognise is shown unadorned rather than guessed at.

`run_all.R` adds one more case: a script that exits before building anything failed in its setup,
not in a panel, so it says so and quotes the error rather than the last lines of output.

### If a panel skips

Two reasons, and the message says which:

- **`CANNOT BE RUN BY THE PUBLIC`**: the panel reads individual-level data that is not in this
  repository and may not be published. ED1A, ED3D and ED3E. See *What cannot be run or released*.
- **`needs a fit that has not been run`**: run the `*_fit.R` script the message names, then
  rebuild. See below.

Either way the other panels of that figure still build.

### The four fits

The four `*_fit.R` scripts hold the computations a panel is later drawn from rather than
anything that draws. Each sits in the folder of the figure that reads it. They are separate
scripts because a fit's output is numbered and several panels select on those numbers:
refitting per panel would let two panels draw two different LV 33s.

**Run the fit before the figure that reads it.** Nothing runs it for you, and a panel whose fit
is missing skips rather than fails.

| Run this | Fits | Before you can build | Configured by |
|---|---|---|
| `extended_data_3/sex_da_fit.R` | sex-stratified differential analysis | ED3A, and `ST2g`, `ST2h` through it | `extended_data_3/sex_da.env` |
| `extended_data_3/celltype_da_fit.R` | blood DA re-fit with cell-type covariates | ED3E | `extended_data_3/celltype_da.env` |
| `figure_3/clinical_omics_fit.R` | baseline clinical x omics regression | FIG3EFG, ED4G, `ST3f`, `ST3g`, `ST3h`, `EDT1` | `figure_3/clinical_omics.env` |
| `figure_4/plier_fit.R` | the two cross-tissue PLIER models, RNA and metabolomics | FIG4G, FIG4H, every ED6 panel | `figure_4/plier.env` |

So to build Extended Data 6 from a clean checkout:

```bash
Rscript figures/landscape/figure_4/plier_fit.R      # once, ~2 h
Rscript figures/landscape/figure_4/ED6.R            # minutes, and repeatable
```

Each fit reads its parameters from the `.env` beside it. Those files hold what was chosen by
hand, so changing a fit is an edit to one file. `figure_4/plier_lvs.env` is separate and holds
the latent-variable selections the *panels* draw, not the fit's own parameters.

After refitting PLIER, the latent-variable numbers move, and this reports the new ones:

```bash
Rscript figures/landscape/figure_4/plier_fit.R report rna     # or metab
```

**These are expensive and are cached.** The PLIER fit takes a bit over two hours and about
8.5 GB. The sex and cell-type differential analyses are each hours of mixed models, and blood
transcriptomics alone leaves a ~258 MB and a ~330 MB table behind. A fit already on disk is
reused rather than recomputed. `extended_data_3/celltype_da_fit.R` cannot be run at all without
consortium access to the CIBERSORTx input, as above, so ED3E is not buildable here.

### Layout

```
figures/landscape/
├── figure_1/ … figure_7/                   a main figure, its Extended Data supplement,
│                                           their tables, their helper and their sources
├── extended_data_3/                        ED3 and the two fits it reads
├── supplementary_figure_1/                 SF1
├── run_all.R                               build every figure and table
├── single_feature_plots.R                  the sixteen single-feature panels, defined once
├── assembled/                              the figures and tables as submitted (§ below)
├── lib/                                    panel and table export, highlights, shared helpers
├── config/                                 panel_map.json, table_map.json, highlights.json,
│                                           landscape.env, package manifests
└── outputs/panels/, outputs/tables/        what a run writes
```

A main figure and its Extended Data figure are usually one computation split across two figure
numbers, which is why they share a folder and a `*_helpers.R` file in it. A helper read by one
figure only is folded into that figure's script.

Curated selections, meaning the pathways, features and genes the manuscript highlights, are in
`config/highlights.json` rather than in a script. Those selections are made after multiple-testing
correction: BH adjustment runs across all tested pathways first.
