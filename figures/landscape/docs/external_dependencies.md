# Inputs

Everything a panel needs that the two data packages do not carry. Paths here are relative to
`figures/landscape/`, which `LANDSCAPE_ROOT` names. Almost all of it ships with the repo, under
the folder of the figure that reads it:

```
sources/<figure_dir>/
```

Most files take an environment-variable override, which is how you point a panel at a newer copy
without editing it:

```bash
ED3D_CIBERSORTX_CSV=/path/to/newer.csv Rscript figures/landscape/ED3.R ED3D
```

| Panel | File | Override |
|---|---|---|
| ED1A | `1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS` | `ED1A_PCA_RDS` |
| ED1A | `1kg_samples.tsv` | `ED1A_1KG_SAMPLES` |
| ED1A | `1kg_population_colors.csv` | `ED1A_1KG_POP_COLORS` |
| ED3D, `analysis/02_celltype_da.R` | `CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv` | `ED3D_CIBERSORTX_CSV`, `CELLTYPE_DA_CIBERSORTX_CSV` |
| ED2E | `landscape_figureED2_PTMSigDB_ORA_results_v2_Aug2026.csv` | `ED2E_PTMSIGDB_ORA_CSV` |
| FIG5B-FIG5H, ED7A, ED7B, ST5a | `Precovid_DEG_HOMER_KnownTF_Results/` (22 files) | — |
| FIG5C-FIG5H, ED7A, ED7B, ST5a | `tfproanno.RDS` | — |
| FIG6D, ED8A, ED8E, ST6c | `Gambardella-et-al-2020-table-s2.xlsx`, `TFEB_TARGET_GENES.v2024.1.Hs.tsv`, `GSM2354032_TFEB_e5_peaks.bed` | — |
| FIG6D, FIG6E, ED8A, ST6c | `TFEB_targets_v2.xlsx` | — |
| FIG6E, ST6a | `SCION_supplemental_table_edges.csv` | `FIG6_SCION_EDGES_CSV` |
| ST6b | `Merged_Network_combined_NMS_v2.csv` | `FIG6_SCION_NMS_CSV` |
| FIG7A, ED9A, ED9B, ST7a | `compartments_extracellular.tsv` | `ACUTE_COMPARTMENTS_TSV` |

**Three of the files above are not in this repository.** The ED1A genotype PCA
(`1kg.motrpac.merged.maf_0_05.ld_pruned.pca.RDS`) and the CIBERSORTx deconvolution result are
individual-level MoTrPAC data: they are absent here and may not be published. Drawing ED1A,
ED3D or ED3E, or running `analysis/02_celltype_da.R`, needs an approved consortium data-access
request and a local copy at the override named above. Every other file is tracked, and a fresh
clone builds every other panel in scope. `1kg_samples.tsv` is individual-level too, but it is the
public 1000 Genomes release and carries no MoTrPAC participants, so it ships.

Figure 5's two inputs take no override, and `tables/ST5.R` reads them from `sources/figure_5/`
rather than carrying a second copy. Both are the *record of a computation that happened once*
rather than a released resource that gets newer copies: pointing a panel at a different HOMER run or
a different motif table would not update the figure, it would make a different figure. They live
under `sources/figure_5/` and ED7 reads them from there.

## What is still fetched at run time

Six things. All but the last come from the consortium buckets rather than from a sibling repo:

| What | Who reads it | Cached to |
|---|---|---|
| The released ancestry PCs | ED1A | downloaded at run time; not kept in the repo |
| ATAC / methylCap QC matrices | FIG1C, ED1B, ED1D, SF1A–SF1F | `staging/raw-files/` |
| Raw RSEM gene counts | `analysis/01_sex_da.R`, `analysis/02_celltype_da.R` | `outputs/fits/*/raw-files/` |
| Epigenomics DA tables | FIG2A, FIG2B, FIG2Bii, FIG2C, FIG2F, ED2A, ED2C, ED2F, ED4A, SF1G, ST2d | `EPIGEN_QC_DIR/data/tmp/` |
| Raw muscle ATAC peak counts | SF1C–SF1F | `EPIGEN_QC_DIR` |
| The clinical x omics fit | ED4G, FIG3EFG, EDT1 | `outputs/fits/clinical_omics/` (`analysis/03_clinical_omics.R`) |

Every one of those is defaulted in `config/landscape.env`, so nothing has to be exported by hand
and no panel is skipped for want of a variable. The epigen cache is `EPIGEN_QC_DIR`; the directory
not existing yet is the normal first-run state. A panel that declares a fit output as `must_exist`
is skipped rather than failed when it is absent, so a figure script still draws everything else.

Every bucket read needs `gsutil` on `PATH` and consortium read access; without them the download
fails where it is attempted rather than up front.

## FIG1C, ED1B and ED1D — the epigen QC matrices

Neither data package ships the ATAC-seq or methylCap-seq QC objects. They are five tab-delimited
matrices plus their sample metadata under the analysis collection's `epigenomics/` prefix, ~6 GB
between them, downloaded on the first cold run and read from the cache after it. ED1B reads all
five, for the vial labels on their columns; ED1D excludes methylCap, so it fetches only the two
~0.5 GB ATAC matrices.

They are loaded through **`load_qc(epigen = TRUE)`** from
`MotrpacHumanPreSuspensionData`, which lists the epigenomics staging prefix, downloads what is
missing through `MotrpacBicQC::dl_read_gcp()`, and caches it under the directory `EPIGEN_QC_DIR`
names. It needs `gsutil` on `PATH` and consortium read access.

`EPIGEN_QC_DIR` is one variable for all three panels and SF1A–SF1F deliberately. They read the same
matrices out of the same collection, and two names for one cache is two caches on the first machine
where only one of them is set. `config/landscape.env` defaults it to `staging/raw-files/`, the same
place precovid-repro puts its bucket downloads, so a plain `Rscript ED1.R` fills it without setup. The
directory not existing yet is the normal first-run state.

## ED1A — ancestry PCA

Nothing in this pipeline computes a principal component. The variant calling, the merge with 1000
Genomes, the LD pruning and the `SNPRelate::snpgdsPCA` that produced these coordinates ran once,
outside any pipeline. The panel consumes the result — which is also why the manifest lists no
stochastic step for it: what fixes the figure is which copy of the PCA is read, not an RNG seed.

**The participants come from the data hub.** Their coordinates are published as
`resources/motrpac_human-precovid_1kg_pca.csv` in the analysis collection, pinned to a collection
version in the panel script and downloaded once. 175 rows: 174 participants keyed on vialLabel, plus
a `var_explained_pct` row carrying the percentage each of the 32 components explains, which is where
the axis labels come from.

**The reference cloud does not.** That resource holds the MoTrPAC projection and nothing else — no
coordinate for any of the 2,504 reference samples. The background points come from the vendored
`snpgdsPCA` object (892 KB, 2,742 samples × 32 components), which is also what
`ED1A_PC_SOURCE=cached_pca` draws the participants from.

`1kg_samples.tsv` (~4,979 rows) and `1kg_population_colors.csv` (27 rows) are public IGSR tables and
are not consortium-gated.

## ED3A — sex-stratified differential analysis

**Fitted by `analysis/01_sex_da.R`** into `outputs/fits/sex_da/`, under
`~ 0 + sex_group_timepoint + <covariates> + (1 | pid)`, with three contrasts per exercise group ×
post-baseline timepoint: the female change from baseline, the same for males, and the difference
between those two. The control arm is not subtracted, so these are not the released acute DA's
`exercise_with_controls` quantity — an effect here carries whatever the control arm also did over
the same interval. `SEX_STRATIFIED_DA_DIR` points there by
default and takes an override: set it to a directory of tables fitted elsewhere and the fit is
bypassed, which is the one way to build this panel without refitting anything.

`analysis/01_sex_da.R` records the model, the inputs and the package versions per table in
`outputs/fits/sex_da/sex_da_manifest.tsv`.

**The fit is the expensive part of this repo.** Blood transcriptomics is hours of mixed models and
the table is ~258 MB. The fit is idempotent — a table recorded in the manifest is reused unless
`SEX_DA_FORCE=TRUE` — and it fits blood, muscle and adipose by default, because ST2g reports
sex-stratified results in those tissues; ED3A plots blood alone. `SEX_DA_TISSUES=all` fits every
tissue that qualifies.

**Only the initial acute bout is fitted.** The metadata is restricted to visitcode `ADU_BAS`, which
for blood transcriptomics is 807 of the 951 samples. The other 144 are the post-training bout, and
they carry the same timepoint labels as the pre-training one — including them would put both bouts
in a single sex × group × timepoint cell and give those participants two rows under one `(1 | pid)`
intercept. This is not configurable.

**The file set is pinned by name.** The legacy code globbed `list.files()` for `.txt`, which also
matches ED3E's cell-type DA table — whose contrasts contain neither `Male` nor `Female`, so they
survive the sex filter and corrupt the contrast parser. The ported script matches on
`_da_dream-sex_differences_.*\.txt$` plus the tissue and assay.

## ED3D — CIBERSORTx deconvolution result

951 mixtures × 25 columns: 22 LM22 leukocyte subsets plus p-value, correlation and RMSE. It is the
download from a run of the CIBERSORTx web service against the LM22 signature matrix. The service has
no R implementation, so this is an archived artifact rather than a regenerable one, which is
exactly why the run parameters are recorded here. Two runs of the service over the same mixtures
do not agree if these differ, and nothing in the CSV itself says which run produced it.

`CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv`:

| Parameter | Value |
|---|---|
| Date | 2026-08-12 10:29:31 |
| Job type | Impute Cell Fractions |
| Signature matrix | `LM22.update-gene-symbols.txt` |
| Mixture file | `bloodrnanorm_081226.txt` |
| Batch correction | disabled |
| Disable quantile normalisation | true |
| Run mode | relative |
| Permutations | 100 |


ED3D plots these fractions and `analysis/02_celltype_da.R` fits ED3E's table on them, and
both read this same file. Two copies is how the model and the figure end up describing different
deconvolutions.

## ED3E — blood DA re-fit with cell-type covariates

~330 MB, tab-delimited: the five z-scored CIBERSORTx cell fractions added to the covariate table and
the standard acute dream model re-run on the blood RSEM gene counts.

**Fitted by `analysis/02_celltype_da.R`** into `outputs/fits/celltype_da/`.
`CELLTYPE_DA_DIR` points there by default; `ED3E_DA_WITH_CELL_TYPES` overrides the lookup with a
single table fitted elsewhere and bypasses the fit, which is the one way to build this panel
without refitting anything.

The five covariates are T cells, B cells, NK cells, monocytes and neutrophils, z-scored, entered as
`Technical`/`numerical` rows appended to the installed `COVARIATES_FILE`. They are not all the
populations the collapse returns: macrophages and mast cells are estimated in whole blood, where
neither is resident in any quantity, and dendritic cells fall below the 1% floor. That selection is
the legacy one and is reproduced rather than revisited.

The fit is blood-only, and that is the scope of the input rather than a default worth changing: the
deconvolution is of bulk RNA-seq against LM22, a leukocyte signature matrix, so no other tissue has
a cell fraction to add to its design.

Tables are written at `v2.0`, not the legacy `v1.3`. The version suffix records the collection
version the table was fitted *from*, so a table fitted here against the 2.0.0 data packages cannot
be mistaken for the one fitted against the v1.3 objects.

## Curated pathways and named features — config/highlights.json

Not an external input: the hand-curated selections live in `config/highlights.json` and are read
through `lib/highlights.R`. No script produces them and no criterion is recorded for what
put a pathway or feature on them.

| Section | Key | Read by | Was |
|---|---|---|---|
| `pathways` | `FIG2Bii`, `FIG2F`, `ED2E` | FIG2Bii, FIG2F, ED2E | `curated_pathways_v2.xlsx`, sheets `2B`, `2F_v3`, `ED2e` |
| `pathways` | `ED4B`, `ED5B` | ED4B and ST3d, ED5B and FIG4D | `curated_pathways.xlsx`, sheets `5`, `3` |
| `pathways` | `ED4C`, `cmeans` | ED4C; FIG4D, FIG4F, ED5C, ED5D, ED5E | a literal in ED4C; `CMEANS_PATHWAY` in `config/cmeans.env` |
| `pathways` | `FIG6D`, `ED8E` | FIG6D, ED8E | `TFEB_ChIP_targets_DA_3.5_4_hr_ORA_oct2025_highlighted.xlsx`, sheet `Heatmap terms`; the autophagy category of `TFEB_ORA_category_heatmaps_v2.R` |
| `features` | one key per panel | the sixteen single-feature plots | literals in `single_feature_plots.R` |
| `gene_lists` | `ED8C` | ED8C | `PP2a_subunits` in `Figure_S7_script.R` |
| `named_features` | — | ED1C and ST1e | `motrpac_named_features.csv` |
| `clinical_omics_features` | `EDT1` | EDT1 | Extended Data Table 1 of `Reorg of PreCAWG Landscape - Working Group Draft.docx` |

**A selection decides which pathways appear, never the statistics drawn.** The panels draw their ORA
or CAMERA results as computed (ED2E from its vendored CSV). A curated identifier those results do not
return stops FIG2Bii and FIG2F, and is logged and left out by ED2E, ED5B and FIG6D.

## ED4G, FIG3EFG, ST3g and ST3h — the clinical x omics fit

Not a file that ships with the repo and not one you supply: **`analysis/03_clinical_omics.R` fits it**, into
`outputs/fits/clinical_omics/clinical_omics_associations.tsv`, and FIG3EFG, ED4G, EDT1, ST3g and
ST3h read it from `CLINICAL_OMICS_TABLE`. Run it before any of them.

It is listed here because all five declare it as an input with `must_exist`, which means they
are **skipped rather than failed** when it is absent — the same treatment ED3A and ED3E give
their own fits.

`FIG3_CLINICAL_OMICS_TABLE` points the panels at a fit produced elsewhere and makes the fit a
no-op. To draw against a fit you did not run here, point `CLINICAL_OMICS_TABLE` at it, with its
`clinical_omics_manifest.tsv` beside it.

The legacy had no equivalent. It wrote `clinical_by_omics_res_v4_Aug2026.RData` into whatever
directory the analyst was in, and a second script `load()`ed that name back; neither file is in any
repo, so the dumbbell panels could not be rebuilt from the code that drew them.


## FIG4G, FIG4H and every ED6 panel — the cross-tissue PLIER fits

Not files that ship with the repo and not ones you supply: **`analysis/04_plier.R` fits them**, into
`outputs/fits/plier/`, and the twelve panels read them from `PLIER_DIR`. Run it before them.

Six objects, two arms:

| File | What it is |
|---|---|
| `crosstissue_rna_plier.rds` | PLIER over ~14,400 transcripts x ~1,550 samples of three tissues, against 8,482 gene sets |
| `crosstissue_rna_input.rds` | the matrix that fit was run on, and its sample annotation |
| `crosstissue_rna_response.rds` | the LV x comparison response grid FIG4G draws |
| `crosstissue_metab_*.rds` | the same three over ~460 metabolites x ~1,420 samples against the RefMet classes |

Listed here because the panels declare them as an input with `must_exist`, which means they are
**skipped rather than failed** when they are absent — the same treatment ED3A, ED3E and FIG3EFG give
their own fits.

The RNA arm is the expensive one: a bit over two hours and ~8.5 GB, against about a minute for the
metabolomics arm. Both are idempotent; `PLIER_FORCE=TRUE` refits.

**A refit renumbers the latent variables, and twelve panels select on those numbers.** That is what
`config/plier_lvs.env` exists for and what its header explains: the LV numbers there were read off the
legacy's fit, they are not derived from anything, and they hold because `analysis/04_plier.R`
reproduces that fit. A new fit, with other pins, k, frac or seed, has to be
re-checked before any of these panels is believed. Every PLIER panel except ED6Aii, ED6D and
ED6I prints the LVs it drew and how many features each carried, which is the cheap version of that check.

The legacy had no equivalent to this fit. Both scripts fitted PLIER inline and `saveRDS()`ed the
result into the analyst's working directory (`crosstissuernacombofin_plierresult_081026.rds`,
`metab_plier_081426.rds`); neither file is in any repo, so twelve panels could not be rebuilt from the
code that drew them without refitting.


## Supplementary Figure 1 — the muscle ATAC-seq raw counts

One file beyond what the loaders above fetch: `*_counts_*.txt.gz`, the raw peak counts with one
column per library, under `QUANT_BUCKET/*/epigenomics/t06-muscle/epigen-atac-seq/`. It is
resolved by listing that prefix with `gsutil ls` and matching on the tissue and ome, and a match
count other than one is an error; `SF1.R` holds the resolver. Read by SF1C–SF1F.

It lands in `EPIGEN_QC_DIR` through `MotrpacBicQC::dl_read_gcp()`, next to the qc-norm matrix
`load_qc(epigen = TRUE)` puts there. The counts are voom-transformed at run time; the qc-norm
matrix supplies which peaks and which libraries are kept.

## The epigenomics differential-analysis tables

Five files, ~7.7 GB, read by every panel that draws a cross-tissue overlap. They are **not** an
unauthenticated HTTPS download, whatever an older note may have said:
`load_differential_analysis(epigen = TRUE)` lists

```
gs://pre-cawg/staging_20260806
```

with `gsutil ls -R`, filters to `*_da_*.txt`, and pulls each match through
`MotrpacBicQC::dl_read_gcp()`. It needs `gsutil` on `PATH` and consortium read access. The bucket
is the package constant `.STAGING_BUCKET`, which mirrors `STAGING_BUCKET` in precovid-repro's
`config/pipeline.env`; a released prefix can be passed as `bucket=` instead.

| Tissue | Assay | Version |
|---|---|---|
| muscle, PBMC | `epigen-atac-seq` | v2.1 |
| muscle, EDTA, adipose | `epigen-methylcap-seq` | v1.2 |

**The cache does not make this offline-capable.** Files land under `EPIGEN_QC_DIR/data/tmp/` and a
second build reuses them, but the cache is consulted *per file, after the bucket listing names it*.
With an expired credential the listing returns nothing, the loader reports "No epigenetic
ome/tissue combination was found" and **returns the five package-shipped omes as though nothing
were missing** — a panel then builds cleanly from five omes instead of seven. Measured: the
`Muscle_only` ORA background halves, 30,168 genes to 15,237.

`cross_tissue_da()` in `lib/da_overlap_helpers.R` therefore checks that both epigen assays
came back and stops if they did not. Refresh with `gcloud auth login`.

Only the versions above are in the bucket. A cache holding a superseded copy — a `v2.0` ATAC table,
say — is dead weight rather than a hazard, because the bucket listing is what selects.


## ED2E — the PTMsigDB over-representation result

`landscape_figureED2_PTMSigDB_ORA_results_v2_Aug2026.csv`, 837 rows over three phosphosite sets.
ED2E draws from it and `ST2f` carries it; it is the published result.

The panel still **recomputes** the pass on every build and logs the comparison, but nothing from
that recomputation reaches the figure or the table. It exists so a divergence shows up in the
build log rather than being discovered later. The check reads site confidence **per tissue**, from
`*_PROT_PH_QC$feature_metadata` rather than from `HUMAN_FEATURE_TO_GENE`, because the two disagree
by construction:

`HUMAN_FEATURE_TO_GENE` is keyed on `(assay, feature_id)` with no tissue column, so its `prot-ph`
rows are the **union** of the two tissues — 18,548 muscle plus 21,022 adipose over 7,865 shared,
exactly 31,705. Localization confidence is a per-tissue measurement and **859 of those shared sites
disagree** (782 muscle-only, 77 adipose-only); with one row per site that table collapses every
disagreement to `FALSE`. `flanking_sequence` does not diverge — all 7,865 agree — because it is a
property of the protein, not the assay.

Where the check lands:

| | recomputed | published |
|---|---|---|
| `Muscle_only` input | 3,850 | 3,836 |
| `Muscle_only` background | 16,390 | 15,818 |
| `Adipose_only` input | 206 | 225 |
| `Adipose_Muscle` background | 4,592 (adipose ∩ muscle) | 15,818 (the muscle background) |

Nine of the 48 drawn cells cross the adjusted-p cut between the two. The recomputation differs from
the published result in the confidence definition it applies and in testing the adipose-and-muscle
set against the two-tissue intersection rather than the muscle background.

## Figure 6, ED8 and ST6 — the SC-ION network and the TFEB ChIP evidence

Six files under `sources/figure_6/`, all from precovid-analyses PR #104
(`nmclark2/Figure7_v2.0`, `figures/landscape/figure_7/`). Two kinds.

**The ChIP evidence is a released resource read verbatim.** Three published TFEB ChIP datasets —
Gambardella et al. 2020 (peak table, bound = a promoter peak within 1 kb), the MSigDB GTRD set
`TFEB_TARGET_GENES` (membership), and Doronzo et al. 2019 (hg19 peaks, annotated with `annotatr`
against `hg19_basicgenes`, bound = a promoter or 1-5 kb peak) — plus Settembre et al. 2013, which
is two genes by ChIP-qPCR and is coded in `helpers/FIG6_ED8.R`. `tfeb_chip_evidence()` scores any
vector of gene symbols against all four; FIG6D calls it over every muscle transcript and writes
the result as `ST6c`, ED8A and ED8E call it again on their own gene lists. The `annotatr`
annotation is built from the installed `TxDb.Hsapiens.UCSC.hg19.knownGene` and `org.Hs.eg.db`,
so nothing is downloaded.

**The network tables are the record of a computation that is not reproduced here.** The SC-ION
inference (`run_SCION()` in the retired `MotrpacHumanPreSuspension` package, driven by
`scion_figures.Rmd` on a cluster with 100 permutations per network), the trim to edge weight
> 0.1, the merge of the four networks (muscle-only and blood-to-muscle, EE and RE) in Cytoscape,
and the Network Motif Score per node all happened once, outside this repository.
None of that code is in the PR beyond the Rmd, whose inputs and outputs are not versioned, and the
four per-network SC-ION output TSVs it produced sit on a lab share this repo cannot reach. What the
PR does carry is the merged result: the edge table (`SCION_supplemental_table_edges.csv`,
116,419 edges, which `ST6a` re-emits and FIG6E reads for the TFEB and RXRA target sets), the node
scores (`Merged_Network_combined_NMS_v2.csv`, which `ST6b` re-emits) and the 705 TFEB targets read
off the merged network (`TFEB_targets_v2.xlsx`, which `ST6c` flags and ED8A draws). FIG6E checks
the target list against the edge table on every build and reports any difference. Pointing a
panel at a different network through the override would make a different figure, which is why
the two overrides exist for a rerun and not for a newer copy.

The curated FIG6D term list and the ED8E gene set are in `config/highlights.json`; the PP2A subunit
list ED8C draws is there too, under `gene_lists`.

## Figure 5, ED7 and ST5a — the HOMER runs

`findMotifs.pl` is not in this pipeline and is not run by it. Twenty-two differentially expressed
gene lists — one per tissue x arm x timepoint, thresholded on nominal p — were written out of the
transcript DA tables, HOMER was run on each against its human background, and its `knownResults.txt`
is what every Figure 5 panel reads. They are the enrichment; nothing here recomputes it.

Each file scores the same 471 motifs and carries, in its target-sequence column header, the number
of genes that run was given. FIG5A does not read these files: it counts the nominally significant
transcripts straight from the DA tables. The lists themselves were not kept, and the
files are named for a timepoint slug (`early`, `mid`, `late`, `immpost`, `during20`, `during40`)
that no longer appears anywhere else; `TF_HOMER_TIMEPOINT_SLUGS` in
`helpers/FIG5_ED7.R` is the only link between a file on disk
and the contrast it came from.

## Figure 5, ED7 and ST5a — the motif annotation

`tfproanno.RDS` is a curated HOMER motif -> gene table carried over from PASS1B: 402 of the 471
motifs, keyed by full motif name. It is what maps `ETS:RUNX(ETS,Runt)/Jurkat-RUNX1-ChIP-Seq...` to
RUNX1 and `TATA-Box(TBP)/Promoter/Homer` to TBP, and it raises the number of motifs resolving to a
gene this study measured from 299 to 416.

**Only its `Gene.Name` column is read.** Its `Ensembl` column holds rat ids (`ENSRNOG...`) and its
six `*.Pro.ID` columns are rat tissue protein ids; every human id a panel uses is resolved from the
symbol against `HUMAN_FEATURE_TO_GENE`. The legacy applied fifteen hand corrections to this file
by row number; they are not reproduced.

## FIG7A, ED9A, ED9B and ST7a — the COMPARTMENTS extracellular scores

`compartments_extracellular.tsv` is the rows of the COMPARTMENTS integrated channel
(`https://download.jensenlab.org/human_compartment_integrated_full.tsv`, downloaded 2026-09-14,
4,856,985 rows, md5 `9889080cc2675c13c10f13a30cb25929`) whose location is Extracellular region,
Extracellular space, Extracellular exosome or Extracellular vesicle: 58,460 rows, with a header
added. The helper takes the best score over the four locations per gene symbol. PR #106 fetched the
full file at run time; the channel is rebuilt without a version, so the figure is pinned to this copy.
