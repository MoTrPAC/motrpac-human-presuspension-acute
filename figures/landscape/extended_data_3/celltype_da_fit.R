#!/usr/bin/env Rscript
# analysis/02_celltype_da.R — fit the cell-type sensitivity analysis ED3E is
# drawn from.
#
# ---------------------------------------------------------------------------
# THIS FIT CANNOT BE RUN FROM THIS REPOSITORY, AND ITS INPUT CANNOT BE RELEASED.
#
# The five covariates are rebuilt from a CIBERSORTx deconvolution of MoTrPAC
# blood RNA-seq: 951 rows keyed by vial label, 22 cell-type proportions per
# sample. That is individual-level data: it is not distributed with this
# repository and may not be published. Running this fit requires an approved
# MoTrPAC consortium data-access request and a local copy of the file at
# extended_data_3/sources/CIBERSORTx_PrecovidBlood_RNADecon_Results_081226.csv,
# or CELLTYPE_DA_CIBERSORTX_CSV pointing at one.
#
# ED3D plots the same file and is gated for the same reason. ED3E reads this
# fit's output and says so rather than failing.
# ---------------------------------------------------------------------------
#
# One table per tissue x ome, under the freeze naming convention, into
# CELLTYPE_DA_DIR. The model is the released acute model
#
#   ~ 0 + group_timepoint + <covariates> + (1 | pid)
#
# with five deconvolution-derived cell fractions — T cells, B cells, NK cells,
# monocytes, neutrophils — added to the covariate table as z-scored numeric
# terms. ED3E plots this fit against the released one, so everything except those
# five terms has to be the same fit.
#
# THE REGENERATION CHAIN THIS SCRIPT CLOSES. The legacy fit
# (revisions/landscape/cellular_deconvolution/DA_with_cell_types.R in this
# repository's history) was broken in two places and its output was an
# operator-supplied 318 MB file that had to be treated as frozen:
#
#   1. It read a `cibersort_with_metadata.csv` intermediate that nothing produced
#      any more — it was written by a chunk of deconvolution_diagnosis.Rmd that
#      had to be run by hand first. This script rebuilds that intermediate from
#      the CIBERSORTx result and writes it out as celltype_covariates_<tissue>.tsv,
#      so the covariate table has the same provenance as the fit.
#   2. It source()d every file in a data-raw directory of a SOURCE CHECKOUT of
#      MotrpacHumanPreSuspensionAnalysis, at whatever commit was present. That
#      directory is now retired. The engine is vendored below.
#
# Blood only, and that is the scope of the input rather than a default: the
# deconvolution is of bulk RNA-seq against LM22, a leukocyte signature matrix, so
# no other tissue has a cell fraction to add to its design.
#
# The engine is vendored separately from analysis/01_sex_da.R's, which vendors
# the same upstream code for a different model. That one dropped
# `custom_covariates`; this fit is the case that argument exists for, since the
# five cell fractions are exactly a second covariate table. Each fit's engine is
# pinned to the one model it fits, so neither can acquire the other's behaviour.
#
# Configuration comes from extended_data_3/celltype_da.env, through the environment.
# Nothing here is a literal that a run might want to change.
#
#   Rscript figures/landscape/analysis/02_celltype_da.R
#
# Writes:
#   <CELLTYPE_DA_DIR>/human-precovid-sed-adu_<tissue_code>_<ome>_da_dream-acute-cell_types_v<version>.txt
#   <CELLTYPE_DA_DIR>/celltype_covariates_<tissue>.tsv   the covariate table, per tissue
#   <CELLTYPE_DA_DIR>/celltype_da_manifest.tsv           one row per table: the formula,
#                                      which five covariates entered it, the feature, contrast
#                                      and sample counts, the counts source, the deconvolution
#                                      it was built from, the two data-package versions and the
#                                      eBayes estimator. A table listed there is reused, not refitted.

suppressPackageStartupMessages({
  library(dplyr)
  # Both data packages have to be ATTACHED, not merely namespace-qualified.
  # load_qc() resolves its lazy-loaded objects with eval(parse(text = "BLOOD_TRNSCRPT_QC"))
  # and similar, which reaches the package's lazydata environment only through
  # the search path — so a qualified call fails with *object 'BLOOD_TRNSCRPT_QC'
  # not found*.
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

# One level below the landscape root. panel_export.R is sourced for
# landscape_root() and for config/landscape.env, which it reads on load.
here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))

# squish_proportions() lives with the ED3 panels, the ones that draw the
# collapsed fractions, and is sourced across rather than copied. That coupling is
# the point: the covariates this script fits on have to be the same numbers ED3D
# plots, and two copies of the LM22 -> coarse-population map is exactly how they
# would stop being.
source(file.path(landscape_root(), "extended_data_3", "ED3_helpers.R"))

# ---- extended_data_3/celltype_da.env -------------------------------------------------

#' Read extended_data_3/celltype_da.env into the environment, once per run.
#'
#' Shell syntax, because it is a list of settings and reads like one; parsed
#' here rather than run through a shell so an R session needs no subprocess.
#' Values already exported WIN, which is how a one-off override works:
#'
#'   CELLTYPE_DA_FORCE=TRUE Rscript figures/landscape/analysis/02_celltype_da.R
.celltype_da_load_env <- function() {
  root <- landscape_root()
  path <- file.path(root, "extended_data_3", "celltype_da.env")
  if (!file.exists(path)) {
    stop("extended_data_3/celltype_da.env is missing from ", root, call. = FALSE)
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export ")) next
    body <- sub("^export ", "", line)
    name <- sub("=.*$", "", body)
    if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", name)) next
    if (nzchar(Sys.getenv(name, unset = ""))) next          # already set: leave it
    value <- sub("^[^=]*=", "", body)
    value <- sub('^"', "", sub('"$', "", value))
    # "${NAME:-default}" -> default
    value <- sub("^\\$\\{[A-Za-z0-9_]+:-", "", value)
    value <- sub("\\}$", "", value)
    value <- gsub("\\$\\{?LANDSCAPE_ROOT\\}?", root, value)
    # This is not a shell. A value still carrying a substitution cannot be
    # resolved here, and setting the literal is worse than leaving it unset.
    if (grepl("[$`]", value)) next
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}

.celltype_da_load_env()

# ---- the model -------------------------------------------------------------
#
# The engine ED3E's input is fitted with. It is unexported internal code in a
# SOURCE CHECKOUT of MotrpacHumanPreSuspensionAnalysis
# (`data-raw/generate_differential_analysis/`), and that directory is now
# retired — every file in it is a tombstone pointing at precovid-repro. It is
# vendored here.
#
# Provenance, in order:
#   1. MotrpacHumanPreSuspensionAnalysis, data-raw/generate_differential_analysis/
#      generate_differential_modeling_functions.R and generate_DA_inputs.R, at the
#      commit that produced the published tables. Now retired in place.
#   2. precovid-repro, scripts/10_build_data/09_build_da/da_common.R, which
#      vendored the same functions for the acute model and is the maintained
#      copy. The bodies below match it, including the two corrections it
#      documents (the topTable keying, the explicit eBayes estimator).
#   3. revisions/landscape/cellular_deconvolution/DA_with_cell_types.R in this
#      repository's history, which is what called the engine for this particular
#      fit and which supplied the five extra covariate rows.

#' Construct covariate metadata and the model formulas.
#'
#' Vendored from the upstream `process_covariates()`, with `custom_covariates`
#' kept: it is how the five cell-type fractions enter the design. When it is
#' NULL the installed `COVARIATES_FILE` is used and the result is the released
#' acute model — which is what makes the with/without comparison ED3E draws a
#' comparison of one thing.
#'
#' Everything else — which covariates are scaled, which are factors, the three
#' interaction terms, the order the formula string is assembled in — is
#' unchanged, and the order matters: it decides which columns
#' makeContrastsDream() drops.
#'
#' @param meta Sample metadata for the tissue x ome being fitted, one row per
#'   sample, rownames set to vialLabel.
#' @param selected_ome Assay code, e.g. "transcript-rna-seq".
#' @param tissue_input Tissue name, e.g. "blood".
#' @param include_technical Whether technical covariates enter the model.
#' @param custom_covariates A covariate table replacing the installed
#'   COVARIATES_FILE. Same four columns: ome, tech_or_design, data_type,
#'   covariate, tissue.
#'
#' @return A list carrying the processed metadata, the formulas, and the original
#'   metadata (which .convert_dream_output() does not use but the upstream
#'   signature passes through).
process_covariates <- function(meta,
                               selected_ome,
                               tissue_input,
                               include_technical = TRUE,
                               custom_covariates = NULL) {
  covariates_return <- list()
  covariates_return[["original_meta"]] <- meta

  input_covariates <- if (!is.null(custom_covariates)) {
    custom_covariates
  } else {
    MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE
  }

  covariates <- input_covariates %>%
    as.data.frame() %>%
    dplyr::filter(ome == selected_ome) %>%
    dplyr::filter(tissue == "all" | tissue == tissue_input)

  num_cov <- covariates %>% dplyr::filter(data_type == "numerical")
  factor_cov <- covariates %>% dplyr::filter(data_type == "factor")

  sel_meta <- meta %>%
    dplyr::select(dplyr::all_of(covariates$covariate)) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(num_cov$covariate),
                                ~ scale(.) %>% as.numeric())) %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(factor_cov$covariate),
                                ~ as.factor(.) %>% droplevels())) %>%
    dplyr::mutate(group_timepoint = droplevels(interaction(randomGroupCode, Timepoint))) %>%
    dplyr::mutate(visit_group_timepoint = droplevels(interaction(visitcode, randomGroupCode, Timepoint))) %>%
    dplyr::mutate(sex_group_timepoint = droplevels(interaction(Sex, randomGroupCode, Timepoint)))

  technical_covs <- covariates %>% dplyr::filter(tech_or_design == "Technical")
  full_formula <- names(sel_meta)[!names(sel_meta) %in%
                                    c("randomGroupCode", "Timepoint", "visitcode", "pid",
                                      "group_timepoint", "visit_group_timepoint",
                                      "sex_group_timepoint")]
  design_covs <- c(full_formula[!full_formula %in% as.character(technical_covs$covariate)],
                   "group_timepoint")
  if (!include_technical) {
    full_formula <- full_formula[!full_formula %in% technical_covs$covariate]
  }

  sex_diff_covs <- full_formula[!full_formula %in% c("Sex", "codedsiteid")]
  sex_formula_string <- paste(sex_diff_covs, collapse = " + ")
  formula_string_sex_differences <- paste("~ 0 + sex_group_timepoint + ",
                                          sex_formula_string, "+ (1 | pid)")

  formula_string <- paste(full_formula, collapse = " + ")

  # group_timepoint is re-added first deliberately. The order of the terms in the
  # string changes which columns makeContrastsDream() drops when a factor level
  # interacts with another term, so it is fixed here rather than left to however
  # the covariate table happens to be sorted.
  formula_string_full <- paste("~ 0 + group_timepoint + ", formula_string, "+ (1 | pid)")
  formula_string_training <- paste("~ 0 + visit_group_timepoint + ", formula_string,
                                   "+ (visitcode | pid)")
  non_mixed_model <- paste("~ 0 + group_timepoint + ", formula_string)

  covariates_return[["technical_cov"]] <- technical_covs
  covariates_return[["design_cov"]] <- design_covs
  covariates_return[["full_formula"]] <- formula_string_full
  covariates_return[["training_formula"]] <- formula_string_training
  covariates_return[["sex_differences_formula"]] <- formula_string_sex_differences
  covariates_return[["metadata"]] <- sel_meta
  covariates_return[["non_mixed_model"]] <- non_mixed_model

  return(covariates_return)
}

#' The acute contrast set, verbatim from upstream.
#'
#' Baseline-anchored rather than fully pairwise: for each post-baseline
#' timepoint, the endurance and resistance changes from pre-exercise both on
#' their own and relative to the control change, endurance against resistance,
#' and the control change itself. Resistance has no during-exercise blood draw,
#' so those three cells are skipped rather than requested and dropped. The three
#' pre-exercise between-group contrasts are appended last.
#'
#' @param metadata The processed metadata matrix.
.generate_contrasts_acute <- function(metadata) {
  pre_contrast_expressions <- c()
  timepoints <- unique(metadata$Timepoint)
  for (tp in timepoints) {
    if (tp == "pre_exercise") next

    contrast_Endur_Cntrl <- sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
    contrast_Endur <- sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise", tp)
    contrast_Resist_Cntrl <- sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise - group_timepointADUControl.%s + group_timepointADUControl.pre_exercise", tp, tp)
    contrast_Resist <- sprintf("group_timepointADUResist.%s - group_timepointADUResist.pre_exercise", tp)
    contrast_Endur_Resist <- sprintf("group_timepointADUEndur.%s - group_timepointADUEndur.pre_exercise - group_timepointADUResist.%s + group_timepointADUResist.pre_exercise", tp, tp)
    contrast_Cntrls <- sprintf("group_timepointADUControl.%s - group_timepointADUControl.pre_exercise", tp)

    if (tp == "during_20_min" || tp == "during_40_min") {
      contrast_Resist_Cntrl <- NULL
      contrast_Endur_Resist <- NULL
      contrast_Resist <- NULL
    }

    pre_contrast_expressions <- c(pre_contrast_expressions,
                                  contrast_Endur_Cntrl,
                                  contrast_Endur,
                                  contrast_Resist_Cntrl,
                                  contrast_Resist,
                                  contrast_Endur_Resist,
                                  contrast_Cntrls)
  }

  pre_ex_endur_res <- "group_timepointADUEndur.pre_exercise - group_timepointADUResist.pre_exercise"
  pre_ex_res_cntrl <- "group_timepointADUResist.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_ex_endur_cntrl <- "group_timepointADUEndur.pre_exercise - group_timepointADUControl.pre_exercise"
  pre_contrast_expressions <- c(pre_contrast_expressions, pre_ex_endur_res,
                                pre_ex_res_cntrl, pre_ex_endur_cntrl)
  return(pre_contrast_expressions)
}

#' The empirical-Bayes estimator eBayes uses: limma's moment estimator.
#'
#' Passed explicitly, never left to limma's default. squeezeVar() resolves a NULL
#' `legacy` as identical(min(df), max(df)), and dream gives every feature its own
#' fractional residual df, so that is always FALSE and limma >= 3.62.0 silently
#' takes the newer estimator. The choice moves df.prior/s2.prior, hence every
#' feature's moderated variance, t and p — so it is stated, and it matches the
#' released acute DA this table is compared against. ED3E is that comparison, so
#' a mismatch here would show up as a difference the cell-type covariates get
#' blamed for.
.ebayes_legacy <- function() {
  if (utils::packageVersion("limma") < "3.62.0") {
    stop("limma ", utils::packageVersion("limma"), " has no `legacy` argument to ",
         "eBayes(); this fit requires >= 3.62.0 so the estimator is chosen ",
         "explicitly rather than by the data", call. = FALSE)
  }
  TRUE
}

#' The parallel backend, or NULL when parallelism is off.
#'
#' MulticoreParam, not SnowParam("SOCK"). A SOCK worker is a cold R session that
#' inherits the manager's options(warn = 2) through exportglobals, so the first
#' findbars() it evaluates turns lme4's once-per-session reformulas deprecation
#' warning into an error and bptry drops that feature from the fit — silently,
#' one feature per worker. MulticoreParam forks, so each worker inherits an
#' already-armed warning cache and never fires it.
.celltype_da_bpparam <- function() {
  cores <- suppressWarnings(as.integer(Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", "")))
  if (is.na(cores) || cores < 2) return(NULL)
  BiocParallel::MulticoreParam(cores, progressbar = TRUE)
}

#' Fit the acute mixed model.
#'
#' The `model_type == "acute"` branch of the upstream `run_dream()`, with the
#' training and sex branches dropped: this script fits one model and a
#' `model_type` argument here could only ever be wrong.
#'
#' @param expression_object A model-ready expression object — for RNA-seq, an
#'   edgeR DGEList of raw counts.
#' @param process_metadata The list returned by process_covariates().
#' @param voom Whether voom precision weights are applied first. TRUE for
#'   RNA-seq.
#'
#' @return The eBayes-moderated dream fit.
run_dream_acute <- function(expression_object,
                            process_metadata,
                            voom = FALSE) {
  param <- .celltype_da_bpparam()

  meta_matrix <- process_metadata$metadata
  # Reorder the covariate rows to follow the expression columns. They should
  # already agree; the fit is what breaks silently if they ever do not, so the
  # match is done rather than assumed. Upstream does the same.
  meta_matrix <- meta_matrix[match(colnames(expression_object), rownames(meta_matrix)), ]

  contrast_expressions <- .generate_contrasts_acute(meta_matrix)
  formula <- stats::as.formula(process_metadata$full_formula)

  L <- variancePartition::makeContrastsDream(formula, meta_matrix,
                                             contrasts = contrast_expressions)

  if (voom) {
    if (!is.null(param)) {
      suppressWarnings({
        expression_object <- variancePartition::voomWithDreamWeights(
          expression_object, formula, meta_matrix, BPPARAM = param)
      })
    } else {
      expression_object <- variancePartition::voomWithDreamWeights(
        expression_object, formula, meta_matrix)
    }
  }

  if (!is.null(param)) {
    suppressWarnings({
      fit <- variancePartition::dream(expression_object, formula, meta_matrix,
                                      L = L, BPPARAM = param)
    })
  } else {
    fit <- variancePartition::dream(expression_object, formula, meta_matrix, L = L)
  }
  fit <- variancePartition::eBayes(fit, legacy = .ebayes_legacy())
  return(fit)
}

# ---- fit -> table ----------------------------------------------------------

#' Convert a dream fit into the long table the freeze ships.
#'
#' One block of rows per contrast, concatenated and sorted by adjusted p value.
#' Only coefficients naming an actual contrast (those carrying a "-") are kept,
#' which drops the covariate coefficients — including the five cell-type ones.
#'
#' The column set matches the released acute DA, which is what ED3E joins this
#' table against. CI.L/CI.R are not among them: the freeze does not ship them for
#' the acute model, and variancePartition's topTable computes them wrongly for a
#' dream fit anyway (every contrast's interval uses the first contrast's df,
#' because df.total is a features x contrasts matrix indexed linearly).
#' analysis/01_sex_da.R rebuilds them because ED3A draws error bars; ED3E does
#' not, so they are left out rather than shipped wrong.
#'
#' @param fit A fit from run_dream_acute().
#' @param formula The model formula, recorded in the full_model column.
#' @param tissue,ome What was fitted, recorded per row.
.convert_dream_output <- function(fit,
                                  formula = NULL,
                                  tissue = NULL,
                                  ome = NULL) {
  full_contrasts <- colnames(stats::coef(fit))
  comparison_subset <- full_contrasts[grep("-", full_contrasts)]

  res_tissue <- data.frame()
  for (contrast in comparison_subset) {
    res <- variancePartition::topTable(fit, coef = contrast, number = Inf,
                                       p.value = 1, confint = FALSE)

    # Keyed by name, not by position. topTable defaults to sort.by = "p", so
    # `res` comes back permuted relative to the fit, while fit$rdf and fit$logLik
    # are per-feature NAMED vectors in the fit's original order — rdf varies per
    # feature because variancePartition derives an approximate residual df from
    # each feature's variance components. Assigning them positionally gives every
    # row another feature's df and logLik, and does it silently, because
    # feature_id comes from rownames() and stays correct.
    res_single <- res %>%
      dplyr::mutate(degrees_of_freedom = fit$rdf[rownames(.)]) %>%
      dplyr::mutate(logLik = fit$logLik[rownames(.)]) %>%
      dplyr::mutate(feature_id = rownames(.)) %>%
      dplyr::mutate(assay = ome) %>%
      dplyr::mutate(contrast = contrast) %>%
      dplyr::mutate(full_model = formula) %>%
      dplyr::mutate(tissue = tissue) %>%
      dplyr::rename(p_value = P.Value) %>%
      dplyr::rename(adj_p_value = adj.P.Val) %>%
      dplyr::select(assay, feature_id, z.std, logFC,
                    degrees_of_freedom, logLik, t, AveExpr,
                    p_value, adj_p_value, contrast, full_model)

    res_tissue <- rbind(res_tissue, res_single)
  }
  res_tissue <- res_tissue %>% dplyr::arrange(adj_p_value)
  return(res_tissue)
}

# ---- naming ----------------------------------------------------------------

# The freeze stem prefix. Every released table under the analysis bucket carries
# it, and ED3E finds its input by matching the rest of the name, so the table
# this script writes uses the same convention.
CELLTYPE_DA_STEM_PREFIX <- "human-precovid-sed-adu"

# The model name in the file stem. It is `dream-acute-cell_types` and not
# `dream-acute`, which is what keeps this table distinguishable from the released
# acute DA it is compared against — including inside a directory holding both.
CELLTYPE_DA_MODEL <- "dream-acute-cell_types"

#' The BIC tissue code for a tissue x ome, from the installed OME_TISSUE_CODE.
#'
#' The code is not derivable from the tissue name — blood transcriptomics is
#' t04-blood-rna while blood proteomics is not — so it is looked up rather than
#' composed.
celltype_da_tissue_code <- function(ome, tissue) {
  hit <- MotrpacHumanPreSuspensionAnalysis::OME_TISSUE_CODE %>%
    dplyr::filter(ome == !!ome, tissue == !!tissue)
  if (nrow(hit) == 0) {
    stop("no tissue code for ", tissue, " / ", ome, " in OME_TISSUE_CODE",
         call. = FALSE)
  }
  hit[["tissue_code"]][1]
}

#' The file name for one written table.
#'
#' `<prefix>_<tissue_code>_<ome>_da_dream-acute-cell_types_v<version>.txt`, which
#' is what the upstream write_with_path_name() composed and what ED3E matches on:
#' the panel pins its file set by the model name plus the tissue and assay, so a
#' table written under any other name is invisible to it.
#'
#' @param ome Assay code.
#' @param tissue Tissue name.
#' @param version Version label. Not a bucket release — see extended_data_3/celltype_da.env.
celltype_da_file_name <- function(ome, tissue, version) {
  paste0(paste(CELLTYPE_DA_STEM_PREFIX, celltype_da_tissue_code(ome, tissue), ome,
               "da", CELLTYPE_DA_MODEL, sep = "_"),
         "_v", version, ".txt")
}

#' Whether a table on disk is one this script finished writing.
#'
#' Existence is not the test. A fit killed partway through write.table() leaves a
#' file that opens, parses, and is missing however many features the process did
#' not get to — and the run would then skip it forever as "already fitted".
#' A manifest row is written only after the table is closed, so the row is the
#' completion record and the file alone is not.
#'
#' @param manifest_tsv Path to celltype_da_manifest.tsv.
#' @param file_name Base name of the table.
.celltype_da_table_is_complete <- function(manifest_tsv, file_name) {
  if (!file.exists(manifest_tsv)) return(FALSE)
  manifest <- tryCatch(
    read.csv(manifest_tsv, sep = "\t", check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(manifest) || !"file" %in% names(manifest)) return(FALSE)
  file_name %in% manifest$file
}

# ---- configuration ---------------------------------------------------------

env_or_stop <- function(name) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) stop(name, " is not set — it is defined in ",
                       "extended_data_3/celltype_da.env or config/landscape.env",
                       call. = FALSE)
  v
}
env_flag <- function(name, default = FALSE) {
  v <- toupper(Sys.getenv(name, unset = ""))
  if (!nzchar(v)) return(default)
  v %in% c("TRUE", "1", "YES")
}

out_dir        <- env_or_stop("CELLTYPE_DA_DIR")
manifest_tsv   <- env_or_stop("CELLTYPE_DA_MANIFEST_TSV")
version        <- env_or_stop("CELLTYPE_DA_VERSION")
quant_bucket   <- env_or_stop("QUANT_BUCKET")
cibersortx_csv <- env_or_stop("CELLTYPE_DA_CIBERSORTX_CSV")
gsutil         <- Sys.getenv("GSUTIL", unset = "gsutil")
force          <- env_flag("CELLTYPE_DA_FORCE")

requested_tissues <- trimws(strsplit(env_or_stop("CELLTYPE_DA_TISSUES"), ",")[[1]])
requested_tissues <- requested_tissues[nzchar(requested_tissues)]

# The raw counts are downloaded once and reused. They are large and the download
# is the only part of this script that needs the network.
cache_dir <- file.path(out_dir, "raw-files")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# Transcriptomics, and only transcriptomics. The deconvolution is of bulk RNA-seq
# against a leukocyte signature matrix, so there is no cell fraction to add to any
# other assay's design.
CELLTYPE_DA_OME <- "transcript-rna-seq"

# The five coarse populations that enter the design, in the make.names() form the
# formula needs. Not all of the populations squish_proportions() returns:
# macrophages and mast cells are estimated in whole blood — where neither is
# resident in any quantity — and dendritic cells fall below the 1% floor. The
# legacy fit used these five and the choice is reproduced rather than revisited.
# A deconvolution that does not carry all five stops the run rather than
# quietly fitting a smaller model.
CELLTYPE_DA_COVARIATES <- c("T.cells", "B.cells", "NK.cells", "Monocytes", "Neutrophils")

# What the stage driver wrote to its report TSV goes to the console instead. The
# status word is kept so a run still reads as a list of checks.
record <- function(status, check, detail) {
  message(sprintf("[celltype-da] %-4s %s — %s", status, check, detail))
}
manifest_rows <- list()

# ---- is there anything to fit? ---------------------------------------------
# ED3E_DA_WITH_CELL_TYPES names a table fitted elsewhere — the frozen v1.3
# artifact, say — and bypasses this fit entirely.

override_table <- Sys.getenv("ED3E_DA_WITH_CELL_TYPES", unset = "")
if (nzchar(override_table)) {
  record("SKIP", "celltype_da:override",
         paste0("ED3E_DA_WITH_CELL_TYPES names a table fitted elsewhere (",
                override_table,
                ") — ED3E will read that file and nothing is fitted here"))
  quit(status = 0, save = "no")
}

# ---- what this fit needs on the machine ------------------------------------
# The blood RSEM counts come from QUANT_BUCKET, and the five covariates from the
# CIBERSORTx result, which is not in this repository — see the note at the top.

if (!nzchar(Sys.which(gsutil)) && !file.exists(gsutil)) {
  stop(gsutil, " is not on PATH — the fit reads raw RSEM counts from ",
       quant_bucket, call. = FALSE)
}
bucket_readable <- suppressWarnings(
  system2(gsutil, c("ls", shQuote(quant_bucket)),
          stdout = FALSE, stderr = FALSE) == 0
)
if (!bucket_readable) {
  stop("cannot read ", quant_bucket,
       " — the fit needs the raw counts; run gcloud auth login", call. = FALSE)
}
record("PASS", "celltype_da:bucket", paste0("readable: ", quant_bucket))

if (!file.exists(cibersortx_csv)) {
  record("FAIL", "celltype_da:cibersortx",
         paste0("no deconvolution result at ", cibersortx_csv,
                " — the five cell-fraction covariates come from it. That file ",
                "is individual-level data, is not distributed with this ",
                "repository and may not be published; a copy needs an approved ",
                "MoTrPAC consortium data-access request."))
  quit(status = 1, save = "no")
}
record("PASS", "celltype_da:cibersortx",
       paste0("deconvolution result: ", basename(cibersortx_csv)))

# ---- the covariate table ---------------------------------------------------

#' Rebuild the cell-type covariate table for one tissue.
#'
#' This is `cibersort_with_metadata.csv`, produced rather than supplied. The
#' chain is the legacy one: drop the three CIBERSORTx run-quality columns, collapse
#' the 22 LM22 subsets into coarse populations and rescale to percentages
#' (squish_proportions()), then join to the sample metadata of the tissue x ome
#' being fitted.
#'
#' The legacy joined with `metadata %>% right_join(cell_types)`, which was safe
#' only because the CSV it read had already been inner-joined against the same
#' metadata one chunk earlier. Reading the raw CIBERSORTx result instead, a
#' right join would carry every deconvolved mixture that is not in the fitted
#' sample set through as a row of NA covariates and hand those to the model. The
#' join is an inner one here and what it drops is counted.
#'
#' @param metadata Sample metadata for the tissue x ome, one row per sample.
#' @param tissue Tissue name, for the messages.
#' @return The metadata with one numeric column per entry of
#'   CELLTYPE_DA_COVARIATES, restricted to samples the deconvolution covers.
build_celltype_covariates <- function(metadata, tissue) {
  cibersortx <- read.csv(cibersortx_csv, check.names = FALSE)

  quality_cols <- c("P-value", "Correlation", "RMSE")
  missing_quality <- setdiff(quality_cols, colnames(cibersortx))
  if (length(missing_quality) > 0) {
    stop(basename(cibersortx_csv), " is missing the CIBERSORTx run-quality ",
         "column(s) ", paste(missing_quality, collapse = ", "),
         " — this does not look like a CIBERSORTx Fractions result",
         call. = FALSE)
  }

  collapsed <- cibersortx %>%
    dplyr::select(-dplyr::all_of(quality_cols)) %>%
    squish_proportions() %>%
    dplyr::rename(vialLabel = Mixture) %>%
    dplyr::mutate(vialLabel = as.character(vialLabel))

  # squish_proportions() returns the population names with spaces ("T cells").
  # The formula needs syntactic names, and the legacy got them for free by
  # round-tripping the table through read.csv()'s default name repair. That
  # round trip is gone, so the rename is done explicitly.
  names(collapsed) <- make.names(names(collapsed))

  absent <- setdiff(CELLTYPE_DA_COVARIATES, colnames(collapsed))
  if (length(absent) > 0) {
    stop("the deconvolution does not carry ", paste(absent, collapse = ", "),
         " — squish_proportions() drops a population averaging 1% or less across ",
         "samples, so a CIBERSORTx run whose estimates moved can silently lose a ",
         "covariate the model expects. Populations present: ",
         paste(setdiff(colnames(collapsed), "vialLabel"), collapse = ", "),
         call. = FALSE)
  }

  collapsed <- collapsed %>%
    dplyr::select(vialLabel, dplyr::all_of(CELLTYPE_DA_COVARIATES))

  fitted_samples <- as.character(metadata$vialLabel)
  covered <- intersect(fitted_samples, collapsed$vialLabel)
  dropped <- setdiff(fitted_samples, collapsed$vialLabel)

  if (length(covered) == 0) {
    stop("no sample of ", tissue, " / ", CELLTYPE_DA_OME, " is in ",
         basename(cibersortx_csv), " — the deconvolution and the qc-norm object ",
         "do not describe the same samples", call. = FALSE)
  }
  if (length(dropped) > 0) {
    # Not silent. A sample the deconvolution does not cover leaves the fit, so
    # this fit and the released one are over different sample sets — which is a
    # difference ED3E would otherwise attribute to the covariates.
    message(sprintf("[celltype-da] %s: %d of %d fitted samples are not in %s and leave the fit",
                    tissue, length(dropped), length(fitted_samples),
                    basename(cibersortx_csv)))
    record("WARN", paste0("celltype_da:", tissue, ":coverage"),
           paste0(length(dropped), " of ", length(fitted_samples),
                  " sample(s) absent from ", basename(cibersortx_csv),
                  " and dropped from the fit"))
  }

  full <- metadata %>%
    dplyr::mutate(vialLabel = as.character(vialLabel)) %>%
    dplyr::inner_join(collapsed, by = "vialLabel")

  # z-scored before they reach the design, as the legacy did. process_covariates()
  # scales numerical covariates again, which is a no-op on an already-standardised
  # column; doing it here as well keeps the written covariate table and the fitted
  # design the same numbers.
  full <- full %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(CELLTYPE_DA_COVARIATES),
                                ~ as.numeric(scale(.))))

  rownames(full) <- full$vialLabel
  full
}

# ---- inputs ----------------------------------------------------------------

#' The gs:// path of the raw RSEM gene counts for one tissue.
#'
#' Resolved by listing the bucket rather than by composing a path, the way the
#' legacy hardcoded one. More than one match is an error rather than a first-hit,
#' because picking silently between two versions of the counts is how a table ends
#' up fitted from data nobody can name.
raw_counts_path <- function(tissue) {
  tissue_code <- celltype_da_tissue_code(CELLTYPE_DA_OME, tissue)
  pattern <- file.path(quant_bucket, "*", "transcriptomics", tissue_code,
                       CELLTYPE_DA_OME, "*rsem-genes-count*.txt")
  hits <- suppressWarnings(
    system2(gsutil, c("ls", shQuote(pattern)), stdout = TRUE, stderr = FALSE)
  )
  hits <- hits[nzchar(hits)]
  if (length(hits) == 0) {
    stop("no rsem-genes-count file under ", pattern,
         " — the fit needs raw counts, which neither data package carries",
         call. = FALSE)
  }
  if (length(hits) > 1) {
    stop("several rsem-genes-count files match ", pattern, ": ",
         paste(basename(hits), collapse = ", "),
         " — pin one before fitting", call. = FALSE)
  }
  hits[[1]]
}

# ---- the fit ---------------------------------------------------------------

#' Fit and write one tissue's cell-type sensitivity table.
#'
#' RNA-seq is modelled on raw counts, not on the qc-norm matrix: the qc-norm
#' object supplies which genes and which samples survived the log-CPM cutoffs,
#' an edgeR DGEList is built from the raw counts for exactly that selection, and
#' voom precision weights are fitted inside the model. That is the upstream
#' arrangement and changing it would change every logFC.
#'
#' @param tissue One of the tissues the deconvolution covers. In practice blood.
#' @return The written table, invisibly.
fit_one_tissue <- function(tissue) {
  counts_path <- raw_counts_path(tissue)

  raw_counts <- MotrpacBicQC::dl_read_gcp(counts_path, sep = "\t",
                                          tmpdir = cache_dir,
                                          gsutil_path = gsutil,
                                          check_first = TRUE)

  # load_acute_only = TRUE, which is the released acute model's sample set and
  # what the legacy fit used. analysis/01_sex_da.R reaches the same restriction
  # differently: it loads both bouts and filters the metadata to ADU_BAS, because
  # its gene set is decided over the full sample set.
  qc <- MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = CELLTYPE_DA_OME,
    selected_tissues = tissue,
    load_acute_only = TRUE,
    verbose = FALSE
  )
  if (length(qc) == 0 || is.null(qc[[tissue]][[CELLTYPE_DA_OME]])) {
    stop("no qc-norm object for ", tissue, " / ", CELLTYPE_DA_OME, call. = FALSE)
  }

  metadata <- qc[[tissue]][[CELLTYPE_DA_OME]][["sample_metadata"]]
  full_metadata <- build_celltype_covariates(metadata, tissue)

  covariate_tsv <- file.path(out_dir, paste0("celltype_covariates_", tissue, ".tsv"))
  write.table(full_metadata, file = covariate_tsv, sep = "\t", quote = FALSE,
              row.names = FALSE)

  raw_counts <- raw_counts %>%
    dplyr::filter(gene_id %in% rownames(qc[[tissue]][[CELLTYPE_DA_OME]][["qc_norm"]])) %>%
    tibble::column_to_rownames("gene_id") %>%
    dplyr::select(dplyr::all_of(as.character(full_metadata$vialLabel)))

  dge <- edgeR::DGEList(counts = raw_counts)
  dge <- edgeR::calcNormFactors(dge)

  # The five cell fractions are appended to the installed covariate table as
  # Technical numerical terms for this ome x tissue, which is how they enter the
  # formula. "Technical" rather than "Design" matters: process_covariates() uses
  # that column to decide which covariates go into design_covs, and the released
  # acute model puts everything that is not a design variable there.
  added_covariates <- MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE %>%
    as.data.frame()
  for (cell_type in CELLTYPE_DA_COVARIATES) {
    added_covariates <- rbind(
      added_covariates,
      data.frame(ome = CELLTYPE_DA_OME, tech_or_design = "Technical",
                 data_type = "numerical", covariate = cell_type, tissue = tissue)
    )
  }

  process_metadata <- process_covariates(meta = full_metadata,
                                         selected_ome = CELLTYPE_DA_OME,
                                         tissue_input = tissue,
                                         custom_covariates = added_covariates)

  fit <- run_dream_acute(expression_object = dge,
                         process_metadata = process_metadata,
                         voom = TRUE)

  out <- .convert_dream_output(fit,
                               formula = process_metadata[["full_formula"]],
                               tissue = tissue,
                               ome = CELLTYPE_DA_OME)
  out$tissue <- tissue

  path <- file.path(out_dir, celltype_da_file_name(CELLTYPE_DA_OME, tissue, version))
  utils::write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE)

  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    file = basename(path),
    tissue = tissue,
    assay = CELLTYPE_DA_OME,
    model = CELLTYPE_DA_MODEL,
    formula = process_metadata[["full_formula"]],
    cell_type_covariates = paste(CELLTYPE_DA_COVARIATES, collapse = ";"),
    n_features = dplyr::n_distinct(out$feature_id),
    n_contrasts = dplyr::n_distinct(out$contrast),
    n_samples = ncol(dge),
    ebayes_legacy = .ebayes_legacy(),
    counts_source = counts_path,
    cibersortx_source = basename(cibersortx_csv),
    covariate_table = basename(covariate_tsv),
    version = version,
    data_pkg_version = as.character(utils::packageVersion("MotrpacHumanPreSuspensionData")),
    analysis_pkg_version = as.character(utils::packageVersion("MotrpacHumanPreSuspensionAnalysis")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    fitted_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  message(sprintf("[celltype-da] %s: %d features x %d contrasts over %d samples -> %s (%.0f MB)",
                  tissue, dplyr::n_distinct(out$feature_id),
                  dplyr::n_distinct(out$contrast), ncol(dge), basename(path),
                  file.size(path) / 1024^2))
  invisible(out)
}

# ---- drive -----------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[celltype-da] fitting into %s — tissues: %s, version: v%s",
                out_dir, paste(requested_tissues, collapse = ","), version))
message(sprintf("[celltype-da]   cores=%s  deconvolution=%s",
                Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", unset = "1"),
                basename(cibersortx_csv)))
message("[celltype-da]   a table that is not already in the manifest takes hours to fit")

failed <- 0L
for (tissue in requested_tissues) {
  path <- file.path(out_dir, celltype_da_file_name(CELLTYPE_DA_OME, tissue, version))

  # A complete table is reused. The fit is hours and the table is hundreds of
  # megabytes, so re-running has to be cheap or nobody will re-run. Completeness
  # is judged on the manifest, not on the file existing: a table killed halfway
  # through its write is on disk and is not a table.
  if (!force && file.exists(path) &&
        .celltype_da_table_is_complete(manifest_tsv, basename(path))) {
    record("SKIP", paste0("celltype_da:", tissue),
           paste0(basename(path), " already fitted — CELLTYPE_DA_FORCE=TRUE to refit"))
    next
  }

  message(sprintf("[celltype-da] fitting %s / %s — this is not quick",
                  tissue, CELLTYPE_DA_OME))
  result <- tryCatch(fit_one_tissue(tissue), error = function(e) {
    record("FAIL", paste0("celltype_da:", tissue), conditionMessage(e))
    NULL
  })
  if (is.null(result)) {
    failed <- failed + 1L
    next
  }
  record("PASS", paste0("celltype_da:", tissue),
         paste0(basename(path), " (", round(file.size(path) / 1024^2), " MB)"))
}

# ---- manifest --------------------------------------------------------------
# One row per table, replaced in place on a refit rather than appended, so the
# manifest always describes what is on disk now.
if (length(manifest_rows) > 0) {
  rows <- do.call(rbind, manifest_rows)
  if (file.exists(manifest_tsv)) {
    prior <- read.csv(manifest_tsv, sep = "\t", check.names = FALSE)
    prior <- prior[!prior$file %in% rows$file, , drop = FALSE]
    rows <- rbind(prior, rows)
  }
  rows <- rows[order(rows$tissue), , drop = FALSE]
  dir.create(dirname(manifest_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(rows, file = manifest_tsv, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  record("PASS", "celltype_da:manifest",
         paste0(nrow(rows), " table(s) recorded in ", basename(manifest_tsv)))
}

quit(status = if (failed > 0) 1 else 0, save = "no")
