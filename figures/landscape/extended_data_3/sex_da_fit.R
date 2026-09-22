#!/usr/bin/env Rscript
# analysis/01_sex_da.R — fit the sex-stratified differential analysis ED3A and
# ST2g are drawn from.
#
# One table per tissue x ome, under the freeze naming convention, into
# SEX_DA_DIR. The model is
#
#   ~ 0 + sex_group_timepoint + <covariates> + (1 | pid)
#
# and the contrasts are, per exercise group x post-baseline timepoint: the female
# change from baseline, the same for males, and the difference between those two.
# The control arm is not subtracted, so these are not the released acute DA's
# `exercise_with_controls` quantity — which is the set ED3A joins against for its
# significance filter.
#
# SEX_DA_TISSUES defaults to blood,muscle,adipose, the tissues ST2g reports; ED3A
# plots blood only. SEX_DA_TISSUES=all fits everything that qualifies. Only the
# initial acute bout is fitted, and that is not configurable.
#
# Configuration comes from extended_data_3/sex_da.env, through the environment. Nothing
# here is a literal that a run might want to change.
#
# Needs consortium data access, and gsutil authenticated against QUANT_BUCKET:
# transcriptomics is modelled on the raw RSEM counts, which neither data package
# carries. The fit is hours and the tables are hundreds of megabytes, so a table
# already recorded in the manifest is reused rather than refitted;
# SEX_DA_FORCE=TRUE overrides that.
#
#   Rscript figures/landscape/analysis/01_sex_da.R
#
# Writes:
#   <SEX_DA_DIR>/human-precovid-sed-adu_<tissue_code>_<ome>_da_dream-sex_differences_v<version>.txt
#   <SEX_DA_DIR>/sex_da_manifest.tsv   one row per table: the formula, the feature
#                                      and contrast counts, the sample count, the
#                                      counts source, the two data-package
#                                      versions and the eBayes estimator. A table
#                                      listed there is reused, not refitted.

suppressPackageStartupMessages({
  library(dplyr)
  # Both data packages have to be ATTACHED, not merely namespace-qualified.
  # load_qc() resolves its lazy-loaded objects with eval(parse(text = "ADIPOSE_TRNSCRPT_QC"))
  # and similar, which reaches the package's lazydata environment only through
  # the search path — so a qualified call fails with *object 'ADIPOSE_TRNSCRPT_QC'
  # not found*. attach_data_packages() in lib/panel_export.R does the same thing
  # for the panels, driven from the manifest; this fit is not a panel and has no
  # manifest entry, so it attaches them here.
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

# One level below the landscape root. panel_export.R is sourced for
# landscape_root() and for config/landscape.env, which it reads on load.
here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "panel_export.R"))

# The written tables carry CI.L and CI.R, which the released DA tables do not:
# precovid-repro forbids those two columns on every DA schema. Here they are
# produced by .add_confidence_interval() below — recomputed per contrast rather
# than taken from topTable, whose version uses the first contrast's degrees of
# freedom for every contrast — and ED3A draws its error bars from them.

# ---- extended_data_3/sex_da.env ------------------------------------------------------

#' Read extended_data_3/sex_da.env into the environment, once per run.
#'
#' Shell syntax, because it is a list of settings and reads like one; parsed
#' here rather than run through a shell so an R session needs no subprocess.
#' Values already exported WIN, which is how a one-off override works:
#'
#'   SEX_DA_TISSUES=blood Rscript figures/landscape/analysis/01_sex_da.R
.sex_da_load_env <- function() {
  root <- landscape_root()
  path <- file.path(root, "extended_data_3", "sex_da.env")
  if (!file.exists(path)) {
    stop("extended_data_3/sex_da.env is missing from ", root, call. = FALSE)
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

.sex_da_load_env()

# ---- the model -------------------------------------------------------------
#
# The engine ED3A's input is fitted with. It is unexported internal code in a
# SOURCE CHECKOUT of MotrpacHumanPreSuspensionAnalysis
# (`data-raw/generate_differential_analysis/`), and as of the current release
# that directory is retired: every file in it is a tombstone pointing at
# precovid-repro, which reimplements the engine for the acute model and never
# fits the sex one. So there is no installed, no released and no runnable copy
# of the code that produced ED3A's input, and it is vendored here.
#
# Provenance, in order:
#   1. MotrpacHumanPreSuspensionAnalysis, data-raw/generate_differential_analysis/
#      generate_differential_modeling_functions.R and generate_DA_inputs.R, at the
#      commit that produced the published tables. Now retired in place.
#   2. precovid-repro, scripts/10_build_data/09_build_da/da_common.R, which
#      vendored the same functions for the acute model. The bodies below match
#      that copy, which is the maintained one; where it documents a correction to
#      the upstream (the topTable keying, the eBayes estimator), the correction is
#      carried here too and marked.
#
# One body departs from both: .generate_sex_contrasts() does not subtract the
# control arm. Those two contrast sets are different quantities, not two
# spellings of one, and the uncontrolled set is what ED3A shows.

#' Construct covariate metadata and the model formulas.
#'
#' Vendored from the upstream `process_covariates()`. The only change from the
#' upstream body is that `custom_covariates` is gone: this fit has exactly one
#' covariate table, the installed `COVARIATES_FILE`, and a second one reaching
#' this function would silently refit the panel's input against a different
#' design. Everything else - which covariates are scaled, which are factors, the
#' three interaction terms, the order the formula string is assembled in - is
#' unchanged, and the order matters: it decides which columns makeContrastsDream()
#' drops.
#'
#' @param meta Sample metadata for the tissue x ome being fitted, one row per
#'   sample, rownames set to vialLabel.
#' @param selected_ome Assay code, e.g. "transcript-rna-seq".
#' @param tissue_input Tissue name, e.g. "blood".
#' @param include_technical Whether technical covariates enter the model.
#'
#' @return A list carrying the processed metadata, the three formulas, and the
#'   original metadata (which .convert_dream_output() does not use but the
#'   upstream signature passes through).
process_covariates <- function(meta,
                               selected_ome,
                               tissue_input,
                               include_technical = TRUE) {
  covariates_return <- list()
  covariates_return[["original_meta"]] <- meta

  input_covariates <- MotrpacHumanPreSuspensionAnalysis::COVARIATES_FILE

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

  # Sex and codedsiteid come out of the sex-stratified formula, and for different
  # reasons. Sex is inside sex_group_timepoint, so leaving it in makes the design
  # rank-deficient. codedsiteid is dropped because within a sex x group x
  # timepoint cell the site is close to collinear with the participant.
  sex_diff_covs <- full_formula[!full_formula %in% c("Sex", "codedsiteid")]
  sex_formula_string <- paste(sex_diff_covs, collapse = " + ")
  formula_string_sex_differences <- paste("~ 0 + sex_group_timepoint + ",
                                          sex_formula_string, "+ (1 | pid)")

  formula_string <- paste(full_formula, collapse = " + ")
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

#' The sex-stratified contrast set.
#'
#' Three contrasts per exercise group x post-baseline timepoint:
#'   the female change from baseline;
#'   the same for males;
#'   the difference of those two, positive where the female change is larger.
#'
#' The changes are within-sex and within-group, and the control arm is NOT
#' subtracted from them. That makes these contrasts a different quantity from the
#' released acute DA's `exercise_with_controls` set, which does subtract it: an
#' effect here carries whatever the control arm also did over the same interval —
#' circadian drift, the effect of the visit itself — rather than the
#' exercise-attributable part alone. It is what ED3A is drawn from.
#'
#' Resistance exercise has no during-exercise sample, so those two cells are
#' skipped rather than requested and dropped.
#'
#' @param metadata The processed metadata matrix.
.generate_sex_contrasts <- function(metadata) {
  contrast_expressions <- c()
  timepoints <- unique(metadata$Timepoint)
  for (tp in timepoints) {
    if (tp == "pre_exercise") next
    for (group in c("ADUEndur", "ADUResist")) {
      if ((tp == "during_20_min" | tp == "during_40_min") && group == "ADUResist") next

      female_change <- paste0(
        "(sex_group_timepointFemale.", group, ".", tp, " - ",
        "sex_group_timepointFemale.", group, ".pre_exercise)"
      )
      male_change <- gsub("Female", "Male", female_change)
      female_change_vs_male_change <- paste0(female_change, " - ", male_change)

      contrast_expressions <- c(contrast_expressions,
                                female_change,
                                male_change,
                                female_change_vs_male_change)
    }
  }
  return(contrast_expressions)
}

#' The empirical-Bayes estimator eBayes uses: limma's moment estimator.
#'
#' Passed explicitly, never left to limma's default. squeezeVar() resolves a NULL
#' `legacy` as identical(min(df), max(df)), and dream gives every feature its own
#' fractional residual df, so that is always FALSE and limma >= 3.62.0 silently
#' takes the newer estimator. The choice moves df.prior/s2.prior, hence every
#' feature's moderated variance, t and p — so it is stated, and it matches the
#' released acute DA these tables are compared against.
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
#' warning into an error and bptry drops that feature from the fit - silently,
#' one feature per worker. MulticoreParam forks, so each worker inherits an
#' already-armed warning cache and never fires it.
.sex_da_bpparam <- function() {
  cores <- suppressWarnings(as.integer(Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", "")))
  if (is.na(cores) || cores < 2) return(NULL)
  BiocParallel::MulticoreParam(cores, progressbar = TRUE)
}

#' Fit the sex-stratified mixed model.
#'
#' The `model_type == "sex_differences"` branch of the upstream `run_dream()`,
#' with the acute and training branches dropped: this script fits one model and a
#' `model_type` argument here could only ever be wrong.
#'
#' @param expression_object A model-ready expression object - for RNA-seq, an
#'   edgeR DGEList of raw counts.
#' @param process_metadata The list returned by process_covariates().
#' @param voom Whether voom precision weights are applied first. TRUE for
#'   RNA-seq.
#'
#' @return The eBayes-moderated dream fit.
run_dream_sex <- function(expression_object,
                          process_metadata,
                          voom = FALSE) {
  param <- .sex_da_bpparam()

  meta_matrix <- process_metadata$metadata
  # Reorder the covariate rows to follow the expression columns. They should
  # already agree; the fit is what breaks silently if they ever do not, so the
  # match is done rather than assumed. Upstream does the same.
  meta_matrix <- meta_matrix[match(colnames(expression_object), rownames(meta_matrix)), ]

  contrast_expressions <- .generate_sex_contrasts(meta_matrix)
  formula <- stats::as.formula(process_metadata$sex_differences_formula)

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
#' which drops the covariate coefficients.
#'
#' Every column comes from topTable EXCEPT three. degrees_of_freedom and logLik
#' are read off the fit by feature name, and **CI.L and CI.R are computed by
#' .add_confidence_interval() below** — topTable is called with confint = FALSE
#' and never supplies them. So a CI.L or CI.R in one of these tables is this
#' repository's own number, per contrast, and not variancePartition's.
#'
#' @param fit A fit from run_dream_sex().
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
    # confint = FALSE. topTable supplies no interval here at all: CI.L and CI.R
    # are added by .add_confidence_interval() a few lines down.
    #
    # topTable's own interval is unusable for a dream fit. It builds it as
    # `se * qt(alpha, df = eb$df.total[top])`, and df.total is a features x
    # contrasts MATRIX while `top` is a linear index, so every contrast's
    # interval comes out computed with the FIRST contrast's degrees of freedom.
    # Still true as of variancePartition 1.41.5.
    #
    # The released DA tables answer this by shipping no interval — precovid-repro
    # forbids CI.L/CI.R on every DA schema. These tables keep them, recomputed
    # correctly, because ED3A draws the error bars on both of its axes from them.
    res <- variancePartition::topTable(fit, coef = contrast, number = Inf,
                                       p.value = 1, confint = FALSE)

    # Keyed by name, not by position. topTable defaults to sort.by = "p", so
    # `res` comes back permuted relative to the fit, while fit$rdf and fit$logLik
    # are per-feature NAMED vectors in the fit's original order - rdf varies per
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
      dplyr::rename(adj_p_value = adj.P.Val)

    res_single <- .add_confidence_interval(res_single, fit, contrast)

    res_single <- res_single %>%
      dplyr::select(assay, feature_id, z.std, logFC, CI.L, CI.R,
                    degrees_of_freedom, logLik, t, AveExpr,
                    p_value, adj_p_value, contrast, full_model)

    res_tissue <- rbind(res_tissue, res_single)
  }
  res_tissue <- res_tissue %>% dplyr::arrange(adj_p_value)
  return(res_tissue)
}

#' Add the 95% confidence interval for one contrast, keyed by feature.
#'
#' **This function is where CI.L and CI.R in the written tables come from.**
#' Nothing else produces them: topTable is called with confint = FALSE.
#'
#' CI.L = logFC - se * qt(0.975, df.total[, contrast])
#' CI.R = logFC + se * qt(0.975, df.total[, contrast])
#'   with se = logFC / t
#'
#' `logFC / t` is the standard error. dream's t is
#' `coefficients / stdev.unscaled / sqrt(s2.post)`, so the ratio returns the
#' MODERATED standard error, the one the t and p_value in this same row are
#' built from. `stdev.unscaled * sigma` is not it: eBayes leaves sigma
#' unmoderated, so that puts the interval and the test on different variance
#' estimates, and an interval can then exclude zero on a feature the test calls
#' null.
#'
#' df.total is read with TWO subscripts. It is features x contrasts for a dream
#' fit, and a single subscript is a linear index into it — the upstream bug this
#' function exists to avoid.
#'
#' Note df.total is not the same as the degrees_of_freedom column, which is
#' fit$rdf. The interval cannot be reconstructed from the written table for that
#' reason, which is why these two columns are written rather than derived later.
#'
#' @param res A per-contrast result frame, feature ids in the rownames.
#' @param fit The dream fit.
#' @param contrast The contrast column being extracted.
.add_confidence_interval <- function(res, fit, contrast) {
  df <- fit$df.total[rownames(res), contrast]
  margin <- (res$logFC / res$t) * stats::qt(0.975, df = df)

  res$CI.L <- res$logFC - margin
  res$CI.R <- res$logFC + margin
  res
}

# ---- naming ----------------------------------------------------------------

# The freeze stem prefix. Every released table under the analysis bucket carries
# it, and ED3A finds its input by matching the rest of the name, so the tables
# this script writes use the same convention.
SEX_DA_STEM_PREFIX <- "human-precovid-sed-adu"

#' The BIC tissue code for a tissue x ome, from the installed OME_TISSUE_CODE.
#'
#' The code is not derivable from the tissue name - blood transcriptomics is
#' t04-blood-rna while blood proteomics is not - so it is looked up rather than
#' composed.
sex_da_tissue_code <- function(ome, tissue) {
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
#' `<prefix>_<tissue_code>_<ome>_da_dream-sex_differences_v<version>.txt`, which
#' is what the upstream write_with_path_name() composes and what ED3A matches on:
#' the panel pins its file set by the model name plus the tissue and assay, so a
#' table written under any other name is invisible to it.
#'
#' @param ome Assay code.
#' @param tissue Tissue name.
#' @param version Version label. Not a bucket release - see extended_data_3/sex_da.env.
sex_da_file_name <- function(ome, tissue, version) {
  paste0(paste(SEX_DA_STEM_PREFIX, sex_da_tissue_code(ome, tissue), ome,
               "da", "dream-sex_differences", sep = "_"),
         "_v", version, ".txt")
}

#' Whether a table on disk is one this script finished writing.
#'
#' Existence is not the test. A fit killed partway through write.table() leaves a
#' file that opens, parses, and is missing however many features the process did
#' not get to - and the run would then skip it forever as "already fitted".
#' A manifest row is written only after the table is closed, so the row is the
#' completion record and the file alone is not.
#'
#' @param manifest_tsv Path to sex_da_manifest.tsv.
#' @param file_name Base name of the table.
.sex_da_table_is_complete <- function(manifest_tsv, file_name) {
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
  if (!nzchar(v)) stop(name, " is not set — it is defined in extended_data_3/sex_da.env",
                       call. = FALSE)
  v
}
env_flag <- function(name, default = FALSE) {
  v <- toupper(Sys.getenv(name, unset = ""))
  if (!nzchar(v)) return(default)
  v %in% c("TRUE", "1", "YES")
}

out_dir      <- env_or_stop("SEX_DA_DIR")
manifest_tsv <- env_or_stop("SEX_DA_MANIFEST_TSV")
version      <- env_or_stop("SEX_DA_VERSION")
quant_bucket <- env_or_stop("QUANT_BUCKET")
gsutil       <- Sys.getenv("GSUTIL", unset = "gsutil")
force        <- env_flag("SEX_DA_FORCE")

split_tissues <- function(name) {
  v <- trimws(strsplit(Sys.getenv(name, unset = ""), ",")[[1]])
  v[nzchar(v)]
}
requested_tissues <- split_tissues("SEX_DA_TISSUES")
# Metabolomics is scoped by its own tissue list rather than sharing
# SEX_DA_TISSUES, because the two modalities are not wanted for the same
# tissues: ST2g reports transcriptomics for all three and metabolomics for
# blood and muscle only. A single list could not say that, and adipose
# metabolomics does qualify on the >= 3 rule, so leaving it to the rule would
# fit twelve platforms nothing reports.
requested_metab_tissues <- split_tissues("SEX_DA_METAB_TISSUES")

# Either list may be empty — fitting only metabolomics, or only transcriptomics,
# is a legitimate scope. Both empty is not, and is almost always an unset
# environment rather than an intent.
if (length(requested_tissues) == 0 && length(requested_metab_tissues) == 0) {
  stop("neither SEX_DA_TISSUES nor SEX_DA_METAB_TISSUES is set — both are ",
       "defined in extended_data_3/sex_da.env", call. = FALSE)
}

# The raw counts are downloaded once and reused. They are large and the download
# is the only part of this script that needs the network.
cache_dir <- file.path(out_dir, "raw-files")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# The two modalities this script fits. The upstream driver
# (MotrpacHumanPreSuspensionAnalysis, data-raw/generate_differential_analysis/
# sex_differential_analysis.R) computed a qualifying set across every ome and then
# skipped metabolomics and proteomics inside the loop with two `next` statements.
# Metabolomics is fitted here because ST2g reports it; proteomics still is not,
# because nothing in this repository reports it.
#
# METAB_SELECTOR is what load_qc() is asked for. It expands to one entry per
# platform, and each platform is fitted separately — the platform IS the ome for
# metabolomics, which is what .convert_dream_output() records as the assay.
TRANSCRIPT_OME <- "transcript-rna-seq"
METAB_SELECTOR <- "metab"

#' TRUE for a metabolomics platform code.
is_metab_ome <- function(ome) grepl("^metab", ome)

# What the stage driver wrote to its report TSV goes to the console instead. The
# status word is kept so a run still reads as a list of checks.
record <- function(status, check, detail) {
  message(sprintf("[sex-da] %-4s %s — %s", status, check, detail))
}
manifest_rows <- list()

# ---- is there anything to fit? ---------------------------------------------
# ED3A reads SEX_STRATIFIED_DA_DIR. Pointed anywhere but this script's output
# directory, it names tables fitted elsewhere and there is nothing to do here.

sex_stratified_dir <- Sys.getenv("SEX_STRATIFIED_DA_DIR", unset = "")
if (nzchar(sex_stratified_dir) && !identical(sex_stratified_dir, out_dir)) {
  record("SKIP", "sex_da:override",
         paste0("SEX_STRATIFIED_DA_DIR points outside this fit (",
                sex_stratified_dir,
                ") — ED3A will read that directory and nothing is fitted here"))
  quit(status = 0, save = "no")
}

# ---- what this fit needs on the machine ------------------------------------
# Transcriptomics is modelled on the raw RSEM counts, which neither data package
# carries; they are read from QUANT_BUCKET.

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
record("PASS", "sex_da:bucket", paste0("readable: ", quant_bucket))

# ---- which tissues qualify -------------------------------------------------

#' Tissues whose sex x group x timepoint cells are all large enough to fit.
#'
#' Vendored from the qualifying-set block of the upstream driver: count samples
#' per group x timepoint x sex, drop controls, and keep a tissue x ome only when
#' no remaining cell holds fewer than three samples.
#'
#' Two things about that rule are worth stating rather than leaving to be
#' rediscovered. Controls are dropped BEFORE the count, which is what lets
#' adipose qualify at all; and the contrasts then reference the control cells
#' anyway, so a tissue can qualify on cells the model does not fit and be fitted
#' on cells nobody counted. That is the published behaviour and it is reproduced,
#' not corrected.
#'
#' The flag is grouped by tissue x assay, which is what the upstream did over
#' every ome. Grouping by tissue alone would be equivalent only while a single
#' ome is loaded, and two are now.
#'
#' @return A data frame of qualifying tissue x assay pairs.
qualifying_tissue_omes <- function() {
  qc <- MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = c(TRANSCRIPT_OME, METAB_SELECTOR),
    selected_tissues = "all",
    epigen = FALSE,
    remove_redundant_metab = TRUE,
    verbose = FALSE
  )
  pheno_data <- MotrpacHumanPreSuspensionData::pheno$data

  blocks <- list()
  for (tissue_name in names(qc)) {
    for (ome_name in names(qc[[tissue_name]])) {
      mat <- qc[[tissue_name]][[ome_name]][["qc_norm"]]
      if (is.null(mat) || nrow(mat) == 0) next
      blocks[[paste(tissue_name, ome_name)]] <- pheno_data %>%
        dplyr::filter(vialLabel %in% colnames(mat)) %>%
        dplyr::group_by(randomGroupCode, Timepoint, Sex) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(tissue = tissue_name, assay = ome_name)
    }
  }
  counts <- dplyr::bind_rows(blocks)

  counts %>%
    dplyr::filter(randomGroupCode != "ADUControl") %>%
    dplyr::group_by(tissue, assay) %>%
    dplyr::mutate(flag_low_n = any(n < 3)) %>%
    dplyr::filter(!flag_low_n) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(tissue, assay) %>%
    as.data.frame()
}

# ---- inputs ----------------------------------------------------------------

#' The gs:// path of the raw RSEM gene counts for one tissue.
#'
#' Transcriptomics only. There is no count matrix behind a metabolomics platform,
#' which is why fit_one() takes the qc-norm matrix for those and never calls this.
#'
#' Resolved by listing the bucket rather than by composing a path, the way the
#' upstream `.find_path_name()` did: the collection version and the file version
#' are both in the path and neither is this repository's to assert. More than one
#' match is an error rather than a first-hit, because picking silently between two
#' versions of the counts is how a table ends up fitted from data nobody can name.
raw_counts_path <- function(tissue) {
  tissue_code <- sex_da_tissue_code(TRANSCRIPT_OME, tissue)
  pattern <- file.path(quant_bucket, "*", "transcriptomics", tissue_code,
                       TRANSCRIPT_OME, "*rsem-genes-count*.txt")
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

#' Fit and write one tissue's sex-stratified table.
#'
#' RNA-seq is modelled on raw counts, not on the qc-norm matrix: the qc-norm
#' object supplies which genes and which samples survived the log-CPM cutoffs,
#' an edgeR DGEList is built from the raw counts for exactly that selection, and
#' voom precision weights are fitted inside the model. That is the upstream
#' arrangement and changing it would change every logFC.
#'
#' Metabolomics is modelled on the qc-norm matrix itself, with no voom step and
#' no raw-count download: there is no count matrix behind a metabolomics platform,
#' and the abundances are already on the scale the model expects. That is the
#' arrangement precovid-repro's generate_metab_da.R uses for the acute model, with
#' `model_type` swapped for the sex-stratified one.
#'
#' @param tissue One of blood, muscle, adipose.
#' @param ome    "transcript-rna-seq", or one metabolomics platform code.
#' @return The written table, invisibly.
fit_one <- function(tissue, ome) {
  metab <- is_metab_ome(ome)

  # load_acute_only = FALSE, and the bout filter is applied to the metadata
  # below instead. For transcriptomics the qc-norm matrix is used for its gene
  # set, which is decided over the full sample set; restricting the load would
  # decide it over a subset.
  qc <- MotrpacHumanPreSuspensionData::load_qc(
    selected_omes = ome,
    selected_tissues = tissue,
    load_acute_only = FALSE,
    verbose = FALSE
  )
  if (length(qc) == 0 || is.null(qc[[tissue]][[ome]])) {
    stop("no qc-norm object for ", tissue, " / ", ome, call. = FALSE)
  }

  metadata <- qc[[tissue]][[ome]][["sample_metadata"]]

  # The initial acute bout only. ADU_PAS is the post-training bout, and it
  # carries the same timepoint labels as the pre-training one, so including it
  # would put both bouts in one sex x group x timepoint cell and give those
  # participants two rows under a single (1 | pid) intercept.
  metadata <- metadata %>% dplyr::filter(visitcode == "ADU_BAS")
  rownames(metadata) <- metadata$vialLabel

  if (metab) {
    counts_path <- NA_character_
    expression_object <- qc[[tissue]][[ome]][["qc_norm"]] %>%
      dplyr::select(dplyr::all_of(as.character(metadata$vialLabel)))
    n_samples <- ncol(expression_object)
  } else {
    counts_path <- raw_counts_path(tissue)
    raw_counts <- MotrpacBicQC::dl_read_gcp(counts_path, sep = "\t",
                                            tmpdir = cache_dir,
                                            gsutil_path = gsutil,
                                            check_first = TRUE)
    raw_counts <- raw_counts %>%
      dplyr::filter(gene_id %in% rownames(qc[[tissue]][[ome]][["qc_norm"]])) %>%
      tibble::column_to_rownames("gene_id") %>%
      dplyr::select(dplyr::all_of(as.character(metadata$vialLabel)))
    expression_object <- edgeR::DGEList(counts = raw_counts)
    expression_object <- edgeR::calcNormFactors(expression_object)
    n_samples <- ncol(expression_object)
  }

  process_metadata <- process_covariates(meta = metadata,
                                         selected_ome = ome,
                                         tissue_input = tissue)

  fit <- run_dream_sex(expression_object = expression_object,
                       process_metadata = process_metadata,
                       voom = !metab)

  out <- .convert_dream_output(fit,
                               formula = process_metadata[["sex_differences_formula"]],
                               tissue = tissue,
                               ome = ome)
  out$tissue <- tissue

  path <- file.path(out_dir, sex_da_file_name(ome, tissue, version))
  utils::write.table(out, file = path, sep = "\t", quote = FALSE, row.names = FALSE)

  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    file = basename(path),
    tissue = tissue,
    assay = ome,
    model = "dream-sex_differences",
    formula = process_metadata[["sex_differences_formula"]],
    n_features = dplyr::n_distinct(out$feature_id),
    n_contrasts = dplyr::n_distinct(out$contrast),
    n_samples = n_samples,
    ebayes_legacy = .ebayes_legacy(),
    # Metabolomics has no raw-count source: the qc-norm matrix in the data
    # package IS the input, so the version of that package is the provenance.
    counts_source = if (metab) {
      paste0("qc_norm from MotrpacHumanPreSuspensionData ",
             utils::packageVersion("MotrpacHumanPreSuspensionData"))
    } else {
      counts_path
    },
    version = version,
    data_pkg_version = as.character(utils::packageVersion("MotrpacHumanPreSuspensionData")),
    analysis_pkg_version = as.character(utils::packageVersion("MotrpacHumanPreSuspensionAnalysis")),
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    fitted_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  )

  message(sprintf("[sex-da] %s / %s: %d features x %d contrasts -> %s (%.1f MB)",
                  tissue, ome, dplyr::n_distinct(out$feature_id),
                  dplyr::n_distinct(out$contrast), basename(path),
                  file.size(path) / 1024^2))
  invisible(out)
}

# ---- drive -----------------------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[sex-da] fitting into %s — tissues: %s, metab tissues: %s, version: v%s",
                out_dir, paste(requested_tissues, collapse = ","),
                paste(requested_metab_tissues, collapse = ","), version))
message(sprintf("[sex-da]   cores=%s",
                Sys.getenv("VARIANCEPARTITION_PARALLEL_CORES", unset = "1")))
message("[sex-da]   a table that is not already in the manifest takes hours to fit")

qualifying <- tryCatch(qualifying_tissue_omes(), error = function(e) {
  record("FAIL", "sex_da:qualifying",
         paste0("cannot decide which tissue x ome pairs qualify: ",
                conditionMessage(e)))
  NULL
})

if (is.null(qualifying)) {
  quit(status = 1, save = "no")
}

record("PASS", "sex_da:qualifying",
       paste0(nrow(qualifying), " tissue x ome pair(s) qualify: ",
              paste(unique(qualifying$tissue), collapse = ", ")))

# ---- what was asked for ----------------------------------------------------
# One tissue list per modality, because the two are not wanted for the same
# tissues. See SEX_DA_METAB_TISSUES in extended_data_3/sex_da.env for why adipose
# metabolomics is not among them.
expand_scope <- function(tissues, ome_filter) {
  if (length(tissues) == 0) return(qualifying[0, , drop = FALSE])
  pool <- qualifying[ome_filter(qualifying$assay), , drop = FALSE]
  if (identical(tolower(tissues), "all")) return(pool)
  pool[pool$tissue %in% tissues, , drop = FALSE]
}

requested <- rbind(
  expand_scope(requested_tissues, function(a) a == TRANSCRIPT_OME),
  expand_scope(requested_metab_tissues, is_metab_ome)
)

# A tissue named for a modality it cannot be fitted for is worth saying out
# loud, rather than silently producing fewer tables than were asked for.
for (tissue in setdiff(requested_tissues, c("all", qualifying$tissue[qualifying$assay == TRANSCRIPT_OME]))) {
  record("WARN", paste0("sex_da:", tissue, ":", TRANSCRIPT_OME),
         "requested but does not qualify — some sex x group x timepoint cell holds fewer than 3 samples")
}
metab_qualifying <- unique(qualifying$tissue[is_metab_ome(qualifying$assay)])
for (tissue in setdiff(requested_metab_tissues, c("all", metab_qualifying))) {
  record("WARN", paste0("sex_da:", tissue, ":metab"),
         "requested but no metabolomics platform qualifies")
}

if (nrow(requested) == 0) {
  record("WARN", "sex_da:scope",
         paste0("nothing to fit. SEX_DA_TISSUES='",
                paste(requested_tissues, collapse = ","),
                "', SEX_DA_METAB_TISSUES='",
                paste(requested_metab_tissues, collapse = ","), "'"))
  quit(status = 0, save = "no")
}

requested <- requested[order(requested$tissue, requested$assay), , drop = FALSE]
record("PASS", "sex_da:scope",
       paste0(nrow(requested), " table(s) in scope: ",
              paste(sprintf("%s/%s", requested$tissue, requested$assay),
                    collapse = ", ")))

# ---- fit -------------------------------------------------------------------
failed <- 0L
for (i in seq_len(nrow(requested))) {
  tissue <- requested$tissue[i]
  ome <- requested$assay[i]
  tag <- paste0("sex_da:", tissue, ":", ome)
  path <- file.path(out_dir, sex_da_file_name(ome, tissue, version))

  # A complete table is reused. A transcriptomics fit is hours and its table is
  # hundreds of megabytes, so re-running has to be cheap or nobody will re-run.
  # Completeness is judged on the manifest, not on the file existing: a table
  # killed halfway through its write is on disk and is not a table.
  if (!force && file.exists(path) && .sex_da_table_is_complete(manifest_tsv, basename(path))) {
    record("SKIP", tag,
           paste0(basename(path), " already fitted — SEX_DA_FORCE=TRUE to refit"))
    next
  }

  if (is_metab_ome(ome)) {
    message(sprintf("[sex-da] fitting %s / %s", tissue, ome))
  } else {
    message(sprintf("[sex-da] fitting %s / %s — this is not quick", tissue, ome))
  }
  result <- tryCatch(fit_one(tissue, ome), error = function(e) {
    record("FAIL", tag, conditionMessage(e))
    NULL
  })
  if (is.null(result)) {
    failed <- failed + 1L
    next
  }
  record("PASS", tag,
         paste0(basename(path), " (", round(file.size(path) / 1024^2, 1), " MB)"))
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
  rows <- rows[order(rows$tissue, rows$assay), , drop = FALSE]
  dir.create(dirname(manifest_tsv), recursive = TRUE, showWarnings = FALSE)
  write.table(rows, file = manifest_tsv, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  record("PASS", "sex_da:manifest",
         paste0(nrow(rows), " table(s) recorded in ", basename(manifest_tsv)))
}

quit(status = if (failed > 0) 1 else 0, save = "no")
