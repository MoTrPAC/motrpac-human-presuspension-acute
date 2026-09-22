#!/usr/bin/env Rscript
# Baseline clinical x omics associations — the fit FIG3EFG, ED4G, ST3g, ST3h
# and EDT1 read.
#
# One linear model per feature per trait per tissue per ome:
#
#   feature ~ trait + BMI + Sex + age  (+ TOTAL_LEAN, optionally)
#
# For the Lactate trait, randomGroupCode joins the covariates: peak lactate is
# a post-exercise quantity and the two bouts produce different amounts of it.
#
# A feature is reported only if its fit has at least CLINICAL_OMICS_MIN_RESIDUAL_DF
# residual degrees of freedom.
#
# Site is not a covariate. The legacy fitted it as a ten-level factor, which
# in the small omes (adipose proteomics, n = 22 over 7 sites; blood Olink,
# n = 43 over 5) spent most of the residual degrees of freedom on site.
#
# Each model is fitted on the pre-exercise samples only, and reported as the
# trait coefficient. The output is one long table, and ED4G, FIG3EFG, ST3g,
# ST3h and EDT1 are five readings of it.
#
# A fit rather than panel code because the provenance of a figure's input is
# part of the provenance of the figure, and because five scripts would otherwise
# refit it five ways. The legacy wrote it to
# `clinical_by_omics_res_v4_Aug2026.RData` in the analyst's working directory
# and load()ed that name back in a second script; neither file is in any repo,
# so the dumbbell panels could not be rebuilt from the code that produced them.
#
# Needs consortium data access: it reads individual-level QC matrices and the
# clinical tables out of MotrpacHumanPreSuspensionData. Tens of minutes over the
# 108-combination grid. A table already recorded in the manifest is reused
# rather than refitted.
#
# Parameters: figure_3/clinical_omics.env, plus CLINICAL_OMICS_TABLE and
# CLINICAL_OMICS_TRAITS from config/landscape.env.
#
#   Rscript figures/landscape/analysis/03_clinical_omics.R
#   CLINICAL_OMICS_FORCE=TRUE Rscript figures/landscape/analysis/03_clinical_omics.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(MotrpacHumanPreSuspensionData)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R for landscape_root(), and for the config/landscape.env read it
# performs on load — CLINICAL_OMICS_TABLE and CLINICAL_OMICS_TRAITS come from
# there.
source(file.path(here, "..", "lib", "panel_export.R"))

# ---- figure_3/clinical_omics.env ---------------------------------------------

# The same line parser load_landscape_env() uses, over this fit's own file.
# Values already exported WIN, which is how a one-off override works. Anything
# still carrying a substitution after ${LANDSCAPE_ROOT} is resolved cannot be
# evaluated here, and setting the literal is worse than leaving the variable
# unset for the code default to handle.
load_clinical_omics_env <- function() {
  root <- landscape_root()
  path <- file.path(root, "figure_3", "clinical_omics.env")
  if (!file.exists(path)) {
    stop("figure_3/clinical_omics.env is missing from ", root, call. = FALSE)
  }
  for (line in readLines(path, warn = FALSE)) {
    line <- trimws(line)
    if (!startsWith(line, "export CLINICAL_OMICS_")) next
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
    if (grepl("[$`]", value)) next
    do.call(Sys.setenv, stats::setNames(list(value), name))
  }
  invisible(NULL)
}

load_clinical_omics_env()

# Progress and outcome go to the console. These rows were once written to a TSV
# that a shell stage folded into its own accounting; there is no stage here, and
# a run's report is what the reader sees on screen.
record <- function(status, check, detail) {
  message(sprintf("[%s] %s — %s", status, check, detail))
}

env_list <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) stop(name, " is unset", call. = FALSE)
  trimws(strsplit(value, ",", fixed = TRUE)[[1]])
}

TISSUES <- env_list("CLINICAL_OMICS_TISSUES")
OMES <- env_list("CLINICAL_OMICS_OMES")
TRAITS <- env_list("CLINICAL_OMICS_TRAITS")
ADJUST_LEAN <- toupper(Sys.getenv("CLINICAL_OMICS_ADJUST_LEAN", "FALSE")) == "TRUE"
MIN_RESIDUAL_DF <- as.integer(Sys.getenv("CLINICAL_OMICS_MIN_RESIDUAL_DF", "5"))
OUT_TABLE <- Sys.getenv("CLINICAL_OMICS_TABLE")
MANIFEST <- Sys.getenv("CLINICAL_OMICS_MANIFEST_TSV")
FORCE <- toupper(Sys.getenv("CLINICAL_OMICS_FORCE", "FALSE")) == "TRUE"

if (!nzchar(OUT_TABLE)) {
  stop("CLINICAL_OMICS_TABLE is unset; config/landscape.env defaults it",
       call. = FALSE)
}
if (!nzchar(MANIFEST)) {
  stop("CLINICAL_OMICS_MANIFEST_TSV is unset; figure_3/clinical_omics.env ",
       "defaults it", call. = FALSE)
}

# ---- is there anything to fit? ----------------------------------------------

# FIG3_CLINICAL_OMICS_TABLE names a fit produced elsewhere, which the panels and
# tables read in place of this one. Nothing is fitted here when it is set.
override <- Sys.getenv("FIG3_CLINICAL_OMICS_TABLE", unset = "")
if (nzchar(override)) {
  record("SKIP", "clinical_omics:override",
         paste0("FIG3_CLINICAL_OMICS_TABLE names a fit produced elsewhere (",
                override, ") — nothing is fitted here"))
  quit(save = "no", status = 0)
}

# Which traits come from the screening table and which from blood chemistry.
# Stated rather than inferred from the name, which is what the legacy did with
# an `if (trait %in% c(...))` listing seven of the nine and defaulting the rest.
TRAIT_SOURCE <- c(
  VO2max_L               = "clinical",
  cp_rest_corrected_rmat = "clinical",
  lp_rest_corrected_rmat = "clinical",
  le_rest_corrected_rmat = "clinical",
  peaktorq_iske          = "clinical",
  peak_handgrip          = "clinical",
  HOMA_IR                = "blood_analytes",
  Lactate                = "blood_analytes",
  NEFA                   = "blood_analytes"
)

unknown <- setdiff(TRAITS, names(TRAIT_SOURCE))
if (length(unknown) > 0) {
  stop("trait with no source: ", paste(unknown, collapse = ", "),
       "\n  Add it to TRAIT_SOURCE in analysis/03_clinical_omics.R.",
       call. = FALSE)
}

# ---- reuse ------------------------------------------------------------------

if (!FORCE && file.exists(OUT_TABLE) && file.exists(MANIFEST)) {
  record("SKIP", "clinical_omics:reuse",
         paste0("already fitted: ", basename(OUT_TABLE),
                " — set CLINICAL_OMICS_FORCE=TRUE to refit"))
  quit(save = "no", status = 0)
}

message(sprintf("Fitting into %s — version: v%s", dirname(OUT_TABLE),
                Sys.getenv("CLINICAL_OMICS_VERSION")))
message("  tissues=", paste(TISSUES, collapse = ","))
message("  omes=", paste(OMES, collapse = ","))
message("  traits=", paste(TRAITS, collapse = ","))
message(sprintf("  adjust_lean=%s  min_residual_df=%d", ADJUST_LEAN,
                MIN_RESIDUAL_DF))

# ---- the clinical side, loaded once -----------------------------------------

# The legacy called load_pheno() three times per combination — 324 times over
# the grid — and reloaded the clinical tables with it. They do not vary by
# tissue or ome, so they are loaded once here.
pheno <- MotrpacHumanPreSuspensionData::load_pheno()$data

baseline_covariates <- pheno |>
  filter(.data$visitcode == "ADU_BAS", .data$Timepoint == "pre_exercise") |>
  distinct(.data$pid, .keep_all = TRUE) |>
  transmute(
    pid = as.character(.data$pid),
    BMI = .data$BMI,
    Sex = factor(.data$Sex),
    calculatedAge = .data$calculatedAge
  )

screening <- MotrpacHumanPreSuspensionData::load_clinical_data(
)$curated$cln_curated_screening$data |>
  filter(.data$visit_code == "ADU_SCP") |>
  mutate(pid = as.character(.data$pid))

# One row per participant. The legacy joined these tables without checking
# cardinality and then took distinct(pid, .keep_all = TRUE) afterwards, which
# silently keeps whichever repeat measurement sorted first.
assert_one_row_per_pid <- function(df, what) {
  dupes <- df$pid[duplicated(df$pid)]
  if (length(dupes) > 0) {
    stop(what, ": ", length(unique(dupes)), " participant(s) have more than one ",
         "row, e.g. ", paste(utils::head(sort(unique(dupes)), 3), collapse = ", "),
         "\n  The join would fan out and an arbitrary value would enter the model.",
         call. = FALSE)
  }
  df
}

screening <- assert_one_row_per_pid(screening, "screening table")

dxa <- MotrpacHumanPreSuspensionData::cln_raw_dxa_analysis$data |>
  filter(.data$visitcode == "ADU_SCP") |>
  transmute(pid = as.character(.data$pid), TOTAL_LEAN = .data$TOTAL_LEAN) |>
  assert_one_row_per_pid("DXA table")

# Blood chemistry, wide. Lactate is taken at its peak rather than at baseline —
# post 10 min for both exercise groups — because a resting lactate is not the
# quantity the figure is about.
chemistry_wide <- bind_rows(
  MotrpacHumanPreSuspensionData::load_clinical_data()$chemistry
) |>
  column_to_rownames(var = "analyte_name") |>
  t() |>
  as.data.frame(check.names = FALSE) |>
  rownames_to_column(var = "vialLabel") |>
  inner_join(
    pheno |>
      transmute(vialLabel = .data$vialLabel, pid = as.character(.data$pid),
                randomGroupCode = .data$randomGroupCode,
                Timepoint = .data$Timepoint),
    by = "vialLabel"
  )

resting_chemistry <- chemistry_wide |>
  filter(.data$Timepoint == "pre_exercise") |>
  mutate(
    Insulin = as.numeric(.data$Insulin),
    Glucose = as.numeric(.data$Glucose),
    HOMA_IR = (.data$Insulin * .data$Glucose) / (405 * 34.8),
    NEFA = as.numeric(.data$NEFA)
  ) |>
  select(-"vialLabel", -"randomGroupCode", -"Timepoint") |>
  distinct(.data$pid, .keep_all = TRUE) |>
  assert_one_row_per_pid("resting blood chemistry")

# randomGroupCode travels with it: it is a covariate of the Lactate fit.
peak_lactate <- chemistry_wide |>
  filter(.data$randomGroupCode %in% c("ADUEndur", "ADUResist"),
         .data$Timepoint == "post_10_min") |>
  transmute(pid = .data$pid, Lactate = as.numeric(.data$Lactate),
            randomGroupCode = factor(.data$randomGroupCode)) |>
  distinct(.data$pid, .keep_all = TRUE) |>
  assert_one_row_per_pid("peak lactate")

# Covariates a trait brings with it, beyond BMI, Sex and age.
TRAIT_COVARIATES <- list(Lactate = "randomGroupCode")

trait_values <- function(trait) {
  if (identical(unname(TRAIT_SOURCE[trait]), "clinical")) {
    source_df <- screening
  } else if (identical(trait, "Lactate")) {
    source_df <- peak_lactate
  } else {
    source_df <- resting_chemistry
  }
  if (!trait %in% names(source_df)) {
    stop("trait '", trait, "' is not a column of its source table", call. = FALSE)
  }
  source_df[, c("pid", trait, TRAIT_COVARIATES[[trait]])]
}

# ---- the omics side ---------------------------------------------------------

# Baseline abundances, one row per participant, one column per feature.
baseline_matrix <- function(tissue, ome) {
  if (identical(ome, "metab")) {
    qc <- MotrpacHumanPreSuspensionData::load_qc(selected_tissues = tissue,
                                                 selected_omes = ome)
    matrix_df <- MotrpacHumanPreSuspensionData::combine_qc_matrixes(
      qc, scale = FALSE, make_metab_rownames = TRUE
    )
    rownames(matrix_df) <- gsub(paste0("..", ome, "..", tissue), "",
                                rownames(matrix_df))
    # Metabolomics columns are "<pid>..<timepoint>"; the other omes are keyed by
    # vial label and carry their own sample metadata.
    keep <- grepl("\\.\\.pre_exercise$", colnames(matrix_df))
    baseline <- matrix_df[, keep, drop = FALSE]
    colnames(baseline) <- sub("\\.\\.pre_exercise$", "", colnames(baseline))
    sample_pid <- colnames(baseline)
  } else {
    qc <- MotrpacHumanPreSuspensionData::load_qc(selected_tissues = tissue,
                                                 selected_omes = ome)[[tissue]][[ome]]
    metadata <- qc$sample_metadata |>
      filter(.data$Timepoint == "pre_exercise")
    keep <- colnames(qc$qc_norm) %in% as.character(metadata$vialLabel)
    baseline <- qc$qc_norm[, keep, drop = FALSE]
    sample_pid <- as.character(
      metadata$pid[match(colnames(baseline), as.character(metadata$vialLabel))]
    )
  }

  if (ncol(baseline) == 0) return(NULL)

  transposed <- as.data.frame(t(baseline), check.names = FALSE)
  transposed$pid <- sample_pid
  # A participant with two baseline vials would otherwise contribute two rows to
  # a model that treats each row as one person.
  transposed <- transposed[!duplicated(transposed$pid), , drop = FALSE]
  transposed
}

# ---- the fit ----------------------------------------------------------------

# expression_df is passed in rather than loaded here: it depends on tissue and
# ome but not on the trait, and loading it per trait meant 108 load_qc passes
# and 108 transposes of a matrix that only changes 12 times.
fit_combination <- function(expression_df, tissue, ome, trait) {
  if (is.null(expression_df)) return(NULL)

  # A feature named after the trait would be renamed <trait>.x / <trait>.y by
  # the join below, and the legacy's select(all_of(trait)) then failed for the
  # whole combination inside a tryCatch that reported it as "skipping failed
  # combo". Blood metabolomics carries a feature named "NEFA"; no matrix
  # carries one named "Lactate". The feature is dropped from this combination's
  # input, and the rest are fitted.
  collisions <- intersect(names(expression_df), trait)
  if (length(collisions) > 0) {
    record("WARN", "clinical_omics:collision",
           paste0(tissue, " | ", ome, " | ", trait, " — dropped feature(s) named ",
                  "after the trait from this combination's input: ",
                  paste(collisions, collapse = ", ")))
    expression_df <- expression_df[, setdiff(names(expression_df), collisions),
                                   drop = FALSE]
  }

  covariates <- baseline_covariates
  if (ADJUST_LEAN) covariates <- inner_join(covariates, dxa, by = "pid")

  model_df <- expression_df |>
    inner_join(covariates, by = "pid") |>
    inner_join(trait_values(trait), by = "pid")

  trait_covariates <- TRAIT_COVARIATES[[trait]]
  feature_cols <- setdiff(
    names(expression_df),
    c("pid", "BMI", "Sex", "calculatedAge", "TOTAL_LEAN", trait, trait_covariates)
  )
  rhs <- c("trait_value", "BMI", "Sex", "calculatedAge",
           if (ADJUST_LEAN) "TOTAL_LEAN", trait_covariates)

  results <- lapply(feature_cols, function(feature) {
    frame <- data.frame(
      expression = model_df[[feature]],
      trait_value = model_df[[trait]],
      BMI = model_df$BMI,
      Sex = model_df$Sex,
      calculatedAge = model_df$calculatedAge,
      TOTAL_LEAN = if (ADJUST_LEAN) model_df$TOTAL_LEAN else NA_real_
    )
    for (covariate in trait_covariates) frame[[covariate]] <- model_df[[covariate]]
    frame <- droplevels(stats::na.omit(
      frame[, c("expression", rhs), drop = FALSE]
    ))

    n_obs <- nrow(frame)
    empty <- tibble(feature_id = feature, beta = NA_real_, se = NA_real_,
                    t_stat = NA_real_, p_value = NA_real_, n_obs = n_obs)

    if (n_obs == 0) return(empty)
    if (dplyr::n_distinct(frame$trait_value) < 2) return(empty)
    if (dplyr::n_distinct(frame$Sex) < 2) return(empty)
    for (covariate in trait_covariates) {
      if (dplyr::n_distinct(frame[[covariate]]) < 2) return(empty)
    }

    fit <- tryCatch(
      stats::lm(stats::as.formula(paste("expression ~", paste(rhs, collapse = " + "))),
                data = frame),
      error = function(e) NULL
    )
    if (is.null(fit)) return(empty)

    # A fit with fewer than MIN_RESIDUAL_DF residual degrees of freedom is
    # reported as missing and stays out of the BH correction. So is a feature
    # whose trait term is aliased and has no row in the coefficient table.
    if (fit$df.residual < MIN_RESIDUAL_DF) return(empty)

    coefficients <- summary(fit)$coefficients
    if (!"trait_value" %in% rownames(coefficients)) return(empty)

    tibble(
      feature_id = feature,
      beta = coefficients["trait_value", "Estimate"],
      se = coefficients["trait_value", "Std. Error"],
      t_stat = coefficients["trait_value", "t value"],
      p_value = coefficients["trait_value", "Pr(>|t|)"],
      n_obs = n_obs
    )
  })

  bind_rows(results) |>
    mutate(adj_p = stats::p.adjust(.data$p_value, method = "BH"))
}

# ---- the grid ---------------------------------------------------------------

dir.create(dirname(OUT_TABLE), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(MANIFEST), recursive = TRUE, showWarnings = FALSE)

fitted <- list()
n_skipped <- 0L
for (tissue in TISSUES) {
  for (ome in OMES) {
    # Loaded once per tissue x ome, then reused across every trait.
    expression_df <- tryCatch(baseline_matrix(tissue, ome), error = function(e) {
      record("WARN", "clinical_omics:matrix",
             paste0(tissue, " | ", ome, " — ", conditionMessage(e)))
      NULL
    })
    if (is.null(expression_df)) {
      n_skipped <- n_skipped + length(TRAITS)
      next
    }
    message(sprintf("  %s | %s: %d features x %d participants",
                    tissue, ome, ncol(expression_df) - 1L, nrow(expression_df)))

    for (trait in TRAITS) {
      label <- paste(trait, tissue, ome, sep = " | ")
      result <- tryCatch(
        fit_combination(expression_df, tissue, ome, trait),
        error = function(e) {
          record("WARN", "clinical_omics:combination",
                 paste0(label, " — ", conditionMessage(e)))
          NULL
        }
      )
      if (is.null(result) || nrow(result) == 0) {
        n_skipped <- n_skipped + 1L
        next
      }
      fitted[[label]] <- result |>
        mutate(clinical_trait = trait, tissue = tissue, omics_platform = ome,
               .before = 1)
    }
  }
}

if (length(fitted) == 0) {
  record("FAIL", "clinical_omics:fit", "no trait x tissue x ome combination fitted")
  quit(save = "no", status = 1)
}

associations <- bind_rows(fitted)

write.table(associations, file = OUT_TABLE, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")

manifest <- data.frame(
  table = basename(OUT_TABLE),
  version = Sys.getenv("CLINICAL_OMICS_VERSION"),
  combinations = length(fitted),
  skipped = n_skipped,
  rows = nrow(associations),
  covariates = paste(c("BMI", "Sex", "calculatedAge",
                       if (ADJUST_LEAN) "TOTAL_LEAN"), collapse = "+"),
  trait_covariates = paste(sprintf("%s:%s", names(TRAIT_COVARIATES),
                                   vapply(TRAIT_COVARIATES, paste, "", collapse = "+")),
                           collapse = ";"),
  adjust_lean = ADJUST_LEAN,
  min_residual_df = MIN_RESIDUAL_DF,
  data_pkg_version = as.character(utils::packageVersion("MotrpacHumanPreSuspensionData")),
  fitted_utc = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
  stringsAsFactors = FALSE
)
write.table(manifest, file = MANIFEST, sep = "\t", row.names = FALSE,
            quote = FALSE, na = "")

record("PASS", "clinical_omics:fit",
       sprintf("%d combination(s) fitted, %d skipped, %d rows -> %s",
               length(fitted), n_skipped, nrow(associations),
               basename(OUT_TABLE)))
