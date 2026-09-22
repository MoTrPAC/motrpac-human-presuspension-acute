#!/usr/bin/env Rscript
# Extended Data Table 1 — Selected molecular features associated with baseline
# clinical traits and acute exercise
#
# The features the manuscript selected from the baseline clinical x omics fit,
# one row each: the baseline association with its clinical trait, and the
# feature's EE and RE response at the timepoint the manuscript reports. The
# selection is config/highlights.json (clinical_omics_features, EDT1); every
# statistic is read from the current fit and differential analysis, so a refit
# changes the numbers and not the rows.
#
# ee_vs_re_diff is YES when the EE-RE contrast at that timepoint is FDR < 0.05,
# NO otherwise, and empty where the EE-RE contrast was not fitted (blood has no
# RE samples during exercise).
#
# Standalone: no panel draws these rows together. FIG3EFG draws VO2peak
# transcripts from the same fit, chosen by |t| rather than by hand.
#
# Needs consortium data access, and the fit, which is not in this repository:
#
#   Rscript figures/landscape/analysis/03_clinical_omics.R
#
#   Rscript figures/landscape/tables/EDT1.R

suppressPackageStartupMessages({
  library(dplyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R for landscape_root(), which highlights.R needs, and for the
# config/landscape.env read it performs on load; table_export.R after it.
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "..", "lib", "highlights.R"))
source(file.path(here, "FIG3_ED4_helpers.R"))

# ---- shared ----------------------------------------------------------------

CONTRASTS <- c("EE-CON", "RE-CON", "EE-RE")

# The timepoint shorthand the manuscript tables use, applied on export; pick()
# below joins on the differential analysis names. Same labels as FIG3C and
# lib/plier_helpers.R.
TIMEPOINT_LABELS <- c(
  pre_exercise      = "Pre",
  during_20_min     = "D20M",
  during_40_min     = "D40M",
  post_10_min       = "P10M",
  post_15_30_45_min = "P15-45M",
  post_3.5_4_hr     = "P3.5/4H",
  post_24_hr        = "P24H"
)

# ---- EDT1 — the selection, its baseline association and its acute response --

edt1 <- function() {
  skip <- skip_if_no_fit("EDT1")
  if (!is.null(skip)) return(skip)
  spec <- table_init("EDT1")

  selected <- highlight_clinical_omics_features("EDT1")

  baseline <- clinical_omics_associations() |>
    select("clinical_trait", "tissue", "omics_platform", "feature_id",
           "beta", "t_stat", "adj_p")

  rows <- selected |>
    left_join(baseline,
              by = c(trait = "clinical_trait", "tissue", assay = "omics_platform",
                     "feature_id"),
              relationship = "one-to-one")

  absent <- is.na(rows$t_stat)
  if (any(absent)) {
    stop("selected feature(s) not in the clinical x omics fit: ",
         paste(sprintf("%s (%s/%s/%s)", rows$gene, rows$trait, rows$tissue,
                       rows$assay)[absent], collapse = ", "),
         "\n  Refit with:  CLINICAL_OMICS_FORCE=TRUE Rscript ",
         "figures/landscape/analysis/03_clinical_omics.R", call. = FALSE)
  }

  not_significant <- is.na(rows$adj_p) | rows$adj_p >= 0.05
  if (any(not_significant)) {
    message(sprintf("        EDT1: %d selected feature(s) not FDR < 0.05 in the current fit: %s",
                    sum(not_significant),
                    paste(sprintf("%s/%s", rows$trait, rows$gene)[not_significant],
                          collapse = ", ")))
  }

  # ---- the acute response ---------------------------------------------------

  # One load per tissue x assay, filtered to the selected features.
  slices <- distinct(rows, .data$tissue, .data$assay)
  responses <- lapply(seq_len(nrow(slices)), function(i) {
    tissue <- slices$tissue[i]
    assay <- slices$assay[i]
    ids <- rows$feature_id[rows$tissue == tissue & rows$assay == assay]
    MotrpacHumanPreSuspensionAnalysis::load_differential_analysis(
      selected_omes = assay,
      selected_tissues = tissue,
      single_matrix = TRUE,
      combine_with_featgene = FALSE,
      epigen = FALSE,
      verbose = FALSE
    ) |>
      as.data.frame() |>
      filter(.data$contrast_type %in% c("exercise_with_controls", "Endur_vs_Resist"),
             .data$contrast_category %in% CONTRASTS,
             .data$feature_id %in% ids) |>
      transmute(tissue = tissue, assay = assay,
                feature_id = as.character(.data$feature_id),
                timepoint = as.character(.data$Timepoint),
                contrast = .data$contrast_category,
                log2fc = .data$logFC, adj_p_value = .data$adj_p_value)
  }) |>
    bind_rows()

  pick <- function(contrast, column) {
    hit <- responses[responses$contrast == contrast, , drop = FALSE]
    key <- paste(hit$tissue, hit$assay, hit$feature_id, hit$timepoint)
    hit[[column]][match(paste(rows$tissue, rows$assay, rows$feature_id,
                              rows$timepoint), key)]
  }

  rows$ee_log2fc <- pick("EE-CON", "log2fc")
  rows$ee_adj_p <- pick("EE-CON", "adj_p_value")
  rows$re_log2fc <- pick("RE-CON", "log2fc")
  rows$re_adj_p <- pick("RE-CON", "adj_p_value")
  ee_re_adj_p <- pick("EE-RE", "adj_p_value")
  rows$ee_vs_re_diff <- ifelse(is.na(ee_re_adj_p), NA_character_,
                               ifelse(ee_re_adj_p < 0.05, "YES", "NO"))

  no_response <- is.na(rows$ee_log2fc) & is.na(rows$re_log2fc)
  if (any(no_response)) {
    stop("no differential analysis row at the reported timepoint for: ",
         paste(sprintf("%s %s/%s @ %s", rows$gene, rows$tissue, rows$assay,
                       rows$timepoint)[no_response], collapse = ", "),
         call. = FALSE)
  }

  unlabelled <- setdiff(unique(rows$timepoint), names(TIMEPOINT_LABELS))
  if (length(unlabelled) > 0) {
    stop("timepoint with no shorthand: ", paste(unlabelled, collapse = ", "),
         "\n  Add it to TIMEPOINT_LABELS.", call. = FALSE)
  }

  # ---- export ---------------------------------------------------------------

  out <- rows |>
    transmute(
      clinical_trait = .data$trait,
      tissue = .data$tissue,
      omics_platform = .data$assay,
      feature_id = .data$feature_id,
      gene_symbol = .data$gene,
      t_stat = .data$t_stat,
      adj_p = .data$adj_p,
      beta = .data$beta,
      ee_log2fc = .data$ee_log2fc,
      ee_adj_p = .data$ee_adj_p,
      re_log2fc = .data$re_log2fc,
      re_adj_p = .data$re_adj_p,
      timepoint = unname(TIMEPOINT_LABELS[.data$timepoint]),
      ee_vs_re_diff = .data$ee_vs_re_diff,
      molecule_info = .data$notes
    ) |>
    as.data.frame()

  export_table(out[, spec$columns, drop = FALSE], "EDT1")
}

run_tables(list(EDT1 = edt1))
