#!/usr/bin/env Rscript
# Supplementary Table 1 — Cohort composition, sampling profile and outliers
#
# Sub-tables: ST1a  participants per randomization group
#             ST1b  participants per group and sex
#             ST1d  number of samples excluded as outliers, per tissue and assay
#
# Not built here, see each sub-table's `script` in config/table_map.json:
#   ST1c  participants per tissue and omics assay   written by FIG1.R (FIG1C)
#   ST1e  named features                            written by ED1.R  (ED1C)
#
# ST1a and ST1b need consortium data access; ST1d does not, it counts the
# OUTLIERS object shipped in MotrpacHumanPreSuspensionAnalysis.
#
#   Rscript figures/landscape/tables/ST1.R          every sub-table
#   Rscript figures/landscape/tables/ST1.R ST1d     one sub-table

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(MotrpacHumanPreSuspensionData)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(landscape_root(), "lib", "supp_table_helpers.R"))

# ---- shared ----------------------------------------------------------------

PROFILE_LEVELS <- c("Early", "Middle", "Late", "All")
GROUP_LABELS <- c(ADUControl = "CON", ADUEndur = "EE", ADUResist = "RE")

# ST1a and ST1b are the same counts with and without the sex split, off one
# participant-level frame. Lazy and memoised rather than loaded at the top of
# the script: pheno is read once for the two of them, and not at all when only
# ST1d is selected.
#
# Unlike every other sub-table here these describe the cohort rather than a
# measurement, so nothing about them moves when the omics data is reprocessed.
#
# The temporal sampling profile is not derived. pheno carries it per participant
# in tempSampProfile: Early, Middle and Late mean sampled ONLY in that
# post-exercise window, All means sampled in all three. The categories are
# mutually exclusive and every participant falls in exactly one, which is
# asserted below rather than assumed.
participants <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)

    # ADU_BAS is the acute visit. pheno carries one row per vial, so the
    # distinct() reduces it to one row per participant.
    df <- MotrpacHumanPreSuspensionData::pheno$data %>%
      dplyr::filter(.data$visitcode == "ADU_BAS") %>%
      dplyr::distinct(.data$pid, .data$tempSampProfile, .data$randomGroupCode,
                      .data$Sex)

    if (nrow(df) != dplyr::n_distinct(df$pid)) {
      stop("ST1a/ST1b: a participant carries more than one profile, group or ",
           "sex at ADU_BAS, so the cross-tabulation would count them twice.",
           call. = FALSE)
    }
    unexpected <- setdiff(unique(df$tempSampProfile), PROFILE_LEVELS)
    if (length(unexpected) > 0) {
      stop("ST1a/ST1b: tempSampProfile values with no row in the table: ",
           paste(sort(unexpected), collapse = ", "), call. = FALSE)
    }

    cache <<- df %>%
      dplyr::mutate(
        profile = factor(.data$tempSampProfile, levels = PROFILE_LEVELS),
        group = GROUP_LABELS[as.character(.data$randomGroupCode)]
      )
    cache
  }
})

#' "n (%)", the percentage taken within the column as the reference does.
#'
#' The denominator is the column total, not the row total, so the four profile
#' rows of one column sum to 100% and the columns do not sum across.
n_pct <- function(n, denom) {
  ifelse(is.na(n) | n == 0, "0 (0%)",
         sprintf("%d (%d%%)", n, as.integer(round(100 * n / denom))))
}

#' Cross-tabulate the profile against one or more grouping columns.
#'
#' A Total row is emitted as data. The reference carries the denominators in the
#' header instead ("Overall\nN = 175"), which a single-header TSV cannot hold,
#' and a percentage without its denominator is not interpretable.
tabulate_profiles <- function(df, col_of) {
  counts <- df %>%
    dplyr::count(.data$profile, .data$column, name = "n") %>%
    tidyr::complete(profile = PROFILE_LEVELS, column = col_of, fill = list(n = 0))
  denom <- counts %>%
    dplyr::group_by(.data$column) %>%
    dplyr::summarise(N = sum(.data$n), .groups = "drop")

  body <- counts %>%
    dplyr::left_join(denom, by = "column") %>%
    dplyr::mutate(value = n_pct(.data$n, .data$N)) %>%
    dplyr::select("profile", "column", "value") %>%
    tidyr::pivot_wider(names_from = "column", values_from = "value") %>%
    dplyr::arrange(match(.data$profile, PROFILE_LEVELS)) %>%
    dplyr::rename(Characteristic = "profile")

  total <- denom %>%
    dplyr::mutate(value = as.character(.data$N)) %>%
    dplyr::select("column", "value") %>%
    tidyr::pivot_wider(names_from = "column", values_from = "value") %>%
    dplyr::mutate(Characteristic = "Total")

  as.data.frame(dplyr::bind_rows(body, total))
}

# ---- ST1a — participants per randomization group ---------------------------

st1a <- function() {
  spec <- table_init("ST1a")

  by_group <- dplyr::bind_rows(
    dplyr::mutate(participants(), column = "Overall"),
    dplyr::mutate(participants(), column = .data$group)
  )
  out <- tabulate_profiles(by_group, c("Overall", unname(GROUP_LABELS)))

  export_table(out[, spec$columns, drop = FALSE], "ST1a")
}

# ---- ST1b — participants per group and sex ---------------------------------

st1b <- function() {
  spec <- table_init("ST1b")

  by_group_sex <- dplyr::bind_rows(
    dplyr::mutate(participants(), column = paste("Overall", .data$Sex, sep = "_")),
    dplyr::mutate(participants(), column = paste(.data$group, .data$Sex, sep = "_"))
  )
  out <- tabulate_profiles(by_group_sex, setdiff(spec$columns, "Characteristic"))

  export_table(out[, spec$columns, drop = FALSE], "ST1b")
}

# ---- ST1d — samples excluded as outliers, per tissue and assay -------------

# Standalone, and cheaply so: the whole table is a count over the OUTLIERS
# object, 160 rows shipped in MotrpacHumanPreSuspensionAnalysis. Nothing here
# loads a QC matrix, which is why it sits here rather than in FIG1.R beside
# ST1c — the two look like siblings and share no computation.
#
# OUTLIERS is the complement of what ST1c counts. A sample measured on a tissue
# x assay either passed QC and is in ST1c's counts, or failed and is in this
# one.
st1d <- function() {
  spec <- table_init("ST1d")

  declared <- spec$columns
  tissue_cols <- setdiff(declared, "assay")

  # One row per excluded sample. distinct() is defensive rather than corrective:
  # OUTLIERS carries one row per sample with a `reason` string, so a sample
  # failing two QC criteria is one row today. If that ever became two, this
  # table would report the sample twice without it.
  counts <- OUTLIERS %>%
    distinct(.data$vialLabel, .data$tissue, .data$ome) %>%
    count(.data$ome, .data$tissue, name = "n")

  unexpected_tissue <- setdiff(unique(counts$tissue), tissue_cols)
  if (length(unexpected_tissue) > 0) {
    stop("ST1d: outliers recorded for tissue(s) the manifest does not declare: ",
         paste(sort(unexpected_tissue), collapse = ", "),
         "\n  Add them to `columns` for ST1d in config/table_map.json.",
         call. = FALSE)
  }

  # Every assay in the study gets a row, not only those with an outlier. A blank
  # cell then means "none excluded", which is a result; an absent row would mean
  # "not asked", which is not the same claim and is the one a reader would have
  # to guess at.
  out <- data.frame(assay = supp_table_assays(), stringsAsFactors = FALSE) %>%
    left_join(counts, by = c("assay" = "ome")) %>%
    pivot_wider(names_from = "tissue", values_from = "n") %>%
    as.data.frame()

  # pivot_wider() emits a column per tissue it actually saw, and adds an "NA"
  # column for the assays that matched nothing in the join.
  out[["NA"]] <- NULL
  for (missing_col in setdiff(tissue_cols, names(out))) {
    out[[missing_col]] <- NA_integer_
  }

  out <- order_by_assay(out, "assay")[, declared]

  export_table(out, "ST1d")
}

run_tables(list(ST1a = st1a, ST1b = st1b, ST1d = st1d))
