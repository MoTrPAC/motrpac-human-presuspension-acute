#!/usr/bin/env Rscript
# Supplementary Table 5 — the standalone sub-table
#
# Sub-tables:  ST5a  transcription factor enrichment and expression
#
# For every HOMER motif: its enrichment for each tissue's differentially
# expressed genes, and the exercise response of the transcription factor that
# motif maps to, in each ome that measured it. Four measures over one
# tissue x timepoint grid, one row per motif x measure x arm.
#
#   enrichment_q       HOMER q for the motif, in that tissue, arm and timepoint
#   transcript_adj_p   the factor's own transcript, FDR
#   protein_adj_p      the factor's protein — mass spec in muscle and adipose,
#                      Olink in blood — FDR
#   phospho_adj_p      the strongest of the factor's phosphosites, FDR
#
# A response is the exercise-with-controls contrast, and a factor carrying
# several features of an ome is reported by its smallest adj p, through
# tf_min_adj_p_by_motif() — the same function FIG5H reads, so the panel and the
# table cannot report different numbers for the same factor.
#
# Standalone rather than written by a panel: no panel computes this grid. FIG5B
# and FIG5C read the enrichment block, FIG5H the muscle enrichment and phospho
# blocks filtered to 15 motifs, and nothing draws the protein or transcript
# blocks at all. Every number here still comes from helpers/FIG5_ED7.R.
#
# Blank means not measured, not "not significant". The legacy initialised the
# whole table to 1 and left that value wherever a factor had no feature, so a
# reader could not tell an unmeasured factor from one that did not respond.
#
# Needs no consortium data access: the HOMER results and tfproanno.RDS ship
# under figure_5/sources/, and the DA tables are in the Analysis package.
#
#   Rscript figures/landscape/tables/ST5.R          every sub-table here
#   Rscript figures/landscape/tables/ST5.R ST5a     one of them

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
# panel_export.R first, for the panel_source() that resolves the vendored HOMER
# files; table_export.R after it, then the helper the figures read.
source(file.path(here, "..", "lib", "panel_export.R"))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(here, "FIG5_ED7_helpers.R"))

# ---- shared ----------------------------------------------------------------

# The ome each measure reads, per tissue. Blood has no phosphoproteomics and
# adipose measured protein and phosphoprotein at 3.5/4 hr only; both are absences
# in the data, and the grid below is intersected with what each DA table holds
# rather than being asserted here.
MEASURE_ASSAYS <- list(
  transcript_adj_p = c(Muscle = "transcript-rna-seq",
                       Adipose = "transcript-rna-seq",
                       Blood = "transcript-rna-seq"),
  protein_adj_p    = c(Muscle = "prot-pr", Adipose = "prot-pr", Blood = "prot-ol"),
  phospho_adj_p    = c(Muscle = "prot-ph", Adipose = "prot-ph")
)

MEASURE_DA <- list(
  transcript_adj_p = c(Muscle = "MUSCLE_TRNSCRPT_DA",
                       Adipose = "ADIPOSE_TRNSCRPT_DA",
                       Blood = "BLOOD_TRNSCRPT_DA"),
  protein_adj_p    = c(Muscle = "MUSCLE_PROT_PR_DA",
                       Adipose = "ADIPOSE_PROT_PR_DA",
                       Blood = "BLOOD_PROT_OL_DA"),
  phospho_adj_p    = c(Muscle = "MUSCLE_PROT_PH_DA", Adipose = "ADIPOSE_PROT_PH_DA")
)

MEASURE_ORDER <- c("enrichment_q", "transcript_adj_p", "protein_adj_p",
                   "phospho_adj_p")

# ---- ST5a — transcription factor enrichment and expression -----------------

# One sub-table, so the HOMER parse is read straight through tf_enrichment()
# rather than through the memoised accessors in helpers/FIG5_ED7.R, which are
# there for the figures' several panels.

st5a <- function() {
  table_init("ST5a")

  enrichment <- tf_enrichment()
  annotation <- tf_annotation(enrichment$motifs)
  comparisons <- enrichment$comparisons

  #' The comparisons of one tissue that a DA table actually holds.
  #'
  #' Adipose proteomics and phosphoproteomics were run at 3.5/4 hr only. Asking
  #' tf_da_stats() for the other two timepoints is an error, correctly, so the ask
  #' is narrowed to what is there instead of the absence being hardcoded.
  available <- function(da, tissue) {
    rows <- comparisons[comparisons$Tissue == tissue, ]
    da <- as.data.frame(da)
    da <- da[da$contrast_type == "exercise_with_controls", ]
    have <- unique(paste(da$randomGroupCode, da$Timepoint))
    rows[paste(rows$randomGroupCode, rows$Timepoint) %in% have, ]
  }

  #' One measure, as a motif x comparison matrix over the full 22-column grid.
  measure_matrix <- function(measure) {
    out <- matrix(NA_real_, nrow = nrow(annotation), ncol = nrow(comparisons),
                  dimnames = list(rownames(annotation), comparisons$comparison))
    if (identical(measure, "enrichment_q")) {
      # Unfloored. The panels replace HOMER's reported 0 with 1e-5 so a -log10 is
      # finite; a table that did the same would assert a q value HOMER never
      # estimated. 0 here means "below the 1e-4 HOMER reports to".
      out[, ] <- enrichment$q_raw[rownames(out), colnames(out)]
      return(out)
    }

    for (tissue in names(MEASURE_ASSAYS[[measure]])) {
      da <- get(MEASURE_DA[[measure]][[tissue]],
                envir = asNamespace("MotrpacHumanPreSuspensionAnalysis"))
      rows <- available(da, tissue)
      if (nrow(rows) == 0) next
      m <- tf_min_adj_p_by_motif(annotation, MEASURE_ASSAYS[[measure]][[tissue]],
                                 da, rows)
      out[rownames(m), colnames(m)] <- m
    }
    out
  }

  # One row per motif x measure x arm; one column per tissue x timepoint, which
  # is the shape ST1c and ST2c use for the same grid.
  long <- bind_rows(lapply(MEASURE_ORDER, function(measure) {
    m <- measure_matrix(measure)
    data.frame(
      motif = rep(rownames(m), times = ncol(m)),
      comparison = rep(colnames(m), each = nrow(m)),
      measure = measure,
      value = as.vector(m),
      stringsAsFactors = FALSE
    )
  }))

  long <- long |>
    left_join(comparisons[, c("comparison", "Tissue", "randomGroupCode", "Timepoint")],
              by = "comparison") |>
    mutate(column = paste(.data$Tissue, .data$Timepoint, sep = "_"))

  value_columns <- comparisons |>
    distinct(.data$Tissue, .data$Timepoint) |>
    arrange(.data$Tissue, match(.data$Timepoint, names(TF_TIMEPOINT_CODES))) |>
    mutate(column = paste(.data$Tissue, .data$Timepoint, sep = "_")) |>
    pull("column")

  wide <- long |>
    select("motif", "measure", "randomGroupCode", "column", "value") |>
    pivot_wider(names_from = "column", values_from = "value") |>
    as.data.frame()

  absent <- setdiff(value_columns, names(wide))
  for (column in absent) wide[[column]] <- NA_real_

  # A motif with no feature in an ome contributes a row of blanks in every
  # tissue. Dropping those is the difference between a 3,768-row table and one a
  # reader can scan; a motif absent from a measure entirely is absent from that
  # measure's rows.
  values <- as.matrix(wide[, value_columns, drop = FALSE])
  keep <- rowSums(!is.na(values)) > 0
  message(sprintf("        %d of %d motif x measure x arm rows carry a value",
                  sum(keep), nrow(wide)))
  wide <- wide[keep, ]

  wide$gene_symbol <- annotation$gene_symbol[match(wide$motif, rownames(annotation))]

  st5a_table <- wide |>
    mutate(
      measure = factor(.data$measure, levels = MEASURE_ORDER),
      motif_order = match(.data$motif, rownames(annotation))
    ) |>
    arrange(.data$measure, .data$randomGroupCode, .data$motif_order) |>
    mutate(measure = as.character(.data$measure)) |>
    select("motif", "gene_symbol", "measure", "randomGroupCode",
           all_of(value_columns))

  export_table(st5a_table, "ST5a")
}

run_tables(list(ST5a = st5a))
