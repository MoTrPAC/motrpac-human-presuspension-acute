#!/usr/bin/env Rscript
# Supplementary Table 4 — the standalone sub-table
#
# Sub-tables:  ST4a  ORA of the fuzzy c-means clusters, filtered
#
# Not built here, see each sub-table's `script` in config/table_map.json:
#   ST4b  cross-tissue RNA-seq PLIER     unassigned
#   ST4c  transcript affinity            unassigned
#   ST4d  cross-tissue metabolomics PLIER  unassigned
#   ST4e  metabolite affinity            unassigned
# All four are blocks of the fits analysis/04_plier.R writes; no script owns
# them yet.
#
# ST4a needs no consortium data access: FCM_ORA ships in
# MotrpacHumanPreSuspensionAnalysis.
#
#   Rscript figures/landscape/tables/ST4.R          every sub-table here
#   Rscript figures/landscape/tables/ST4.R ST4a     one of them

suppressPackageStartupMessages({
  library(dplyr)
  library(MotrpacHumanPreSuspensionAnalysis)
})

here <- dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "..", "lib", "table_export.R"))
source(file.path(landscape_root(), "lib", "supp_table_helpers.R"))

# ---- ST4a — ORA of the fuzzy c-means clusters, filtered --------------------

# Standalone. FCM_ORA in MotrpacHumanPreSuspensionAnalysis is the whole ORA over
# every cluster of every tissue x assay; this table is the significant part of
# it. FIG4D, ED5B and ED5C read FCM_ORA through helpers/FIG4_ED5.R at
# FDR < 0.05; no panel draws this nominal-p cut, so there is nothing to
# co-generate with.
#
# What the table reports is decided entirely by which version of the Analysis
# package is installed: the clustering is in the object, not computed here. The
# pin in config/required_packages.tsv is what fixes it, and the export manifest
# records the version every run was built from.

st4a <- function() {
  spec <- table_init("ST4a")

  # Nominal p, not adjusted. The adjusted value is reported alongside, so a reader
  # can apply a stricter cut, but the row set is the nominally significant one.
  out <- FCM_ORA %>%
    dplyr::filter(.data$p_value < 0.05)

  missing <- setdiff(spec$columns, names(out))
  if (length(missing) > 0) {
    stop("ST4a: FCM_ORA has no column(s) ", paste(missing, collapse = ", "),
         "\n  It carries: ", paste(names(FCM_ORA), collapse = ", "),
         "\n  Reconcile `columns` for ST4a in config/table_map.json with the object.",
         call. = FALSE)
  }

  # Tissue, then assay in the shared modality order, then cluster. Sorting is
  # stable, so rows within one cluster keep the order the object stores them in,
  # which is by ascending p value.
  out <- out[order(out$tissue,
                   match(out$assay, ASSAY_DISPLAY_ORDER),
                   out$cluster), , drop = FALSE]

  unplaced <- setdiff(unique(out$assay), ASSAY_DISPLAY_ORDER)
  if (length(unplaced) > 0) {
    stop("ST4a: assay(s) with no place in ASSAY_DISPLAY_ORDER: ",
         paste(sort(unplaced), collapse = ", "),
         "\n  Add them to lib/supp_table_helpers.R.", call. = FALSE)
  }

  # set_id is dropped: it is MSigDB's numeric handle for the set, and `set` and
  # `set_short` already name it.
  out <- as.data.frame(out)[, spec$columns, drop = FALSE]

  export_table(out, "ST4a")
}

run_tables(list(ST4a = st4a))
